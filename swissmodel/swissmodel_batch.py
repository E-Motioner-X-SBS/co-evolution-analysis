#!/usr/bin/env python3
"""
Batch SWISS-MODEL automodel pipeline for the Omicron Spike dataset.

Ships 299 UNIQUE sequences (dedup of the 1,299 in Spike_protein.aln-fasta),
polls each project to COMPLETED/FAILED, downloads model PDB/CIF files into a
proper folder structure, and produces a verification manifest covering ALL
1,299 original sequences.

Rate limiting (verified from https://swissmodel.expasy.org/docs/help#modelling_api):
  - rapid submission 200/min, prolonged 10000/6h, excess -> 429
  - status checks share the same budget (we cap total at ~90/min and back off on 429)
Server-side caching: re-submitting a seen sequence returns 200 + same project_id.

State (crash-safe, resumable) lives in batch/state/*.json. Run phases:
  python swissmodel_batch.py run          # submit -> poll -> download until all done
  python swissmodel_batch.py submit       # submit phase only
  python swissmodel_batch.py poll         # poll + download loop (one pass per invocation)
  python swissmodel_batch.py verify       # build manifest.json for all 1299 sequences
"""

import argparse
import json
import os
import re
import sys
import time

import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
FASTA = os.path.join(REPO_ROOT, "Spike_protein.aln-fasta")
BATCH = os.path.join(SCRIPT_DIR, "batch")
STATE = os.path.join(BATCH, "state")
PDB_DIR = os.path.join(BATCH, "pdbs")
LOG = os.path.join(BATCH, "orchestrator.log")

API = "https://swissmodel.expasy.org"
TOKEN = None

RAPID_LIMIT = 200  # per minute (docs)
SAFE_SUBMIT_PER_MIN = 30  # submissions (polite)
SAFE_POLL_PER_MIN = 90  # status checks (leave headroom)
MAX_429_BACKOFF = 300  # seconds
MAX_FAILED_RETRIES = 1


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def load_token():
    global TOKEN
    if TOKEN:
        return TOKEN
    t = os.environ.get("SWISSMODEL_API_TOKEN")
    if not t:
        env_path = os.path.join(SCRIPT_DIR, ".env")
        if os.path.exists(env_path):
            for line in open(env_path):
                line = line.strip()
                if line.startswith("SWISSMODEL_API_TOKEN="):
                    t = line.split("=", 1)[1].strip()
    TOKEN = t
    return t


def state_path(name):
    return os.path.join(STATE, name)


def load_state(name, default):
    p = state_path(name)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return default


def save_state(name, data):
    p = state_path(name)
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=1)
    os.replace(tmp, p)  # atomic


# ──────────────────────────── sequences ────────────────────────────


def parse_fasta(path):
    seqs = []
    h = None
    s = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if h:
                seqs.append((h, "".join(s)))
            h = line[1:]
            s = []
        elif line:
            s.append(line.upper())
    if h:
        seqs.append((h, "".join(s)))
    return seqs


def build_uniques():
    """Return list of unique-sequence records + per-original-index mapping."""
    raw = parse_fasta(FASTA)
    by_md5 = {}
    order = []
    for idx, (hdr, s0) in enumerate(raw):
        clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s0)
        md5 = __import__("hashlib").md5(clean.encode()).hexdigest()
        if md5 not in by_md5:
            by_md5[md5] = {
                "md5": md5,
                "seq": clean,
                "len": len(clean),
                "headers": [hdr.split()[0]],
                "indices": [],
            }
            order.append(md5)
        by_md5[md5]["indices"].append(idx)
    for u in by_md5.values():
        u["uid"] = f"u{order.index(u['md5']):04d}"
    return [by_md5[m] for m in order]


def seed_state():
    """(Re)create state from the FASTA if not present. Safe to re-run."""
    if not os.path.exists(state_path("uniques.json")):
        uniques = build_uniques()
        save_state("uniques.json", uniques)
        save_state(
            "orig_map.json", {str(i): u["uid"] for u in uniques for i in u["indices"]}
        )
        log(
            f"seeded {len(uniques)} unique sequences (covering "
            f"{len(open(FASTA).readlines())} raw lines)"
        )
    return load_state("uniques.json", [])


# ──────────────────────────── submission ────────────────────────────


def api_post(url, payload, retries=5):
    for attempt in range(retries):
        try:
            r = requests.post(
                url,
                headers={"Authorization": f"Token {load_token()}"},
                json=payload,
                timeout=90,
            )
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", "60"))
                log(f"429 on POST: backing off {wait}s (attempt {attempt + 1})")
                time.sleep(min(wait, MAX_429_BACKOFF))
                continue
            return r
        except requests.RequestException as e:
            log(f"POST error {e} (attempt {attempt + 1})")
            time.sleep(10 * (attempt + 1))
    return None


def submit_all():
    uniques = load_state("uniques.json", [])
    submitted = load_state("submitted.json", {})
    done = load_state("done.json", {})
    failed = load_state("failed.json", {})
    # retry policy: transient failures get MAX_FAILED_RETRIES re-submissions
    retryable = {
        md5: f for md5, f in failed.items() if f.get("retries", 0) < MAX_FAILED_RETRIES
    }
    for md5, f in retryable.items():
        log(
            f"retrying failed {f.get('uid')} {f.get('header')} (retries={f.get('retries', 0)})"
        )
        del failed[md5]
        if md5 in submitted:
            del submitted[md5]
    save_state("failed.json", failed)
    save_state("submitted.json", submitted)
    pending = [
        u
        for u in uniques
        if u["md5"] not in submitted and u["md5"] not in done and u["md5"] not in failed
    ]
    log(f"submit: {len(pending)} to submit, {len(submitted)} already submitted")
    for u in pending:
        r = api_post(
            f"{API}/automodel",
            {
                "target_sequences": u["seq"],
                "project_title": f"E-motioner {u['uid']} ({u['headers'][0]})",
            },
        )
        if r is None:
            log(f"submit FAILED (network) for {u['uid']} — will retry next run")
            break
        if r.status_code in (200, 201, 202):
            submitted[u["md5"]] = r.json().get("project_id")
            save_state("submitted.json", submitted)
            code = "cached" if r.status_code == 200 else f"accepted({r.status_code})"
            log(
                f"submitted {u['uid']} {u['headers'][0]} -> {submitted[u['md5']]} ({code})"
            )
        else:
            log(f"submit ERROR {r.status_code} for {u['uid']}: {r.text[:200]}")
            failed[u["md5"]] = {"error": r.text[:500]}
            save_state("failed.json", failed)
            if r.status_code in (400, 401, 403):
                log("fatal auth/validation error — aborting submit loop")
                break
        time.sleep(60.0 / SAFE_SUBMIT_PER_MIN)  # pace: 30/min
    save_state("submitted.json", submitted)
    log(f"submit done: {len(submitted)} submitted")


# ──────────────────────────── polling + download ────────────────────────────


def api_get(url, retries=5):
    for attempt in range(retries):
        try:
            r = requests.get(
                url, headers={"Authorization": f"Token {load_token()}"}, timeout=120
            )
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", "60"))
                log(f"429 on GET: backing off {wait}s")
                time.sleep(min(wait, MAX_429_BACKOFF))
                continue
            return r
        except requests.RequestException as e:
            log(f"GET error {e} (attempt {attempt + 1})")
            time.sleep(10 * (attempt + 1))
    return None


def download_project(u, payload):
    """Download all model coordinates for a completed project."""
    outdir = os.path.join(PDB_DIR, u["uid"])
    os.makedirs(outdir, exist_ok=True)
    files = []
    models = payload.get("models") or []
    for m in models:
        for key in ("coordinates_url", "modelcif_url"):
            url = m.get(key)
            if not url:
                continue
            fname = os.path.basename(url)
            dest = os.path.join(outdir, fname)
            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                files.append(dest)
                continue
            try:
                r = requests.get(url, timeout=300)
                if r.ok:
                    with open(dest, "wb") as f:
                        f.write(r.content)
                    files.append(dest)
                    log(f"  downloaded {dest} ({len(r.content)} bytes)")
                else:
                    log(f"  download FAILED {r.status_code} {url}")
            except requests.RequestException as e:
                log(f"  download EXC {e} {url}")
        # keep model metadata
        meta = {
            k: v for k, v in m.items() if k not in ("coordinates_url", "modelcif_url")
        }
        with open(os.path.join(outdir, f"model_{m.get('model_id')}.json"), "w") as f:
            json.dump(meta, f, indent=1)
    # store raw payload
    with open(os.path.join(outdir, "full_details.json"), "w") as f:
        json.dump(payload, f, indent=1)
    # per-model match summary (the "proper Omicron match data")
    summary = []
    for m in models:
        tpl = None
        for tg in m.get("targets") or []:
            for ch in tg.get("model_chains") or []:
                for t in ch.get("templates") or []:
                    tpl = t
                    break
        qmean = m.get("qmean_global") or {}
        summary.append(
            {
                "model_id": m.get("model_id"),
                "gmqe": m.get("gmqe"),
                "qsqe": m.get("qsqe"),
                "qmean4_z": qmean.get("qmean4_z_score"),
                "qmean_global": qmean.get("qmean_global"),
                "oligo_state": m.get("oligo_state"),
                "coverage": m.get("coverage"),
                "sequence_identity": m.get("sequence_identity"),
                "sequence_similarity": m.get("sequence_similarity"),
                "template_qualified_name": tpl.get("qualified_name") if tpl else None,
                "template_found_by": tpl.get("found_by") if tpl else None,
                "template_method": tpl.get("method") if tpl else None,
                "template_gmqe": tpl.get("gmqe") if tpl else None,
                "template_title": (tpl.get("title") if tpl else None),
            }
        )
    with open(os.path.join(outdir, "templates_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    return files


def poll_once():
    uniques = load_state("uniques.json", [])
    submitted = load_state("submitted.json", {})
    done = load_state("done.json", {})
    failed = load_state("failed.json", {})
    unfinished = [
        u
        for u in uniques
        if u["md5"] in submitted and u["md5"] not in done and u["md5"] not in failed
    ]
    if not unfinished:
        log("poll: nothing unfinished")
        return 0
    log(f"poll: {len(unfinished)} unfinished projects")
    interval = max(1.0, 60.0 / SAFE_POLL_PER_MIN * len(unfinished) / 299.0)
    n = 0
    for u in unfinished:
        r = api_get(f"{API}/project/{submitted[u['md5']]}/models/full-details/")
        if r is None:
            log(f"poll network fail for {u['uid']} — continue")
            continue
        if r.status_code != 200:
            log(f"poll HTTP {r.status_code} for {u['uid']} {r.text[:150]}")
            continue
        payload = r.json()
        status = payload.get("status", "UNKNOWN")
        n += 1
        if status in ("COMPLETED", "FAILED"):
            if status == "COMPLETED":
                files = download_project(u, payload)
                done[u["md5"]] = {
                    "uid": u["uid"],
                    "header": u["headers"][0],
                    "project_id": submitted[u["md5"]],
                    "n_models": len(payload.get("models") or []),
                    "files": files,
                    "completed_at": time.time(),
                }
                save_state("done.json", done)
                log(
                    f"DONE {u['uid']} {u['headers'][0]} models={len(payload.get('models') or [])}"
                )
            else:
                failed[u["md5"]] = {
                    "uid": u["uid"],
                    "header": u["headers"][0],
                    "project_id": submitted[u["md5"]],
                    "message": payload.get("message", ""),
                    "retries": 0,
                }
                save_state("failed.json", failed)
                log(
                    f"FAILED {u['uid']} {u['headers'][0]}: {payload.get('message', '')[:200]}"
                )
        time.sleep(interval)
    return n


# ──────────────────────────── verification ────────────────────────────


def verify():
    uniques = load_state("uniques.json", [])
    orig_map = load_state("orig_map.json", {})
    done = load_state("done.json", {})
    failed = load_state("failed.json", {})
    submitted = load_state("submitted.json", {})
    manifest = {
        "n_original": len(orig_map),
        "n_unique": len(uniques),
        "summary": {
            "done": len(done),
            "failed": len(failed),
            "submitted": len(submitted),
        },
        "sequences": {},
    }
    for i in sorted(orig_map, key=int):
        uid = orig_map[i]
        u = next((x for x in uniques if x["uid"] == uid), None)
        entry = {"index": int(i), "uid": uid, "header": u["headers"][0] if u else None}
        md5 = u["md5"] if u else None
        if md5 in done:
            d = done[md5]
            entry["status"] = "done"
            entry["project_id"] = d["project_id"]
            entry["files"] = [f for f in d.get("files", []) if os.path.exists(f)]
            entry["pdb_present"] = any(
                f.endswith(".pdb") or f.endswith(".pdb.gz") for f in entry["files"]
            )
            entry["model_pdb"] = next(
                (
                    f
                    for f in d.get("files", [])
                    if f.endswith(".pdb") or f.endswith(".pdb.gz")
                ),
                None,
            )
        elif md5 in failed:
            entry["status"] = "failed"
            entry["message"] = failed[md5].get("message", "")
        elif md5 in submitted:
            entry["status"] = "pending"
            entry["project_id"] = submitted[md5]
        else:
            entry["status"] = "not_submitted"
        manifest["sequences"][str(i)] = entry
    out = os.path.join(BATCH, "manifest.json")
    with open(out, "w") as f:
        json.dump(manifest, f, indent=1)
    ok = sum(1 for e in manifest["sequences"].values() if e["status"] == "done")
    have_pdb = sum(1 for e in manifest["sequences"].values() if e.get("pdb_present"))
    log(
        f"verify: {ok}/{len(orig_map)} original sequences done; "
        f"{have_pdb} with PDB on disk; failed={len(failed)}; "
        f"pending_unique={len(submitted) - len(done) - len(failed)}"
    )
    return manifest


# ──────────────────────────── main ────────────────────────────


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("phase", choices=["run", "submit", "poll", "verify"])
    ap.add_argument("--max-poll-rounds", type=int, default=0, help="0 = unlimited")
    ap.add_argument("--max-hours", type=float, default=168)
    args = ap.parse_args()

    os.makedirs(STATE, exist_ok=True)
    os.makedirs(PDB_DIR, exist_ok=True)
    if not load_token():
        sys.exit("no token found (swissmodel/.env or SWISSMODEL_API_TOKEN)")

    uniques = seed_state()
    log(f"phase={args.phase} unique_sequences={len(uniques)}")

    if args.phase in ("run", "submit"):
        submit_all()
    if args.phase == "run":
        t0 = time.time()
        rounds = 0
        while time.time() - t0 < args.max_hours * 3600:
            verify()
            submitted = load_state("submitted.json", {})
            done = load_state("done.json", {})
            failed = load_state("failed.json", {})
            remaining = [
                u
                for u in uniques
                if u["md5"] in submitted
                and u["md5"] not in done
                and u["md5"] not in failed
            ]
            if not remaining:
                log("ALL UNIQUE SEQUENCES COMPLETE")
                break
            poll_once()
            rounds += 1
            if args.max_poll_rounds and rounds >= args.max_poll_rounds:
                break
        verify()
    elif args.phase == "poll":
        poll_once()
    elif args.phase == "verify":
        verify()
    log("done")


if __name__ == "__main__":
    main()
