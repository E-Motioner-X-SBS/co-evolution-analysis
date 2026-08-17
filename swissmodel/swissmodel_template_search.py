#!/usr/bin/env python3
"""
SWISS-MODEL template search for sequences from the Omicron Spike FASTA.

Two modes:
  --mode coreapi   Official Modelling API (token REQUIRED; runs full automodel
                   pipeline incl. model building). Token from
                   https://swissmodel.expasy.org/account  → env SWISSMODEL_API_TOKEN
                   or --token.
  --mode guest     Anonymous web-form path (no account needed). POST /interactive
                   with target=sequence, poll /interactive/{project}/templates/
                   until the template search finishes, parse embedded results.

Usage:
  python swissmodel/swissmodel_template_search.py --fasta Spike_protein.aln-fasta --seq-index 0 --mode guest
  python swissmodel/swissmodel_template_search.py --fasta X.fa --mode coreapi   # token from swissmodel/.env
"""

import argparse
import json
import os
import re
import sys
import time

import requests

UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Chrome/126 Safari/537.36"}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)  # co-evolution-analysis/


def load_token():
    """Token from env SWISSMODEL_API_TOKEN, else from swissmodel/.env (gitignored)."""
    token = os.environ.get("SWISSMODEL_API_TOKEN")
    if token:
        return token
    env_path = os.path.join(SCRIPT_DIR, ".env")
    if os.path.exists(env_path):
        for line in open(env_path):
            line = line.strip()
            if line.startswith("SWISSMODEL_API_TOKEN="):
                return line.split("=", 1)[1].strip()
    return None


def default_fasta():
    return os.path.join(REPO_ROOT, "Spike_protein.aln-fasta")


def read_fasta(fasta_path, seq_index=0):
    """Parse FASTA, return (header, gap-stripped protein sequence) at index."""
    seqs = []
    h = None
    s = []
    for line in open(fasta_path):
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
    h0, s0 = seqs[seq_index]
    clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s0)
    return h0, clean


# ──────────────────────────── GUEST MODE ────────────────────────────


def guest_template_search(sequence, project_title, poll_interval=20, max_wait_s=1800):
    """Anonymous template search via the /interactive form. Returns project id."""
    s = requests.Session()
    s.headers.update(UA)
    r0 = s.get("https://swissmodel.expasy.org/interactive/", timeout=90)
    m = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', r0.text)
    if not m:
        raise RuntimeError("could not get CSRF token")
    csrf = m.group(1)
    data = {
        "csrfmiddlewaretoken": csrf,
        "target": sequence,
        "project_title": project_title,
        "email": "",
        "automodel": "false",
        "is_alignment": "false",
    }
    r = s.post("https://swissmodel.expasy.org/interactive", data=data, timeout=300)
    m = re.search(r"/interactive/([A-Za-z0-9]+)/templates/", r.url)
    if not m:
        errs = re.findall(r'errorlist">([^<]+)<', r.text)
        raise RuntimeError(f"submission failed: {errs or r.status_code}")
    pid = m.group(1)
    print(f"[guest] project created: {pid}")
    t0 = time.time()
    while time.time() - t0 < max_wait_s:
        page = s.get(
            f"https://swissmodel.expasy.org/interactive/{pid}/templates/", timeout=120
        ).text
        text = re.sub(r"<script.*?</script>", " ", page, flags=re.S)
        text = re.sub(r"<[^>]+>", " ", text)
        text = re.sub(r"\s+", " ", text)
        m_st = re.search(r"currently ([\w ]+?)\.", text)
        st = m_st.group(1).lower() if m_st else ""
        if "queueing" not in st and "running" not in st:
            return pid, page
        time.sleep(poll_interval)
    raise TimeoutError(f"template search not finished after {max_wait_s}s")


def parse_guest_templates(html):
    """Extract the full template list from `globalTemplateData` embedded in the page."""
    m = re.search(r"const globalTemplateData = (\[.*?\]);", html, re.S)
    if not m:
        return []
    try:
        return json.loads(m.group(1))
    except json.JSONDecodeError:
        return []


# ──────────────────────────── COREAPI MODE ────────────────────────────


def coreapi_automodel(sequence, token, project_title):
    """Submit automodel job; returns project_id."""
    r = requests.post(
        "https://swissmodel.expasy.org/automodel",
        headers={"Authorization": f"Token {token}"},
        json={"target_sequences": sequence, "project_title": project_title},
        timeout=90,
    )
    if r.status_code not in (200, 201, 202):
        raise RuntimeError(f"automodel failed: {r.status_code} {r.text[:500]}")
    pid = r.json()["project_id"]
    print(f"[coreapi] project created: {pid} (http {r.status_code})")
    return pid


def coreapi_poll(pid, token, poll_interval=15, max_wait_s=3600):
    """Poll full-details until COMPLETED/FAILED; returns JSON payload."""
    t0 = time.time()
    while time.time() - t0 < max_wait_s:
        r = requests.get(
            f"https://swissmodel.expasy.org/project/{pid}/models/full-details/",
            headers={"Authorization": f"Token {token}"},
            timeout=120,
        )
        r.raise_for_status()
        data = r.json()
        st = data.get("status", "UNKNOWN")
        print(f"[coreapi] status: {st}")
        if st in ("COMPLETED", "FAILED"):
            return data
        time.sleep(poll_interval)
    raise TimeoutError(f"job not finished after {max_wait_s}s")


# ──────────────────────────────── MAIN ────────────────────────────────


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--fasta",
        default=default_fasta(),
        help="Omicron Spike aligned FASTA (default: repo Spike_protein.aln-fasta)",
    )
    ap.add_argument("--seq-index", type=int, default=0)
    ap.add_argument("--mode", choices=["guest", "coreapi"], default="guest")
    ap.add_argument("--token", default=load_token())
    ap.add_argument("--out", default=None, help="JSON output path")
    ap.add_argument("--max-wait-s", type=int, default=3600)
    args = ap.parse_args()

    header, seq = read_fasta(args.fasta, args.seq_index)
    print(f"sequence {args.seq_index}: {header.split()[0]}, {len(seq)} aa")
    title = f"E-motioner co-evolution seq{args.seq_index} ({header.split()[0]})"

    if args.mode == "guest":
        pid, html = guest_template_search(seq, title, max_wait_s=args.max_wait_s)
        open(f"/tmp/opencode/guest_{pid}.html", "w").write(html)
        templates = parse_guest_templates(html)
        summary = [
            {
                "rank": t.get("rank"),
                "qname": t.get("qname"),
                "pdb_id": t.get("pdb_id"),
                "chain": t.get("chain"),
                "method": t.get("short_method") or t.get("method"),
                "seq_identity": t.get("seq_id"),
                "coverage": t.get("coverage"),
                "seq_similarity": t.get("seq_sim"),
                "found_by": t.get("found_by"),
                "resolution_A": t.get("resolution"),
                "oligo_state": t.get("oligo_state"),
                "qsqe": t.get("qsqe"),
                "gmqe_afdb_plddt": t.get("pred_lddt"),
                "title": t.get("title"),
            }
            for t in templates
        ]
        result = {
            "mode": "guest",
            "project_id": pid,
            "sequence": header.split()[0],
            "n_templates": len(summary),
            "templates": summary,
        }
    else:
        if not args.token:
            sys.exit(
                "coreapi mode needs --token or SWISSMODEL_API_TOKEN (see /account)"
            )
        pid = coreapi_automodel(seq, args.token, title)
        data = coreapi_poll(pid, args.token, max_wait_s=args.max_wait_s)
        result = {
            "mode": "coreapi",
            "project_id": pid,
            "sequence": header.split()[0],
            "status": data.get("status"),
            "raw": data,
        }
    if args.out:
        json.dump(result, open(args.out, "w"), indent=2)
        print(f"saved -> {args.out}")
    else:
        print(json.dumps(result, indent=2)[:4000])


if __name__ == "__main__":
    main()
