#!/usr/bin/env python3
"""
Materialize ONE folder per FASTA ENTRY (1,299) under batch/pdbs_by_entry/.

The models themselves are stored once per UNIQUE sequence (batch/pdbs/u####/,
299 folders — the server builds one model per distinct sequence and caches
duplicates, so identical FASTA entries share a model). This script gives every
original FASTA entry its own folder containing hardlinks to its model files,
so downstream pipelines can process entry-by-entry without dedup logic.

Usage: python materialize_entries.py [--fasta X] [--pdbs DIR] [--out DIR] [--state DIR]
Exit 0 = all 1299 entry folders created.
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--fasta",
        default=os.path.join(os.path.dirname(SCRIPT_DIR), "Spike_protein.aln-fasta"),
    )
    ap.add_argument("--pdbs", default=os.path.join(SCRIPT_DIR, "batch", "pdbs"))
    ap.add_argument("--state", default=os.path.join(SCRIPT_DIR, "batch", "state"))
    ap.add_argument("--out", default=os.path.join(SCRIPT_DIR, "batch", "pdbs_by_entry"))
    args = ap.parse_args()

    uniques = json.load(open(os.path.join(args.state, "uniques.json")))
    uid_of_md5 = {u["md5"]: u["uid"] for u in uniques}
    submitted = json.load(open(os.path.join(args.state, "submitted.json")))
    done = json.load(open(os.path.join(args.state, "done.json")))

    seqs = parse_fasta(args.fasta)
    os.makedirs(args.out, exist_ok=True)

    n_ok = 0
    problems = []
    for idx, (hdr, raw) in enumerate(seqs):
        acc = hdr.split()[0]
        md5 = hashlib.md5(
            re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", raw).encode()
        ).hexdigest()
        uid = uid_of_md5.get(md5)
        if uid is None:
            problems.append(f"{idx} {acc}: no uid")
            continue
        src = os.path.join(args.pdbs, uid)
        dest = os.path.join(args.out, f"entry_{idx:04d}_{acc}")
        os.makedirs(dest, exist_ok=True)
        linked = []
        for f in sorted(os.listdir(src)):
            if f.endswith((".pdb.gz", ".cif.gz")):
                sp, dp = os.path.join(src, f), os.path.join(dest, f)
                if not os.path.exists(dp):
                    try:
                        os.link(sp, dp)  # hardlink: zero extra disk
                    except OSError:
                        shutil.copy2(sp, dp)  # cross-device fallback
                linked.append(f)
        # copy the match summary too
        ss = os.path.join(src, "templates_summary.json")
        if os.path.exists(ss) and not os.path.exists(
            os.path.join(dest, "templates_summary.json")
        ):
            try:
                os.link(ss, os.path.join(dest, "templates_summary.json"))
            except OSError:
                shutil.copy2(ss, os.path.join(dest, "templates_summary.json"))
        meta = {
            "index": idx,
            "accession": acc,
            "full_header": hdr,
            "uid": uid,
            "md5": md5,
            "project_id": submitted.get(md5),
            "n_models": len(done.get(md5, {}).get("files", [])),
            "pdb_files": [f for f in linked if f.endswith(".pdb.gz")],
        }
        json.dump(meta, open(os.path.join(dest, "entry.json"), "w"), indent=1)
        if linked:
            n_ok += 1
        else:
            problems.append(f"{idx} {acc}: no files linked")

    print(f"entry folders with files: {n_ok}/{len(seqs)}")
    if problems:
        print(f"PROBLEMS ({len(problems)}):")
        for p in problems[:20]:
            print("  -", p)
        return 1
    print(f"ALL {len(seqs)} ENTRY FOLDERS MATERIALIZED under {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
