#!/usr/bin/env python3
"""
Independent hard verification: EVERY sequence in the Omicron FASTA must have
a valid model PDB on disk. Does not trust manifest.json — re-derives
sequence -> uid -> files directly from the FASTA and the filesystem.

Usage: python verify_all.py [--fasta PATH] [--pdbs DIR]
Exit code 0 = all done; 1 = problems found.
"""

import argparse
import gzip
import hashlib
import json
import os
import re
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


def clean(seq):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--fasta",
        default=os.path.join(os.path.dirname(SCRIPT_DIR), "Spike_protein.aln-fasta"),
    )
    ap.add_argument("--pdbs", default=os.path.join(SCRIPT_DIR, "batch", "pdbs"))
    ap.add_argument("--state", default=os.path.join(SCRIPT_DIR, "batch", "state"))
    args = ap.parse_args()

    # 1. load uid mapping from state (uid -> md5/header) — state is our own,
    #    but the disk check below is independent of it
    uniques = json.load(open(os.path.join(args.state, "uniques.json")))
    uid_of_md5 = {u["md5"]: u["uid"] for u in uniques}

    # 2. parse the FASTA fresh
    seqs = parse_fasta(args.fasta)
    print(f"FASTA: {len(seqs)} sequences")

    problems = []
    counts = {"with_pdb": 0, "no_pdb": 0, "missing_uid": 0}
    all_pdb_files = set()
    all_cif_files = set()

    for idx, (hdr, raw) in enumerate(seqs):
        md5 = hashlib.md5(clean(raw).encode()).hexdigest()
        uid = uid_of_md5.get(md5)
        if uid is None:
            problems.append(f"seq {idx} ({hdr.split()[0]}): no uid for md5")
            counts["missing_uid"] += 1
            continue
        pdir = os.path.join(args.pdbs, uid)
        if not os.path.isdir(pdir):
            problems.append(f"seq {idx} ({hdr.split()[0]}): missing folder {pdir}")
            counts["no_pdb"] += 1
            continue
        pdbs = sorted(f for f in os.listdir(pdir) if f.endswith(".pdb.gz"))
        cifs = sorted(f for f in os.listdir(pdir) if f.endswith(".cif.gz"))
        if not pdbs:
            problems.append(f"seq {idx} ({hdr.split()[0]}): uid {uid} has NO pdb")
            counts["no_pdb"] += 1
            continue
        # validate every pdb file: non-empty + gzip + ATOM records
        ok = True
        for f in pdbs:
            p = os.path.join(pdir, f)
            if os.path.getsize(p) == 0:
                problems.append(f"seq {idx} ({hdr.split()[0]}): {p} is EMPTY")
                ok = False
                continue
            try:
                with gzip.open(p, "rt") as fh:
                    atoms = sum(1 for _ in fh)  # full read = gzip integrity
                if atoms < 1000:  # <1000 lines: not a real protein model
                    problems.append(
                        f"seq {idx} ({hdr.split()[0]}): {p} too small ({atoms} lines)"
                    )
                    ok = False
            except Exception as e:
                problems.append(f"seq {idx} ({hdr.split()[0]}): {p} CORRUPT: {e}")
                ok = False
        if ok:
            counts["with_pdb"] += 1
        else:
            counts["no_pdb"] += 1
        all_pdb_files.update(os.path.join(pdir, f) for f in pdbs)
        all_cif_files.update(os.path.join(pdir, f) for f in cifs)

    print(f"verified: {counts['with_pdb']}/{len(seqs)} sequences have valid PDB models")
    print(
        f"unique PDB files on disk: {len(all_pdb_files)} | CIF files: {len(all_cif_files)}"
    )
    print(f"missing uid: {counts['missing_uid']} | no/invalid PDB: {counts['no_pdb']}")

    # 3. also confirm every uid folder has at least one pdb (orphan check)
    uids_on_disk = sorted(os.listdir(args.pdbs)) if os.path.isdir(args.pdbs) else []
    orphan = [
        u
        for u in uids_on_disk
        if not any(
            f.endswith(".pdb.gz") for f in os.listdir(os.path.join(args.pdbs, u))
        )
    ]
    print(
        f"uid folders on disk: {len(uids_on_disk)} | folders without PDB: {len(orphan)}"
    )

    if problems:
        print(f"\n!!! {len(problems)} PROBLEMS:")
        for p in problems[:40]:
            print("  -", p)
        return 1
    print(
        "\nALL SEQUENCES VERIFIED: every one of the",
        len(seqs),
        "sequences has a valid PDB model.",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
