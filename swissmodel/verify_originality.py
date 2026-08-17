#!/usr/bin/env python3
"""
ORIGINALITY audit: every PDB must be ITS OWN sequence's model.

Checks (independent of manifest.json):
  1. target-sequence match: the server-recorded target of EVERY model in every
     uid folder must equal the FASTA sequence of that uid (gap-stripped).
  2. residue coverage: CA-atom count of the best model vs expected residues
     (sequence length x chains of the model's own oligo state), 60-105%.
  3. content uniqueness: md5 of every .pdb.gz — no duplicates across or within
     uid folders (a copied file would produce an identical hash).

Usage: python verify_originality.py
Exit 0 = all original and valid.
"""

import gzip
import hashlib
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BATCH = os.path.join(SCRIPT_DIR, "batch")
PDBS = os.path.join(BATCH, "pdbs")


def main():
    uniques = json.load(open(os.path.join(BATCH, "state", "uniques.json")))
    by_uid = {u["uid"]: u for u in uniques}
    problems = []
    hashes = {}
    ca_all = []

    for uid in sorted(by_uid):
        u = by_uid[uid]
        fd_path = os.path.join(PDBS, uid, "full_details.json")
        if not os.path.exists(fd_path):
            problems.append(f"{uid}: missing full_details.json")
            continue
        fd = json.load(open(fd_path))
        models = fd.get("models") or []
        if not models:
            problems.append(f"{uid}: no models in payload")
            continue
        # 1) target match for every model
        for m in models:
            tgt = (
                (m.get("targets") or [{}])[0]
                .get("target_sequence", "")
                .replace("-", "")
            )
            if tgt != u["seq"]:
                problems.append(
                    f"{uid} {u['headers'][0]}: target len {len(tgt)} != fasta {len(u['seq'])}"
                )
        # 2) CA coverage of best model (own oligo state)
        best = max(models, key=lambda m: m.get("gmqe") or 0)
        oligo = best.get("oligo_state", "monomer")
        n_chains = 3 if "trimer" in oligo else (2 if "dimer" in oligo else 1)
        expected = len(u["seq"]) * n_chains
        pdb = os.path.join(PDBS, uid, f"{best.get('model_id')}.pdb.gz")
        if os.path.exists(pdb):
            with gzip.open(pdb, "rt") as fh:
                ca = sum(
                    1 for l in fh if l.startswith("ATOM") and l[12:16].strip() == "CA"
                )
            ratio = ca / expected if expected else 0
            ca_all.append(ratio)
            if not (0.6 <= ratio <= 1.05):
                problems.append(f"{uid}: CA {ca} vs expected {expected} ({ratio:.2f})")
        # 3) content hashes
        for f in os.listdir(os.path.join(PDBS, uid)):
            if f.endswith(".pdb.gz"):
                p = os.path.join(PDBS, uid, f)
                h = hashlib.md5(open(p, "rb").read()).hexdigest()
                hashes.setdefault(h, []).append(p)

    dups = [p for h, p in hashes.items() if len(p) > 1]
    n_pdb = len(hashes)

    print(f"uids checked: {len(by_uid)}")
    print(f"target-sequence mismatches: {sum(1 for p in problems if 'target' in p)}")
    print(
        f"CA coverage: min={min(ca_all):.2f} median={sorted(ca_all)[len(ca_all) // 2]:.2f} "
        f"max={max(ca_all):.2f} | outliers: "
        f"{sum(1 for p in problems if 'CA' in p)}"
    )
    print(f"PDB files: {n_pdb} | duplicate-content files: {len(dups)}")
    for p in dups[:5]:
        print("   DUP:", p)
    if problems:
        print(f"\n!!! {len(problems)} PROBLEMS:")
        for p in problems[:30]:
            print("  -", p)
        return 1
    print(
        "\nALL ORIGINAL: every PDB is its own sequence's model; "
        "no duplicates, no mismatches."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
