#!/usr/bin/env python3
"""Fetch representative native PDB structures for contact-map validation.

PSICOV ships its own native PDBs (data/psicov150/raw/pdb/). For the other
datasets, fetch from RCSB: Pfam reps (6TNE, 6AK2, 1HE1), GPCR (7E2X, 2RH1),
unit tests (4P3R, 1WHZ, 1ATZ). DeepMSA natives ship with the dataset.
Resumable + retried; results logged to data/structures/fetch.log.
"""

import os
import sys
import time
from pathlib import Path

import requests

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = DATA / "structures"
IDS = ["6TNE", "6AK2", "1HE1", "7E2X", "2RH1", "4P3R", "1WHZ", "1ATZ"]


def main():
    os.makedirs(OUT, exist_ok=True)
    logp = OUT / "fetch.log"
    for pdb in IDS:
        dest = OUT / f"{pdb.lower()}.pdb"
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"[cached] {pdb}")
            continue
        url = f"https://files.rcsb.org/download/{pdb}.pdb"
        ok = False
        for attempt in range(3):
            try:
                r = requests.get(url, timeout=120)
                if r.ok and len(r.text) > 1000:
                    dest.write_text(r.text)
                    print(f"[OK] {pdb}: {len(r.text)} bytes")
                    ok = True
                    break
                print(f"[retry {attempt}] {pdb}: http {r.status_code}")
            except Exception as e:
                print(f"[retry {attempt}] {pdb}: {e}")
            time.sleep(5)
        with open(logp, "a") as f:
            f.write(f"{time.strftime('%H:%M:%S')} {pdb} {'OK' if ok else 'FAILED'}\n")
        if not ok:
            print(f"[FAILED] {pdb}")
    # also count PSICOV natives
    ps = DATA / "psicov150" / "raw" / "pdb"
    if ps.exists():
        print(f"psicov native pdbs: {len(list(ps.glob('*.pdb')))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
