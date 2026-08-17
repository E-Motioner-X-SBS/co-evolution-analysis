#!/usr/bin/env python3
"""
PHASE 7: Contact-map validation — do co-evolving pairs/rules map to 3D contacts?

For each target with a native structure:
  1. contact map from PDB chain A: Cb-Cb < 8 A (Ca for Gly), |i-j| >= 6
  2. map MSA columns -> PDB residues via Needleman-Wunsch alignment of the
     MSA consensus/reference sequence against the PDB chain sequence
  3. score top co-evolving pairs (from master_boolean / full_length results):
       precision@K = fraction of top-K pairs that are contacts
       enrichment  = precision vs random-pair contact rate (same sep dist)
  4. rules: are the positions in the top Boolean rules in contact?

Output: results/contacts/<dataset>/<target>.json + consolidated table.
Usage: python scripts/validate_contacts.py [--dataset X --target Y] [--all]
"""

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
RESULTS = REPO / "results"
CONTACT_CUTOFF = 8.0
MIN_SEP = 6


def parse_pdb_chain(pdb_path, chain="A"):
    """Return (seq, coords) for the first chain: list of (aa, (x,y,z))."""
    residues = []
    last = None
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
            resi = int(line[22:26])
            name = line[12:16].strip()
            if name == "CA":
                aa3 = line[17:20]
                aa = {
                    "ALA": "A",
                    "CYS": "C",
                    "ASP": "D",
                    "GLU": "E",
                    "PHE": "F",
                    "GLY": "G",
                    "HIS": "H",
                    "ILE": "I",
                    "LYS": "K",
                    "LEU": "L",
                    "MET": "M",
                    "ASN": "N",
                    "PRO": "P",
                    "GLN": "Q",
                    "ARG": "R",
                    "SER": "S",
                    "THR": "T",
                    "VAL": "V",
                    "TRP": "W",
                    "TYR": "Y",
                }.get(aa3, "X")
                if resi != last:
                    residues.append(
                        [
                            aa,
                            [
                                float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54]),
                            ],
                        ]
                    )
                    last = resi
    return residues


def cb_coords(residues):
    """Cb coords: real Cb if present, else Ca (Gly). Rebuild atom list."""
    out = []
    for line_aa, _ in residues:
        pass
    # simpler: re-parse atoms properly
    return None


def contact_map(pdb_path, chain="A"):
    """Return (seq, contact set as frozenset of (i,j) 0-indexed pairs)."""
    # collect per-residue CA + CB
    res_atoms = {}  # resi -> {'CA': (x,y,z), 'CB': ...}
    for line in open(pdb_path):
        if line.startswith("ATOM") and (
            line[21] == chain or (chain == "A" and line[21].strip() == "")
        ):
            resi = int(line[22:26])
            name = line[12:16].strip()
            if name in ("CA", "CB"):
                res_atoms.setdefault(resi, {})[name] = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
    # residue order + sequence
    resis = sorted(res_atoms)
    aa3 = {}
    for line in open(pdb_path):
        if (
            line.startswith("ATOM")
            and (line[21] == chain or (chain == "A" and line[21].strip() == ""))
            and line[12:16].strip() == "CA"
        ):
            aa3[int(line[22:26])] = line[17:20]
    AA3 = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y",
    }
    seq = "".join(AA3.get(aa3.get(r, "XXX"), "X") for r in resis)
    contacts = set()
    n = len(resis)
    for a in range(n):
        for b in range(a + MIN_SEP, n):
            ca_a = res_atoms[resis[a]].get("CA")
            cb_a = res_atoms[resis[a]].get("CB", ca_a)
            ca_b = res_atoms[resis[b]].get("CA")
            cb_b = res_atoms[resis[b]].get("CB", ca_b)
            if cb_a is None or cb_b is None:
                continue
            d = sum((u - v) ** 2 for u, v in zip(cb_a, cb_b)) ** 0.5
            if d < CONTACT_CUTOFF:
                contacts.add((a, b))
    return seq, contacts


def nw_align(q, t):
    """Global alignment (Needleman-Wunsch, simple scoring) -> q_col->t_col map."""
    import numpy as np

    G = -2
    M = 2
    S = -1
    nq, nt = len(q), len(t)
    score = np.zeros((nq + 1, nt + 1))
    for i in range(1, nq + 1):
        score[i, 0] = i * G
    for j in range(1, nt + 1):
        score[0, j] = j * G
    for i in range(1, nq + 1):
        for j in range(1, nt + 1):
            diag = score[i - 1, j - 1] + (M if q[i - 1] == t[j - 1] else S)
            up = score[i - 1, j] + G
            left = score[i, j - 1] + G
            score[i, j] = max(diag, up, left)
    mapping = {}
    i, j = nq, nt
    while i > 0 and j > 0:
        if q[i - 1] == t[j - 1] and score[i, j] == score[i - 1, j - 1] + M:
            mapping[i - 1] = j - 1
            i -= 1
            j -= 1
        elif score[i, j] == score[i - 1, j - 1] + S:
            mapping[i - 1] = j - 1
            i -= 1
            j -= 1
        elif score[i, j] == score[i - 1, j] + G:
            i -= 1
        else:
            j -= 1
    return mapping


def a2m_reference_map(a2m_path, pdb_seq):
    """Map MSA columns -> PDB residues using the A2M reference frame.

    In A2M, the first (query) sequence is uppercase at match-state columns
    and lowercase at insert columns. We strip lowercase to get the reference
    protein, align it to the PDB chain, and return {msa_col: pdb_res}.
    Insert columns get NO mapping (correctly excluded from structure checks).
    """
    first = []
    for line in open(a2m_path):
        if line.startswith(">"):
            if first:
                break
            continue
        first.append(line.strip())
    s = "".join(first)
    ref_prot = "".join(ch.upper() for ch in s if ch.isalpha())
    aln = nw_align(ref_prot, pdb_seq)
    # build col map: walk match-state columns in order
    mapping = {}
    ref_idx = 0
    for col, ch in enumerate(s):
        if ch.isalpha():
            if ch.isupper():
                # match-state column: maps through the alignment
                if ref_idx in aln:
                    mapping[col] = aln[ref_idx]
                ref_idx += 1
            # lowercase = insert column -> no map
        # '-' = gap in query -> no map
    return mapping


def reference_seq(fasta_path):
    """First sequence of the MSA = reference anchor (A2M: query frame)."""
    cur = []
    for line in open(fasta_path):
        if line.startswith(">"):
            if cur:
                return "".join(cur).upper()
            cur = []
        else:
            cur.append(line.strip())
    return "".join(cur).upper() if cur else ""


def consensus_seq(fasta_path, n=200):
    """Most common residue per column (gaps excluded) over first n sequences."""
    import numpy as np
    from collections import Counter

    seqs = []
    cur = []
    for line in open(fasta_path):
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    seqs = seqs[:n]
    L = max(len(s) for s in seqs)
    cons = []
    for c in range(L):
        cnt = Counter(s[c] for s in seqs if c < len(s) and s[c] != "-")
        if cnt:
            cons.append(cnt.most_common(1)[0][0])
        else:
            cons.append("-")
    return "".join(cons)


def load_pairs_from_results(target_dir):
    """Top co-evolving pairs from master_boolean_summary.json or full_length summary."""
    pairs = []
    mb = Path(target_dir) / "master_boolean" / "master_boolean_summary.json"
    if mb.exists():
        d = json.load(open(mb))
        for p in d.get("top_co_evolving_pairs", []):
            pairs.append((p["pos_i"], p["pos_j"], p.get("mi", 0)))
        return pairs
    fl = Path(target_dir) / "full_length_analysis" / "full_length_summary.json"
    if fl.exists():
        d = json.load(open(fl))
        for p in d.get("top_pairs", []):
            pairs.append((p["pos_i"], p["pos_j"], p.get("mi", 0)))
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset")
    ap.add_argument("--target")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    # structure lookup per target
    def find_pdb(ds, tgt):
        cands = [
            RESULTS / ds / tgt,
        ]
        # explicit structure dirs
        if ds == "unit_tests":
            m = {"DHFR": "4p3r", "1whzA": "1whz", "1atzA": "1atz"}
            return REPO / "data" / "structures" / f"{m.get(tgt, '')}.pdb"
        if ds == "gpcrdb":
            m = {"b2ar": "2rh1", "5ht1a": "7e2x"}
            return REPO / "data" / "structures" / f"{m.get(tgt, '')}.pdb"
        if ds == "psicov150":
            return REPO / "data" / "psicov150" / "raw" / "pdb" / f"{tgt}.pdb"
        return None

    out_dir = RESULTS / "contacts"
    out_dir.mkdir(parents=True, exist_ok=True)
    table = []

    targets = []
    for ds_dir in sorted(RESULTS.iterdir()):
        if not ds_dir.is_dir() or ds_dir.name in ("consistency", "contacts"):
            continue
        for tgt_dir in sorted(ds_dir.iterdir()):
            if not tgt_dir.is_dir():
                continue
            if args.dataset and ds_dir.name != args.dataset:
                continue
            if args.target and tgt_dir.name != args.target:
                continue
            targets.append((ds_dir.name, tgt_dir.name))

    for ds, tgt in targets:
        pdb = find_pdb(ds, tgt)
        fasta = REPO / "data" / ds / "converted" / f"{tgt}.fasta"
        if not pdb or not pdb.exists() or not fasta.exists():
            continue
        seq_pdb, contacts = contact_map(str(pdb))
        if not seq_pdb:
            continue
        # A2M-based targets: use the original a2m reference frame for mapping
        a2m_src = None
        if ds == "unit_tests" and tgt == "DHFR":
            a2m_src = REPO / "data" / "unit_tests" / "DHFR.a2m"
        elif ds == "evmutation":
            a2m_src = REPO / "data" / "evmutation" / "raw" / f"{tgt}.a2m"
        if a2m_src is not None and a2m_src.exists():
            mapping = a2m_reference_map(str(a2m_src), seq_pdb)
            print(f"  [{tgt}] a2m reference-frame mapping: {len(mapping)} cols", flush=True)
        else:
            ref = reference_seq(str(fasta))
            mapping = nw_align(ref, seq_pdb)
        if len(mapping) < 10 and ref:
            # fallback: consensus-based mapping (e.g., reference too gappy)
            cons = consensus_seq(str(fasta))
            mapping = nw_align(cons, seq_pdb)
        pairs = load_pairs_from_results(str(RESULTS / ds / tgt))
        if not pairs:
            continue
        # map to structure indices
        mapped = []
        for i, j, mi in pairs:
            ti, tj = mapping.get(i), mapping.get(j)
            if ti is not None and tj is not None:
                mapped.append((i, j, ti, tj, mi))
        if not mapped:
            continue
        n_contact = sum(1 for _, _, ti, tj, _ in mapped if (ti, tj) in contacts)
        prec = n_contact / len(mapped)
        # random expectation: same separation distribution
        seps = [abs(ti - tj) for _, _, ti, tj, _ in mapped]
        random.seed(42)
        n_rand = 2000
        hits = 0
        for _ in range(n_rand):
            a = random.randrange(len(seq_pdb))
            b = a + random.choice(seps)
            if b < len(seq_pdb) and (a, b) in contacts:
                hits += 1
        exp = hits / n_rand
        enrich = prec / exp if exp > 0 else float("nan")
        row = {
            "dataset": ds,
            "target": tgt,
            "pairs_mapped": len(mapped),
            "pairs_on_contact": n_contact,
            "precision": round(prec, 3),
            "random_expectation": round(exp, 3),
            "enrichment": round(enrich, 2),
            "seq_len_pdb": len(seq_pdb),
            "contacts": len(contacts),
        }
        table.append(row)
        json.dump(row, open(out_dir / f"{ds}_{tgt}.json", "w"), indent=1)
        print(
            f"{ds}/{tgt}: prec={prec:.3f} exp={exp:.3f} enrich={enrich:.1f}x "
            f"({n_contact}/{len(mapped)} pairs on contact)",
            flush=True,
        )

    json.dump(table, open(out_dir / "contact_validation.json", "w"), indent=1)
    print(f"table -> {out_dir / 'contact_validation.json'} ({len(table)} targets)")


if __name__ == "__main__":
    main()
