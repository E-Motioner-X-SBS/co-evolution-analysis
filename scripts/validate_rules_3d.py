#!/usr/bin/env python3
"""
RULE → 3D VALIDATION for the OMICRON models.

Question: the Boolean prime implicants / co-evolving position pairs were
deduced by pure mathematical logic (MI + K-map + Quine-McCluskey) from the
FASTA alignment. Are those same position pairs CLOSE IN 3D in the
SwissModel-built models of the Omicron Spike?

Method (per unique sequence u0000..u0298):
  1. alignment column -> model residue map: col - n_gaps_before(col) + 1
  2. best model PDB (highest GMQE from templates_summary.json)
  3. contact: Cb-Cb < 8 A (Ca for Gly), residue separation >= 6
  4. aggregate contact fraction over all models; enrichment vs random pairs
     with the same separation distribution (2,000 draws per model)

Output: results/contacts/rules_3d_omicron.json + printed table.
"""
import gzip
import json
import os
import random
import sys
from collections import Counter
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
SWISS = REPO / "swissmodel"
FASTA = REPO / "Spike_protein.aln-fasta"
MB = "/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/master_boolean/master_boolean_summary.json"
CORR = REPO / "analysis" / "corrected_results.json"
OUT = REPO / "results" / "contacts" / "rules_3d_omicron.json"
CONTACT = 8.0
MIN_SEP = 6


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


def col_to_residue_map(aligned_seq):
    """aligned col (0-based) -> model residue (1-based); None for gap cols."""
    m = {}
    res = 1
    for c, ch in enumerate(aligned_seq):
        if ch == "-":
            continue
        m[c] = res
        res += 1
    return m


def model_coords(pdb_gz_path):
    """Return {residue: (aa, cb_coords)} — Cb (Ca for Gly)."""
    atoms = {}
    with gzip.open(pdb_gz_path, "rt") as f:
        for line in f:
            if line.startswith("ATOM"):
                name = line[12:16].strip()
                if name not in ("CA", "CB"):
                    continue
                resi = int(line[22:26])
                aa3 = line[17:20].strip()
                aa1 = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G",
                       "HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N",
                       "PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V",
                       "TRP":"W","TYR":"Y"}.get(aa3, "X")
                atoms.setdefault(resi, {})[name] = (
                    float(line[30:38]), float(line[38:46]), float(line[46:54]))
    coords = {}
    for resi, at in atoms.items():
        aa = at.get("CA") and "X"  # placeholder
        cb = at.get("CB", at.get("CA"))
        if cb is not None:
            coords[resi] = cb
    return coords


def dist(a, b):
    return sum((u - v) ** 2 for u, v in zip(a, b)) ** 0.5


def main():
    seqs = parse_fasta(str(FASTA))
    uniques = json.load(open(SWISS / "batch" / "state" / "uniques.json"))
    mb = json.load(open(MB))
    corr = json.load(open(CORR))

    # rule position pairs (alignment coords, 0-based)
    pair_sets = {
        "master_boolean_top10": [(p["pos_i"], p["pos_j"])
                                 for p in mb["top_co_evolving_pairs"]],
        "essential_rules": [(r["pos_i"], r["pos_j"]) for r in mb["inferences"]],
        "corrected_top12": [(p["pos_i"], p["pos_j"]) for p in corr["top_mi_pairs"]],
    }
    all_pairs = sorted({p for ps in pair_sets.values() for p in ps})

    # per-uid aligned sequence (first occurrence)
    uid_seq = {}
    for u in uniques:
        uid_seq[u["uid"]] = seqs[u["indices"][0]][1]

    per_pair = {p: {"models_checked": 0, "models_in_contact": 0,
                    "distances": []} for p in all_pairs}
    random_contact_total = 0
    random_contact_n = 0
    models_used = 0

    for u in uniques:
        uid = u["uid"]
        aligned = uid_seq[uid]
        cmap = col_to_residue_map(aligned)
        # best model (highest GMQE)
        tsum = os.path.join(SWISS, "batch", "pdbs", uid, "templates_summary.json")
        if not os.path.exists(tsum):
            continue
        best = max(json.load(open(tsum)), key=lambda m: m.get("gmqe") or 0)
        pdb = os.path.join(SWISS, "batch", "pdbs", uid, f"{best.get('model_id')}.pdb.gz")
        if not os.path.exists(pdb):
            continue
        coords = model_coords(pdb)
        if not coords:
            continue
        models_used += 1
        residues = sorted(coords)
        # rule pairs (sep >= 3 — only exclude backbone-adjacent)
        for (i, j) in all_pairs:
            ri, rj = cmap.get(i), cmap.get(j)
            if ri is None or rj is None or ri not in coords or rj not in coords:
                continue
            if abs(ri - rj) < 3:
                continue
            per_pair[(i, j)]["models_checked"] += 1
            d = dist(coords[ri], coords[rj])
            per_pair[(i, j)]["distances"].append(d)
            if d < CONTACT:
                per_pair[(i, j)]["models_in_contact"] += 1
        # random control: same separation distribution AS THE RULES
        seps = [abs(cmap.get(i) - cmap.get(j))
                for (i, j) in all_pairs
                if cmap.get(i) and cmap.get(j) and abs(cmap.get(i) - cmap.get(j)) >= 3]
        if not seps:
            continue
        random.seed(42 + models_used)
        n_rand = 2000
        for _ in range(n_rand):
            a = random.choice(residues)
            s = random.choice(seps)
            b = a + s
            if b in coords and a in coords:
                if dist(coords[a], coords[b]) < CONTACT:
                    random_contact_total += 1
            random_contact_n += 1

    exp = random_contact_total / random_contact_n if random_contact_n else 0
    print(f"models used: {models_used}/299")
    print(f"random contact expectation: {exp:.4f}")
    print(f"\n{'pair':>10} {'checked':>8} {'contacts':>9} {'frac':>7} {'enrich':>7} {'median_d':>9}")
    rows = []
    for p in all_pairs:
        r = per_pair[p]
        n, c = r["models_checked"], r["models_in_contact"]
        if n == 0:
            continue
        frac = c / n
        enrich = frac / exp if exp > 0 else float("nan")
        md = sorted(r["distances"])[len(r["distances"]) // 2]
        rows.append({"pair": p, "models_checked": n, "in_contact": c,
                     "fraction": round(frac, 3), "enrichment": round(enrich, 2),
                     "median_dist_A": round(md, 2)})
        print(f"{str(p):>10} {n:>8} {c:>9} {frac:>7.3f} {enrich:>7.2f} {md:>9.2f}")
    # essential rules specifically
    for name, ps in pair_sets.items():
        tot_n = sum(per_pair[p]["models_checked"] for p in ps)
        tot_c = sum(per_pair[p]["models_in_contact"] for p in ps)
        if tot_n:
            print(f"\n{name}: {tot_c}/{tot_n} models with >=1 contact "
                  f"({tot_c/tot_n:.3f}, enrich {tot_c/tot_n/exp:.2f}x)")
    json.dump({"models_used": models_used, "random_expectation": exp,
               "per_pair": rows,
               "sets": {name: {"checked": sum(per_pair[p]["models_checked"] for p in ps),
                               "contacts": sum(per_pair[p]["models_in_contact"] for p in ps)}
                        for name, ps in pair_sets.items()}},
              open(OUT, "w"), indent=1)
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
