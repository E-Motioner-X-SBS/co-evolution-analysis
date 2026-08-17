#!/usr/bin/env python3
"""
validate_residue_rules_3d.py — RESIDUE-LEVEL validation of the 36 prime
implicants (and the 2 essential rules).

A prime implicant is a rule of the form
    IF pos_i = aa_i AND pos_j = aa_j THEN co-evolutionary
deduced purely from sequence logic (K-map + Quine-McCluskey).

Test: in the 3D models, for every sequence/model where the rule's residues
ACTUALLY OCCUR at those positions (i.e., the sequence satisfies the rule's
left-hand side), are the two residues in 3D contact?

This is the strongest test of the rule's structural meaning: it validates the
rule on the instances where it fires.

Method:
  - for each of the 36 residue rules (pos_i, aa_i, pos_j, aa_j):
      for each of the 299 models:
        map alignment cols -> model residues (canonical-only, verified)
        if model residue at pos_i == aa_i AND at pos_j == aa_j (the rule fires)
          and sep >= 3:
            measure Cb-Cb distance; contact if < 8 A
  - control: random residue pairs (same separation distribution) that "fire"
    a random rule-like condition -> contact rate
  - report per-rule: n fired, n contacts, fraction, enrichment, binomial p

Output: results/contacts/validate_residue_rules_3d.json
"""
import json
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO / "scripts"))
from rules3d_common import (OmicronModels, CONTACT, MIN_SEP, N_RAND,
                            dist)  # noqa: E402

INV = REPO / "results" / "contacts" / "rules_inventory.json"
OUT = REPO / "results" / "contacts" / "validate_residue_rules_3d.json"


def main():
    inv = json.load(open(INV))
    rules = inv["residue_rules"]          # 36 on-set cells
    essentials = inv["essential_rules"]   # 2
    models = OmicronModels()

    # per-rule stats
    stats = {tuple(sorted((r["pos_i"], r["pos_j"]))): {
        "rule": r, "fired": 0, "contacts": 0, "dists": []}
        for r in rules}
    # also track by rule identity (pos_i,aa_i,pos_j,aa_j)
    fired_map = {}
    for r in rules:
        fired_map[(r["pos_i"], r["aa_i"], r["pos_j"], r["aa_j"])] = {
            "fired": 0, "contacts": 0, "dists": []}
    for r in essentials:
        fired_map[(r["pos_i"], r["aa_i"], r["pos_j"], r["aa_j"])] = {
            "fired": 0, "contacts": 0, "dists": []}

    rnd_contacts = 0
    rnd_n = 0
    model_count = 0
    import random
    for uid, aligned, cmap, coords in models.iter_models():
        model_count += 1
        chains = sorted(coords)
        if not chains:
            continue
        chain = chains[0]
        # evaluate every residue rule on this model
        for (pi, ai, pj, aj), st in fired_map.items():
            ri, rj = cmap.get(pi), cmap.get(pj)
            if ri is None or rj is None:
                continue
            if ri not in coords[chain] or rj not in coords[chain]:
                continue
            if abs(ri - rj) < MIN_SEP:
                continue
            # does the rule fire? need residue identity at those positions
            # (recover aa from the model: we only stored coords... use the
            # model's target sequence instead — the PDB parse drops aa info)
            # -> we reconstruct via the aligned sequence column char
            fired = False
            # aa at alignment col pi/pj from the model's target sequence
            # (accessible via cmap: residue ri in target == aligned char)
            # We need the target sequence; recover from the aligned seq:
            # the residue at column pi is aligned[pi] if canonical.
            # NOTE: aligned[pi] may be '-'? no — cmap maps only canonical
            # columns, and the char at a mapped column is canonical by
            # construction, but it may differ from the RULE's aa_i (the rule
            # aa is from the MSA consensus of that column).
            if aligned[pi] == ai and aligned[pj] == aj:
                fired = True
            if not fired:
                continue
            st["fired"] += 1
            d = dist(coords[chain][ri], coords[chain][rj])
            st["dists"].append(d)
            if d < CONTACT:
                st["contacts"] += 1
        # control: random pairs with matched separations, "firing" at random
        residues = sorted(coords[chain])
        if len(residues) < 10:
            continue
        seps = [abs(ri - rj) for (pi, ai, pj, aj), st in fired_map.items()
                if (ri := cmap.get(pi)) and (rj := cmap.get(pj))
                and abs(ri - rj) >= MIN_SEP]
        if not seps:
            continue
        rng = random.Random(7 + model_count)
        for _ in range(N_RAND):
            a = rng.choice(residues)
            s = rng.choice(seps)
            b = a + s
            if b not in coords[chain]:
                continue
            if dist(coords[chain][a], coords[chain][b]) < CONTACT:
                rnd_contacts += 1
            rnd_n += 1

    control = rnd_contacts / rnd_n if rnd_n else 0.0
    print(f"models used: {model_count} | random control rate: {control:.4f}")

    from scipy.stats import binomtest
    rows = []
    for key, st in fired_map.items():
        if st["fired"] == 0:
            continue
        frac = st["contacts"] / st["fired"]
        enrich = frac / control if control else float("nan")
        pval = binomtest(st["contacts"], st["fired"], control,
                         alternative="greater").pvalue
        rows.append({"rule": list(key), "fired": st["fired"],
                     "contacts": st["contacts"], "fraction": round(frac, 3),
                     "enrichment": round(enrich, 2), "p_value": pval,
                     "median_dist": round(sorted(st["dists"])[len(st["dists"]) // 2], 2)})
    rows.sort(key=lambda r: -r["enrichment"])
    print(f"\n{'rule':>22} {'fired':>6} {'ct':>4} {'frac':>7} {'enr':>7} {'p':>9} {'medA':>6}")
    for r in rows:
        print(f"{str(r['rule']):>22} {r['fired']:>6} {r['contacts']:>4} "
              f"{r['fraction']:>7.3f} {r['enrichment']:>7.2f} "
              f"{r['p_value']:>9.2g} {r['median_dist']:>6.1f}")

    # pooled over all 36 PIs
    n = sum(st["fired"] for st in fired_map.values())
    c = sum(st["contacts"] for st in fired_map.values())
    frac = c / n if n else 0
    pv = binomtest(c, n, control, alternative="greater").pvalue if n else 1.0
    print(f"\nALL 36 PIs pooled: {c}/{n} ({frac:.3f}, "
          f"enrich {frac/control:.2f}x, p={pv:.3g})")
    # essentials only
    ne = sum(fired_map[(r["pos_i"], r["aa_i"], r["pos_j"], r["aa_j"])]["fired"]
             for r in essentials)
    ce = sum(fired_map[(r["pos_i"], r["aa_i"], r["pos_j"], r["aa_j"])]["contacts"]
             for r in essentials)
    fe = ce / ne if ne else 0
    pe = binomtest(ce, ne, control, alternative="greater").pvalue if ne else 1.0
    print(f"2 ESSENTIAL rules pooled: {ce}/{ne} ({fe:.3f}, "
          f"enrich {fe/control:.2f}x, p={pe:.3g})")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"models_used": model_count, "control_rate": control,
               "per_rule": rows,
               "pooled_all_36": {"n": n, "contacts": c, "fraction": frac,
                                 "p_value": pv},
               "pooled_essential": {"n": ne, "contacts": ce,
                                    "fraction": fe, "p_value": pe}},
              open(OUT, "w"), indent=1)
    print(f"saved -> {OUT}")


if __name__ == "__main__":
    main()
