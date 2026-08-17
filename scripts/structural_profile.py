#!/usr/bin/env python3
"""
structural_profile.py — what does the 3D structure tell us about the
sequence-space features (entropy / conservation)?

Analyses:
  1. BURIAL profile: per-residue local density (number of Cb atoms within
     12 A) in the model = burial proxy. Classify buried vs exposed.
  2. ENTROPY vs BURIAL: are high-entropy (variable) positions surface-exposed
     and conserved positions buried? (Spearman rho; the "housekeeping"
     hypothesis: conserved = structural core, variable = immune-exposed)
  3. SECONDARY STRUCTURE context of the co-evolving positions (from the
     model's dssp annotation when present, else inferred from backbone
     geometry): are co-evolving pairs in the same SS element?
  4. The 3D context of the co-evolving pairs: are the contact-supported
     pairs (407,410), (212,215), (210,215) in structured (sheet/helix) vs
     loop regions?

Output: results/contacts/structural_profile.json
"""
import json
import math
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO / "scripts"))
from rules3d_common import (OmicronModels, parse_model_pdb, CONTACT,
                            MIN_SEP)  # noqa: E402

OUT = REPO / "results" / "contacts" / "structural_profile.json"
BURIAL_RADIUS = 12.0


def local_density(coords, residue):
    """Count Cb atoms within 12 A of residue (excl. self)."""
    c0 = coords[residue]
    n = 0
    for r, c in coords.items():
        if r == residue:
            continue
        if math.sqrt(sum((u - v) ** 2 for u, v in zip(c0, c))) < BURIAL_RADIUS:
            n += 1
    return n


def infer_ss(phi_psi_scores):
    """Not used — dssp from the payload is preferred; fallback below."""
    return None


def main():
    models = OmicronModels()
    # entropy per alignment column (from the canonical analysis)
    full_len = json.load(open(
        "/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/"
        "full_length_results/full_length_summary.json"))
    # entropy vector: recompute properly from the alignment (the summary has
    # only the top-20)
    import numpy as np
    from collections import Counter
    from rules3d_common import parse_fasta
    seqs = parse_fasta(str(REPO / "Spike_protein.aln-fasta"))
    L = len(seqs[0][1])
    H = np.zeros(L)
    for c in range(L):
        cnt = Counter(s[c] for _, s in seqs if len(s) > c and s[c] in
                      "ACDEFGHIKLMNPQRSTVWY")
        t = sum(cnt.values())
        if t == 0:
            continue
        H[c] = -sum((v / t) * math.log2(v / t) for v in cnt.values())

    burial_by_col = {}
    ss_by_col = {}
    n_models = 0
    for uid, aligned, cmap, coords in models.iter_models():
        n_models += 1
        chain = sorted(coords)[0]
        residues = sorted(coords[chain])
        dens = {r: local_density(coords[chain], r) for r in residues}
        # burial threshold: median density over the protein
        med = sorted(dens.values())[len(dens) // 2]
        for col, res in cmap.items():
            if res in dens:
                burial_by_col.setdefault(col, []).append(
                    1 if dens[res] > med else 0)
        # secondary structure from dssp annotation if available
        fd = REPO / "swissmodel" / "batch" / "pdbs" / uid / "full_details.json"
        try:
            details = json.load(open(fd))
            for m in details.get("models", []):
                ann = (m.get("targets") or [{}])[0].get(
                    "target_sequence_annotations") or {}
                dssp = ann.get("dssp")
                if dssp:
                    for col, res in cmap.items():
                        if res - 1 < len(dssp):
                            ch = dssp[res - 1]
                            ss_by_col.setdefault(col, []).append(ch)
        except Exception:
            pass

    # collapse per column
    cols = sorted(burial_by_col)
    frac_buried = {c: sum(v) / len(v) for c, v in burial_by_col.items()}
    # entropy quantiles
    H_vals = [H[c] for c in cols]
    hs = sorted(H_vals)
    q = {0.2: hs[int(0.2 * len(hs))], 0.8: hs[int(0.8 * len(hs))]}
    groups = {"conserved": [c for c in cols if H[c] < 1e-9],
              "intermediate": [c for c in cols if 1e-9 <= H[c] <= q[0.8]],
              "variable": [c for c in cols if H[c] > q[0.8]]}
    for g, cl in groups.items():
        fb = [frac_buried[c] for c in cl]
        print(f"{g}: n={len(cl)} mean_frac_buried={sum(fb)/len(fb) if fb else 0:.3f}")
    # Spearman entropy vs burial
    xs = [H[c] for c in cols]
    ys = [frac_buried[c] for c in cols]
    from scipy.stats import spearmanr
    rho = spearmanr(xs, ys)
    print(f"Spearman(entropy, burial): rho={rho.statistic:.3f} p={rho.pvalue:.3g}")
    # variable positions buried fraction vs conserved
    fb_var = [frac_buried[c] for c in groups["variable"]]
    fb_cons = [frac_buried[c] for c in groups["conserved"]]
    print(f"variable positions mean buried frac: {sum(fb_var)/len(fb_var):.3f} "
          f"| conserved: {sum(fb_cons)/len(fb_cons):.3f}")

    # SS context of the KEY co-evolving positions
    print("\nsecondary-structure context of co-evolving positions:")
    ss_rows = {}
    for col in [212, 213, 214, 215, 216, 210, 373, 378, 407, 410, 442, 448,
                454, 488, 495, 498, 18, 26, 66, 94]:
        if col in ss_by_col:
            from collections import Counter as C2
            ss_rows[col] = dict(C2(ss_by_col[col]))
            print(f"  pos {col}: {ss_rows[col]}")
        elif col in burial_by_col:
            print(f"  pos {col}: burial frac {frac_buried[col]:.3f} "
                  f"(no dssp), H={H[col]:.2f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"n_models": n_models, "spearman_entropy_burial": rho.statistic,
               "spearman_p": rho.pvalue,
               "burial_by_group": {g: round(sum(
                   [frac_buried[c] for c in cl]) / len(cl), 3) if cl else None
                   for g, cl in groups.items()},
               "frac_buried": frac_buried, "ss_context": ss_rows},
              open(OUT, "w"), indent=1)
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
