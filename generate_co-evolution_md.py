#!/usr/bin/env python3
"""
Generate comprehensive co-evolution markdown with K-maps,
Boolean functions, coupling constants, and inference rules.
"""

import sys
import json
import numpy as np
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")

from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012
from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap


def parse_fasta(filepath):
    sequences = []
    current_header = None
    current_seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line.upper())
    if current_header is not None:
        sequences.append((current_header, "".join(current_seq)))
    return sequences


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")

    # Load data
    with open(base_dir / "master_boolean" / "master_boolean_summary.json") as f:
        summary = json.load(f)

    aa_list = list(AMINO_HE_2012)
    encoder = Base20AminoEncoder(version=1)
    sequences = parse_fasta(base_dir / "Spike_protein.aln-fasta")
    n_all = len(sequences)

    # Build position arrays
    # BUG FIX: Original code used clean[:80] (first 80 positions only), but
    # master_boolean.py computes co-evolving pairs across the FULL length
    # (e.g., positions 413, 427, ...). This mismatch caused the K-map to be
    # all zeros → 0 rules. Fixed to use the full sequence length.
    pos_arrays = []
    for _, seq in sequences:
        # CORRECTED: aligned columns, gap = 20 (was gap-stripped, misaligned)
        arr = np.array([encoder.encode.get(aa, 20) for aa in seq], dtype=np.int32)
        pos_arrays.append(arr)

    def get_majority(pos, n):
        return Counter(
            int(a[pos]) for a in pos_arrays[:n] if pos < len(a) and 0 <= a[pos] < 20
        ).most_common(1)[0][0]

    def compute_coupling(pos_i, pos_j, n):
        joint = np.zeros((20, 20), dtype=np.float64)
        for arr in pos_arrays[:n]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if 0 <= ci < 20 and 0 <= cj < 20:
                    joint[ci, cj] += 1
        total = joint.sum()
        if total == 0:
            return []
        joint /= total
        marg_i = joint.sum(axis=1)
        marg_j = joint.sum(axis=0)
        epsilon = 1e-10
        J = np.log(
            (joint + epsilon)
            / ((marg_i[:, None] + epsilon) * (marg_j[None, :] + epsilon))
        )
        couplings = []
        for ai in range(20):
            for aj in range(20):
                if joint[ai, aj] > 0.005:
                    couplings.append(
                        (
                            aa_list[ai],
                            aa_list[aj],
                            float(J[ai, aj]),
                            float(joint[ai, aj]),
                        )
                    )
        couplings.sort(key=lambda x: abs(x[2]), reverse=True)
        return couplings

    # Generate markdown
    md = []
    md.append("# Co-evolution Boolean Functions for SARS-CoV-2 Spike Protein\n")
    md.append("**Dataset:** {} sequences".format(n_all))
    md.append("**Encoding:** Base-20 (He 2012 order)")
    md.append("**Method:** Variable-position K-map with don't-care conditions")
    md.append("**Quine-McCluskey minimization**\n")
    md.append("---\n")

    md.append("## Variables\n")
    md.append("| Variable | Meaning |")
    md.append("|----------|---------|")
    md.append("| s4, s3, s2, s1, s0 | Binary code for residue at position i (0-19, 5 bits) |")
    md.append("| t4, t3, t2, t1, t0 | Binary code for residue at position j (0-19, 5 bits) |")
    md.append("| ~s4 | NOT s4 (bit is 0) |")
    md.append("| s4 | bit is 1 |")
    md.append("| s4.s3.s2.s1.s0 | AND of bits |")
    md.append("")
    md.append("---\n")

    rule_num = 0

    for idx, pair in enumerate(summary["top_co_evolving_pairs"][:15]):
        pos_i = int(pair["pos_i"])
        pos_j = int(pair["pos_j"])
        mi = float(pair["mi"])
        n_muts = int(pair["n"])
        ref_i = pair["ref_i"]
        ref_j = pair["ref_j"]

        ref_aa_i = ref_i
        ref_aa_j = ref_j
        ref_i_code = aa_list.index(ref_i) if ref_i in aa_list else 0
        ref_j_code = aa_list.index(ref_j) if ref_j in aa_list else 0

        # Build K-map
        # CORRECTED (FIX A2): 32x32 padded, rows/cols 20-31 don't-care.
        # 20x20 region default 0 (off-set, never observed); observed
        # mutation -> 1; reference -> -1 (don't-care).
        kmap = np.full((32, 32), -1, dtype=np.int32)
        kmap[:20, :20] = 0
        for arr in pos_arrays[:n_all]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if 0 <= ci < 20 and 0 <= cj < 20:
                    kmap[ci, cj] = 1 if (ci != ref_i_code or cj != ref_j_code) else -1

        # Run QM
        bool_flat = kmap.flatten().astype(int)
        result = boolean_minimize_kmap(bool_flat, algorithm="qm")

        # Decode QM cubes to residue pairs; deduplicate by (aa_i, aa_j),
        # OR-ing the essential flag across cubes that decode to the same
        # pair (they differ only in which don't-care cells they cover).
        pi_map = {}  # (aa_i, aa_j) -> is_essential
        pi_order = []
        for pi in result["prime_implicants"]:
            values = list(pi["values"])
            mask = list(pi["mask"])
            # CORRECTED (FIX A2): kmap is 32x32 -> 5 bits per axis (10 bits)
            while len(values) < 10:
                values.append(0)
                mask.append(False)
            row_code = sum(
                values[j] * (2 ** (4 - j)) for j in range(5) if not mask[j]
            )
            col_code = sum(
                values[j + 5] * (2 ** (4 - j)) for j in range(5) if not mask[j + 5]
            )
            if row_code < 20 and col_code < 20:
                aa_i = aa_list[row_code]
                aa_j = aa_list[col_code]
                is_essential = pi in result["essential_prime_implicants"]
                key = (aa_i, aa_j)
                if key not in pi_map:
                    pi_map[key] = is_essential
                    pi_order.append(key)
                else:
                    pi_map[key] = pi_map[key] or is_essential
        pi_entries = [(k[0], k[1], pi_map[k]) for k in pi_order]

        # Count cells
        n_on = int((kmap == 1).sum())
        n_off = int((kmap == 0).sum())
        n_dc = int((kmap == -1).sum())

        # Compute couplings
        couplings = compute_coupling(pos_i, pos_j, n_all)

        md.append("## Position Pair ({}, {})".format(pos_i, pos_j))
        md.append("")
        md.append("| Property | Value |")
        md.append("|----------|-------|")
        md.append("| Mutual Information | {:.4f} |".format(mi))
        md.append("| Total mutations | {} |".format(n_muts))
        md.append("| Reference pos {} | {} |".format(pos_i, ref_aa_i))
        md.append("| Reference pos {} | {} |".format(pos_j, ref_aa_j))
        md.append("| On-set cells | {} |".format(n_on))
        md.append("| Off-set cells | {} |".format(n_off))
        md.append("| Don't-care cells | {} |".format(n_dc))
        md.append("| Prime implicants | {} |".format(len(pi_entries)))
        md.append("| Essential PIs | {} |".format(sum(1 for _, _, e in pi_entries if e)))
        md.append("")

        # K-map compact
        md.append("### K-map (Compact View)")
        md.append("")
        md.append("```")
        md.append(
            "Position pair ({}, {}): Reference = ({}, {})".format(
                pos_i, pos_j, ref_aa_i, ref_aa_j
            )
        )
        md.append("")
        md.append("Co-evolutionary residue pairs (on-set):")
        on_pairs = []
        for ai in range(20):
            for aj in range(20):
                if kmap[ai, aj] == 1:
                    on_pairs.append("{}-{}".format(aa_list[ai], aa_list[aj]))
        md.append("  {}".format(", ".join(on_pairs)))
        md.append("")
        md.append("Don't-care positions (conserved): {} cells".format(n_dc))
        md.append("Never-seen pairs (off-set): {} cells".format(n_off))
        md.append("```")
        md.append("")

        # Boolean function
        md.append("### Boolean Function")
        md.append("")
        md.append("```")
        md.append(
            "f(pos_{}, pos_{}) = 1 if ANY of these residue pairs appear:".format(
                pos_i, pos_j
            )
        )
        md.append("")

        pi_num = 0
        for aa_i, aa_j, is_essential in pi_entries:
            marker = "*" if is_essential else " "
            pi_num += 1
            md.append(
                "  {} PI_{}: pos_{}={} AND pos_{}={}".format(
                    marker, pi_num, pos_i, aa_i, pos_j, aa_j
                )
            )

        md.append("")
        md.append("(* = essential prime implicant)")
        md.append("```")
        md.append("")

        # Coupling constants
        md.append("### Coupling Constants (J_ij)")
        md.append("")
        md.append("| Residue Pair | J_ij | Frequency | Type |")
        md.append("|-------------|------|-----------|------|")
        for aa1, aa2, j_val, freq in couplings[:10]:
            jtype = "co-evolutionary" if j_val > 0 else "anti-correlated"
            md.append(
                "| {}-{} | {:.4f} | {:.4f} | {} |".format(aa1, aa2, j_val, freq, jtype)
            )
        md.append("")

        # Inference rules
        md.append("### Inference Rules (Natural Language)")
        md.append("")
        for aa_i, aa_j, is_essential in pi_entries:
            rule_num += 1
            if is_essential:
                md.append(
                    "**Rule {}:** IF position {} = **{}** AND position {} = **{}** THEN co-evolutionary (MI = {:.3f})".format(
                        rule_num, pos_i, aa_i, pos_j, aa_j, mi
                    )
                )
            else:
                md.append(
                    "Rule {}: IF position {} = {} AND position {} = {} THEN co-evolutionary (MI = {:.3f})".format(
                        rule_num, pos_i, aa_i, pos_j, aa_j, mi
                    )
                )

        md.append("")
        md.append("---")
        md.append("")

    # Summary
    md.append("## Summary\n")
    md.append("| Metric | Value |")
    md.append("|--------|-------|")
    md.append("| Sequences | {} |".format(n_all))
    md.append("| Variable positions | {} |".format(summary["variable_positions"]))
    md.append("| Co-evolutionary pairs | {} |".format(summary["co_evolving_pairs"]))
    md.append("| Total inference rules | {} |".format(rule_num))
    md.append("| Position pairs with rules | {} |".format(len(summary["top_co_evolving_pairs"])))
    md.append("")
    md.append("## How to Apply\n")
    md.append("1. Extract the residue pair at each co-evolving position pair from a new sequence")
    md.append("2. For each position pair, check if the residue pair matches any rule below")
    md.append("3. If YES: that position pair is co-evolutionary (consistent with the observed data)")
    md.append(
        "4. If position i mutates: find which residue at position j satisfies the co-evolutionary constraint"
    )
    md.append("")
    if summary["inferences"]:
        ex = summary["inferences"][0]
        md.append("**Example:** If position {} mutates to {}, check rules for position {}.".format(ex["pos_i"], ex["aa_i"], ex["pos_i"]))
        md.append("Essential rule: IF pos {} = {} AND pos {} = {} THEN co-evolutionary.".format(ex["pos_i"], ex["aa_i"], ex["pos_j"], ex["aa_j"]))
        md.append("So position {} must also show {} to satisfy the co-evolutionary constraint.".format(ex["pos_j"], ex["aa_j"]))
    else:
        md.append("**Example:** No essential inference rules found; see per-pair rule lists above.")

    # Write
    output_path = base_dir / "kmap_boolean_coevolution" / "COEVOLUTION_KMAP_BOOLEAN.md"
    with open(output_path, "w") as f:
        f.write("\n".join(md))

    print("Written {} lines to {}".format(len(md), output_path))
    print("Total rules: {}".format(rule_num))




    # ============================================================
    # COMBINED MI + PERPLEXITY ANALYSIS (all experiments)
    # ============================================================
    print("\n=== Combined MI + Perplexity Analysis ===")
    try:
        from coevolution_shared import (
            combined_pair_scores, compute_entropy_vectorized,
            load_position_arrays as _lpa,
        )
        _pa, _na, _fl = _lpa(max_pos=None, aligned=True)
        _ent = compute_entropy_vectorized(_pa, _na, _fl)
        _var = [p for p in range(_fl) if _ent[p] > 0.3]
        _pairs = [(i, j) for idx, i in enumerate(_var)
                  for j in _var[idx + 1:] if j - i <= 30]
        _scored = combined_pair_scores(_pa, _pairs, _na, _ent)
        print(f"  Variable positions: {len(_var)}")
        print(f"  Pairs scored (MI + perplexity ratio): {len(_scored)}")
        print(f"  Top 5 combined (MI + ratio):")
        for _s in _scored[:5]:
            print(f"    ({_s['pos_i']},{_s['pos_j']}): MI={_s['mi']:.3f} "
                  f"ratio={_s['ratio']:.2f} combined={_s['combined']:.3f}")
        _mi_top = sorted(_scored, key=lambda s: -s['mi'])[:5]
        print(f"  Top 5 by MI alone:")
        for _s in _mi_top:
            print(f"    ({_s['pos_i']},{_s['pos_j']}): MI={_s['mi']:.3f} "
                  f"ratio={_s['ratio']:.2f}")
        # ranking agreement
        _r_mi = {(_s['pos_i'], _s['pos_j']): idx
                 for idx, _s in enumerate(sorted(_scored, key=lambda s: -s['mi']))}
        _r_cb = {(_s['pos_i'], _s['pos_j']): idx
                 for idx, _s in enumerate(_scored)}
        _same = sum(1 for k in _r_mi if _r_mi[k] == _r_cb[k])
        print(f"  Ranking agreement (MI vs combined, top-5 same): "
              f"{len([k for k in _r_mi if k in _r_cb and _r_mi[k] < 5 and _r_cb[k] < 5])}/5")
    except Exception as _e:
        print(f"  Combined analysis skipped: {_e}")


if __name__ == "__main__":
    main()
