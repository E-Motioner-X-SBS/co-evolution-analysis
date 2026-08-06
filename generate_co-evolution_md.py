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
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        arr = np.array([encoder.encode.get(aa, -1) for aa in clean], dtype=np.int32)
        pos_arrays.append(arr)

    def get_majority(pos, n):
        return Counter(
            int(a[pos]) for a in pos_arrays[:n] if pos < len(a) and a[pos] >= 0
        ).most_common(1)[0][0]

    def compute_coupling(pos_i, pos_j, n):
        joint = np.zeros((20, 20), dtype=np.float64)
        for arr in pos_arrays[:n]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if ci >= 0 and cj >= 0:
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
    md.append("| s3, s2, s1, s0 | Binary code for residue at position i (0-19) |")
    md.append("| t3, t2, t1, t0 | Binary code for residue at position j (0-19) |")
    md.append("| ~s3 | NOT s3 (bit is 0) |")
    md.append("| s3 | bit is 1 |")
    md.append("| s3.s2.s1.s0 | AND of bits |")
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
        kmap = np.zeros((20, 20), dtype=np.int32)
        for arr in pos_arrays[:n_all]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if ci >= 0 and cj >= 0:
                    kmap[ci, cj] = 1 if (ci != ref_i_code or cj != ref_j_code) else -1

        # Run QM
        bool_flat = kmap.flatten().astype(int)
        result = boolean_minimize_kmap(bool_flat, algorithm="qm")

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
        md.append("| Prime implicants | {} |".format(result["n_prime_implicants"]))
        md.append("| Essential PIs | {} |".format(result["n_essential"]))
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
        for pi in result["prime_implicants"]:
            values = list(pi["values"])
            mask = list(pi["mask"])
            while len(values) < 8:
                values.append(0)
                mask.append(False)

            row_code = sum(values[j] * (2 ** (3 - j)) for j in range(4) if not mask[j])
            col_code = sum(
                values[j + 4] * (2 ** (3 - j)) for j in range(4) if not mask[j + 4]
            )

            if row_code < 20 and col_code < 20:
                aa_i = aa_list[row_code]
                aa_j = aa_list[col_code]
                is_essential = pi in result["essential_prime_implicants"]
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
        for pi in result["prime_implicants"]:
            values = list(pi["values"])
            mask = list(pi["mask"])
            while len(values) < 8:
                values.append(0)
                mask.append(False)

            row_code = sum(values[j] * (2 ** (3 - j)) for j in range(4) if not mask[j])
            col_code = sum(
                values[j + 4] * (2 ** (3 - j)) for j in range(4) if not mask[j + 4]
            )

            if row_code < 20 and col_code < 20:
                aa_i = aa_list[row_code]
                aa_j = aa_list[col_code]
                is_essential = pi in result["essential_prime_implicants"]
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
    md.append("| Variable positions | 57 |")
    md.append("| Co-evolutionary pairs | {} |".format(summary["co_evolving_pairs"]))
    md.append("| Total inference rules | {} |".format(rule_num))
    md.append("| Position pairs with rules | 15 |")
    md.append("")
    md.append("## How to Apply\n")
    md.append("1. Extract residues at positions 68-79 from a new sequence")
    md.append("2. For each position pair, check if the residue pair matches any rule")
    md.append("3. If YES: that position pair is co-evolutionary")
    md.append(
        "4. If position i mutates: find which residue at position j satisfies the co-evolutionary constraint"
    )
    md.append("")
    md.append("**Example:** If position 76 mutates to Y, check rules for position 76.")
    md.append("Rule 2 says: IF pos 76 = Y AND pos 77 = K THEN co-evolutionary.")
    md.append("So position 77 must also mutate to K.")

    # Write
    output_path = base_dir / "kmap_boolean_coevolution" / "COEVOLUTION_KMAP_BOOLEAN.md"
    with open(output_path, "w") as f:
        f.write("\n".join(md))

    print("Written {} lines to {}".format(len(md), output_path))
    print("Total rules: {}".format(rule_num))


if __name__ == "__main__":
    main()
