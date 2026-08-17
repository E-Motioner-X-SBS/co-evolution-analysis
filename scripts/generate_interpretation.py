#!/usr/bin/env python3
"""
PHASE 9: Per-dataset interpretation docs for COMPLETED datasets.

For each dataset with full results, writes:
  My_Own_Interpretation_Across_New_Datasets/<dataset>/NN_<script>.md
  (what it does, algorithm, key formulas, results table, honest notes)

Run as results complete: python scripts/generate_interpretation.py [--dataset X]
"""

import json
import os
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
RESULTS = REPO / "results"
OUT = REPO / "My_Own_Interpretation_Across_New_Datasets"

SCRIPT_DOC = {
    "master_boolean": "Master Boolean function: per-pair mutation K-maps -> "
    "Quine-McCluskey -> prime implicants + essential rules.",
    "kmap_boolean_coevolution": "K-map Boolean equations (5-bit Gray residue codes) "
    "for each co-evolving pair.",
    "flipped_boolean_coevolution": "Flipped (forbidden) K-maps: residue pairs never "
    "observed -> negative-selection constraints.",
    "nary_kmap_results": "Base-20 (n-ary) K-map analysis: 20x20 frequency maps -> "
    "Boolean minimization in base-20.",
    "perplexity_results": "Perplexity ratio PP(j)/mean PP(j|i): determinism of "
    "co-evolution (ratio>1 = constrained).",
    "full_length_results": "Full-length entropy + mutual information over ALL "
    "positions (window 30).",
    "boolean_results": "Whole-alignment dipeptide Boolean minimization + coupling "
    "constants J = MI.",
    "position_kmap_results": "Per-position-pair K-maps with MI.",
    "full_position_results": "All-sequence MI matrix (window 30).",
    "kmap_results": "Binary Gray K-map pipeline: H1 adjacency enrichment, signatures, "
    "Walsh-Hadamard, mutation analysis.",
    "variable_position_results": "Variable positions only (H>0.3) with don't-care K-maps.",
    "constraint_function_results": "Constraint function C=ln(P/P_exp) train/test "
    "prediction (proportional 62/38 split).",
    "allseq_constraint_results": "LOO-CV constraint-function prediction.",
    "dca_boolean_results": "Local precision matrix -> Boolean (NOT real DCA).",
    "dca_results": "Proper mfDCA: Frobenius norm + APC + direct information.",
    "advanced_analysis_results": "Network, clustering, Walsh-Hadamard, variant signatures.",
    "full_gpu_results": "GPU full analysis: full MI matrix, entropy, refs.",
    "mi_heatmap": "MI heatmap PNG/PDF + summary.",
}


def load(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset")
    args = ap.parse_args()

    for ds_dir in sorted(RESULTS.iterdir()):
        if not ds_dir.is_dir() or ds_dir.name in ("consistency", "contacts"):
            continue
        ds = ds_dir.name
        if args.dataset and ds != args.dataset:
            continue
        ds_out = OUT / ds
        os.makedirs(ds_out, exist_ok=True)
        idx = []
        for tgt_dir in sorted(ds_dir.iterdir()):
            if not tgt_dir.is_dir():
                continue
            tgt = tgt_dir.name
            for dirname, desc in SCRIPT_DOC.items():
                rdir = tgt_dir / dirname
                if not rdir.exists():
                    continue
                runlog = load(rdir / "run_log.json")
                status = runlog.get("status", "?") if runlog else "?"
                secs = runlog.get("seconds", "?") if runlog else "?"
                # collect result files
                files = sorted(
                    f.name for f in rdir.iterdir() if f.name != "run_log.json"
                )
                md = [
                    f"# {ds} / {tgt} — {dirname}",
                    "",
                    f"**Status:** {status} ({secs}s)  ",
                    "",
                    f"**What it does:** {desc}",
                    "",
                    "**Outputs:** " + (", ".join(files[:8]) if files else "—"),
                    "",
                ]
                # add a results summary for the main JSONs
                for jf in (
                    "master_boolean_summary.json",
                    "boolean_functions.json",
                    "perplexity_summary.json",
                    "nary_analysis_summary.json",
                ):
                    jp = rdir / jf
                    if jp.exists():
                        d = load(jp)
                        md.append(f"### {jf}")
                        md.append("```json")
                        md.append(
                            json.dumps(
                                {
                                    k: d[k]
                                    for k in list(d)[:8]
                                    if not isinstance(d.get(k), (dict, list))
                                },
                                indent=1,
                            )[:1500]
                        )
                        md.append("```")
                        md.append("")
                outfile = ds_out / f"{tgt[:40]}_{dirname}.md"
                outfile.write_text("\n".join(md))
                idx.append((tgt, dirname, status, str(secs)))
        if idx:
            readme = [f"# {ds} — Interpretation Index", ""]
            readme.append("| Target | Script | Status | Time(s) |")
            readme.append("|--------|--------|--------|---------|")
            for t, s, st, sec in sorted(idx):
                readme.append(f"| {t[:40]} | {s} | {st} | {sec} |")
            (ds_out / "README.md").write_text("\n".join(readme))
            print(f"{ds}: {len(idx)} docs")
    print("Phase 9 docs done")


if __name__ == "__main__":
    main()
