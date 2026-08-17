#!/usr/bin/env python3
"""
COMPLETENESS CHECK: every dataset x target x script must have a run_log
(status ok, or a documented failure with reason). Reports anything missing.

Usage: python scripts/check_completeness.py
Exit 0 = complete; 1 = gaps found.
"""
import json
import os
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
RESULTS = REPO / "results"

# canonical dir names per script
DM = {
    "master_boolean.py": "master_boolean",
    "boolean_co-evolution.py": "boolean_results",
    "nary_kmap_co-evolution.py": "nary_kmap_results",
    "position_kmap_coevolution.py": "position_kmap_results",
    "run_allseq_analysis.py": "full_position_results",
    "run_kmap_analysis.py": "kmap_results",
    "flipped_boolean_coevolution.py": "flipped_boolean_results",
    "kmap_boolean_coevolution.py": "kmap_boolean_coevolution",
    "variable_position_coevolution.py": "variable_position_results",
    "predictive_constraint_function.py": "constraint_function_results",
    "allseq_constraint_function.py": "allseq_constraint_results",
    "dca_boolean_coevolution.py": "dca_boolean_results",
    "dca_mf_analysis.py": "dca_results",
    "perplexity_coevolution.py": "perplexity_results",
    "advanced_co-evolution_analysis.py": "advanced_analysis_results",
    "full_length_analysis.py": "full_length_results",
    "gpu_full_analysis.py": "full_gpu_results",
    "create_mi_heatmap.py": "mi_heatmap",
}
GEN = ["generate_co-evolution_md.py", "generate_full_analysis_md.py",
       "generate_full_pipeline_doc.py"]

FULL = list(DM) + GEN
CORE = ["full_length_analysis.py", "master_boolean.py",
        "kmap_boolean_coevolution.py", "flipped_boolean_coevolution.py",
        "nary_kmap_co-evolution.py", "perplexity_coevolution.py",
        "run_kmap_analysis.py", "variable_position_coevolution.py",
        "position_kmap_coevolution.py", "run_allseq_analysis.py",
        "advanced_co-evolution_analysis.py"]


def main():
    gaps = []
    ok = 0
    for ds_dir in sorted(RESULTS.iterdir()):
        if not ds_dir.is_dir() or ds_dir.name in ("consistency", "contacts"):
            continue
        ds = ds_dir.name
        targets = sorted(d for d in ds_dir.iterdir() if d.is_dir())
        scripts = CORE if ds == "psicov150" else FULL
        for t in targets:
            for s in scripts:
                dn = DM.get(s, s.replace(".py", ""))
                logf = t / dn / "run_log.json"
                if not logf.exists():
                    gaps.append(f"{ds}/{t.name}/{s}: NO RUN LOG")
                    continue
                d = json.load(open(logf))
                if d.get("status") != "ok":
                    gaps.append(f"{ds}/{t.name}/{s}: {d.get('status')} "
                                f"({d.get('seconds')}s) {str(d.get('tail'))[-80:]}")
                else:
                    ok += 1
    print(f"completeness: {ok} ok, {len(gaps)} gaps/failures")
    for g in gaps:
        print("  -", g)
    return 0 if not gaps else 1


if __name__ == "__main__":
    sys.exit(main())
