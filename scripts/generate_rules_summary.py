#!/usr/bin/env python3
"""
PHASE 6: Consolidate rules, prime implicants, Boolean equations and K-map
tables per dataset into My_Own_Interpretation_Across_New_Datasets/<dataset>/.

For each target with results, gather from the per-script JSON outputs:
  - master_boolean: PIs, essential rules (main K-map)
  - kmap_boolean_coevolution: Boolean equations (expression strings)
  - flipped_boolean_coevolution: forbidden rules
  - nary_kmap_results: n-ary PIs/essential
  - perplexity_results: ratio summaries
Write: rules_and_tables.md + raw JSON dump per dataset.
"""

import json
import os
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
RESULTS = REPO / "results"
OUT = REPO / "My_Own_Interpretation_Across_New_Datasets"


def load(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def main():
    os.makedirs(OUT, exist_ok=True)
    for ds_dir in sorted(RESULTS.iterdir()):
        if not ds_dir.is_dir() or ds_dir.name in ("consistency", "contacts"):
            continue
        ds = ds_dir.name
        ds_out = OUT / ds
        os.makedirs(ds_out, exist_ok=True)
        rows = []
        for tgt_dir in sorted(ds_dir.iterdir()):
            if not tgt_dir.is_dir():
                continue
            tgt = tgt_dir.name
            mb = load(tgt_dir / "master_boolean" / "master_boolean_summary.json")
            kb = load(tgt_dir / "kmap_boolean_coevolution" / "boolean_functions.json")
            fb = load(
                tgt_dir / "flipped_boolean_results" / "flipped_boolean_summary.json"
            )
            ny = load(tgt_dir / "nary_kmap_results" / "nary_analysis_summary.json")
            px = load(tgt_dir / "perplexity_results" / "perplexity_summary.json")
            entry = {"target": tgt}
            if mb:
                entry["master_boolean"] = {
                    "variable_positions": mb.get("variable_positions"),
                    "co_evolving_pairs": mb.get("co_evolving_pairs"),
                    "n_prime_implicants": mb.get("total_prime_implicants"),
                    "n_essential": mb.get("essential_prime_implicants"),
                    "rules": mb.get("inferences", [])[:10],
                    "top_pairs": mb.get("top_co_evolving_pairs", [])[:5],
                }
            if kb:
                entry["boolean_equations"] = [
                    {
                        "aa_i": r.get("aa_i"),
                        "aa_j": r.get("aa_j"),
                        "expr": r.get("expr"),
                        "coverage": r.get("coverage"),
                    }
                    for r in kb.get("rules", [])[:10]
                ]
            if fb:
                entry["flipped"] = {
                    "n_forbidden": fb.get("n_forbidden_rules")
                    or fb.get("total_rules")
                    or len(fb.get("rules", [])),
                    "sample": fb.get("rules", [])[:5],
                }
            if ny:
                ny_min = ny.get("minimization") or {}
                entry["nary"] = {
                    "n_PI": ny_min.get("n_prime_implicants")
                    or ny.get("n_prime_implicants"),
                    "n_essential": ny_min.get("n_essential") or ny.get("n_essential"),
                }
            if px:
                entry["perplexity"] = px.get("top_pairs") or px.get("pairs", [])[:5]
            if mb or kb or fb or ny or px:
                rows.append(entry)
        if not rows:
            continue
        with open(ds_out / "rules_and_tables.json", "w") as f:
            json.dump(rows, f, indent=1)
        # markdown
        md = [f"# {ds} — Rules, Prime Implicants & Boolean Equations", ""]
        for r in rows:
            md.append(f"## {r['target']}")
            if "master_boolean" in r:
                m = r["master_boolean"]
                md.append(
                    f"- variable positions: {m['variable_positions']}, "
                    f"pairs: {m['co_evolving_pairs']}, "
                    f"**PIs: {m['n_prime_implicants']}, essential: {m['n_essential']}**"
                )
                for rule in m["rules"]:
                    md.append(
                        f"  - IF pos {rule.get('pos_i')}={rule.get('aa_i')} "
                        f"AND pos {rule.get('pos_j')}={rule.get('aa_j')} → co-evolutionary"
                    )
            if "boolean_equations" in r:
                md.append("- Boolean equations (kmap_boolean):")
                for e in r["boolean_equations"]:
                    md.append(f"  - `{e['expr']}` (coverage {e.get('coverage')})")
            if "flipped" in r:
                md.append(f"- Flipped/forbidden rules: {r['flipped']['n_forbidden']}")
            if "nary" in r:
                md.append(
                    f"- N-ary: {r['nary']['n_PI']} PIs, {r['nary']['n_essential']} essential"
                )
            md.append("")
        with open(ds_out / "rules_and_tables.md", "w") as f:
            f.write("\n".join(md))
        print(f"{ds}: {len(rows)} targets -> {ds_out}/rules_and_tables.{{json,md}}")
    print("Phase 6 consolidation done")


if __name__ == "__main__":
    main()
