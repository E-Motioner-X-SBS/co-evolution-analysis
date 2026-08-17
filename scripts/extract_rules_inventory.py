#!/usr/bin/env python3
"""
Extract ALL co-evolution rule types from the 23-script analysis into ONE
clean JSON (results/contacts/rules_inventory.json):

  1. position_pairs  — from every script that emits position pairs
     (master_boolean, full_length, gpu_full, run_allseq, position_kmap,
      allseq_constraint, variable_position, perplexity, network edges)
  2. residue_rules   — the 36 prime-implicant cells (pos_i=aa_i ∧ pos_j=aa_j)
     parsed from KMAP_BOOLEAN_TABLES.md on-set cells
  3. essential_rules — the 2-rule minimal irredundant core
  4. flipped_rules   — the 490 forbidden residue pairs
  5. perplexity_pairs— determinism-ranked pairs

Honest note: nary_kmap and boolean_co-evolution MOTIFS are protein-wide
dipeptide preferences (row_aa/col_aa only, no positions) — recorded separately;
they cannot be mapped to 3D positions and are excluded from position-level
structural validation.
"""
import json
import re
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
BASE = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
OUT = REPO / "results" / "contacts" / "rules_inventory.json"


def load(p):
    try:
        return json.load(open(BASE / p))
    except Exception:
        return None


def pairs_from_list(items, key_i="pos_i", key_j="pos_j", score=None):
    out = []
    for it in items or []:
        if key_i in it and key_j in it:
            out.append({"pos_i": it[key_i], "pos_j": it[key_j],
                        **({"score": it[score]} if score and score in it else {})})
    return out


def main():
    inv = {"sources": {}, "position_pairs": {}, "residue_rules": [],
           "essential_rules": [], "flipped_rules": [], "perplexity_pairs": [],
           "nary_motifs": [], "boolean_motifs": [], "network": None}

    # 1. position pairs per source script
    sources = {
        "master_boolean": load("master_boolean/master_boolean_summary.json"),
        "full_length": load("full_length_results/full_length_summary.json"),
        "gpu_full": load("full_gpu_results/gpu_summary.json"),
        "run_allseq": load("full_position_results/full_analysis_summary.json"),
        "position_kmap": load("position_kmap_results/position_kmap_summary.json"),
        "allseq_constraint": load("allseq_constraint_results/allseq_constraint_summary.json"),
        "variable_position": load("variable_position_results/variable_position_summary.json"),
        "perplexity": load("perplexity_results/perplexity_summary.json"),
    }
    p = inv["position_pairs"]
    mb = sources["master_boolean"]
    if mb:
        p["master_boolean"] = pairs_from_list(mb.get("top_co_evolving_pairs"), score="mi")
    fl = sources["full_length"]
    if fl:
        p["full_length"] = pairs_from_list(fl.get("top_20_co_evolving_pairs"))
    ra = sources["run_allseq"]
    if ra:
        p["run_allseq"] = pairs_from_list(ra.get("top_mi_pairs"), score="mi")
    pk = sources["position_kmap"]
    if pk:
        p["position_kmap"] = pairs_from_list(pk.get("mi_results"), score="mi")
    ac = sources["allseq_constraint"]
    if ac:
        p["allseq_constraint"] = pairs_from_list(ac.get("results"), score="accuracy")
    vp = sources["variable_position"]
    if vp:
        ce = vp.get("co_evolving_pairs") or []
        p["variable_position"] = [
            {"pos_i": int(c[0]), "pos_j": int(c[1]), "score": c[2]}
            for c in ce if isinstance(c, (list, tuple)) and len(c) >= 3]
    px = sources["perplexity"]
    if px and px.get("results"):
        inv["perplexity_pairs"] = [
            {"pos_i": r.get("pos_i"), "pos_j": r.get("pos_j"),
             "ratio": r.get("ratio"), "mi": r.get("mi")}
            for r in px.get("results", []) if isinstance(r, dict) and r.get("pos_i") is not None][:20]
        p["perplexity"] = inv["perplexity_pairs"]
    else:
        # The saved perplexity summary is empty (verified) — compute the
        # combined MI+perplexity ranking properly from the alignment via the
        # shared module (identical to what every script's combined section does).
        try:
            import sys as _sys
            _sys.path.insert(0, str(REPO))
            _sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
            _sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
            from coevolution_shared import (load_position_arrays,
                                            compute_entropy_vectorized,
                                            combined_pair_scores)
            _pa, _na, _fl = load_position_arrays(
                fasta_path=str(BASE / "Spike_protein.aln-fasta"),
                max_pos=None, aligned=True, clear_cache=True)
            _ent = compute_entropy_vectorized(_pa, _na, _fl)
            _var = [p for p in range(_fl) if _ent[p] > 0.3]
            _pairs = [(i, j) for a, i in enumerate(_var)
                      for j in _var[a + 1:] if j - i <= 30]
            scored = combined_pair_scores(_pa, _pairs, _na, _ent)
            inv["perplexity_pairs"] = [
                {"pos_i": s["pos_i"], "pos_j": s["pos_j"],
                 "ratio": round(s["ratio"], 3), "mi": round(s["mi"], 3),
                 "combined": round(s["combined"], 3)}
                for s in scored[:20]]
            p["perplexity"] = inv["perplexity_pairs"]
        except Exception as e:
            print(f"  [note] perplexity computation failed: {e}")
    gpu = sources["gpu_full"]
    if gpu:
        cp = gpu.get("couplings")
        if isinstance(cp, dict):
            cp = cp.get("top_couplings")
        if isinstance(cp, list):
            p["gpu_full"] = [
                {"pos_i": c.get("pos_i"), "pos_j": c.get("pos_j"), "score": c.get("mi")}
                for c in cp if isinstance(c, dict) and "pos_i" in c]
    # ── Extended sources: cover ALL pair-emitting scripts (14 total) ──
    # 6. mfDCA direct-information pairs (dca_mf_analysis)
    dca = load("dca_results/dca_mf_summary.json")
    if dca:
        di = dca.get("top_20_di") or dca.get("top_20_apc") or []
        p["dca_mf"] = [
            {"pos_i": r.get("pos_i"), "pos_j": r.get("pos_j"), "score": r.get("score")}
            for r in di if isinstance(r, dict) and r.get("pos_i") is not None]
    # 7. local-precision -> Boolean pairs (dca_boolean_coevolution)
    dcb = load("dca_boolean_results/dca_boolean_summary.json")
    if dcb:
        p["dca_boolean"] = [
            {"pos_i": r.get("pos_i"), "pos_j": r.get("pos_j"),
             "score": r.get("n_prime_implicants")}
            for r in dcb.get("results", []) if isinstance(r, dict) and r.get("pos_i") is not None]
    # 8. predictive constraint-function pairs (predictive_constraint_function)
    pcf = load("constraint_function_results/constraint_function_summary.json")
    if pcf:
        p["predictive_constraint"] = [
            {"pos_i": r.get("pos_i"), "pos_j": r.get("pos_j"), "score": r.get("mi")}
            for r in pcf.get("results", []) if isinstance(r, dict) and r.get("pos_i") is not None]
    # 9. flipped (forbidden) UNIQUE position pairs (flipped_boolean_coevolution)
    fbp = load("flipped_boolean_results/flipped_boolean_summary.json")
    if fbp:
        seen_fp = set()
        fp_pairs = []
        for r in fbp.get("rules", []):
            key = (r.get("pos_i"), r.get("pos_j"))
            if key in seen_fp or key[0] is None:
                continue
            seen_fp.add(key)
            fp_pairs.append({"pos_i": key[0], "pos_j": key[1], "score": r.get("mi")})
        p["flipped_positions"] = fp_pairs
    # 10. K-map Boolean essential-rule pairs (kmap_boolean_coevolution)
    kbr = load("kmap_boolean_coevolution/boolean_functions.json")
    if kbr:
        p["kmap_boolean"] = [
            {"pos_i": r.get("pos_i"), "pos_j": r.get("pos_j"), "score": r.get("mi")}
            for r in kbr.get("rules", []) if isinstance(r, dict) and r.get("pos_i") is not None]
    # 11. co-evolution network edges (advanced_co-evolution_analysis)
    net = load("advanced_analysis_results/coevolution_network.json")
    if net:
        inv["network"] = {"nodes": len(net.get("nodes", [])),
                          "edges": net.get("edges", [])[:20]}
        p["network"] = [
            {"pos_i": e.get("source"), "pos_j": e.get("target"), "score": e.get("mi")}
            for e in net.get("edges", []) if isinstance(e, dict) and e.get("source") is not None]

    # 2. residue rules from the K-map tables (on-set cells)
    tables = (BASE / "kmap_boolean_coevolution" / "KMAP_BOOLEAN_TABLES.md")
    if tables.exists():
        txt = tables.read_text()
        # sections: "## Position Pair (i, j)"
        for m in re.finditer(r"## Position Pair \((\d+), (\d+)\)", txt):
            i, j = int(m.group(1)), int(m.group(2))
            sec_end = txt.find("\n## ", m.end())
            sec = txt[m.end(): sec_end if sec_end != -1 else len(txt)]
            # parse the 20x20 table: find header row then rows
            tbl = re.search(r"\| AA_i \\ AA_j \|(.+?)\n((?:\|.*\n)+)", sec)
            if not tbl:
                continue
            cols = [c.strip() for c in tbl.group(1).split("|") if c.strip()]
            for row in tbl.group(2).strip().split("\n"):
                cells = [c.strip() for c in row.split("|")]
                if len(cells) < 3:
                    continue
                aa_i = cells[1].strip()
                for k, c in enumerate(cells[2:2 + len(cols)]):
                    if c == "**1**" or c == "1":
                        inv["residue_rules"].append(
                            {"pos_i": i, "pos_j": j, "aa_i": aa_i,
                             "aa_j": cols[k], "type": "on-set"})
    # 3. essential rules
    if mb:
        inv["essential_rules"] = [{"pos_i": r["pos_i"], "pos_j": r["pos_j"],
                                   "aa_i": r["aa_i"], "aa_j": r["aa_j"],
                                   "mi": r.get("mi")} for r in mb.get("inferences", [])]
    # 4. flipped (forbidden) rules
    fb = load("flipped_boolean_results/flipped_boolean_summary.json")
    if fb:
        seen = set()
        for r in fb.get("rules", []):
            key = (r.get("pos_i"), r.get("pos_j"), r.get("aa_i"), r.get("aa_j"))
            if key in seen:
                continue
            seen.add(key)
            inv["flipped_rules"].append({"pos_i": r.get("pos_i"), "pos_j": r.get("pos_j"),
                                         "aa_i": r.get("aa_i"), "aa_j": r.get("aa_j"),
                                         "essential": r.get("is_essential", False)})
    # 5. n-ary / boolean motifs (protein-wide — recorded, not position-mappable)
    ny = load("nary_kmap_results/nary_analysis_summary.json")
    if ny:
        inv["nary_motifs"] = [{"row_aa": m.get("row_aa"), "col_aa": m.get("col_aa"),
                               "coverage": m.get("coverage"), "count": m.get("motif_count")}
                              for m in ny.get("motifs", [])[:50]]
    bc = load("boolean_results/boolean_analysis_summary.json")
    if bc:
        inv["boolean_motifs"] = [{"row_aa": m.get("row_aa"), "col_aa": m.get("col_aa"),
                                  "coverage": m.get("coverage"),
                                  "count": m.get("motif_count")}
                                 for m in bc.get("motifs", [])[:50]]

    # summary
    inv["summary"] = {
        "position_pair_sources": {k: len(v) for k, v in inv["position_pairs"].items()},
        "n_residue_rules": len(inv["residue_rules"]),
        "n_essential": len(inv["essential_rules"]),
        "n_flipped": len(inv["flipped_rules"]),
        "n_nary_motifs": len(inv["nary_motifs"]),
        "n_boolean_motifs": len(inv["boolean_motifs"]),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(inv, open(OUT, "w"), indent=1)
    print(json.dumps(inv["summary"], indent=1))
    print("saved ->", OUT)


if __name__ == "__main__":
    main()
