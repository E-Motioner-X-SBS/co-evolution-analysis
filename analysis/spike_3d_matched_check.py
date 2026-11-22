#!/usr/bin/env python3
"""Spike co-evolving pairs in 3D: separation-matched controls, templates as units.

scripts/validate_3d_complete.py reported the Omicron "contact" pairs
(500,503), (503,507), (407,410), (210,214), (212,215), (210,215) at 3.4-10.5x
enrichment, p <= 1e-8.  Two problems with that test:

  * the control is pooled over the separations of ALL rule pairs, while
    these pairs are 3-5 residues apart, where a Cb contact is close to
    automatic (helix i,i+3/i+4, turns);
  * the 299 models are SwissModel homology models built on a handful of
    templates, so their contacts are not independent observations: a binomial
    test over models multiplies one template's geometry by its model count.

For each pair this script compares the pair's contact rate with the rate of
all residue pairs at the SAME sequence separation in the same model (global)
and within +-30 residues of the pair (local, same domain), and then repeats
the test with templates as the unit: a template "has the contact" when most
of its models do; the binomial test is over templates against the mean
separation-matched rate.

Pairs: the 10 master_boolean co-evolving pairs plus the six pairs above.
Columns are 0-based alignment columns; residues come from the same
column -> model-residue map as rules3d_common (contact: Cb-Cb < 8 A, chain A).

Output: analysis/spike_3d_matched_check.json
Run:    python analysis/spike_3d_matched_check.py
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]

import numpy as np  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import rules3d_common as rc  # noqa: E402

LOCAL = 30
EXTRA = [(500, 503), (503, 507), (407, 410), (210, 214), (212, 215), (210, 215)]


def best_template(uid: str) -> str:
    tsum = rc.SWISS / "batch" / "pdbs" / uid / "templates_summary.json"
    best = max(json.load(open(tsum)), key=lambda m: m.get("gmqe") or 0)
    return best.get("template_qualified_name") or "?"


def main() -> int:
    mb = json.load(open(REPO / "master_boolean" / "master_boolean_summary.json"))
    pairs = sorted({(p["pos_i"], p["pos_j"]) for p in mb["top_co_evolving_pairs"]} | set(EXTRA))
    stats = {p: {"by_template": defaultdict(lambda: {"models": 0, "contacts": 0,
                                                      "exp_global": [], "exp_local": []}),
                 "seps": set()} for p in pairs}
    models = rc.OmicronModels()
    n_models = 0
    templates = defaultdict(int)
    for uid, _aligned, cmap, coords in models.iter_models():
        n_models += 1
        tpl = best_template(uid)
        templates[tpl] += 1
        chain = sorted(coords)[0]
        res = coords[chain]
        keys = np.array(sorted(res))
        xyz = np.array([res[k] for k in keys])
        index = {int(k): n for n, k in enumerate(keys)}
        for (ci, cj) in pairs:
            ri, rj = cmap.get(ci), cmap.get(cj)
            if ri is None or rj is None or ri not in index or rj not in index:
                continue
            s = rj - ri
            a = np.array([index[k] for k in keys if k + s in index])
            b = np.array([index[int(keys[n]) + s] for n in a])
            dd = np.linalg.norm(xyz[a] - xyz[b], axis=1)
            near = np.abs(keys[a] - ri) <= LOCAL
            t = stats[(ci, cj)]["by_template"][tpl]
            t["models"] += 1
            t["contacts"] += int(np.linalg.norm(xyz[index[ri]] - xyz[index[rj]]) < rc.CONTACT)
            t["exp_global"].append(float((dd < rc.CONTACT).mean()))
            t["exp_local"].append(float((dd[near] < rc.CONTACT).mean()) if near.any() else np.nan)
            stats[(ci, cj)]["seps"].add(int(s))

    rows = []
    for p in pairs:
        bt = stats[p]["by_template"]
        m = sum(t["models"] for t in bt.values())
        if m == 0:
            rows.append({"pair": list(p), "models": 0})
            continue
        c = sum(t["contacts"] for t in bt.values())
        eg = float(np.mean([x for t in bt.values() for x in t["exp_global"]]))
        el = float(np.nanmean([x for t in bt.values() for x in t["exp_local"]]))
        n_t = len(bt)
        t_contact = sum(t["contacts"] * 2 > t["models"] for t in bt.values())
        e_t = float(np.mean([np.mean(t["exp_global"]) for t in bt.values()]))
        e_tl = float(np.nanmean([np.nanmean(t["exp_local"]) for t in bt.values()]))
        rows.append({
            "pair": list(p), "separation": sorted(stats[p]["seps"]),
            "models": m, "contact_fraction": c / m,
            "expected_same_separation": eg, "expected_same_separation_local": el,
            "enrichment_global": (c / m) / eg if eg else None,
            "enrichment_local": (c / m) / el if el else None,
            "p_models_vs_global": binomtest(c, m, min(max(eg, 1e-12), 1.0), "greater").pvalue,
            "templates": n_t, "templates_with_contact": t_contact,
            "p_templates_vs_global": binomtest(t_contact, n_t, min(max(e_t, 1e-12), 1.0), "greater").pvalue,
            "p_templates_vs_local": binomtest(t_contact, n_t, min(max(e_tl, 1e-12), 1.0), "greater").pvalue,
        })
        r = rows[-1]
        print(f"{str(p):12s} sep {r['separation']} models {m:3d} contact {r['contact_fraction']:.2f} "
              f"exp {eg:.2f}/{el:.2f} enr {r['enrichment_global']:.2f}/{r['enrichment_local']:.2f} "
              f"templates {t_contact}/{n_t} p_tpl={r['p_templates_vs_local']:.2g}", flush=True)
    out = {"n_models": n_models, "skipped": len(models.skipped),
           "templates": dict(sorted(templates.items(), key=lambda kv: -kv[1])),
           "contact": f"Cb-Cb < {rc.CONTACT} A, chain A", "local_window": LOCAL, "pairs": rows}
    (REPO / "analysis" / "spike_3d_matched_check.json").write_text(json.dumps(out, indent=1, default=str))
    print("templates:", out["templates"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
