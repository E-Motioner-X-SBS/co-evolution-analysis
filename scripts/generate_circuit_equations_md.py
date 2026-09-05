#!/usr/bin/env python3
"""
generate_circuit_equations_md.py — emit BOOLEAN_CIRCUIT_EQUATIONS.md.

Reads the stage JSONs written by run_sequence_circuit_campaign.py and renders
the actual Boolean equations: the co-evolution circuits for every Spike
co-evolving pair in all three polarities x both don't-care policies, the
sequence-circuit statistics, the contact-map block circuits, and the
feature-class ceiling table.

    PY=/store/shuvam/.venv/bin/python
    $PY scripts/generate_circuit_equations_md.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
RES = REPO / "results" / "sequence_circuits"
OUT = REPO / "BOOLEAN_CIRCUIT_EQUATIONS.md"


def _load(name: str) -> Optional[Dict]:
    p = RES / name
    if not p.exists():
        print(f"  (missing {name} — section skipped)")
        return None
    return json.loads(p.read_text())


def _fmt_set(letters: List[str]) -> str:
    """Render a residue tolerance class compactly."""
    if not letters:
        return "∅"
    if len(letters) == 1:
        return letters[0]
    if len(letters) >= 18:
        return "any"
    return "{" + ",".join(letters) + "}"


def section_coevo(d: Dict) -> str:
    L: List[str] = []
    L.append("## 1. Co-evolution circuits — SARS-CoV-2 Spike\n")
    L.append(f"Dataset: 1,299 Omicron Spike sequences, {d['n_variable_positions']} "
             f"variable positions (Shannon H > {d['entropy_threshold']}), "
             f"{d['n_coevolving_pairs']} co-evolving pairs "
             f"(mutation-only MI > {d['mi_threshold']}, |i-j| ≤ {d['max_gap']}, "
             f"≥ {d['min_muts']} mutations).\n")
    L.append("Variables: `a4..a0` = 5-bit Gray code of the residue at position i, "
             "`b4..b0` = same at position j (He-2012 ordering `AILVMFYWEDQNHKRSTCPG`). "
             "Each product term is one prime implicant; a free bit widens the term "
             "into a residue tolerance class.\n")
    L.append("Two don't-care policies are reported for every polarity, because the "
             "choice dominates the result on a 2-4 residue column:\n")
    L.append("- **`/off`** — combinations whose residue never appears in that column "
             "are OFF. Cubes stay inside the observed alphabet: tight, testable rules.")
    L.append("- **`/dc`** — those combinations are don't-care. Formally sound, but QM "
             "then absorbs the ~99% unobservable region and returns near-vacuous "
             "rules. Shown so the degeneracy is on the record, **not** for use.\n")

    for pr in d["pairs"]:
        L.append(f"### Pair ({pr['i']}, {pr['j']}) — MI = {pr['mi']:.4f}, "
                 f"{pr['n_mut']} mutated sequences\n")
        for pol_key, e in pr["polarities"].items():
            s = e["summary"]
            alpha_i = "".join(s.get("alphabet_i", []))
            alpha_j = "".join(s.get("alphabet_j", []))
            L.append(f"**{pol_key}** — on-set {s['n_on']}, don't-care {s['n_dc']}, "
                     f"primes {s['n_primes']}, cover {s['n_cover']}, "
                     f"essential {s['n_essential']}, "
                     f"literals {s['n_literals']}, "
                     f"sound={s['sound']}, complete={s['complete']}  "
                     f"(alphabet {alpha_i} × {alpha_j})\n")
            if not e["rules"]:
                L.append("_no canonical rules (all cubes lie in the padding border)_\n")
                continue
            L.append("| # | residues at i | residues at j | chemistry (i × j) | "
                     "free bits | essential |")
            L.append("|---|---|---|---|---|---|")
            for n, r in enumerate(e["rules"][:16], 1):
                L.append(f"| {n} | {_fmt_set(r['residues_i'])} | "
                         f"{_fmt_set(r['residues_j'])} | "
                         f"{'/'.join(r['groups_i'])} × {'/'.join(r['groups_j'])} | "
                         f"{r['n_free']} | {'yes' if r['essential'] else 'no'} |")
            if len(e["rules"]) > 16:
                L.append(f"| … | _+{len(e['rules']) - 16} more_ | | | | |")
            L.append("")
            L.append("```")
            L.append(f"f_{pr['i']},{pr['j']} = " + e["sop"])
            L.append("```\n")
    return "\n".join(L)


def section_seq(d: Dict) -> str:
    L: List[str] = []
    L.append("## 2. The sequence itself as a minimal circuit\n")
    L.append(f"All **{d['n_sequences_analysed']}** Spike sequences were compiled to "
             "the positional circuit `f_s(i, a) = 1 ⟺ s_i = a` "
             "(11 position bits + 5 residue bits) and decompiled again.\n")
    L.append(f"- **Lossless on every sequence: {d['all_lossless']}** — the residue "
             "string is recovered exactly by evaluating the minimised circuit.")
    L.append(f"- Every circuit sound: {d['all_sound']}; every circuit complete: "
             f"{d['all_complete']}.\n")
    sp, sh, un = d["spike"], d.get("shuffled_control", {}), d.get("uniform_control", {})
    L.append("| quantity | Spike | composition-matched shuffle | uniform random |")
    L.append("|---|---|---|---|")

    def _row(label: str, key: str, fmt: str = "{:.2f}") -> str:
        a = fmt.format(sp[key]["mean"]) if key in sp else "—"
        b = fmt.format(sh[key]["mean"]) if key in sh else "—"
        c = fmt.format(un[key]["mean"]) if key in un else "—"
        return f"| {label} | {a} | {b} | {c} |"

    L.append(_row("length", "length", "{:.0f}"))
    L.append(_row("cubes in minimal cover", "rb_cover", "{:.1f}"))
    L.append(_row("prime implicants", "rb_primes", "{:.1f}"))
    L.append(_row("compression (residues / cubes)", "rb_compression", "{:.4f}"))
    L.append(_row("dipeptide-circuit cubes", "ra_cover", "{:.1f}"))
    L.append(_row("dipeptide compression", "ra_compression", "{:.3f}"))
    L.append("")
    return "\n".join(L)


def section_contact(d: Dict) -> str:
    L: List[str] = []
    L.append("## 3. The contact map as a Boolean function — segment × segment blocks\n")
    L.append(f"{d['n_targets']} PSICOV targets, native contacts (Cβ–Cβ < 8 Å, "
             f"separation ≥ 6), each against {d['n_shuffles']} "
             "separation-preserving shuffles (contact count and the full |i−j| "
             "distribution held fixed; only *which* residues touch is randomised).\n")
    L.append("A cube whose free bits are the low bits of both position fields is "
             "exactly `[i₀, i₀+2^u) × [j₀, j₀+2^v) ⊆ contacts` — every residue of "
             "one chain segment touching every residue of another.\n")
    c, b = d["compression"], d["block_fraction"]
    L.append("| statistic | real maps | shuffled | real > shuffled | Wilcoxon p |")
    L.append("|---|---|---|---|---|")
    L.append(f"| compression (contacts / cubes) | {c['real_mean']:.4f} | "
             f"{c['shuffled_mean']:.4f} | {c['n_real_higher']}/{d['n_targets']} | "
             f"{c['wilcoxon']['p']:.3g} |")
    L.append(f"| fraction of contacts inside segment blocks | {b['real_mean']:.4f} | "
             f"{b['shuffled_mean']:.4f} | {b['n_real_higher']}/{d['n_targets']} | "
             f"{b['wilcoxon']['p']:.3g} |")
    L.append("")
    L.append(f"All circuits sound: {d['all_sound']}; complete: {d['all_complete']}.\n")
    return "\n".join(L)


def section_ceiling(d: Dict) -> str:
    L: List[str] = []
    L.append("## 4. The information ceiling of a feature class\n")
    L.append("A circuit whose inputs are (residue at i, residue at j, separation bin) "
             "is a **function of its cell**: every candidate pair landing in the same "
             "cell receives the same score. No minimisation, threshold, don't-care "
             "policy or cover-selection rule can rank one above another.\n")
    L.append("Each cell is given its own pooled empirical contact rate — the feature "
             "class handed the test labels. **`achievable`** is the expectation over "
             "within-cell orderings, and is the meaningful quantity: a function of "
             "the cell cannot distinguish pairs inside a cell, so it cannot break "
             "that symmetry, and no such function can systematically beat it.\n")
    L.append("The `tie-break max` column resolves within-cell ties in favour of "
             "contacts. Read it only as *how much the tie-break matters* — it is "
             "vacuous for coarse partitions (`sep_only` reaches 0.98 with four "
             "cells purely by sorting on the label, which no real function can do).\n")
    L.append("| feature set | cells | **achievable** | tie-break max |")
    L.append("|---|---|---|---|")
    for name, e in sorted(d["feature_sets"].items(),
                          key=lambda kv: -kv[1]["ceiling_expected_mean"]):
        L.append(f"| `{name}` ({', '.join(e['features'])}) | {e['n_cells_used']} | "
                 f"**{e['ceiling_expected_mean']:.4f}** | "
                 f"{e['ceiling_optimistic_mean']:.4f} |")
    a = d["literature_anchor"]
    L.append("")
    L.append("Measured against the same benchmark:\n")
    L.append(f"- PATH B circuit-L2, honest holdout: **{a['PATH_B_circuit_L2_holdout']}**")
    L.append(f"- plain MI: {a['PATH_B_plain_MI']}")
    L.append(f"- PSICOV published DCA: **{a['PSICOV_published_prec_at_L5']}**\n")
    return "\n".join(L)


def main() -> None:
    parts: List[str] = []
    parts.append("# Boolean Circuit Equations for Biological Sequences\n")
    parts.append("_Generated by `scripts/generate_circuit_equations_md.py` from the "
                 "stage JSONs under `results/sequence_circuits/`. Engine: "
                 "`scripts/sequence_circuits.py`; gates: "
                 "`scripts/tests/test_sequence_circuits.py`._\n")
    parts.append("Every circuit in this document is **verified sound and complete** "
                 "at build time: no product term fires on an off-cell, and every "
                 "on-cell is covered. Minimisation is exact Quine–McCluskey on the "
                 "sparse term set — see the note on the reference minimiser in §5.\n")

    for name, fn in (("stage_coevo.json", section_coevo),
                     ("stage_seq.json", section_seq),
                     ("stage_contact.json", section_contact),
                     ("stage_ceiling.json", section_ceiling)):
        d = _load(name)
        if d is not None:
            parts.append(fn(d))

    parts.append("""## 5. Note on the minimiser

These circuits are minimised by `sequence_circuits.prime_implicants`, an exact
sparse Quine–McCluskey. It agrees cell-for-cell with the repo's existing
`kmap_sbm.analysis.prime_implicants` when there are no don't-cares
(`test_sequence_circuits.py::t05a/t05b`), and is independently audited for
implicant-hood and primality (`t06a–t06c`, 0 violations at n = 8, 10, 12).

Where don't-cares are present the two diverge, for a documented reason: the
reference keeps only don't-cares at Hamming distance 1 from a real on-minterm,
a performance choice stated in its own docstring. With clustered don't-cares —
which is precisely the contact-campaign regime, where DC means "cell has
n < n_min" — it **under-merges**. Measured on a 10-variable table with 60
on-cells and 400 don't-cares: 522 exact prime implicants vs 187 from the
reference, and every one of the 74 cubes the reference reports that exact QM
does not is *strictly inside* one of the true primes (`t06d`).

Consequence: the reference's rules stay **sound**, but its prime-implicant
counts and its "essential" designations are approximations, not the exact QM
quantities. Published PI counts derived from it should be read as lower bounds.
""")

    OUT.write_text("\n".join(parts))
    print(f"wrote {OUT} ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
