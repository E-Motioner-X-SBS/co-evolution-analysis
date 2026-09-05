#!/usr/bin/env python3
"""
run_sequence_circuit_campaign.py — can a sequence be written as a minimal
Boolean circuit, and does that circuit carry 3D information?

Stages (each resumable, each writes one JSON under results/sequence_circuits/)

  seq      R-B lossless positional circuits + R-A dipeptide circuits for the
           1,299 Spike sequences, against length/composition-matched random
           controls.  Answers: "is a real protein sequence more circuit-
           compressible than a random one with the same composition?"

  coevo    R-C co-evolution circuits on the Spike co-evolving pairs, in all
           three polarities (standard / flipped / perplexity).  Emits the
           actual Boolean equations.

  contact  R-D — the native contact map AS a Boolean function of position
           bits, over all 150 PSICOV targets, against separation-preserving
           shuffles.  This is PATH F at full scale (it was n=6 before, limited
           by the dense reference minimiser) and with plain-binary position
           coding so that a cube means "segment x segment".

  chem     the physicochemical contact circuit (group_i x group_j x sep)
           rebuilt with EXACT QM, cross-checked against the PATH B numbers.

Usage
-----
    PY=/store/shuvam/.venv/bin/python; export OMP_PROC_BIND=FALSE
    $PY scripts/run_sequence_circuit_campaign.py --stage seq
    $PY scripts/run_sequence_circuit_campaign.py --stage coevo
    $PY scripts/run_sequence_circuit_campaign.py --stage contact
    $PY scripts/run_sequence_circuit_campaign.py --stage chem

Determinism: every RNG is seeded (seed=42 for splits, 20260905 for controls);
re-running a stage reproduces its JSON byte-for-byte.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kmap_structure as ks          # noqa: E402
import sequence_circuits as sc       # noqa: E402

REPO = HERE.parent
OUT = REPO / "results" / "sequence_circuits"
SPIKE = REPO / "Spike_protein.aln-fasta"

CONTROL_SEED = 20260905
SPLIT_SEED = 42

# Spike co-evolution conventions, matching coevolution_shared.py exactly so the
# pair set reproduces the published 21 positions / 10 pairs.
ENTROPY_THRESHOLD = 0.3
MI_THRESHOLD = 0.1
MAX_GAP = 30
MIN_MUTS = 5


def _dump(name: str, obj: Dict) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / name
    p.write_text(json.dumps(obj, indent=1, sort_keys=True, default=str))
    print(f"  -> {p}")
    return p


# ---------------------------------------------------------------------------
# Stage `seq`
# ---------------------------------------------------------------------------

def _circuit_row(codes: Sequence[int], tag: str) -> Dict:
    """R-B and R-A summaries for one sequence, with the losslessness check."""
    cb = sc.seq_circuit_positional(codes, name=tag)
    rec = sc.reconstruct_from_circuit(cb, len(codes))
    ca = sc.seq_circuit_kmer(codes, k=2, name=tag)
    return {
        "tag": tag,
        "length": len(codes),
        "rb_primes": len(cb.primes),
        "rb_cover": len(cb.cover),
        "rb_essential": len(cb.essentials),
        "rb_literals": cb.n_literals,
        "rb_compression": len(codes) / max(1, len(cb.cover)),
        "rb_lossless": rec == list(codes),
        "rb_sound": cb.sound,
        "rb_complete": cb.complete,
        "ra_on": ca.n_on,
        "ra_cover": len(ca.cover),
        "ra_primes": len(ca.primes),
        "ra_compression": ca.n_on / max(1, len(ca.cover)),
        "ra_sound": ca.sound,
        "ra_complete": ca.complete,
    }


def stage_seq(n_seqs: int, n_controls: int) -> Dict:
    print(f"[seq] loading {SPIKE.name}")
    dense, headers = sc.load_fasta_codes(SPIKE)
    n_total, width = dense.shape
    print(f"[seq] MSA {n_total} x {width}")

    rng = np.random.default_rng(CONTROL_SEED)
    rows: List[Dict] = []
    ctrl_rows: List[Dict] = []

    take = min(n_seqs, n_total)
    for r in range(take):
        codes = [int(c) for c in dense[r] if c != ks.GAP_STATE]
        rows.append(_circuit_row(codes, f"spike_{r}"))
        if r < n_controls:
            # Composition-matched control: same multiset of residues, shuffled.
            perm = rng.permutation(len(codes))
            shuffled = [codes[int(k)] for k in perm]
            ctrl_rows.append(_circuit_row(shuffled, f"shuffled_{r}"))
        if (r + 1) % 100 == 0:
            print(f"[seq]   {r + 1}/{take}")

    # Uniform-random control (not composition matched) for a second baseline.
    uni_rows: List[Dict] = []
    med_len = int(np.median([r["length"] for r in rows]))
    for r in range(n_controls):
        codes = [int(x) for x in rng.integers(0, 20, size=med_len)]
        uni_rows.append(_circuit_row(codes, f"uniform_{r}"))

    def _agg(rs: List[Dict], key: str) -> Dict:
        v = np.array([x[key] for x in rs], dtype=float)
        return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1
                else 0.0, "median": float(np.median(v)), "n": len(v)}

    res = {
        "n_sequences_analysed": take,
        "n_controls": n_controls,
        "all_lossless": all(r["rb_lossless"] for r in rows),
        "all_sound": all(r["rb_sound"] and r["ra_sound"] for r in rows),
        "all_complete": all(r["rb_complete"] and r["ra_complete"] for r in rows),
        "spike": {k: _agg(rows, k) for k in
                  ("length", "rb_cover", "rb_primes", "rb_literals",
                   "rb_compression", "ra_on", "ra_cover", "ra_compression")},
        "shuffled_control": {k: _agg(ctrl_rows, k) for k in
                             ("rb_cover", "rb_compression", "ra_cover",
                              "ra_compression")} if ctrl_rows else {},
        "uniform_control": {k: _agg(uni_rows, k) for k in
                            ("rb_cover", "rb_compression", "ra_cover",
                             "ra_compression")} if uni_rows else {},
        "per_sequence": rows[:50],
        "seed": CONTROL_SEED,
    }
    return res


# ---------------------------------------------------------------------------
# Stage `coevo`
# ---------------------------------------------------------------------------

def _entropy_all(dense: np.ndarray) -> np.ndarray:
    return np.array([sc.column_entropy(dense, c) for c in range(dense.shape[1])])


def _mi_mutation_only(dense: np.ndarray, i: int, j: int,
                      ref_i: int, ref_j: int) -> Tuple[float, int]:
    """Mutation-only MI: binary "differs from majority reference" per column.

    Mirrors coevolution_shared.mi_mutation_only — the repo's pair-selection
    statistic (this is what produced the published 10-pair set).
    """
    ci, cj = dense[:, i], dense[:, j]
    ok = (ci != ks.GAP_STATE) & (cj != ks.GAP_STATE) & (ci >= 0) & (cj >= 0)
    a = (ci[ok] != ref_i).astype(np.int8)
    b = (cj[ok] != ref_j).astype(np.int8)
    n = a.size
    if n < 10:
        return 0.0, 0
    joint = np.zeros((2, 2), dtype=np.float64)
    for x in (0, 1):
        for y in (0, 1):
            joint[x, y] = np.sum((a == x) & (b == y))
    n_mut = int(joint[1, :].sum() + joint[:, 1].sum() - joint[1, 1])
    p = joint / n
    pi_ = p.sum(axis=1)
    pj_ = p.sum(axis=0)
    mi = 0.0
    for x in (0, 1):
        for y in (0, 1):
            if p[x, y] > 0 and pi_[x] > 0 and pj_[y] > 0:
                mi += p[x, y] * np.log2(p[x, y] / (pi_[x] * pj_[y]))
    return float(mi), n_mut


def _majority_ref(dense: np.ndarray, col: int) -> int:
    v = dense[:, col]
    v = v[(v != ks.GAP_STATE) & (v >= 0)]
    if v.size == 0:
        return 0
    cnt = Counter(int(x) for x in v)
    best = max(cnt.values())
    return min(c for c, n in cnt.items() if n == best)


def stage_coevo(max_pairs: int) -> Dict:
    print(f"[coevo] loading {SPIKE.name}")
    dense, _ = sc.load_fasta_codes(SPIKE)
    ent = _entropy_all(dense)
    variable = [int(p) for p in np.where(ent > ENTROPY_THRESHOLD)[0]]
    print(f"[coevo] variable positions (H>{ENTROPY_THRESHOLD}): {len(variable)}")

    refs = {p: _majority_ref(dense, p) for p in variable}
    pairs = []
    for idx, pi_ in enumerate(variable):
        for pj in variable[idx + 1:]:
            if abs(pi_ - pj) > MAX_GAP:
                continue
            mi, n_mut = _mi_mutation_only(dense, pi_, pj, refs[pi_], refs[pj])
            if mi > MI_THRESHOLD and n_mut >= MIN_MUTS:
                pairs.append({"i": pi_, "j": pj, "mi": mi, "n_mut": n_mut})
    pairs.sort(key=lambda d: -d["mi"])
    print(f"[coevo] co-evolving pairs (mutation-only MI>{MI_THRESHOLD}): {len(pairs)}")

    out_pairs: List[Dict] = []
    for pr in pairs[:max_pairs]:
        entry = {**pr, "polarities": {}}
        for pol in sc.POLARITIES:
            for unobs in ("off", "dc"):
                circ = sc.coevo_circuit(dense, pr["i"], pr["j"], pol,
                                        unobservable=unobs)
                rules = []
                for imp in sorted(circ.cover, key=lambda p: -p.n_free):
                    dec = sc.decode_residue_pair_implicant(imp)
                    if not dec["slot_i"]["residues"] or \
                            not dec["slot_j"]["residues"]:
                        continue                   # cube lies wholly in padding
                    rules.append({
                        "residues_i": dec["slot_i"]["residues"],
                        "residues_j": dec["slot_j"]["residues"],
                        "groups_i": dec["slot_i"]["groups"],
                        "groups_j": dec["slot_j"]["groups"],
                        "n_free": imp.n_free,
                        "literals": imp.literals(circ.var_names),
                        "essential": (imp.val, imp.mask) in
                                     {(e.val, e.mask) for e in circ.essentials},
                    })
                entry["polarities"][f"{pol}/{unobs}"] = {
                    "summary": circ.summary(),
                    "sop": circ.to_sop(max_terms=12),
                    "rules": rules,
                }
        out_pairs.append(entry)

    return {
        "entropy_threshold": ENTROPY_THRESHOLD,
        "mi_threshold": MI_THRESHOLD,
        "max_gap": MAX_GAP,
        "min_muts": MIN_MUTS,
        "n_variable_positions": len(variable),
        "variable_positions": variable,
        "n_coevolving_pairs": len(pairs),
        "n_pairs_reported": len(out_pairs),
        "pairs": out_pairs,
    }


# ---------------------------------------------------------------------------
# Stage `contact` — the 3D question
# ---------------------------------------------------------------------------

def _shuffle_within_sep_bands(contacts: List[Tuple[int, int]], length: int,
                              rng: np.random.Generator,
                              min_sep: int) -> List[Tuple[int, int]]:
    """Relabel contacts at random WITHIN each separation value.

    Preserves exactly: the number of contacts, and the full distribution of
    |i-j|.  Destroys: which particular residues touch, hence all block/segment
    structure.  This is the control PATH F used, restated for arbitrary L.
    """
    by_sep: Dict[int, int] = defaultdict(int)
    for i, j in contacts:
        by_sep[abs(i - j)] += 1
    out: List[Tuple[int, int]] = []
    for sep, k in by_sep.items():
        starts = np.arange(0, length - sep, dtype=np.int64)
        if starts.size == 0:
            continue
        take = min(k, starts.size)
        pick = rng.choice(starts, size=take, replace=False)
        out.extend((int(s), int(s) + sep) for s in pick)
    return out


def _contact_stats(contacts: List[Tuple[int, int]], length: int) -> Dict:
    circ = sc.contact_circuit_positional(contacts, length, band_as_dc=False)
    segs = [sc.implicant_to_segments(p, circ.meta["pos_bits"], length)
            for p in circ.cover]
    n_block = sum(1 for s in segs
                  if s["i"]["contiguous"] and s["j"]["contiguous"]
                  and s["n_free"] >= 2)
    covered_by_blocks = sum(s["size"] for s in segs
                            if s["i"]["contiguous"] and s["j"]["contiguous"]
                            and s["n_free"] >= 2)
    n_on = circ.n_on
    return {
        "n_contacts": len(set((min(a, b), max(a, b)) for a, b in contacts)),
        "n_on": n_on,
        "n_cover": len(circ.cover),
        "n_primes": len(circ.primes),
        "compression": n_on / max(1, len(circ.cover)),
        "pis_per_contact": len(circ.cover) / max(1, n_on),
        "n_block_cubes": n_block,
        "frac_minterms_in_blocks": covered_by_blocks / max(1, n_on),
        "sound": circ.sound,
        "complete": circ.complete,
    }


def stage_contact(limit: Optional[int], n_shuffles: int) -> Dict:
    targets = ks.list_targets()
    if limit:
        targets = targets[:limit]
    print(f"[contact] {len(targets)} PSICOV targets, {n_shuffles} shuffles each")
    rng = np.random.default_rng(CONTROL_SEED)

    rows: List[Dict] = []
    for n, t in enumerate(targets):
        try:
            pd = ks.load_protein(t, load_dense=False)   # contacts only
        except Exception as exc:                    # noqa: BLE001
            print(f"[contact]   {t}: SKIP ({exc})")
            continue
        L = pd.n_cols
        contacts = [(int(i), int(j)) for i, j in pd.contacts]
        if len(contacts) < 20:
            continue
        real = _contact_stats(contacts, L)
        sh = [_contact_stats(
                  _shuffle_within_sep_bands(contacts, L, rng, ks.MIN_SEPARATION), L)
              for _ in range(n_shuffles)]
        rows.append({
            "target": t, "length": L,
            "real": real,
            "shuffled_mean_compression": float(np.mean([s["compression"] for s in sh])),
            "shuffled_mean_cover": float(np.mean([s["n_cover"] for s in sh])),
            "shuffled_mean_block_frac": float(np.mean([s["frac_minterms_in_blocks"]
                                                       for s in sh])),
        })
        if (n + 1) % 10 == 0:
            print(f"[contact]   {n + 1}/{len(targets)}")

    real_c = np.array([r["real"]["compression"] for r in rows])
    shuf_c = np.array([r["shuffled_mean_compression"] for r in rows])
    real_b = np.array([r["real"]["frac_minterms_in_blocks"] for r in rows])
    shuf_b = np.array([r["shuffled_mean_block_frac"] for r in rows])

    def _wilcoxon(a: np.ndarray, b: np.ndarray) -> Dict:
        """Two-sided Wilcoxon signed-rank via scipy when available."""
        try:
            from scipy.stats import wilcoxon
            st, p = wilcoxon(a, b)
            return {"stat": float(st), "p": float(p)}
        except Exception:                            # noqa: BLE001
            return {"stat": None, "p": None}

    return {
        "n_targets": len(rows),
        "n_shuffles": n_shuffles,
        "compression": {
            "real_mean": float(real_c.mean()), "shuffled_mean": float(shuf_c.mean()),
            "n_real_higher": int((real_c > shuf_c).sum()),
            "wilcoxon": _wilcoxon(real_c, shuf_c),
        },
        "block_fraction": {
            "real_mean": float(real_b.mean()), "shuffled_mean": float(shuf_b.mean()),
            "n_real_higher": int((real_b > shuf_b).sum()),
            "wilcoxon": _wilcoxon(real_b, shuf_b),
        },
        "all_sound": all(r["real"]["sound"] for r in rows),
        "all_complete": all(r["real"]["complete"] for r in rows),
        "per_target": rows,
        "seed": CONTROL_SEED,
    }


# ---------------------------------------------------------------------------
# Stage `chem` — physicochemical contact circuit with EXACT QM
# ---------------------------------------------------------------------------

def stage_chem(limit: Optional[int], t_mult: float, n_min: int) -> Dict:
    """Rebuild the PATH B level-1 circuit using exact QM instead of the
    under-merging reference, and report both rule sets side by side."""
    targets = ks.list_targets()
    if limit:
        targets = targets[:limit]
    folds = ks.protein_folds(targets, n_folds=5, seed=SPLIT_SEED)
    print(f"[chem] {len(targets)} targets, 5 folds")

    def _accumulate(ts: List[str]) -> ks.CellCounts:
        counts1 = ks.new_level1_counts()
        counts2 = ks.new_level2_counts()
        for t in ts:
            try:
                pd = ks.load_protein(t)
            except Exception:                        # noqa: BLE001
                continue
            pairs = ks.pair_universe(pd)
            ks.accumulate_protein(counts1, counts2, pd, pairs)
        return counts1

    fold_rows: List[Dict] = []
    for f, (train, test) in enumerate(folds):
        c1 = _accumulate(train)
        tt = ks.label_truth_table(c1, t_mult=t_mult, n_min=n_min)
        on = [int(i) for i in np.where(tt == 1)[0]]
        dc = [int(i) for i in np.where(tt == -1)[0]]
        exact = sc.prime_implicants(on, dc, 8)
        cover, ess = sc._essential_and_greedy_cover(exact, on)
        off = [int(i) for i in np.where(tt == 0)[0]]
        complete, sound = sc.verify_cover(cover, on, off)
        ref_n = len(sc._ref_pi_keys(on, dc, 8))
        fold_rows.append({
            "fold": f, "n_train": len(train), "n_test": len(test),
            "n_on": len(on), "n_dc": len(dc), "n_off": len(off),
            "exact_primes": len(exact), "reference_primes": ref_n,
            "exact_cover": len(cover), "exact_essential": len(ess),
            "sound": sound, "complete": complete,
        })
        print(f"[chem]   fold {f}: on={len(on)} dc={len(dc)} "
              f"exact_PIs={len(exact)} ref_PIs={ref_n} cover={len(cover)}")

    # Full-data circuit, decoded into human-readable chemistry rules.
    c1 = _accumulate(targets)
    tt = ks.label_truth_table(c1, t_mult=t_mult, n_min=n_min)
    on = [int(i) for i in np.where(tt == 1)[0]]
    dc = [int(i) for i in np.where(tt == -1)[0]]
    off = [int(i) for i in np.where(tt == 0)[0]]
    pis = sc.prime_implicants(on, dc, 8)
    cover, ess = sc._essential_and_greedy_cover(pis, on)
    complete, sound = sc.verify_cover(cover, on, off)

    rules = []
    for imp in sorted(cover, key=lambda p: -p.n_free):
        dec = ks.decode_level1_implicant(
            [(imp.val >> b) & 1 for b in range(7, -1, -1)],
            [(imp.mask >> b) & 1 for b in range(7, -1, -1)])
        gi = dec["groups_i"]
        gj = dec["groups_j"]
        rules.append({
            "groups_i": [ks.GROUP_NAMES[g] for g in gi] if gi is not None else "any",
            "groups_j": [ks.GROUP_NAMES[g] for g in gj] if gj is not None else "any",
            "sep_bins": dec["seps"] if dec["seps"] is not None else "any",
            "n_free": imp.n_free,
            "essential": (imp.val, imp.mask) in {(e.val, e.mask) for e in ess},
        })

    return {
        "t_mult": t_mult, "n_min": n_min,
        "base_rate": float(c1.base_rate()),
        "per_fold": fold_rows,
        "full": {"n_on": len(on), "n_dc": len(dc), "n_off": len(off),
                 "exact_primes": len(pis),
                 "reference_primes": len(sc._ref_pi_keys(on, dc, 8)),
                 "cover": len(cover), "essential": len(ess),
                 "sound": sound, "complete": complete},
        "rules": rules,
        "group_names": ks.GROUP_NAMES,
        "sep_bin_edges": list(ks.SEP_BIN_EDGES),
    }


# ---------------------------------------------------------------------------
# Stage `ceiling` — the information limit of a feature class
# ---------------------------------------------------------------------------
#
# A circuit whose inputs are (residue_i, residue_j, separation-bin) is a
# FUNCTION OF ITS CELL.  Every candidate pair landing in the same cell receives
# the same score, so no minimisation, threshold, or rule-selection scheme can
# rank one above another.  That caps precision at a value determined purely by
# the feature map -- computable by giving the circuit the test labels and
# letting it score each cell by its own empirical contact rate (an ORACLE it
# could never achieve in practice).
#
# The gap between that oracle and PSICOV's published DCA is therefore NOT a
# minimisation failure and cannot be closed by better Boolean machinery; it is
# the feature class running out of information.  Conversely the gap between
# PATH B's 0.179 and the oracle is what better machinery could still win.

FEATURE_SETS = {
    "L1_group_sep":      ("g_i", "g_j", "sbin"),
    "L2_identity_sep":   ("aa_i", "aa_j", "sbin"),
    "L2_identity_only":  ("aa_i", "aa_j"),
    "sep_only":          ("sbin",),
    "gray_h_sep":        ("gray_h", "sbin"),
}


def _cell_keys(feat: Dict[str, np.ndarray], names: Sequence[str]) -> np.ndarray:
    """Pack several small integer features into one flat cell id."""
    key = np.zeros(len(feat["sep"]), dtype=np.int64)
    for nm in names:
        v = feat[nm].astype(np.int64)
        key = key * (int(v.max(initial=0)) + 2) + (v + 1)
    return key


def stage_ceiling(limit: Optional[int], with_mi: bool) -> Dict:
    targets = ks.list_targets()
    if limit:
        targets = targets[:limit]
    print(f"[ceiling] {len(targets)} targets")

    per_protein: List[Dict] = []
    cache: List[Tuple[str, int, Dict[str, np.ndarray]]] = []
    for n, t in enumerate(targets):
        try:
            # dense is required: consensus (hence the pair universe) is derived
            # from the alignment, and is all -1 when the MSA is not loaded.
            pd = ks.load_protein(t, load_dense=True)
        except Exception as exc:                     # noqa: BLE001
            print(f"[ceiling]   {t}: SKIP ({exc})")
            continue
        pairs = ks.pair_universe(pd)
        if pairs.shape[0] < 50:
            continue
        feat = ks.pair_features(pd, pairs)
        if with_mi and pd.dense is not None:
            mi = ks.mi_scores(pd, pairs)
            # Quartile-binned MI, so it enters as a small categorical feature.
            q = np.quantile(mi, [0.25, 0.5, 0.75])
            feat["mi_q"] = np.digitize(mi, q).astype(np.int8)
        cache.append((t, pd.n_cols, feat))
        if (n + 1) % 25 == 0:
            print(f"[ceiling]   loaded {n + 1}/{len(targets)}")

    sets = dict(FEATURE_SETS)
    if with_mi:
        sets["L2_identity_sep_MIq"] = ("aa_i", "aa_j", "sbin", "mi_q")
        sets["MIq_only"] = ("mi_q",)

    results: Dict[str, Dict] = {}
    for set_name, names in sets.items():
        if any(nm not in cache[0][2] for nm in names):
            continue
        # Pooled oracle rate per cell (uses the labels -- deliberately).
        pos: Dict[int, int] = defaultdict(int)
        tot: Dict[int, int] = defaultdict(int)
        keys_per_protein = []
        for t, L, feat in cache:
            k = _cell_keys(feat, names)
            keys_per_protein.append((t, L, k, feat["is_contact"]))
            for kk, lab in zip(k.tolist(), feat["is_contact"].tolist()):
                tot[kk] += 1
                if lab:
                    pos[kk] += 1
        rate = {c: pos[c] / tot[c] for c in tot}

        opt, exp = [], []
        for t, L, k, lab in keys_per_protein:
            kk = max(1, L // 5)
            # A function of the cell assigns ONE score per cell, so the order
            # WITHIN a cell is arbitrary.  Bounding correctly therefore means
            # bounding over tie-breaks, not picking one:
            #   optimistic = contacts first inside the boundary cell  -> a true
            #                upper bound for every function of these features
            #   expected   = random tie-break -> what a real circuit averages
            order = sorted(rate, key=lambda c: -rate[c])
            here_tot: Dict[int, int] = defaultdict(int)
            here_pos: Dict[int, int] = defaultdict(int)
            for c, y in zip(k.tolist(), lab.tolist()):
                here_tot[c] += 1
                if y:
                    here_pos[c] += 1
            rem, hits_o, hits_e = kk, 0.0, 0.0
            for c in order:
                n_c = here_tot.get(c, 0)
                if n_c == 0:
                    continue
                p_c = here_pos.get(c, 0)
                take = min(rem, n_c)
                hits_o += min(take, p_c)                  # contacts first
                hits_e += take * (p_c / n_c)              # random order
                rem -= take
                if rem == 0:
                    break
            opt.append(hits_o / kk)
            exp.append(hits_e / kk)
        results[set_name] = {
            "features": list(names),
            "n_cells_used": len(tot),
            "ceiling_optimistic_mean": float(np.mean(opt)),
            "ceiling_expected_mean": float(np.mean(exp)),
            "ceiling_expected_median": float(np.median(exp)),
            "ceiling_expected_sd": float(np.std(exp, ddof=1)),
        }
        print(f"[ceiling]   {set_name:22s} cells={len(tot):6d} "
              f"upper-bound={np.mean(opt):.4f}  expected={np.mean(exp):.4f}")

    return {
        "n_targets": len(cache),
        "note": ("Each cell is scored by its own pooled empirical contact rate -- "
                 "the feature class handed the test labels. 'optimistic' resolves "
                 "within-cell ties in favour of contacts and is therefore a true "
                 "upper bound for ANY function of these features; 'expected' uses "
                 "random tie-breaking and is what a real circuit averages."),
        "literature_anchor": {"PSICOV_published_prec_at_L5": 0.723,
                              "PATH_B_circuit_L2_holdout": 0.179,
                              "PATH_B_plain_MI": 0.097},
        "feature_sets": results,
    }


# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True,
                    choices=("seq", "coevo", "contact", "chem", "ceiling"))
    ap.add_argument("--with-mi", action="store_true",
                    help="ceiling stage: also bound feature sets that use MI")
    ap.add_argument("--n-seqs", type=int, default=1299)
    ap.add_argument("--n-controls", type=int, default=50)
    ap.add_argument("--max-pairs", type=int, default=10)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--shuffles", type=int, default=5)
    ap.add_argument("--t-mult", type=float, default=2.0)
    ap.add_argument("--n-min", type=int, default=30)
    args = ap.parse_args()

    t0 = time.time()
    if args.stage == "seq":
        res = stage_seq(args.n_seqs, args.n_controls)
        _dump("stage_seq.json", res)
    elif args.stage == "coevo":
        res = stage_coevo(args.max_pairs)
        _dump("stage_coevo.json", res)
    elif args.stage == "contact":
        res = stage_contact(args.limit, args.shuffles)
        _dump("stage_contact.json", res)
    elif args.stage == "ceiling":
        res = stage_ceiling(args.limit, args.with_mi)
        _dump("stage_ceiling.json", res)
    else:
        res = stage_chem(args.limit, args.t_mult, args.n_min)
        _dump("stage_chem.json", res)
    print(f"[{args.stage}] done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
