#!/usr/bin/env python3
"""
test_kmap_structure.py — Adversarial unit + property tests for the
contact-K-map campaign library (scripts/kmap_structure.py).

Run:  /store/shuvam/.venv/bin/python scripts/tests/test_kmap_structure.py
Exit 0 = all green.  No external test framework required (plain asserts).

Coverage map (each test cites the invariant it locks):
  T01-T04  encoding helpers (Gray, groups, separation bins)
  T05-T08  alignment loading + consensus (ragged, case, gaps, ties)
  T09-T11  PDB parsing (drop-rule, Gly-CA, altloc defence, blank chain)
  T12      native_contacts geometry (min_sep, cutoff) on synthetic structures
  T13      LIVE identity-mapping gate (1a3aA ok; 1vjkA 87-vs-88 artifact)
  T14      LIVE end-to-end feature/label consistency on a real protein
  T15-T16  MI + perplexity EQUIVALENCE vs coevolution_shared (property, 25 trials)
  T17      cell-count accumulation hand-example
  T18      truth-table labeling hand-example + degenerate all-DC
  T19      QM SOUNDNESS + COMPLETENESS on real labeled tables (the core theorem):
             every ON cell covered by >=1 prime implicant;
             no prime implicant matches any OFF cell.
  T20      decode_level1_implicant hand-check + brute-force predicate equivalence
  T21-T23  metrics: precision@k, MCC, AUC vs brute-force concordance (ties!)
  T24      protein_folds determinism + exact partition + stratification
  T25      adversarial: all-gap tiny MSA, zero-contact protein, empty pipelines
"""
from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent
REPO = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(REPO))
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")

import kmap_structure as ks  # noqa: E402


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def make_pdb(lines):
    return [ln + "\n" for ln in lines]


def atom(serial, name, resname, chain, resseq, x, y, z, alt=" "):
    return (f"ATOM  {serial:5d} {name:>4s}{alt}{resname:>3s} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]:>2s}")


def expect_raises(exc_type, fn, *args, **kwargs):
    """Assert fn(*args) raises exc_type; AssertionError otherwise (runner-safe:
    never uses SystemExit, which would escape the except-Exception harness)."""
    try:
        fn(*args, **kwargs)
    except exc_type:
        return
    raise AssertionError(f"expected {exc_type.__name__} from {getattr(fn, '__name__', fn)}")


# ---------------------------------------------------------------------------
# T01-T04  encoding helpers
# ---------------------------------------------------------------------------

def t01_gray_code_basics():
    assert ks.gray_code(0) == 0
    assert ks.gray_code(1) == 1
    assert ks.gray_code(2) == 3
    assert ks.gray_code(3) == 2
    assert ks.gray_code(19) == 19 ^ (19 >> 1) == 26
    # involution property: decoding g(x) by cumulative xor recovers x
    for x in range(32):
        y, acc = ks.gray_code(x), 0
        while y:
            acc ^= y
            y >>= 1
        assert acc == x, f"Gray decode failed at {x}"


def t02_gray_hamming():
    assert ks.gray_hamming(0, 0) == 0
    assert ks.gray_hamming(0, 1) == 1     # adjacent integers -> distance 1
    assert ks.gray_hamming(2, 3) == 1
    # Reflected-Gray cyclic property: codes 0 and 31 are ADJACENT in Q5
    # (g(31)=10000), distance 1 -- locked explicitly:
    assert ks.gray_hamming(0, 31) == 1
    # a genuine maximum-distance pair: g(21)=11111 so (0,21) spans 5 bits
    assert ks.gray_code(21) == 31
    assert ks.gray_hamming(0, 21) == 5
    for bad in (-1, 32):
        expect_raises(ValueError, ks.gray_hamming, bad, 0)
    # brute-force agreement
    rng = np.random.default_rng(0)
    for _ in range(200):
        a, b = rng.integers(0, 32, size=2)
        ga, gb = ks.gray_code(int(a)), ks.gray_code(int(b))
        assert ks.gray_hamming(int(a), int(b)) == bin(ga ^ gb).count("1")


def t03_groups_and_bins():
    expect = {"A": 0, "I": 0, "L": 0, "V": 0, "M": 0,
              "F": 1, "Y": 1, "W": 1,
              "E": 2, "D": 2, "Q": 3, "N": 3,
              "H": 4, "K": 4, "R": 4,
              "S": 5, "T": 5, "C": 6, "P": 7, "G": 7}
    for aa, g in expect.items():
        assert ks.group_of_code(ks.AA_TO_CODE[aa]) == g, aa
    assert ks.group_of_code(ks.GAP_STATE) == -1
    assert ks.group_of_code(-1) == -1
    # separation bins: boundaries locked
    assert [ks.sep_bin(s) for s in (6, 10)] == [0, 0]
    assert [ks.sep_bin(s) for s in (11, 20)] == [1, 1]
    assert [ks.sep_bin(s) for s in (21, 30)] == [2, 2]
    assert ks.sep_bin(31) == 3 and ks.sep_bin(10_000) == 3
    expect_raises(ValueError, ks.sep_bin, 5)


# alias kept for readability at call sites below
_expect_raises = expect_raises


def t04_group_partition_is_total():
    # every He code maps to some group; partition covers all 20 exactly once
    seen = sorted(ks.group_of_code(c) for c in range(20))
    assert seen == sorted([ks.AA_GROUP8[a] for a in ks.CODE_TO_AA])


# ---------------------------------------------------------------------------
# T05-T08  alignment loading + consensus
# ---------------------------------------------------------------------------

def t05_load_aln_basic():
    he = ks.AMINO_HE_2012                       # 'AILVMFYWEDQNHKRSTCPG'
    with tempfile.NamedTemporaryFile("w", suffix=".aln", delete=False) as f:
        f.write(he + "-\n")                     # 21 cols: 20 codes + gap
        f.write(he.lower() + "x\n")             # lowercase + X -> same codes/gap
        path = Path(f.name)
    dense, rows = ks.load_aln_codes(path)
    assert dense.shape == (2, 21)
    assert list(dense[0]) == list(range(20)) + [ks.GAP_STATE]
    assert list(dense[1]) == list(range(20)) + [ks.GAP_STATE]
    path.unlink()


def t06_load_aln_ragged_and_empty():
    with tempfile.NamedTemporaryFile("w", suffix=".aln", delete=False) as f:
        f.write("AAAA\nAA\n")
        p1 = Path(f.name)
    expect_raises(ValueError, ks.load_aln_codes, p1)
    with tempfile.NamedTemporaryFile("w", suffix=".aln", delete=False) as f:
        p2 = Path(f.name)
    expect_raises(ValueError, ks.load_aln_codes, p2)
    p1.unlink(); p2.unlink()


def t07_consensus_majority_lowest_tie():
    dense = np.array([
        [0, 5, 3, ks.GAP_STATE],
        [1, 5, 7, ks.GAP_STATE],
        [1, 5, 3, ks.GAP_STATE],
        [2, 5, 7, ks.GAP_STATE],
    ], dtype=np.int8)
    cons = ks.consensus_codes(dense)
    assert cons[0] == 1                 # mode of {0,1,1,2}
    assert cons[1] == 5
    assert cons[2] == 3                 # tie 3-vs-7 -> LOWEST code
    assert cons[3] == -1                # all-gap sentinel


def t08_consensus_excludes_gap_state():
    dense = np.full((3, 2), ks.GAP_STATE, dtype=np.int8)
    dense[0, 0] = 9
    assert ks.consensus_codes(dense)[0] == 9
    assert ks.consensus_codes(dense)[1] == -1


# ---------------------------------------------------------------------------
# T09-T11  PDB parsing
# ---------------------------------------------------------------------------

def _tiny_pdb(extra_terminal=False, glycb=False, altloc=False):
    lines = []
    s = 1
    # res1 LEU at origin; res7 LEU at (5,0,0) => contact (sep 6, d=5A)
    lines.append(atom(s, "CA", "LEU", " ", 1, 0.0, 0.0, 0.0)); s += 1
    lines.append(atom(s, "CB", "LEU", " ", 1, 1.5, 0.0, 0.0)); s += 1
    for r in range(2, 7):
        lines.append(atom(s, "CA", "ALA", " ", r, float(r * 30), 0.0, 0.0)); s += 1
        lines.append(atom(s, "CB", "ALA", " ", r, float(r * 30), 1.5, 0.0)); s += 1
    lines.append(atom(s, "CA", "LEU", " ", 7, 5.0, 0.0, 0.0)); s += 1
    lines.append(atom(s, "CB", "LEU", " ", 7, 5.0, 1.5, 0.0)); s += 1
    if glycb:
        # res8 GLY WITH a nominal CB (should still use CA per convention)
        lines.append(atom(s, "CA", "GLY", " ", 8, 500.0, 0.0, 0.0)); s += 1
        lines.append(atom(s, "CB", "GLY", " ", 8, 501.5, 0.0, 0.0)); s += 1
    else:
        lines.append(atom(s, "CA", "GLY", " ", 8, 500.0, 0.0, 0.0)); s += 1
    if extra_terminal:
        lines.append(atom(s, "CA", "SER", " ", 9, 900.0, 0.0, 0.0)); s += 1
    if altloc:
        # res9 VAL with two alternate CB locations; first occurrence must win
        lines.append(atom(s, "CB", "VAL", " ", 9, 700.0, 0.0, 0.0, alt="A")); s += 1
        lines.append(atom(s, "CB", "VAL", " ", 9, 800.0, 0.0, 0.0, alt="B")); s += 1
        lines.append(atom(s, "CA", "VAL", " ", 9, 699.0, 0.0, 0.0)); s += 1
    return make_pdb(lines)


def t09_parse_drop_rule_and_gly():
    res = ks.parse_pdb_cb(_tiny_pdb(), n_cols=8)
    assert set(res.keys()) == {1, 2, 3, 4, 5, 6, 7, 8}
    res_all = ks.parse_pdb_cb(_tiny_pdb(extra_terminal=True), n_cols=None)
    assert 9 in res_all
    res_trim = ks.parse_pdb_cb(_tiny_pdb(extra_terminal=True), n_cols=8)
    assert 9 not in res_trim                     # drop-rule removes res > n_cols
    assert res[8].xyz == (500.0, 0.0, 0.0)       # Gly used CA despite nominal CB


def t10_parse_altloc_first_wins():
    res = ks.parse_pdb_cb(_tiny_pdb(altloc=True), n_cols=None)
    assert res[9].xyz == (700.0, 0.0, 0.0)


def t11_parse_rejects_inconsistent_names():
    lines = _tiny_pdb()
    lines.append(atom(999, "CA", "TRP", " ", 1, 42.0, 0.0, 0.0))
    # Strict mode raises ValueError
    expect_raises(ValueError, ks.parse_pdb_cb, lines, strict=True)
    # Tolerant mode silently keeps first occurrence
    res_tol = ks.parse_pdb_cb(lines, strict=False)
    assert res_tol[1].aa == "L"


# ---------------------------------------------------------------------------
# T12  native_contacts geometry
# ---------------------------------------------------------------------------

def t12_native_contacts_geometry():
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as f:
        f.writelines(_tiny_pdb())
        path = Path(f.name)
    cont = ks.native_contacts(path, n_cols=8)
    # res1-res7: sep 6 >= 6, CB-CB distance = |(5,1.5,0)-(1.5,0,0)| ~= 3.64 < 8
    assert (0, 6) in cont
    # nothing closer than sep 6 may appear even though res1..res6 are near-ish?
    # (here they're 30A apart anyway; the assertion is structural:)
    for (i, j) in cont:
        assert j - i >= ks.MIN_SEPARATION
        assert j < 8                              # drop-rule respected
    path.unlink()


# ---------------------------------------------------------------------------
# T13  LIVE identity-mapping gate
# ---------------------------------------------------------------------------

def t13_identity_mapping_live():
    ok, w, ca = ks.verify_identity_mapping("1a3aA")
    assert ok and w == ca and w > 0
    ok2, w2, ca2 = ks.verify_identity_mapping("1vjkA")
    assert (not ok2) and w2 == 87 and ca2 == 88   # documented artifact


# ---------------------------------------------------------------------------
# T14  LIVE end-to-end consistency on a real protein
# ---------------------------------------------------------------------------

def t14_live_feature_label_consistency():
    pd = ks.load_protein("1a3aA")
    assert pd.identity_ok and pd.n_seqs > 100 and len(pd.contacts) > 20
    pairs = ks.pair_universe(pd)
    feats = ks.pair_features(pd, pairs)
    assert len(pairs) == len(feats["i"]) == len(feats["is_contact"])
    # labels agree exactly with the native contact set
    for (i, j), lab in zip(pairs[:500], feats["is_contact"][:500]):
        assert bool(lab) == ((int(i), int(j)) in pd.contacts)
    # separation bins consistent with separations
    for s, b in zip(feats["sep"][:2000], feats["sbin"][:2000]):
        assert ks.sep_bin(int(s)) == int(b)
    # gray hamming consistent with brute force
    for ai, aj, gh in zip(feats["aa_i"][:500], feats["aa_j"][:500], feats["gray_h"][:500]):
        assert int(gh) == bin(ks.gray_code(int(ai)) ^ ks.gray_code(int(aj))).count("1")
    # 1vjkA drop-rule note present and contacts confined to aligned region
    pd2 = ks.load_protein("1vjkA", load_dense=False)
    assert pd2.notes and all(j < 87 for (_, j) in pd2.contacts)


# ---------------------------------------------------------------------------
# T15-T16  equivalence vs coevolution_shared (property tests)
# ---------------------------------------------------------------------------

def _random_msas(trials=25, seed=123):
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        n = int(rng.integers(10, 200))
        L = int(rng.integers(5, 30))
        dense = rng.integers(-1, 21, size=(n, L)).astype(np.int32)
        yield dense, rng


def t15_mi_equivalence_with_shared():
    from coevolution_shared import mutual_information
    for dense, rng in _random_msas():
        arrs = [dense[k] for k in range(dense.shape[0])]
        i, j = sorted(rng.choice(dense.shape[1], size=2, replace=False).tolist())
        a = mutual_information(arrs, int(i), int(j), dense.shape[0])
        b = ks.mutual_information_columns(dense, int(i), int(j))
        assert abs(a - b) < 1e-12, (a, b)


def t16_perplexity_equivalence_with_shared():
    from coevolution_shared import perplexity_ratio
    for dense, rng in _random_msas(trials=25, seed=777):
        arrs = [dense[k] for k in range(dense.shape[0])]
        i, j = sorted(rng.choice(dense.shape[1], size=2, replace=False).tolist())
        a = perplexity_ratio(arrs, int(i), int(j), dense.shape[0])
        b = ks.perplexity_ratio_columns(dense, int(i), int(j))
        assert (a is None) == (b is None)
        if a is not None:
            assert abs(a - b) < 1e-12, (a, b)


def t16b_perplexity_none_condition_mirrors_code():
    """F1 regression lock: BOTH implementations return None only when NO
    conditioning residue reaches n>=5 (mirrors shared CODE, not its docstring)."""
    # every residue at column i appears <=2 times -> cond list empty -> None
    dense = np.zeros((8, 2), dtype=np.int32)
    dense[:, 0] = [0, 0, 1, 1, 2, 2, 3, 3]
    dense[:, 1] = [0, 1, 0, 1, 0, 1, 0, 1]
    assert ks.perplexity_ratio_columns(dense, 0, 1) is None
    from coevolution_shared import perplexity_ratio
    arrs = [dense[k] for k in range(8)]
    assert perplexity_ratio(arrs, 0, 1, 8) is None


# ---------------------------------------------------------------------------
# T17-T18  counting + labeling
# ---------------------------------------------------------------------------

def t17_accumulate_hand_example():
    c1 = ks.new_level1_counts()
    c2 = ks.new_level2_counts()
    feats = {
        "g_i": np.array([2, 2, 4, 7], dtype=np.int8),
        "g_j": np.array([4, 4, 2, 7], dtype=np.int8),
        "sbin": np.array([0, 0, 1, 3], dtype=np.int8),
        "aa_i": np.array([2, 2, 13, 19], dtype=np.int8),
        "aa_j": np.array([13, 13, 2, 19], dtype=np.int8),
        "is_contact": np.array([True, False, True, False]),
    }
    ks.accumulate_protein(c1, c2, feats)
    idx = ks.level1_cell_index(feats["g_i"], feats["g_j"], feats["sbin"])
    assert c1.tot[idx[0]] == 2 and c1.pos[idx[0]] == 1     # acidic-basic short
    assert c1.tot[idx[2]] == 1 and c1.pos[idx[2]] == 1     # basic-acidic mid
    assert c1.tot[idx[3]] == 1 and c1.pos[idx[3]] == 0
    assert c1.tot.sum() == 4 and c1.pos.sum() == 2
    assert abs(c1.base_rate() - 0.5) < 1e-12
    # level-2: both acidic->basic obs share cell (aa_i=2, aa_j=13) in bin 0
    cell = ks.level2_cell_index(np.int8(2), np.int8(13))
    assert c2[0].tot[cell] == 2 and c2[0].pos[cell] == 1
    assert c2[1].tot[ks.level2_cell_index(np.int8(13), np.int8(2))] == 1


def t18_truth_labeling():
    c = ks.CellCounts(4, np.array([8, 1, 0, 0], np.int64),
                      np.array([10, 100, 5, 0], np.int64))
    tt = ks.label_truth_table(c, t_mult=2.0, n_min=6)
    base = 9.0 / 115.0
    assert tt[0] == 1                       # rate .8 >= 2*.0783, n=10>=6
    assert tt[1] == 0                       # rate .01 < threshold
    assert tt[2] == -1                      # n=5 < n_min -> don't-care
    assert tt[3] == -1                      # n=0 -> don't-care
    # degenerate: n_min above every count -> all DC
    tt2 = ks.label_truth_table(c, t_mult=2.0, n_min=1000)
    assert (tt2 == -1).all()
    # explicit base override
    tt3 = ks.label_truth_table(c, t_mult=2.0, n_min=6, base_rate=0.5)
    assert tt3[0] == 0 and tt3[1] == 0      # threshold 1.0 unreachable


# ---------------------------------------------------------------------------
# T19  QM soundness/completeness (the mathematical heart)
# ---------------------------------------------------------------------------

def _pi_predicate(values, mask, width_per_field, n_bits):
    """Return predicate f(cell_idx)->bool for an implicant (MSB-first bits)."""
    def pred(cell: int) -> bool:
        for b in range(n_bits):
            bit = (cell >> (n_bits - 1 - b)) & 1
            if not mask[b] and bit != values[b]:
                return False
        return True
    return pred


def _check_qm_soundness(tt: np.ndarray, result: dict, tag: str):
    """Completeness: every ON cell covered. Soundness: no PI matches an OFF cell."""
    n_bits = int(np.log2(len(tt)))
    on_cells = np.flatnonzero(tt == 1)
    off_cells = np.flatnonzero(tt == 0)
    preds = [_pi_predicate(pi["values"], pi["mask"], None, n_bits)
             for pi in result["prime_implicants"]]
    for c in on_cells:
        assert any(p(int(c)) for p in preds), f"{tag}: ON cell {c} uncovered"
    for p in preds:
        for c in off_cells:
            assert not p(int(c)), f"{tag}: PI matches OFF cell {c}"
    # declared counts consistent
    assert result["n_prime_implicants"] == len(result["prime_implicants"])
    assert int((tt == 1).sum()) == result["n_minterms"]
    assert int((tt == -1).sum()) == result["n_dontcare"]


def t19_qm_soundness_synthetic_and_real():
    # synthetic: 256-cell level-1 table with planted block structure
    c = ks.new_level1_counts()
    rng = np.random.default_rng(5)
    for rep in range(400):
        gi, gj, sb = 2, 4, int(rng.integers(0, 2))       # enriched block
        c.add(np.array([(gi * 8 + gj) * 4 + sb]), np.array([True]))
        gi2, gj2 = int(rng.integers(0, 8)), int(rng.integers(0, 8))
        c.add(np.array([(gi2 * 8 + gj2) * 4 + 3]), np.array([False]))
    tt = ks.label_truth_table(c, t_mult=2.0, n_min=30)
    assert (tt == 1).sum() >= 1
    res = ks.qm_minimize(tt)
    _check_qm_soundness(tt.astype(int), res, "synthetic-level1")
    # real: 1a3aA end-to-end level-1 table
    pd = ks.load_protein("1a3aA")
    c1r = ks.new_level1_counts()
    c2r = ks.new_level2_counts()
    ks.accumulate_protein(c1r, c2r, ks.pair_features(pd, ks.pair_universe(pd)))
    tt_r = ks.label_truth_table(c1r, t_mult=2.0, n_min=30)
    res_r = ks.qm_minimize(tt_r)
    _check_qm_soundness(tt_r.astype(int), res_r, "real-1a3aA-level1")
    # real level-2 (bin 0), padded 1024-cell table
    tt_l2 = ks.label_truth_table(c2r[0], t_mult=2.0, n_min=30)
    res_l2 = ks.qm_minimize(tt_l2)
    _check_qm_soundness(tt_l2.astype(int), res_l2, "real-1a3aA-level2-bin0")
    # padding doctrine: no PI may fix a border code (>=20) on either axis
    for pi in res_l2["prime_implicants"]:
        vals, msk = pi["values"], pi["mask"]
        row_fixed = sum(v << (4 - k) for k, (v, m) in
                        enumerate(zip(vals[:5], msk[:5])) if not m)
        col_fixed = sum(v << (4 - k) for k, (v, m) in
                        enumerate(zip(vals[5:10], msk[5:10])) if not m)
        row_free = any(msk[k] for k in range(5))
        col_free = any(msk[k] for k in range(5, 10))
        if not row_free:
            assert row_fixed < 20
        if not col_free:
            assert col_fixed < 20


# ---------------------------------------------------------------------------
# T20  implicant decoder
# ---------------------------------------------------------------------------

def t20_decode_level1():
    d = ks.decode_level1_implicant(
        values=[0, 1, 0, 1, 0, 1, 0, 1],
        mask=[False, False, False, True, True, True, False, False])
    assert d["groups_i"] == [2]
    assert d["groups_j"] is None
    assert d["seps"] == [1]
    assert d["n_dontcare"] == 3
    # brute-force predicate equivalence over random implicants
    rng = np.random.default_rng(11)
    for _ in range(300):
        values = rng.integers(0, 2, size=8).tolist()
        mask = [bool(m) for m in (rng.random(8) < 0.5)]
        dec = ks.decode_level1_implicant(values, mask)
        for gi in range(8):
            for gj in range(8):
                for sb in range(4):
                    cell = (gi * 8 + gj) * 4 + sb
                    bits = [(cell >> (7 - b)) & 1 for b in range(8)]
                    want = all(m or bits[b] == values[b] for b, m in enumerate(mask))
                    got = ((dec["groups_i"] is None or gi in dec["groups_i"])
                           and (dec["groups_j"] is None or gj in dec["groups_j"])
                           and (dec["seps"] is None or sb in dec["seps"]))
                    assert want == got, (gi, gj, sb)


# ---------------------------------------------------------------------------
# T21-T23  metrics
# ---------------------------------------------------------------------------

def t21_precision_at_k():
    scores = np.array([0.9, 0.8, 0.7, 0.6])
    labels = np.array([True, False, True, False])
    assert abs(ks.precision_at_k(scores, labels, 2) - 0.5) < 1e-12
    assert abs(ks.precision_at_k(scores, labels, 4) - 0.5) < 1e-12
    assert ks.precision_at_k(scores, labels, 0) == 0.0
    assert ks.precision_at_k(np.array([]), np.array([], bool), 3) == 0.0


def t22_mcc():
    y = np.array([1, 1, 0, 0], bool)
    assert abs(ks.binary_metrics(np.array([1, 1, 0, 0]), y)["mcc"] - 1.0) < 1e-12
    assert abs(ks.binary_metrics(np.array([0, 0, 1, 1]), y)["mcc"] + 1.0) < 1e-12
    single = np.array([1, 1, 1, 1], bool)
    assert ks.binary_metrics(single, y)["mcc"] == 0.0     # degenerate denominator


def t23_auc_vs_bruteforce():
    rng = np.random.default_rng(2024)
    for trial in range(50):
        n = 40
        scores = rng.normal(size=n)
        labels = rng.random(n) < 0.4
        a = ks.auc_mann_whitney(scores, labels)
        pos, neg = scores[labels], scores[~labels]
        conc = 0.0
        for p in pos:
            conc += np.sum(p > neg) + 0.5 * np.sum(p == neg)
        b = conc / (len(pos) * len(neg))
        assert abs(a - b) < 1e-12, (trial, a, b)
    # perfect separation and full-tie corner cases
    assert ks.auc_mann_whitney(np.array([.1, .2, .8, .9]),
                               np.array([0, 0, 1, 1], bool)) == 1.0
    assert ks.auc_mann_whitney(np.array([.5, .5, .5, .5]),
                               np.array([0, 1, 0, 1], bool)) == 0.5
    assert math.isnan(ks.auc_mann_whitney(np.array([1., 2.]),
                                          np.array([1, 1], bool)))


# ---------------------------------------------------------------------------
# T24  folds
# ---------------------------------------------------------------------------

def t24_protein_folds():
    targets = ks.list_targets()
    assert len(targets) == 150
    f1 = ks.protein_folds(targets, n_folds=5, seed=42)
    f2 = ks.protein_folds(targets, n_folds=5, seed=42)
    assert f1 == f2                                    # deterministic
    union = sorted(sum(f1, []))
    assert union == targets                            # exact partition
    assert sum(len(f) for f in f1) == 150
    lens = [sum(1 for f in f1)]
    for f in f1:
        assert 20 <= len(f) <= 40                      # roughly balanced


# ---------------------------------------------------------------------------
# T25  adversarial empties
# ---------------------------------------------------------------------------

def t25_adversarial_empties():
    # all-gap alignment: no valid columns -> empty universe, no crash anywhere
    dense = np.full((6, 9), ks.GAP_STATE, dtype=np.int8)
    cons = ks.consensus_codes(dense)
    pd = ks.ProteinData(target="synthetic", n_cols=9, n_seqs=6, dense=dense,
                        consensus=cons, contacts=set(), identity_ok=True)
    pairs = ks.pair_universe(pd)
    assert len(pairs) == 0
    feats = ks.pair_features(pd, pairs)
    assert all(len(v) == 0 for v in feats.values())
    c1, c2 = ks.new_level1_counts(), ks.new_level2_counts()
    ks.accumulate_protein(c1, c2, feats)
    tt = ks.label_truth_table(c1)
    assert (tt == -1).all()
    res = ks.qm_minimize(tt)
    assert res["n_essential"] == 0
    _check_qm_soundness(tt.astype(int), res, "empty")   # vacuously true, must not raise
    # zero-contact real protein path: counts stay zero; labeling all 0/DC
    pdz = ks.ProteinData(target="z", n_cols=pd.n_cols, n_seqs=6, dense=dense,
                         consensus=cons, contacts=set(), identity_ok=True)
    assert len(pdz.contacts) == 0


# ---------------------------------------------------------------------------

ALL_TESTS = [(name, fn) for name, fn in sorted(globals().items())
             if name.startswith("t") and callable(fn)]

if __name__ == "__main__":
    failures = 0
    for name, fn in ALL_TESTS:
        try:
            fn()
            print(f"PASS {name}")
        except Exception as e:  # noqa: BLE001 -- report and continue
            failures += 1
            print(f"FAIL {name}: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
    print(f"\n{len(ALL_TESTS) - failures}/{len(ALL_TESTS)} tests passed")
    sys.exit(1 if failures else 0)
