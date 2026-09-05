#!/usr/bin/env python3
"""
test_sequence_circuits.py — correctness gates for the sequence-circuit engine.

Mirrors the style of test_kmap_structure.py: every gate is an assertion with a
printed verdict, exit code 1 on any failure.  The gates that matter most:

  t04-t06  our sparse QM produces the SAME prime-implicant set as the repo's
           reference dense implementation (`kmap_sbm.analysis.prime_implicants`)
  t07-t09  the positional circuit is LOSSLESS — the sequence is recovered
           exactly by evaluating the minimised circuit
  t12-t13  the encoding conventions (binary for position, Gray for residue)
           behave as documented
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import kmap_structure as ks          # noqa: E402
import sequence_circuits as sc       # noqa: E402

FAILS: list[str] = []


def check(tag: str, cond: bool, detail: str = "") -> None:
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {tag}  {detail}")
    if not cond:
        FAILS.append(tag)


# --- §1 implicant algebra ---------------------------------------------------
pi = sc.Implicant(val=0b0100, mask=0b0011, n_vars=4)
check("t01 implicant covers", all(pi.covers(0b0100 | k) for k in range(4)),
      "cube 01** covers its 4 minterms")
check("t02 implicant excludes", not pi.covers(0b1100), "01** excludes 11**")
check("t03 literals", pi.literals(["w", "x", "y", "z"]) == ["~w", "x"],
      f"{pi.literals(['w','x','y','z'])}")

# --- §2 QM agreement with the reference implementation ----------------------
rng = np.random.default_rng(20260905)

# t04: the worked toy from ContactCircuits.lean (ccTableToy / ccPiToy).
toy = sc.prime_implicants([0, 1], [], 2)
check("t04 lean toy", {(p.val, p.mask) for p in toy} == {(0, 1)},
      "On={0,1} -> single PI 0* (matches cc_cover_complete)")

# t05: with NO don't-cares the two implementations must agree exactly.
for tag, n_vars, n_on in (("t05a", 8, 30), ("t05b", 10, 60)):
    cells = 1 << n_vars
    perm = rng.permutation(cells)
    on = sorted(int(x) for x in perm[:n_on])
    rep = sc.qm_agreement_check(on, [], n_vars)
    check(f"{tag} QM vs reference (no DC) n={n_vars}", rep["identical"],
          f"ours={rep['n_ours']} ref={rep['n_ref']}")

# t06: SELF-AUDIT — the real correctness proof, independent of any reference.
# Every returned cube must (a) contain only ON/DC minterms (it is an implicant)
# and (b) fail to grow along any fixed bit (it is PRIME).
def _audit(on: list[int], dc: list[int], n_vars: int) -> tuple[int, int, int]:
    allowed = set(on) | set(dc)
    pis = sc.prime_implicants(on, dc, n_vars)
    bad_impl = bad_prime = 0
    for p in pis:
        free = [b for b in range(n_vars) if (p.mask >> b) & 1]
        fixed = [b for b in range(n_vars) if not (p.mask >> b) & 1]

        def _cube(base: int) -> bool:
            for k in range(1 << len(free)):
                m = base
                for t, b in enumerate(free):
                    if (k >> t) & 1:
                        m |= 1 << b
                if m not in allowed:
                    return False
            return True

        if not _cube(p.val):
            bad_impl += 1
        if any(_cube(p.val ^ (1 << b)) for b in fixed):
            bad_prime += 1
    return len(pis), bad_impl, bad_prime


for tag, n_vars, n_on, n_dc in (("t06a", 8, 30, 20), ("t06b", 10, 60, 400),
                                ("t06c", 12, 80, 600)):
    cells = 1 << n_vars
    perm = rng.permutation(cells)
    on = sorted(int(x) for x in perm[:n_on])
    dc = sorted(int(x) for x in perm[n_on:n_on + n_dc])
    n_pi, bad_i, bad_p = _audit(on, dc, n_vars)
    check(f"{tag} self-audit n={n_vars}", bad_i == 0 and bad_p == 0,
          f"pis={n_pi} non-implicants={bad_i} non-prime={bad_p}")

# t06d: DOCUMENTED DIVERGENCE — characterised, not just observed.
#
# The reference prunes don't-cares to those at Hamming distance 1 from a real
# on-minterm (a performance choice stated in its own docstring).  With DCs
# clustered — which is exactly the regime of the contact campaign, where DC =
# "cell has n < n_min" — cubes that would grow through a chain of don't-cares
# never form.  The consequence is precise and worth pinning down:
#
#   * the reference UNDER-MERGES: it misses true primes, and
#   * every cube it reports that we do not is STRICTLY INSIDE one of ours.
#
# So its rules stay SOUND (each cube still touches no off-cell) but its
# prime-implicant counts, and hence its "essential" designations, are
# approximations rather than the exact QM quantities the docstring claims.
cells = 1 << 10
perm = rng.permutation(cells)
on = sorted(int(x) for x in perm[:60])
dc = sorted(int(x) for x in perm[60:460])
ours_pis = sc.prime_implicants(on, dc, 10)
rep = sc.qm_agreement_check(on, dc, 10)


def _inside(val: int, mask: int, b: sc.Implicant) -> bool:
    """Cube (val, mask) lies inside cube b."""
    return (mask & ~b.mask) == 0 and ((val ^ b.val) & ~b.mask) == 0


# `only_ref` is truncated to 8 for display; recompute the full difference.
ref_only_all = {(v, m) for v, m in
                (sc._ref_pi_keys(on, dc, 10) - {(p.val, p.mask) for p in ours_pis})}
subsumed = sum(1 for (v, m) in ref_only_all
               if any(_inside(v, m, p) for p in ours_pis))
check("t06d reference under-merges (every extra cube is subsumed)",
      rep["n_ref"] < rep["n_ours"] and subsumed == len(ref_only_all),
      f"exact={rep['n_ours']} reference={rep['n_ref']} "
      f"ref-only={len(ref_only_all)} all-subsumed={subsumed == len(ref_only_all)}")

# --- §3 losslessness of the positional circuit ------------------------------
def _roundtrip(tag: str, seq: str) -> None:
    codes = [ks.AA_TO_CODE[c] for c in seq]
    circ = sc.seq_circuit_positional(codes, name=tag)
    rec = sc.reconstruct_from_circuit(circ, len(codes))
    ok = rec[:len(codes)] == codes
    check(f"{tag} lossless round-trip", ok and circ.complete and circ.sound,
          f"L={len(codes)} PIs={len(circ.primes)} cover={len(circ.cover)} "
          f"complete={circ.complete} sound={circ.sound}")


_roundtrip("t07", "ACDEFGHIKLMNPQRSTVWY")
_roundtrip("t08", "AAAAAAAAAAAAAAAA")              # maximally compressible
_roundtrip("t09", "MFVFLVLLPLVSSQCVNLITRTQSYTNSFTRGVYYPDKVFRSS")  # Spike N-term

# t10: a constant sequence must compress to far fewer cubes than a random one.
const_codes = [0] * 64
rand_codes = [int(x) for x in rng.integers(0, 20, size=64)]
c_const = sc.seq_circuit_positional(const_codes, name="const")
c_rand = sc.seq_circuit_positional(rand_codes, name="rand")
check("t10 compression tracks repetition",
      len(c_const.cover) < len(c_rand.cover),
      f"constant cover={len(c_const.cover)} random cover={len(c_rand.cover)}")

# t11: gaps are excluded from the on-set, not encoded as a residue.
gapped = [ks.AA_TO_CODE["A"], ks.GAP_STATE, ks.AA_TO_CODE["C"]]
c_gap = sc.seq_circuit_positional(gapped, name="gap")
check("t11 gaps excluded", c_gap.n_on == 2 and c_gap.meta["n_gaps"] == 1,
      f"n_on={c_gap.n_on} n_gaps={c_gap.meta['n_gaps']}")

# --- §4 encoding conventions ------------------------------------------------
# t12: a plain-binary position subcube with low bits free IS a dyadic interval.
seg = sc.implicant_to_segments(
    sc.Implicant(val=(0b0100 << 4) | 0b1000, mask=(0b0011 << 4) | 0b0011,
                 n_vars=8), pos_bits=4, length=16)
check("t12 binary position subcube = interval",
      seg["i"]["contiguous"] and seg["i"]["start"] == 4 and seg["i"]["end"] == 8
      and seg["j"]["start"] == 8 and seg["j"]["end"] == 12,
      f"i={seg['i']['start']}-{seg['i']['end']} j={seg['j']['start']}-{seg['j']['end']}")

# t13: Gray residue coding puts He-adjacent residues one bit apart.
adj = [sc.ks.gray_hamming(c, c + 1) for c in range(19)]
check("t13 Gray residue adjacency", all(d == 1 for d in adj),
      "consecutive He codes differ in exactly 1 Gray bit")

# t14: decoding a free-bit residue slot yields a multi-residue tolerance class.
dec = sc.decode_residue_pair_implicant(
    sc.Implicant(val=0, mask=0b00001_00000, n_vars=10))
check("t14 free bit = tolerance class", len(dec["slot_i"]["residues"]) == 2,
      f"slot_i admits {dec['slot_i']['residues']}")

# --- §5 co-evolution polarities --------------------------------------------
# Build a tiny synthetic MSA: column 0 in {A,I}, column 1 in {L,V}, with the
# combination (A,V) never observed -> the flipped circuit must fire there.
rows = [("A", "L"), ("A", "L"), ("I", "V"), ("I", "L")]
dense = np.array([[ks.AA_TO_CODE[a], ks.AA_TO_CODE[b]] for a, b in rows],
                 dtype=np.int8)
c_std = sc.coevo_circuit(dense, 0, 1, "standard")
c_flip = sc.coevo_circuit(dense, 0, 1, "flipped")
check("t15 standard polarity", c_std.n_on == 3,
      f"observed combos on={c_std.n_on} (AL, IV, IL)")
check("t16 flipped polarity", c_flip.n_on == 1,
      f"forbidden combos on={c_flip.n_on} (AV only)")
flip_cell = sc.decode_residue_pair_implicant(c_flip.cover[0])
check("t17 flipped decodes to A x V",
      "A" in flip_cell["slot_i"]["residues"] and
      "V" in flip_cell["slot_j"]["residues"],
      f"{flip_cell['slot_i']['residues']} x {flip_cell['slot_j']['residues']}")

# t18: soundness/completeness hold for every circuit built above.
allc = [c_const, c_rand, c_gap, c_std, c_flip]
check("t18 all circuits sound+complete",
      all(c.complete and c.sound for c in allc),
      f"{len(allc)} circuits verified")

# --- §6 contact-map circuit -------------------------------------------------
# A clean two-segment contact block: [0,4) x [8,12), well outside the band.
block = [(i, j) for i in range(0, 4) for j in range(8, 12)]
# band_as_dc=False: with the near-diagonal band left as don't-care the
# minimiser legitimately absorbs it into one huge cube, which is correct but
# tells us nothing about block recovery.  Pin the band OFF to test the claim.
c_ct = sc.contact_circuit_positional(block, length=16, name="t19",
                                     band_as_dc=False)
segs = [sc.implicant_to_segments(p, 4, 16) for p in c_ct.cover]
found = any(s["i"]["contiguous"] and s["j"]["contiguous"] and
            {(s["i"]["start"], s["i"]["end"]),
             (s["j"]["start"], s["j"]["end"])} == {(0, 4), (8, 12)}
            for s in segs)
check("t19 contact block -> one segment-segment cube", found,
      f"cover={len(c_ct.cover)} cubes; block recovered={found}")

print()
if FAILS:
    print(f"FAILED {len(FAILS)}: {', '.join(FAILS)}")
    sys.exit(1)
print("all gates passed")
