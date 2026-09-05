#!/usr/bin/env python3
"""
sequence_circuits.py — Boolean-circuit representations of biological sequences.

Companion to `kmap_structure.py` (the contact campaign).  Where that module asks
"does a K-map circuit PREDICT 3D contacts?", this one asks the prior question:

    can a protein / RNA sequence itself be written as a minimal Boolean circuit,
    and does that circuit's algebraic structure carry 3D information?

Four constructions, in increasing distance from the raw sequence
-----------------------------------------------------------------
R-B  `seq_circuit_positional`  position-indexed indicator  f(i,a) = [s_i = a]
     LOSSLESS: s |-> f_s is injective, so the prime-implicant set of f_s is a
     canonical invariant OF THE SEQUENCE.  Round-trip verified at runtime.

R-A  `seq_circuit_kmer`        k-mer-set indicator  f(a_1..a_k) = [k-mer occurs]
     LOSSY (composition only).  This is the classical Petoukhov/FCGR-style
     K-map that the rest of the repo uses; included for comparison.

R-C  `coevo_circuit`           MSA column-pair circuits over the 1,299 Spike
     sequences, in three polarities:
       standard    ON = observed combination        (positive selection)
       flipped     ON = never-observed combination  (negative selection)
       perplexity  ON = perplexity-ratio enriched   (surprise, not MI)

R-D  `contact_circuit_positional`  the CONTACT MAP as a Boolean function of
     position bits, C(i,j).  See §"Encoding choice" — this is where the 3D
     content lives.

Encoding choice (this matters, and is the one real design departure)
--------------------------------------------------------------------
The repo uses Gray code everywhere.  That is correct for RESIDUE IDENTITY —
He-2012 order + Gray puts chemically similar residues one bit-flip apart, so a
don't-care bit in an implicant means "any residue in this chemical class"
(the algebraic form of mutational tolerance; cf. AminoAcidEncoding.lean).

It is WRONG for POSITION.  A subcube of plain-binary position indices with the
low k bits free is exactly the aligned interval [m*2^k, (m+1)*2^k) — a
contiguous stretch of the chain.  Under Gray coding the same subcube is a
reflected, non-contiguous set.  Since secondary-structure elements are
contiguous stretches, plain binary is what makes an implicant of the contact
map mean "segment A touches segment B".  So:

    positions  -> plain binary   (subcube = dyadic interval = chain segment)
    residues   -> He/Gray        (subcube = chemical class)

Both conventions are unit-tested in tests/test_sequence_circuits.py.

House rules followed from kmap_structure.py: explicit tables, no hidden
globals, every count reproducible, don't-cares carried as -1.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

# Reuse the campaign's verified encoding tables rather than redefining them.
import kmap_structure as ks

# ---------------------------------------------------------------------------
# §1  Implicants and exact Quine-McCluskey (sparse, don't-care aware)
# ---------------------------------------------------------------------------
#
# An implicant over `n` variables is a pair (val, mask) of n-bit integers:
#     mask bit = 1  ->  that variable is a don't-care (free) in the implicant
#     mask bit = 0  ->  that variable is fixed to the corresponding bit of val
# Canonical form keeps val's masked bits at 0 so (val, mask) is a unique key.
#
# The existing `kmap_sbm.analysis.prime_implicants` materialises a dense
# (n_cells x 2k+1) table; that is fine at 1024 cells but not at 2^22.  This
# implementation works from the sparse on-set / don't-care-set directly and is
# cross-validated against the dense one in `qm_agreement_check`.


@dataclass(frozen=True)
class Implicant:
    """A product term (cube). `val`/`mask` follow the convention above."""

    val: int
    mask: int
    n_vars: int

    def covers(self, minterm: int) -> bool:
        """True iff `minterm` lies inside this cube."""
        return (minterm & ~self.mask) == self.val

    @property
    def n_free(self) -> int:
        """Number of don't-care variables (cube dimension)."""
        return bin(self.mask).count("1")

    @property
    def size(self) -> int:
        """Number of minterms the cube covers = 2**n_free."""
        return 1 << self.n_free

    def literals(self, var_names: Sequence[str]) -> List[str]:
        """Human-readable literal list, MSB-first, skipping free variables."""
        out: List[str] = []
        for b in range(self.n_vars - 1, -1, -1):
            if (self.mask >> b) & 1:
                continue
            name = var_names[self.n_vars - 1 - b]
            out.append(name if (self.val >> b) & 1 else f"~{name}")
        return out


def prime_implicants(on_set: Iterable[int], dc_set: Iterable[int],
                     n_vars: int) -> List[Implicant]:
    """Exact Quine-McCluskey prime implicants of a partially-specified function.

    `on_set` minterms must be covered; `dc_set` minterms may be used to grow
    cubes but need not be covered (standard QM don't-care handling).  Runs on
    the sparse term set, so cost scales with |on|+|dc|, not with 2**n_vars.

    Returns the (unique) set of prime implicants.  Uniqueness is what makes
    this object a canonical invariant of the function — see module docstring.
    """
    terms: Set[Tuple[int, int]] = {(m, 0) for m in on_set}
    terms |= {(m, 0) for m in dc_set}
    if not terms:
        return []

    primes: Set[Tuple[int, int]] = set()
    current = terms
    while current:
        # Group by (mask, popcount(val)); only adjacent popcount groups merge.
        buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
        for val, mask in current:
            buckets[(mask, bin(val).count("1"))].append(val)

        merged_any: Set[Tuple[int, int]] = set()
        used: Set[Tuple[int, int]] = set()
        for (mask, pc), vals in buckets.items():
            partners = buckets.get((mask, pc + 1))
            if not partners:
                continue
            partner_set = set(partners)
            for v in vals:
                # A merge flips exactly one 0-bit of v that is not masked.
                free = ~mask & ((1 << n_vars) - 1)
                bits = free & ~v
                b = bits
                while b:
                    low = b & -b
                    w = v | low
                    if w in partner_set:
                        merged_any.add((v & ~low, mask | low))
                        used.add((v, mask))
                        used.add((w, mask))
                    b ^= low
        primes |= (current - used)
        current = merged_any

    return [Implicant(v, m, n_vars) for (v, m) in sorted(primes)]


def _essential_and_greedy_cover(pis: List[Implicant],
                                on_set: Sequence[int]) -> Tuple[List[Implicant],
                                                                List[Implicant]]:
    """Essential prime implicants + a greedy completion of the cover.

    Returns (cover, essentials).  Greedy (max newly-covered per pick) rather
    than Petrick: exact minimum-cover is NP-hard and the on-sets here reach
    10^3-10^4 minterms.  The cover is always SOUND and COMPLETE (asserted by
    `verify_cover`); only its cardinality is not certified minimum.
    """
    if not on_set:
        return [], []
    # Coverage as Python-int bitsets over on-minterm INDICES.  Set algebra on
    # big ints is C-speed, which is what makes 10^3 minterms x 10^4 primes
    # tractable; the naive per-minterm loop is O(|pis| * |on|) per greedy pick.
    pos = {m: k for k, m in enumerate(on_set)}
    masks: List[int] = []
    for pi in pis:
        bits = 0
        for m in on_set:
            if pi.covers(m):
                bits |= 1 << pos[m]
        masks.append(bits)

    # Essential: some on-minterm is covered by exactly one prime.
    cover_count = [0] * len(on_set)
    last_idx = [-1] * len(on_set)
    for idx, bits in enumerate(masks):
        b = bits
        while b:
            low = b & -b
            k = low.bit_length() - 1
            cover_count[k] += 1
            last_idx[k] = idx
            b ^= low
    essential_idx: Set[int] = {last_idx[k] for k in range(len(on_set))
                               if cover_count[k] == 1}

    chosen = set(essential_idx)
    covered = 0
    for idx in essential_idx:
        covered |= masks[idx]
    full = (1 << len(on_set)) - 1
    remaining = full & ~covered
    while remaining:
        best, best_gain = None, 0
        for idx, bits in enumerate(masks):
            if idx in chosen:
                continue
            gain = (bits & remaining).bit_count()
            if gain > best_gain:
                best, best_gain = idx, gain
        if best is None:
            break                      # unreachable if `pis` are true primes
        chosen.add(best)
        remaining &= ~masks[best]

    cover = [pis[i] for i in sorted(chosen)]
    essentials = [pis[i] for i in sorted(essential_idx)]
    return cover, essentials


def cube_minterms(pi: Implicant) -> Iterable[int]:
    """Enumerate every minterm inside a cube (2**n_free of them)."""
    free = [b for b in range(pi.n_vars) if (pi.mask >> b) & 1]
    for k in range(1 << len(free)):
        m = pi.val
        for t, b in enumerate(free):
            if (k >> t) & 1:
                m |= 1 << b
        yield m


def verify_cover(cover: Sequence[Implicant], on_set: Iterable[int],
                 off_set: Optional[Iterable[int]] = None,
                 allowed: Optional[Set[int]] = None) -> Tuple[bool, bool]:
    """(complete, sound) — every ON covered; no cube touches an OFF minterm.

    This is the runtime twin of ContactCircuits.lean's `cc_cover_complete` /
    `cc_off_avoiding` (which prove the property shape on a worked instance).

    Pass `allowed` (= ON u DC) to use the cube-enumeration path, which costs
    sum(2**n_free) instead of |cover| * |off_set| — the difference between
    seconds and hours once the off-set reaches 2**16.
    """
    complete = all(any(pi.covers(m) for pi in cover) for m in on_set)
    if allowed is not None:
        sound = all(m in allowed for pi in cover for m in cube_minterms(pi))
    else:
        sound = not any(pi.covers(m) for pi in cover for m in (off_set or ()))
    return complete, sound


@dataclass
class Circuit:
    """A minimised two-level (AND-OR) circuit plus its provenance."""

    name: str
    n_vars: int
    var_names: List[str]
    cover: List[Implicant]
    primes: List[Implicant]
    essentials: List[Implicant]
    n_on: int
    n_off: int
    n_dc: int
    complete: bool
    sound: bool
    meta: Dict = field(default_factory=dict)

    @property
    def n_literals(self) -> int:
        """Total literal count of the SOP — the standard two-level cost."""
        return sum(self.n_vars - pi.n_free for pi in self.cover)

    def to_sop(self, max_terms: Optional[int] = None) -> str:
        """Sum-of-products string, largest cubes first."""
        terms = sorted(self.cover, key=lambda p: (-p.n_free, p.val))
        shown = terms if max_terms is None else terms[:max_terms]
        parts = []
        for pi in shown:
            lits = pi.literals(self.var_names)
            parts.append("(" + " & ".join(lits) + ")" if lits else "1")
        s = " | ".join(parts)
        if max_terms is not None and len(terms) > max_terms:
            s += f"  | ... (+{len(terms) - max_terms} more terms)"
        return s

    def summary(self) -> Dict:
        return {
            "name": self.name,
            "n_vars": self.n_vars,
            "n_on": self.n_on,
            "n_off": self.n_off,
            "n_dc": self.n_dc,
            "n_primes": len(self.primes),
            "n_cover": len(self.cover),
            "n_essential": len(self.essentials),
            "n_literals": self.n_literals,
            "complete": self.complete,
            "sound": self.sound,
            **self.meta,
        }


def build_circuit(name: str, n_vars: int, var_names: Sequence[str],
                  on_set: Sequence[int], dc_set: Sequence[int],
                  meta: Optional[Dict] = None) -> Circuit:
    """Minimise a partially-specified Boolean function and verify the result."""
    on_list = sorted(set(on_set))
    dc_list = sorted(set(dc_set) - set(on_list))
    pis = prime_implicants(on_list, dc_list, n_vars)
    cover, essentials = _essential_and_greedy_cover(pis, on_list)
    # Soundness is checked by enumerating each cube and asserting every minterm
    # it contains is ON or DC -- exact, and independent of the off-set's size.
    total = 1 << n_vars
    allowed = set(on_list) | set(dc_list)
    n_off = total - len(allowed)
    complete, sound = verify_cover(cover, on_list, allowed=allowed)
    return Circuit(name=name, n_vars=n_vars, var_names=list(var_names),
                   cover=cover, primes=pis, essentials=essentials,
                   n_on=len(on_list), n_off=n_off, n_dc=len(dc_list),
                   complete=complete, sound=sound, meta=meta or {})


def _ref_pi_keys(on_set: Sequence[int], dc_set: Sequence[int],
                 n_vars: int) -> Set[Tuple[int, int]]:
    """Prime implicants of the reference implementation as (val, mask) keys.

    NOTE the reference is NOT exact QM: it keeps only don't-cares at Hamming
    distance 1 from a real on-minterm (its own documented performance choice).
    With clustered don't-cares it therefore under-merges — see
    tests/test_sequence_circuits.py::t06d, which pins down the exact
    consequence.  Kept here for cross-validation, not used in production.
    """
    if n_vars % 2 != 0:
        raise ValueError("reference implementation needs an even variable count")
    flat = np.zeros(1 << n_vars, dtype=int)
    if len(list(dc_set)):
        flat[list(dc_set)] = -1
    flat[list(on_set)] = 1
    ref = ks.qm_minimize(flat)

    keys: Set[Tuple[int, int]] = set()
    for entry in ref.get("prime_implicants", []):
        val = msk = 0
        for k, (v, m) in enumerate(zip(entry["values"], entry["mask"])):
            bit = n_vars - 1 - k
            if m:
                msk |= 1 << bit
            elif v:
                val |= 1 << bit
        keys.add((val, msk))
    return keys


def qm_agreement_check(on_set: Sequence[int], dc_set: Sequence[int],
                       n_vars: int) -> Dict:
    """Cross-validate this QM against `kmap_sbm.analysis.prime_implicants`.

    The set of prime implicants of a partially-specified function is canonical,
    so an exact implementation must reproduce it exactly.  Divergence is
    therefore diagnostic, not cosmetic — see `_ref_pi_keys`.
    """
    ref_pis = _ref_pi_keys(on_set, dc_set, n_vars)
    our_pis = {(p.val, p.mask) for p in prime_implicants(on_set, dc_set, n_vars)}
    return {
        "n_ref": len(ref_pis),
        "n_ours": len(our_pis),
        "identical": ref_pis == our_pis,
        "only_ref": sorted(ref_pis - our_pis)[:8],
        "only_ours": sorted(our_pis - ref_pis)[:8],
    }


# ---------------------------------------------------------------------------
# §2  R-B — the lossless positional sequence circuit
# ---------------------------------------------------------------------------

def _pos_bits(length: int) -> int:
    """Bits needed to index positions 0..length-1 (>=1)."""
    return max(1, (length - 1).bit_length())


def seq_circuit_positional(codes: Sequence[int], name: str = "seq",
                           aa_bits: int = 5, gray_residue: bool = True,
                           pad_as_dc: bool = False) -> Circuit:
    """R-B: minimal circuit for  f(i, a) = 1  iff  s[i] = a.

    Variable layout (MSB -> LSB):  i[p-1..0]  then  a[aa_bits-1..0].
    Positions use PLAIN BINARY (subcube = dyadic chain segment); residues use
    the He-2012 Gray image when `gray_residue` (subcube = chemical class).

    `pad_as_dc=False` (default) marks the padding region OFF, which keeps the
    function exactly equal to the sequence indicator on the whole cube and
    keeps QM cheap.  With `pad_as_dc=True` the padding becomes don't-care —
    smaller covers, but the circuit is then only *consistent with* the
    sequence rather than equal to it.  The lossless claim is made for the
    default; `reconstruct_from_circuit` asserts it either way.
    """
    L = len(codes)
    p = _pos_bits(L)
    n_vars = p + aa_bits
    on: List[int] = []
    for i, c in enumerate(codes):
        if c is None or c < 0 or c == ks.GAP_STATE:
            continue                      # gaps are not part of the on-set
        a = ks.gray_code(int(c)) if gray_residue else int(c)
        on.append((i << aa_bits) | a)

    dc: List[int] = []
    if pad_as_dc:
        total = 1 << n_vars
        onset = set(on)
        for m in range(total):
            i = m >> aa_bits
            a = m & ((1 << aa_bits) - 1)
            if i >= L or a >= (1 << aa_bits):
                dc.append(m)
        dc = [m for m in dc if m not in onset]

    var_names = [f"i{b}" for b in range(p - 1, -1, -1)] + \
                [f"a{b}" for b in range(aa_bits - 1, -1, -1)]
    return build_circuit(f"{name}:positional", n_vars, var_names, on, dc,
                         meta={"length": L, "pos_bits": p, "aa_bits": aa_bits,
                               "gray_residue": gray_residue,
                               "pad_as_dc": pad_as_dc,
                               "n_gaps": L - len(on)})


def reconstruct_from_circuit(circ: Circuit, length: int,
                             aa_bits: int = 5,
                             gray_residue: bool = True) -> List[Optional[int]]:
    """Invert R-B: evaluate the circuit over all (i,a) and read the sequence.

    Returns He codes, or None where the circuit fires on zero / multiple
    residues at a position.  Used to *prove by execution* that the circuit is
    a faithful representation (the losslessness claim).
    """
    inv_gray = {ks.gray_code(c): c for c in range(1 << aa_bits)}
    # Walk the cubes and mark what they fire on, rather than probing all
    # (position, residue) cells: cost becomes sum(2**n_free) instead of
    # length * 2**aa_bits * |cover|.
    hits: Dict[int, Set[int]] = defaultdict(set)
    aa_mask = (1 << aa_bits) - 1
    for pi in circ.cover:
        for m in cube_minterms(pi):
            i = m >> aa_bits
            if i < length:
                hits[i].add(m & aa_mask)
    out: List[Optional[int]] = []
    for i in range(length):
        h = hits.get(i, set())
        if len(h) != 1:
            out.append(None)
            continue
        a = next(iter(h))
        out.append(inv_gray[a] if gray_residue else a)
    return out


# ---------------------------------------------------------------------------
# §3  R-A — the k-mer (dipeptide) composition circuit
# ---------------------------------------------------------------------------

def seq_circuit_kmer(codes: Sequence[int], k: int = 2, name: str = "seq",
                     aa_bits: int = 5, min_count: int = 1) -> Circuit:
    """R-A: minimal circuit for  f(a_1..a_k) = 1  iff the k-mer occurs.

    This is the classical Petoukhov / FCGR-style K-map used elsewhere in the
    repo (32x32 for k=2).  Padding cells (any residue code >= 20) are
    don't-care, matching the FIX-A2 convention.
    """
    n_vars = k * aa_bits
    counts: Dict[int, int] = defaultdict(int)
    L = len(codes)
    for i in range(L - k + 1):
        window = [int(c) for c in codes[i:i + k]]
        if any(c < 0 or c == ks.GAP_STATE for c in window):
            continue
        idx = 0
        for c in window:
            idx = (idx << aa_bits) | ks.gray_code(c)
        counts[idx] += 1

    on = [m for m, n in counts.items() if n >= min_count]
    # Padding: any field whose Gray-decoded residue is >= 20.
    inv_gray = {ks.gray_code(c): c for c in range(1 << aa_bits)}
    dc = []
    for m in range(1 << n_vars):
        fields = [(m >> (aa_bits * (k - 1 - t))) & ((1 << aa_bits) - 1)
                  for t in range(k)]
        if any(inv_gray[f] >= 20 for f in fields):
            dc.append(m)
    onset = set(on)
    dc = [m for m in dc if m not in onset]

    var_names = [f"x{t}_{b}" for t in range(k)
                 for b in range(aa_bits - 1, -1, -1)]
    return build_circuit(f"{name}:{k}mer", n_vars, var_names, on, dc,
                         meta={"k": k, "n_distinct_kmers": len(counts),
                               "min_count": min_count,
                               "total_kmers": sum(counts.values())})


# ---------------------------------------------------------------------------
# §4  R-D — the contact map as a Boolean function of position bits
# ---------------------------------------------------------------------------

def contact_circuit_positional(contacts: Iterable[Tuple[int, int]], length: int,
                               name: str = "contact",
                               min_separation: int = ks.MIN_SEPARATION,
                               band_as_dc: bool = True,
                               pad_as_dc: bool = False) -> Circuit:
    """R-D: minimal circuit for  C(i, j) = 1  iff residues i, j are in contact.

    Positions in PLAIN BINARY, so a prime implicant with the low b_i bits of i
    and the low b_j bits of j free is exactly the statement

        "chain segment [i0, i0+2^b_i)  contacts  chain segment [j0, j0+2^b_j)"

    i.e. a SEGMENT-SEGMENT contact block.  This is the geometric content of the
    minimisation and the reason plain binary is used for positions.

    The near-diagonal band |i-j| < min_separation is excluded as don't-care
    (`band_as_dc`), matching the Morcos/Jones contact convention: trivially
    close in sequence, uninformative about the fold.

    `pad_as_dc=False` (default) marks index pairs beyond the chain end OFF, not
    don't-care.  This is both semantics and survival: a cube reaching past
    residue L would assert a contact between residues that do not exist, and
    for a 130-residue chain padded to 256 the don't-care region is 48,636 of
    65,536 cells — enough for QM to sprawl into cubes that describe the padding
    rather than the protein.  Pinning it OFF keeps every cube inside the chain.
    """
    p = _pos_bits(length)
    n_vars = 2 * p
    cset = {(min(i, j), max(i, j)) for i, j in contacts}
    on: List[int] = []
    for i, j in cset:
        if i >= length or j >= length:
            continue
        on.append((i << p) | j)
        on.append((j << p) | i)          # symmetric map

    dc: List[int] = []
    onset = set(on)
    for i in range(1 << p):
        for j in range(1 << p):
            m = (i << p) | j
            if m in onset:
                continue
            if i >= length or j >= length:
                if pad_as_dc:
                    dc.append(m)         # padding beyond the chain
            elif band_as_dc and abs(i - j) < min_separation:
                dc.append(m)             # excluded near-diagonal band

    var_names = [f"i{b}" for b in range(p - 1, -1, -1)] + \
                [f"j{b}" for b in range(p - 1, -1, -1)]
    return build_circuit(f"{name}:contactmap", n_vars, var_names, on, dc,
                         meta={"length": length, "pos_bits": p,
                               "n_contacts_unique": len(cset),
                               "min_separation": min_separation,
                               "band_as_dc": band_as_dc,
                               "pad_as_dc": pad_as_dc})


def implicant_to_segments(pi: Implicant, pos_bits: int,
                          length: int) -> Dict:
    """Decode an R-D implicant into the two chain segments it relates.

    Returns the (start, end) half-open ranges for i and j, clipped to the
    chain.  A cube whose i-field has f_i free bits spans 2**f_i consecutive
    positions ONLY IF the free bits are the low ones; otherwise the span is a
    union of strided blocks and `contiguous` is False.
    """
    def _field(shift: int) -> Dict:
        m = (pi.mask >> shift) & ((1 << pos_bits) - 1)
        v = (pi.val >> shift) & ((1 << pos_bits) - 1)
        free = bin(m).count("1")
        low_free = (m & (m + 1)) == 0            # mask is 0b0..011..1
        start = v & ~m & ((1 << pos_bits) - 1)
        return {"start": start,
                "end": min(start + (1 << free), length) if low_free else None,
                "n_free": free, "contiguous": bool(low_free),
                "span": 1 << free}
    return {"i": _field(pos_bits), "j": _field(0),
            "n_free": pi.n_free, "size": pi.size}


# ---------------------------------------------------------------------------
# §5  R-C — co-evolution circuits over an MSA (three polarities)
# ---------------------------------------------------------------------------

POLARITIES = ("standard", "flipped", "perplexity")


def coevo_pair_table(dense: np.ndarray, i: int, j: int,
                     aa_bits: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """Joint count table for MSA columns (i, j) on the padded 32x32 grid.

    Returns (counts[32,32], observed_mask[32,32]).  Rows/columns 20-31 are the
    FIX-A2 padding and always stay unobserved.  Gaps excluded.
    """
    n_states = 1 << aa_bits
    counts = np.zeros((n_states, n_states), dtype=np.int64)
    col_i = dense[:, i]
    col_j = dense[:, j]
    ok = (col_i != ks.GAP_STATE) & (col_j != ks.GAP_STATE) & \
         (col_i >= 0) & (col_j >= 0)
    for a, b in zip(col_i[ok], col_j[ok]):
        counts[ks.gray_code(int(a)), ks.gray_code(int(b))] += 1
    return counts, counts > 0


def coevo_circuit(dense: np.ndarray, i: int, j: int, polarity: str,
                  aa_bits: int = 5, min_count: int = 1,
                  perplexity_threshold: float = 1.5,
                  unobservable: str = "off",
                  name: Optional[str] = None) -> Circuit:
    """R-C: circuit over one MSA column pair, in one of three polarities.

    standard    ON = combination observed at least `min_count` times.
                Reads as: "which residue pairs does evolution ALLOW here?"
    flipped     ON = combination never observed, both residues individually
                present in their column.  Reads as: "which pairs are
                FORBIDDEN?"  (negative selection; the repo's 490-rule set.)
                Combinations involving a residue absent from its own column
                are don't-care — they are unobservable, not forbidden.  That
                distinction is what makes the flipped map informative rather
                than a restatement of the marginals.
    perplexity  ON = joint perplexity ratio  PP(a)*PP(b)/PP(a,b)  above
                `perplexity_threshold`.  Surprise-based rather than
                MI-based: it fires on pairs whose JOINT occurrence is far
                more predictable than the two marginals separately.

    Padding cells (Gray-decoded residue >= 20) are always don't-care.

    `unobservable` decides what happens to a canonical cell whose residue is
    absent from its own column — and it matters a great deal, because a Spike
    column carries only 2-4 distinct residues, leaving ~99% of the 32x32 grid
    unobservable:

      "off"  (default) cubes cannot grow into never-seen residues.  Rules stay
             inside the observed alphabet: tight, testable, but silent about
             residues the alignment never sampled.
      "dc"   cubes may absorb the unobservable region.  Formally sound, but on
             a 2-residue column QM then returns near-vacuous rules such as
             "DEFMNQWY x CDEGHKNPQRST is forbidden" — an artefact of the
             don't-care mass, not a finding.  This is the same degeneracy that
             sank Phase 1 of 10_Boolean_Contact_Circuits.md.

    The campaign reports both so the sensitivity is on the record.
    """
    if polarity not in POLARITIES:
        raise ValueError(f"polarity must be one of {POLARITIES}")
    if unobservable not in ("off", "dc"):
        raise ValueError("unobservable must be 'off' or 'dc'")
    n_states = 1 << aa_bits
    n_vars = 2 * aa_bits
    counts, observed = coevo_pair_table(dense, i, j, aa_bits)
    total = int(counts.sum())
    inv_gray = {ks.gray_code(c): c for c in range(n_states)}
    canonical = np.array([inv_gray[g] < 20 for g in range(n_states)])

    row_present = (counts.sum(axis=1) > 0)
    col_present = (counts.sum(axis=0) > 0)

    on: List[int] = []
    dc: List[int] = []
    for ga in range(n_states):
        for gb in range(n_states):
            m = (ga << aa_bits) | gb
            if not (canonical[ga] and canonical[gb]):
                dc.append(m)
                continue
            n_ab = int(counts[ga, gb])
            if polarity == "standard":
                if n_ab >= min_count:
                    on.append(m)          # observed -> allowed; else off-set
            elif polarity == "flipped":
                if not (row_present[ga] and col_present[gb]):
                    if unobservable == "dc":
                        dc.append(m)      # unobservable, not forbidden
                    # else: off-set, cubes stay inside the observed alphabet
                elif n_ab == 0:
                    on.append(m)
            else:  # perplexity
                if not (row_present[ga] and col_present[gb]):
                    if unobservable == "dc":
                        dc.append(m)
                elif n_ab == 0:
                    pass                  # off-set
                else:
                    p_a = counts[ga, :].sum() / total
                    p_b = counts[:, gb].sum() / total
                    p_ab = n_ab / total
                    # PP(x) = 1/p(x); ratio = PP(a)*PP(b)/PP(a,b) = p_ab/(p_a*p_b)
                    ratio = p_ab / (p_a * p_b) if p_a > 0 and p_b > 0 else 0.0
                    if ratio >= perplexity_threshold:
                        on.append(m)

    onset = set(on)
    dc = sorted(set(dc) - onset)
    var_names = [f"a{b}" for b in range(aa_bits - 1, -1, -1)] + \
                [f"b{b}" for b in range(aa_bits - 1, -1, -1)]
    nm = name or f"pair_{i}_{j}:{polarity}"
    return build_circuit(nm, n_vars, var_names, on, dc,
                         meta={"col_i": i, "col_j": j, "polarity": polarity,
                               "n_sequences": total,
                               "n_distinct_observed": int(observed.sum()),
                               "min_count": min_count,
                               "unobservable": unobservable,
                               "alphabet_i": sorted(
                                   ks.CODE_TO_AA[inv_gray[g]]
                                   for g in range(n_states)
                                   if row_present[g] and canonical[g]),
                               "alphabet_j": sorted(
                                   ks.CODE_TO_AA[inv_gray[g]]
                                   for g in range(n_states)
                                   if col_present[g] and canonical[g]),
                               "perplexity_threshold": perplexity_threshold})


def decode_residue_pair_implicant(pi: Implicant, aa_bits: int = 5) -> Dict:
    """Decode an R-C / R-A implicant into the residue sets it matches.

    Returns the canonical (code < 20) residues admitted at each of the two
    slots, plus their physicochemical groups.  A cube with free bits admits a
    SET of residues — that set is the rule's tolerance class.
    """
    n_states = 1 << aa_bits
    inv_gray = {ks.gray_code(c): c for c in range(n_states)}

    def _slot(shift: int) -> Dict:
        m = (pi.mask >> shift) & (n_states - 1)
        v = (pi.val >> shift) & (n_states - 1)
        letters, groups = [], set()
        for g in range(n_states):
            if (g & ~m) != v:
                continue
            code = inv_gray[g]
            if code >= 20:
                continue
            letters.append(ks.CODE_TO_AA[code])
            groups.add(ks.GROUP_NAMES[ks.AA_GROUP8[ks.CODE_TO_AA[code]]])
        return {"residues": sorted(letters), "n_free": bin(m).count("1"),
                "groups": sorted(groups)}

    return {"slot_i": _slot(aa_bits), "slot_j": _slot(0),
            "n_free": pi.n_free}


# ---------------------------------------------------------------------------
# §6  FASTA loading for the Spike alignment
# ---------------------------------------------------------------------------

def load_fasta_codes(path, max_seqs: Optional[int] = None
                     ) -> Tuple[np.ndarray, List[str]]:
    """Load an aligned FASTA -> (dense He codes [N,L] int8, headers).

    Non-canonical characters (gaps, X, B, Z, ...) map to GAP_STATE, matching
    `coevolution_shared.py`'s `aligned=True` convention (the A1 fix: alignment
    columns are preserved, gaps are a state that is excluded from statistics
    rather than stripped).
    """
    headers: List[str] = []
    seqs: List[str] = []
    cur: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
                if max_seqs is not None and len(headers) >= max_seqs:
                    break
                headers.append(line[1:])
            else:
                cur.append(line.strip())
    if cur and (max_seqs is None or len(seqs) < len(headers)):
        seqs.append("".join(cur))
    seqs = seqs[:len(headers)]

    width = max(len(s) for s in seqs)
    dense = np.full((len(seqs), width), ks.GAP_STATE, dtype=np.int8)
    for r, s in enumerate(seqs):
        for c, ch in enumerate(s):
            code = ks.AA_TO_CODE.get(ch)
            if code is not None:
                dense[r, c] = code
    return dense, headers


def column_entropy(dense: np.ndarray, col: int) -> float:
    """Shannon entropy (bits) of one alignment column, gaps excluded."""
    v = dense[:, col]
    v = v[(v != ks.GAP_STATE) & (v >= 0)]
    if v.size == 0:
        return 0.0
    _, cnt = np.unique(v, return_counts=True)
    p = cnt / cnt.sum()
    return float(-(p * np.log2(p)).sum())
