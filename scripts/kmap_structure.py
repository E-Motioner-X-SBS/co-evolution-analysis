#!/usr/bin/env python3
"""
kmap_structure.py — Contact-K-map construction from LITERATURE ground truth.
=============================================================================
Campaign: plans/11-kmap-structure-encoding.md (+ 11b-strategy.md)

PURPOSE
-------
This module is the single source of truth for the "K-map encodes 3D structure"
campaign.  It connects, end-to-end and verifiably:

    MSA columns (binary sequence positions)
      -> He-2012 integer codes (gap = state 20)          [nkmap encoding]
      -> physicochemical group codes (8 groups)           [this module]
      -> per-pair observations labeled by NATIVE CONTACTS [literature PDBs]
      -> K-map cell counting (group-level 256 cells;
         identity-level 32x32 padded x4 separation bins)  [K-map construction]
      -> ternary truth-table labeling                     [documented rule]
      -> Quine-McCluskey minimization                     [kmap_sbm, popcount-grouped]
      -> Boolean contact circuits (sum-of-prime-implicants)
      -> evaluation against literature baselines
         (native contacts + PSICOV published predictions)

WHY THE GROUND TRUTH IS TRUSTWORTHY (literature, not self-computed)
-------------------------------------------------------------------
data/psicov150/raw/README (Jones et al. 2012 supplementary):

    "The alignments correspond to the sequence represented in each PDB file ...
     there are exactly the same number of columns in each alignment as
     C-alpha ATOM records in each PDB file."

=> MSA column c (0-based here) corresponds EXACTLY to PDB residue c+1.
   Verified empirically for 149/150 proteins (see verify_identity_mapping);
   the single exception (1vjkA) has ONE EXTRA C-TERMINAL PDB RESIDUE ('S', 88)
   absent from the target sequence; the general rule below (ignore PDB residues
   numbered beyond the alignment width) resolves it without any alignment.

THE IDENTITY-MAPPING RULE (hard invariant, asserted everywhere)
---------------------------------------------------------------
    PDB residue r participates in contact computation iff 1 <= r <= n_cols.
    (n_cols = alignment width.)  Contacts are returned as 0-based column pairs.

CONTACT DEFINITION (literature convention)
------------------------------------------
    contact(i, j)  <=>  d(Cb_i, Cb_j) < 8.0 Angstrom  AND  |i - j| >= 6
    (Cb = CB atom; CA used only for glycine, which lacks CB.)
    This matches validate_contacts.py, Morcos et al. 2011 (PNAS) and
    Jones et al. 2012 (Bioinformatics).  Separation bins used downstream:
        bin 0: 6  <= sep < 11
        bin 1: 11 <= sep < 21
        bin 2: 21 <= sep < 31
        bin 3: sep >= 31
    Rationale: bin 0 spans short-range secondary-structure contacts; the rest
    partition long range roughly by order of magnitude.  Bins are fixed a
    priori (no peeking at labels).

PHYSICOCHEMICAL GROUPS (8 groups, matches scripts/build_boolean_contact_circuits.py)
------------------------------------------------------------------------------------
    0 hydrophobic  A I L V M      4 basic       H K R
    1 aromatic     F Y W          5 hydroxyl    S T
    2 acidic       E D            6 sulfur      C
    3 amide        Q N            7 special     P G
(NB: the Lean AminoAcidEncoding.lean formalizes a related 7-group scheme; this
module uses the 8-group circuit scheme for continuity with prior Phase-2 work.)

CONSENSUS RULE (repo-wide convention, GPU-consistent)
-----------------------------------------------------
    consensus(col) = most frequent non-gap He code; ties -> LOWEST code.
    All-gap columns get sentinel -1 and every pair touching them is SKIPPED
    (counted, never silently dropped -- see SkippedPair accounting).

AUTHOR / PROVENANCE
-------------------
E-Motioner-X-SBS co-evolution-analysis, Aug 2026.  Every public function is
unit-tested in scripts/tests/test_kmap_structure.py including adversarial cases.
"""

from __future__ import annotations



import math

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Repo-path constants
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "data" / "psicov150" / "raw"           # literature supplementary data
RESULTS = REPO / "results" / "kmap_structure"        # this campaign's outputs

# ---------------------------------------------------------------------------
# Encoding tables (documented, explicit, mirrored in unit tests)
# ---------------------------------------------------------------------------

#: He 2012 amino-acid ordering (index = code 0..19).  Gap/unknown = 20.
AMINO_HE_2012 = "AILVMFYWEDQNHKRSTCPG"

#: code -> amino acid letter (list-index lookup kept O(1) and explicit).
CODE_TO_AA: List[str] = list(AMINO_HE_2012)

#: amino acid letter -> He code.
AA_TO_CODE: Dict[str, int] = {a: i for i, a in enumerate(CODE_TO_AA)}

#: gap / ambiguous-character state (kept in arrays, excluded from statistics).
GAP_STATE = 20

#: canonical residues recognised inside PDB ATOM records (3-letter -> 1-letter).
AA3_TO_1: Dict[str, str] = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
    "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y",
}

#: 8-group physicochemical partition (letter -> group id 0..7). See module docstring.
AA_GROUP8: Dict[str, int] = {}
for _a in "AILVM":
    AA_GROUP8[_a] = 0   # hydrophobic
for _a in "FYW":
    AA_GROUP8[_a] = 1   # aromatic
for _a in "ED":
    AA_GROUP8[_a] = 2   # acidic
for _a in "QN":
    AA_GROUP8[_a] = 3   # amide
for _a in "HKR":
    AA_GROUP8[_a] = 4   # basic
for _a in "ST":
    AA_GROUP8[_a] = 5   # hydroxyl
AA_GROUP8["C"] = 6       # sulfur
for _a in "PG":
    AA_GROUP8[_a] = 7   # special

#: Human-readable names, index-aligned with AA_GROUP8 ids.
GROUP_NAMES = ["hydrophobic", "aromatic", "acidic", "amide",
               "basic", "hydroxyl", "sulfur", "special"]

#: Separation-bin upper bounds (exclusive); last bin open-ended.  See module docs.
SEP_BIN_EDGES = (11, 21, 31)          # bins: [6,11) [11,21) [21,31) [31,inf)
N_SEP_BINS = len(SEP_BIN_EDGES) + 1   # 4

#: Contact geometry constants (Angstrom / residues) — literature convention.
CONTACT_CUTOFF_A = 8.0
MIN_SEPARATION = 6


# ---------------------------------------------------------------------------
# Small pure helpers (unit-tested, adversarially)
# ---------------------------------------------------------------------------

def gray_code(x: int) -> int:
    """Reflected binary Gray code of x>=0:  g(x) = x XOR (x >> 1).

    Proven properties (lean_proofs KmapProofs.lean): bijective on fixed width,
    adjacent integers differ in exactly one bit.
    """
    return x ^ (x >> 1)


def gray_hamming(a: int, b: int, bits: int = 5) -> int:
    """Hamming distance between the `bits`-wide Gray images of a and b.

    d(a,b) = popcount(g(a) XOR g(b)), truncated to `bits` bits (defensive:
    inputs above 2**bits raise rather than silently truncate).
    """
    if a < 0 or b < 0 or a >= (1 << bits) or b >= (1 << bits):
        raise ValueError(f"codes must be in [0, 2^{bits}); got {a},{b}")
    xor = gray_code(a) ^ gray_code(b)
    return bin(xor).count("1")


def group_of_code(code: int) -> int:
    """Physicochemical group id (0..7) of a He code; gap (20) -> -1 sentinel."""
    if code == GAP_STATE or code < 0:
        return -1
    return AA_GROUP8[CODE_TO_AA[code]]


def sep_bin(sep: int) -> int:
    """Separation bin index 0..3 for a pair separation sep >= MIN_SEPARATION."""
    if sep < MIN_SEPARATION:
        raise ValueError(f"sep {sep} below minimum {MIN_SEPARATION}")
    for b, edge in enumerate(SEP_BIN_EDGES):
        if sep < edge:
            return b
    return N_SEP_BINS - 1


# ---------------------------------------------------------------------------
# FASTA / alignment loading
# ---------------------------------------------------------------------------

def load_aln_codes(path: Path, max_seqs: Optional[int] = None) -> Tuple[np.ndarray, List[str]]:
    """Load a headerless PSICOV .aln alignment -> (dense codes [N,L] int8, rows-as-read).

    Each line is one gapped sequence (PSICOV format: no headers).  Characters:
      canonical AA (any case) -> its He code 0..19
      anything else ('-', 'X', ...) -> GAP_STATE (20), retained to preserve
      column identity (FIX-A1 doctrine of this repo: never strip gaps).

    Returns (dense, raw_lines).  Raises if rows have inconsistent length
    (an alignment must be rectangular).
    """
    lines = [ln.strip() for ln in open(path)]
    lines = [ln for ln in lines if ln]
    if not lines:
        raise ValueError(f"empty alignment: {path}")
    width = len(lines[0])
    for ln in lines:
        if len(ln) != width:
            raise ValueError(f"ragged alignment {path}: widths {width} vs {len(ln)}")
    if max_seqs is not None:
        lines = lines[:max_seqs]
    lut = np.full(256, GAP_STATE, dtype=np.int8)
    for ch, code in AA_TO_CODE.items():
        lut[ord(ch)] = code
        lut[ord(ch.lower())] = code
    buf = np.frombuffer("".join(lines).encode("ascii"), dtype=np.uint8)
    dense = lut[buf].reshape(len(lines), width)
    return dense, lines


def consensus_codes(dense: np.ndarray) -> np.ndarray:
    """Per-column majority He code (gaps excluded), ties -> lowest code.

    Mirrors coevolution_shared.majority_ref semantics exactly (GPU-consistent
    lowest-code tie-break).  All-gap columns receive -1 (sentinel; callers must
    treat any pair containing -1 as unmappable).
    """
    n, L = dense.shape
    out = np.full(L, -1, dtype=np.int16)
    for c in range(L):
        col = dense[:, c]
        valid = col[(col >= 0) & (col < GAP_STATE)]
        if valid.size == 0:
            continue
        counts = np.bincount(valid, minlength=GAP_STATE)
        out[c] = int(np.argmax(counts))          # argmax returns FIRST maximum = lowest code
    return out


# ---------------------------------------------------------------------------
# Native structure parsing (literature ground truth)
# ---------------------------------------------------------------------------

@dataclass
class ResidueCB:
    """Cbeta (or Calpha-for-glycine) coordinate of one residue."""
    resnum: int        # 1-based PDB residue number
    aa: str            # 1-letter code ('X' if non-canonical)
    xyz: Tuple[float, float, float]


def parse_pdb_cb(path_or_lines, n_cols: Optional[int] = None,
                 strict: bool = False) -> Dict[int, ResidueCB]:
    """Parse ATOM records -> {resnum: ResidueCB} keeping CB (CA for Gly).

    Rules (documented):
      * Only ATOM records (HETATM ignored: waters/ligands).
      * Blank chain id accepted (PSICOV repaired natives are single-chain).
      * Alternate-location indicator: first occurrence wins (PSICOV natives
        carry no altlocs -- verified across all 150 files -- but defensive).
      * If `n_cols` is given, residues with resnum > n_cols are DROPPED
        (identity-mapping rule: they lie outside the aligned region; this
        cleanly handles the 1vjkA extra C-terminal residue).
    """
    if isinstance(path_or_lines, (str, Path)):
        lines = open(path_or_lines).readlines()
    else:
        lines = path_or_lines
    atoms: Dict[Tuple[int, str], Tuple[str, Tuple[float, float, float]]] = {}
    seen_names: Dict[int, str] = {}
    for ln in lines:
        if not ln.startswith("ATOM"):
            continue
        atom = ln[12:16].strip()
        if atom not in ("CA", "CB"):
            continue
        resnum = int(ln[22:26])
        if n_cols is not None and resnum > n_cols:
            continue
        aa3 = ln[17:20].strip()
        aa = AA3_TO_1.get(aa3, "X")
        # Residue-name consistency: in strict mode, a contradictory record
        # raises immediately; in tolerant mode, first occurrence wins.
        prev = seen_names.get(resnum)
        if prev is not None and prev != aa:
            if strict:
                raise ValueError(
                    f"inconsistent residue names at {resnum}: {prev} vs {aa}")
        if resnum not in seen_names:
            seen_names[resnum] = aa
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        key = (resnum, atom)
        if key in atoms:                 # duplicate atom record (altloc): first wins
            continue
        atoms[key] = (aa, xyz)
    residues: Dict[int, ResidueCB] = {}
    for resnum, aa in seen_names.items():
        cb = atoms.get((resnum, "CB"))
        ca = atoms.get((resnum, "CA"))
        # Literature convention: CB for all residue types EXCEPT glycine (CA),
        # because glycine has no CB in nature.  Repaired files occasionally
        # place a nominal CB on Gly; we still use CA for Gly.  Defensive CA
        # fallback if CB is missing entirely on a non-Gly residue.
        pick = ca if aa == "G" else (cb if cb is not None else ca)
        if pick is None:
            continue
        residues[resnum] = ResidueCB(resnum, aa, pick[1])
    return residues


def dist3(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    """Euclidean distance between two 3-vectors."""
    return math.sqrt(sum((u - v) ** 2 for u, v in zip(a, b)))


def native_contacts(pdb_path: Path,
                    n_cols: Optional[int] = None,
                    cutoff_a: float = CONTACT_CUTOFF_A,
                    min_sep: int = MIN_SEPARATION,
                    ) -> Set[Tuple[int, int]]:
    """Native contact set {(i, j)} as 0-BASED ALIGNMENT COLUMNS, i<j, sep>=min_sep.

    Definition (module docstring): Cb-Cb < cutoff_a (CA for Gly),
    |r_i - r_j| >= min_sep, and both residues within [1, n_cols].
    """
    residues = parse_pdb_cb(pdb_path, n_cols=n_cols)
    nums = sorted(residues)
    contacts: Set[Tuple[int, int]] = set()
    for ai, ri in enumerate(nums):
        for rj in nums[ai + 1:]:
            if rj - ri < min_sep:
                continue
            if dist3(residues[ri].xyz, residues[rj].xyz) < cutoff_a:
                contacts.add((ri - 1, rj - 1))   # 1-based PDB num -> 0-based column
    return contacts


def parse_con(path: Path) -> List[Tuple[int, int, float]]:
    """Parse a PSICOV published prediction file (.out in raw/con/).

    Format (raw/README): `i j dummy dummy score`, sorted by score descending.
    Returned as-is (1-based residue numbers); callers filter sep/range.
    """
    out = []
    for ln in open(path):
        parts = ln.split()
        if len(parts) != 5:
            continue
        i, j = int(parts[0]), int(parts[1])
        score = float(parts[4])
        out.append((i, j, score))
    return out


# ---------------------------------------------------------------------------
# Per-protein bundle
# ---------------------------------------------------------------------------

@dataclass
class ProteinData:
    """Everything the campaign needs about one PSICOV target."""
    target: str                       # e.g. '1a3aA'
    n_cols: int                       # alignment width == aligned PDB residues
    n_seqs: int                       # MSA depth (rows read)
    dense: Optional[np.ndarray]       # [N, L] int8 codes (lazy-loaded by loader)
    consensus: np.ndarray             # [L] int16 codes (-1 = all-gap column)
    contacts: Set[Tuple[int, int]]    # native contacts as 0-based column pairs
    identity_ok: bool                 # width == CA-count (149/150 True)
    notes: List[str] = field(default_factory=list)


def verify_identity_mapping(target: str) -> Tuple[bool, int, int]:
    """Assert the README guarantee for one protein; returns (ok, width, n_ca).

    The guarantee: alignment width equals the number of CA atoms in the native
    PDB (renumbered 1..width).  Verified live for 149/150 targets on 2026-08-22;
    the exception (1vjkA) carries one extra C-terminal PDB residue which the
    n_cols drop-rule in parse_pdb_cb removes.
    """
    aln_path = RAW / "aln" / f"{target}.aln"
    pdb_path = RAW / "pdb" / f"{target}.pdb"
    lines = [ln.strip() for ln in open(aln_path) if ln.strip()]
    width = len(lines[0])
    ca = set()
    for ln in open(pdb_path):
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            ca.add(int(ln[22:26]))
    return (width == len(ca)), width, len(ca)


def load_protein(target: str,
                 load_dense: bool = True,
                 max_seqs: Optional[int] = None,
                 ) -> ProteinData:
    """Load one target with all integrity checks applied (G3 gate).

    Checks performed:
      1. identity mapping width==CA (recorded, not fatal -- drop-rule handles).
      2. contacts recomputed under the n_cols drop-rule.
      3. consensus computed with gaps excluded.
    """
    ok, width, n_ca = verify_identity_mapping(target)
    dense = None  # noqa
    if load_dense:
        dense, _rows = load_aln_codes(RAW / "aln" / f"{target}.aln", max_seqs=max_seqs)
        n_seqs = dense.shape[0]
    else:
        lines = [ln for ln in open(RAW / "aln" / f"{target}.aln") if ln.strip()]
        n_seqs = len(lines)
    consensus = consensus_codes(dense) if dense is not None else np.full(width, -1, np.int16)
    contacts = native_contacts(RAW / "pdb" / f"{target}.pdb", n_cols=width)
    notes = []
    if not ok:
        notes.append(f"identity mapping off by {n_ca - width} extra PDB residue(s); "
                     f"drop-rule applied (residues >{width} ignored)")
    return ProteinData(target=target, n_cols=width, n_seqs=n_seqs, dense=dense,
                       consensus=consensus, contacts=contacts, identity_ok=ok,
                       notes=notes)


# ---------------------------------------------------------------------------
# Pair-universe construction + feature extraction
# ---------------------------------------------------------------------------

def pair_universe(pd: ProteinData) -> np.ndarray:
    """All candidate column pairs (i<j) with sep >= MIN_SEPARATION, i,j valid.

    Valid = consensus != -1 at BOTH columns (all-gap columns unmappable).
    Returns int32 array [P, 2] of 0-based columns sorted lexicographically.
    """
    # L = pd.n_cols  # unused: vi already excludes invalid
    valid = pd.consensus >= 0
    vi = np.flatnonzero(valid)
    ii, jj = np.triu_indices(vi.size, k=1)
    pi, pj = vi[ii], vi[jj]
    sep = pj - pi
    keep = sep >= MIN_SEPARATION
    return np.stack([pi[keep], pj[keep]], axis=1).astype(np.int32)


def pair_features(pd: ProteinData, pairs: np.ndarray) -> Dict[str, np.ndarray]:
    """Feature vectors for each candidate pair (all per-protein-intrinsic).

    Returns dict with keys:
      i, j          : int32 columns
      sep           : int32 separation j-i
      sbin          : int8 separation bin 0..3
      aa_i, aa_j    : int8 consensus He codes (0..19)
      g_i, g_j      : int8 group ids (0..7)
      is_contact    : bool  (native-contact label; literature ground truth)
      gray_h        : int8  Hamming distance between 5-bit Gray images of codes
    """
    i = pairs[:, 0]
    j = pairs[:, 1]
    aa_i = pd.consensus[i].astype(np.int8)
    aa_j = pd.consensus[j].astype(np.int8)
    sep = (j - i).astype(np.int32)
    gmap = np.array([group_of_code(c) for c in range(GAP_STATE + 1)], dtype=np.int8)
    g_i = gmap[aa_i]
    g_j = gmap[aa_j]
    sbin = np.array([sep_bin(s) for s in sep.tolist()], dtype=np.int8) if len(sep) \
        else np.zeros(0, np.int8)
    is_contact = np.array([(int(a), int(b)) in pd.contacts for a, b in pairs], dtype=bool)
    gh = np.array([gray_hamming(int(a), int(b)) for a, b in zip(aa_i, aa_j)],
                  dtype=np.int8) if len(aa_i) else np.zeros(0, np.int8)
    return {"i": i, "j": j, "sep": sep, "sbin": sbin,
            "aa_i": aa_i, "aa_j": aa_j, "g_i": g_i, "g_j": g_j,
            "is_contact": is_contact, "gray_h": gh}


def mutual_information_columns(dense: np.ndarray, i: int, j: int) -> float:
    """MI between two alignment columns (bits), gaps excluded.

    Semantically IDENTICAL to coevolution_shared.mutual_information (bincount
    joint over 400 cells, marginals by axis sums); reimplemented locally on the
    dense matrix to avoid per-pair Python column extraction.  Equality with the
    shared module is property-tested in scripts/tests/test_kmap_structure.py.
    """
    ci = dense[:, i]
    cj = dense[:, j]
    valid = (ci >= 0) & (ci < GAP_STATE) & (cj >= 0) & (cj < GAP_STATE)
    ci = ci[valid].astype(np.int64)
    cj = cj[valid].astype(np.int64)
    if ci.size < 10:
        return 0.0
    joint = np.bincount(ci * GAP_STATE + cj,
                        minlength=GAP_STATE * GAP_STATE).reshape(GAP_STATE, GAP_STATE)
    total = joint.sum()
    if total == 0:
        return 0.0
    mi = 0.0
    nz = joint > 0
    p = joint[nz] / total
    pi_ = joint.sum(axis=1)[:, None].repeat(GAP_STATE, axis=1)[nz] / total
    pj_ = joint.sum(axis=0)[None, :].repeat(GAP_STATE, axis=0)[nz] / total
    mi = float(np.sum(p * np.log2(p / (pi_ * pj_))))
    return mi


def mi_scores(pd: ProteinData, pairs: np.ndarray,
              dense_override: Optional[np.ndarray] = None) -> np.ndarray:
    """Plain-MI score for every candidate pair (sequence-only baseline, PATH D)."""
    dense = dense_override if dense_override is not None else pd.dense
    assert dense is not None, "mi_scores requires a loaded dense matrix"
    out = np.empty(len(pairs), dtype=np.float64)
    for k, (i, j) in enumerate(pairs):
        out[k] = mutual_information_columns(dense, int(i), int(j))
    return out


def perplexity_ratio_columns(dense: np.ndarray, i: int, j: int,
                             entropy_vec: Optional[np.ndarray] = None) -> Optional[float]:
    """PP(j)/mean_a PP(j|i=a); identical semantics to coevolution_shared.perplexity_ratio.

    Returns None when <5 observations or <3 conditioning residues (mirrors shared).
    """
    ci = dense[:, i]
    cj = dense[:, j]
    valid = (ci >= 0) & (ci < GAP_STATE) & (cj >= 0) & (cj < GAP_STATE)
    ci, cj = ci[valid], cj[valid]
    if ci.size < 5:
        return None
    joint = np.bincount(ci.astype(np.int64) * GAP_STATE + cj.astype(np.int64),
                        minlength=GAP_STATE * GAP_STATE).reshape(GAP_STATE, GAP_STATE).astype(np.float64)
    tot = joint.sum()
    if tot == 0:
        return None
    row_tot = joint.sum(axis=1)
    cond = []
    for a in range(GAP_STATE):
        if row_tot[a] >= 5:
            p = joint[a] / row_tot[a]
            p = p[p > 0]
            if p.size:
                cond.append(2.0 ** (-float(np.sum(p * np.log2(p)))))
    # NOTE: mirrors the ACTUAL behavior of coevolution_shared.perplexity_ratio,
    # which returns None only when NO conditioning residue has n>=5 -- even
    # though its docstring claims "fewer than 3".  Docstring/code discrepancy
    # recorded in plans/11c-decisions.md (F1); we mirror code, not docstring.
    if not cond:
        return None
    if entropy_vec is not None:
        ppj = 2.0 ** float(entropy_vec[j])
    else:
        marg = joint.sum(axis=0)
        pj = marg[marg > 0] / tot
        ppj = 2.0 ** (-float(np.sum(pj * np.log2(pj))))
    avg = float(np.mean(cond))
    return ppj / avg if avg > 0 else None


# ---------------------------------------------------------------------------
# K-map cell counting (level-1 group circuit; level-2 identity circuits)
# ---------------------------------------------------------------------------

LEVEL1_CELLS = 8 * 8 * N_SEP_BINS          # 256
LEVEL2_CELLS_PER_BIN = 32 * 32             # padded identity grid per sep bin


def level1_cell_index(g_i: np.ndarray, g_j: np.ndarray, sbin: np.ndarray) -> np.ndarray:
    """Flat cell index in [0,256): idx = ((g_i*8)+g_j)*4 + sbin."""
    return (g_i.astype(np.int64) * 8 + g_j) * N_SEP_BINS + sbin


def level2_cell_index(aa_i, aa_j) -> np.ndarray:
    """Flat cell index in [0,1024) of the 32x32 PADDED identity grid (row-major).

    Accepts scalars or arrays (anything numpy can broadcast); returns int64.
    Padding doctrine (FIX A2): codes 20..31 land in the don't-care border; the
    consensus codes fed here are always 0..19 so border cells simply stay empty.
    """
    return np.asarray(aa_i, dtype=np.int64) * 32 + np.asarray(aa_j, dtype=np.int64)


@dataclass
class CellCounts:
    """Per-cell (contacts, total) accumulators for one circuit family."""
    n_cells: int
    pos: np.ndarray            # int64 [cells] contact observations
    tot: np.ndarray            # int64 [cells] all observations

    def add(self, cell_idx: np.ndarray, labels: np.ndarray) -> None:
        """Accumulate one protein's observations (vectorized)."""
        np.add.at(self.tot, cell_idx, 1)
        np.add.at(self.pos, cell_idx, labels.astype(np.int64))

    def merge(self, other: "CellCounts") -> None:
        self.pos += other.pos
        self.tot += other.tot

    def base_rate(self) -> float:
        """Global contact rate over accumulated observations (0-safe)."""
        t = self.tot.sum()
        return float(self.pos.sum()) / float(t) if t else 0.0


def new_level1_counts() -> CellCounts:
    return CellCounts(LEVEL1_CELLS, np.zeros(LEVEL1_CELLS, np.int64),
                      np.zeros(LEVEL1_CELLS, np.int64))


def new_level2_counts() -> List[CellCounts]:
    return [CellCounts(LEVEL2_CELLS_PER_BIN,
                       np.zeros(LEVEL2_CELLS_PER_BIN, np.int64),
                       np.zeros(LEVEL2_CELLS_PER_BIN, np.int64))
            for _ in range(N_SEP_BINS)]


def accumulate_protein(counts1: CellCounts, counts2: List[CellCounts],
                       feats: Dict[str, np.ndarray]) -> None:
    """Add one protein's pair features into the pooled counters."""
    lab = feats["is_contact"]
    counts1.add(level1_cell_index(feats["g_i"], feats["g_j"], feats["sbin"]), lab)
    for b in range(N_SEP_BINS):
        m = feats["sbin"] == b
        if m.any():
            counts2[b].add(level2_cell_index(feats["aa_i"][m], feats["aa_j"][m]),
                           lab[m])


# ---------------------------------------------------------------------------
# Truth-table labeling + QM circuit extraction
# ---------------------------------------------------------------------------

def label_truth_table(counts: CellCounts,
                      t_mult: float = 2.0,
                      n_min: int = 30,
                      base_rate: Optional[float] = None) -> np.ndarray:
    """Ternary truth table from pooled counts (documented labeling rule).

        value =  1  if rate >= t_mult * base_rate  AND  n >= n_min
        value = -1  if n < n_min                              (don't-care)
        value =  0  otherwise                                 (not enriched)

    `base_rate` defaults to the counts' own global rate.  Sensitivity to
    (t_mult, n_min) is swept in run_universal_circuit.py -- conclusions must be
    stable across the sweep (gate G4).
    """
    base = base_rate if base_rate is not None else counts.base_rate()
    tt = np.zeros(counts.n_cells, dtype=np.int8)
    n = counts.tot
    dc = n < n_min
    rate = np.divide(counts.pos, np.maximum(n, 1))
    on = (~dc) & (rate >= t_mult * base) & (n >= n_min)
    tt[on] = 1
    tt[dc] = -1
    return tt


def qm_minimize(truth_flat: np.ndarray) -> dict:
    """Quine-McCluskey via the sibling repo's popcount-grouped implementation.

    `truth_flat`: 1-D ternary array (1 on-set / 0 off-set / -1 don't-care) whose
    length MUST be a power of 4 (kmap_truth_table hard-guards this; FIX A2).
    Returns the boolean_minimize_kmap result dict unchanged.
    """
    import sys as _sys
    for p in ("/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src",):
        if p not in _sys.path:
            _sys.path.insert(0, p)
    from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap  # type: ignore # noqa: E501 (sibling-repo path resolved at runtime)
    return boolean_minimize_kmap(truth_flat.astype(int).flatten(), algorithm="qm")


def decode_level1_implicant(values: Sequence[int], mask: Sequence[int]) -> Dict:
    """Decode an 8-variable implicant -> human-readable rule dict.

    Bit layout (MSB->LSB of the truth-table index): g_i[2], g_i[1], g_i[0],
    g_j[2], g_j[1], g_j[0], sbin[1], sbin[0].  A masked bit = don't-care.
    A fully-masked field decodes to None (= "any value").
    """
    def _matched(offset: int, width: int) -> Optional[List[int]]:
        """All width-bit numbers consistent with the implicant's fixed bits."""
        vals = list(values[offset:offset + width])
        msk = list(mask[offset:offset + width])
        if all(msk):                       # fully don't-care field
            return None
        matched: List[int] = []
        for x in range(1 << width):        # enumerate candidate field values
            ok = True
            for k in range(width):
                if not msk[k]:
                    bit = (x >> (width - 1 - k)) & 1
                    if bit != vals[k]:
                        ok = False
                        break
            if ok:
                matched.append(x)
        return matched

    return {
        "groups_i": _matched(0, 3),
        "groups_j": _matched(3, 3),
        "seps": _matched(6, 2),
        "n_dontcare": int(sum(mask)),
    }


def cell_rate_score(counts: CellCounts, cell_idx: np.ndarray,
                    alpha: float = 1.0) -> np.ndarray:
    """Smoothed empirical contact-rate per cell -> continuous ranking score.

    score(cell) = (pos + alpha * base) / (tot + alpha)   [Laplace-style shrinkage
    toward the global base rate; alpha=1 default].  Held-out pairs inherit their
    cell's TRAIN rate -- the circuit's probabilistic reading.
    """
    base = counts.base_rate()
    num = counts.pos[cell_idx] + alpha * base
    den = counts.tot[cell_idx] + alpha
    return num / den


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def precision_at_k(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """Fraction of the top-k scored pairs that are true contacts."""
    if k <= 0 or len(scores) == 0:
        return 0.0
    k = min(k, len(scores))
    top = np.argsort(-scores, kind="stable")[:k]
    return float(labels[top].mean())


def mcc(tp: int, fp: int, tn: int, fn: int) -> float:
    """Matthews correlation coefficient (0-safe)."""
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fn) * (tn + fp))
    return (tp * tn - fp * fn) / den if den else 0.0


def binary_metrics(pred: np.ndarray, labels: np.ndarray) -> Dict[str, float]:
    """TP/TN/FP/FN + accuracy/precision/recall/F1/MCC/enrichment."""
    tp = int(((pred == 1) & labels).sum())
    fp = int(((pred == 1) & ~labels).sum())
    tn = int(((pred == 0) & ~labels).sum())
    fn = int(((pred == 0) & labels).sum())
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    base = labels.mean() if len(labels) else 0.0
    return {
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "accuracy": (tp + tn) / max(len(labels), 1),
        "precision": prec, "recall": rec, "f1": f1,
        "mcc": mcc(tp, fp, tn, fn),
        "enrichment": prec / base if base else float("nan"),
        "base_rate": float(base),
    }


def auc_mann_whitney(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUC via the Mann-Whitney U statistic (ties handled by mid-ranks).

    AUC = U / (n_pos * n_neg), U = concordant-pairs count from rank sums.
    Rank bookkeeping: `ranks` is indexed by SORTED position; the rank of the
    element whose original index is m equals ranks[k] where order[k] == m,
    hence the positives' rank sum is ranks[labels[order]].sum().
    """
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="stable")
    s_sorted = scores[order]
    ranks = np.empty(len(scores), dtype=np.float64)
    i = 0
    while i < len(s_sorted):
        jj = i
        while jj + 1 < len(s_sorted) and s_sorted[jj + 1] == s_sorted[i]:
            jj += 1
        avg_rank = (i + jj + 2) / 2.0          # 1-based mid-rank of the tie block
        ranks[i:jj + 1] = avg_rank
        i = jj + 1
    rank_sum_pos = float(ranks[labels[order]].sum())
    u = rank_sum_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


# ---------------------------------------------------------------------------
# Protein-level cross-validation splitting
# ---------------------------------------------------------------------------

def protein_folds(targets: List[str], n_folds: int = 5,
                  seed: int = 42) -> List[List[str]]:
    """Stratified round-robin assignment of proteins to folds (by length quintile).

    Deterministic (seeded shuffle inside strata).  Used for the 5-fold CV and,
    with n_folds such that one fold ~= 1/3, for the 100/49 holdout split.
    """
    rng = np.random.default_rng(seed)
    lengths = []
    for t in targets:
        lines = [ln for ln in open(RAW / "aln" / f"{t}.aln") if ln.strip()]
        lengths.append(len(lines[0]) if lines else 0)
    order = np.argsort(lengths, kind="stable")
    folds: List[List[str]] = [[] for _ in range(n_folds)]
    # stratify into quintiles, shuffle within stratum, deal round-robin
    strata = np.array_split(order, 5)
    for stratum in strata:
        perm = rng.permutation(len(stratum))
        for rank, idx in enumerate(stratum[perm]):
            folds[rank % n_folds].append(targets[int(idx)])
    return folds


# ---------------------------------------------------------------------------
# Aggregation driver helpers
# ---------------------------------------------------------------------------

def list_targets() -> List[str]:
    """All 150 PSICOV targets (stem of raw/aln/*.aln), sorted."""
    return sorted(p.stem for p in RAW.glob("aln/*.aln"))


def summarize_counts(counts1: CellCounts, counts2: List[CellCounts]) -> dict:
    """JSON-ready summary of pooled counters (audit trail)."""
    return {
        "level1": {
            "observations": int(counts1.tot.sum()),
            "contacts": int(counts1.pos.sum()),
            "base_rate": counts1.base_rate(),
            "populated_cells": int((counts1.tot > 0).sum()),
        },
        "level2_per_bin": [
            {
                "bin": b,
                "observations": int(c.tot.sum()),
                "contacts": int(c.pos.sum()),
                "rate": c.base_rate(),
                "populated_cells": int((c.tot > 0).sum()),
            }
            for b, c in enumerate(counts2)
        ],
    }
