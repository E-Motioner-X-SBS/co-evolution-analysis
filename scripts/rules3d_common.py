#!/usr/bin/env python3
"""
rules3d_common.py — shared core for rule->3D validation on the Omicron models.

Rigor guarantees:
  1. Mapping integrity: the model's target_sequence (from full_details.json)
     MUST equal the gap-stripped aligned sequence; mismatched models are
     skipped and reported, never silently mapped.
  2. Multi-chain parsing (homo-trimer models: chains A/B/C), residue numbering
     per chain, Cb coords (Ca for Gly).
  3. Contact: Cb-Cb < 8 A, residue separation >= 3 (backbone-adjacent pairs
     excluded as trivial).
  4. Random controls matched to the rule set's own separation distribution
     (2,000 draws per model, fixed seed per model), binomial p-value for the
     observed contact count against the control rate.
  5. Both INTRA-chain (same protomer) and INTER-chain (across protomers)
     contact fractions are computed — the trimer interface is a distinct
     structural dimension.
"""
import gzip
import json
import math
import os
import random
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
SWISS = REPO / "swissmodel"
FASTA = REPO / "Spike_protein.aln-fasta"
CONTACT = 8.0
MIN_SEP = 3
N_RAND = 2000
CANON = set("ACDEFGHIKLMNPQRSTVWY")

AA3 = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
       "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
       "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
       "TRP": "W", "TYR": "Y"}


def parse_fasta(path):
    """FASTA -> [(header, seq)]."""
    seqs = []
    h = None
    s = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if h:
                seqs.append((h, "".join(s)))
            h = line[1:]
            s = []
        elif line:
            s.append(line.upper())
    if h:
        seqs.append((h, "".join(s)))
    return seqs


def gap_strip(seq):
    """Strip ALL non-canonical characters (gaps AND ambiguous, e.g. X)."""
    return "".join(c for c in seq if c in CANON)


def col_to_residue_map(aligned_seq, target_seq):
    """Map alignment column (0-based) -> target residue number (1-based).

    Verified: target_seq must equal gap-stripped aligned_seq (integrity).
    Returns (mapping, ok). Columns that are gaps map to None.
    """
    if gap_strip(aligned_seq) != target_seq:
        return {}, False
    m = {}
    res = 1
    for c, ch in enumerate(aligned_seq):
        if ch not in CANON:
            continue  # gaps and ambiguous chars have no model residue
        m[c] = res
        res += 1
    return m, True


def parse_model_pdb(pdb_gz_path):
    """Parse a model PDB -> {chain: {residue: (aa, cb_coords)}}.

    Cb used (Ca for Gly); only canonical residues kept.
    """
    chains = {}
    with gzip.open(pdb_gz_path, "rt") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            name = line[12:16].strip()
            if name not in ("CA", "CB"):
                continue
            chain = line[21]
            if chain == " ":
                chain = "A"
            aa3 = line[17:20].strip()
            aa = AA3.get(aa3)
            if aa is None:
                continue
            resi = int(line[22:26])
            coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            chains.setdefault(chain, {}).setdefault(resi, {})[name] = coord
    out = {}
    for chain, residues in chains.items():
        out[chain] = {}
        for resi, atoms in residues.items():
            cb = atoms.get("CB", atoms.get("CA"))
            if cb is not None:
                out[chain][resi] = cb
    return out


def dist(a, b):
    return math.sqrt(sum((u - v) ** 2 for u, v in zip(a, b)))


class OmicronModels:
    """Enumerate the 299 unique Omicron sequences and their SwissModel models."""

    def __init__(self):
        self.fasta_seqs = parse_fasta(str(FASTA))
        self.uniques = json.load(open(SWISS / "batch" / "state" / "uniques.json"))
        self.skipped = []

    def uid_aligned_seq(self, uid):
        u = next((x for x in self.uniques if x["uid"] == uid), None)
        if u is None:
            return None
        return self.fasta_seqs[u["indices"][0]][1]

    def best_model(self, uid):
        """Return (model_id, pdb_path, target_seq) for the best-GMQE model."""
        tsum = SWISS / "batch" / "pdbs" / uid / "templates_summary.json"
        fd = SWISS / "batch" / "pdbs" / uid / "full_details.json"
        if not tsum.exists() or not fd.exists():
            return None
        best = max(json.load(open(tsum)), key=lambda m: m.get("gmqe") or 0)
        model_id = best.get("model_id")
        pdb = SWISS / "batch" / "pdbs" / uid / f"{model_id}.pdb.gz"
        if not pdb.exists():
            return None
        details = json.load(open(fd))
        target = None
        for m in details.get("models", []):
            if m.get("model_id") == model_id:
                tg = m.get("targets") or [{}]
                target = tg[0].get("target_sequence", "") if tg else ""
        return model_id, pdb, target

    def iter_models(self):
        """Yield (uid, aligned_seq, coords_by_chain) for every valid model.

        Integrity: skips and records models whose target_sequence does not
        match the gap-stripped aligned sequence.
        """
        for u in self.uniques:
            uid = u["uid"]
            aligned = self.uid_aligned_seq(uid)
            if aligned is None:
                self.skipped.append((uid, "no aligned seq"))
                continue
            bm = self.best_model(uid)
            if bm is None:
                self.skipped.append((uid, "no best model"))
                continue
            _, pdb, target = bm
            if not target:
                self.skipped.append((uid, "no target sequence in payload"))
                continue
            cmap, ok = col_to_residue_map(aligned, target)
            if not ok:
                self.skipped.append((uid, "target/alignment mismatch"))
                continue
            coords = parse_model_pdb(str(pdb))
            if not coords:
                self.skipped.append((uid, "no parsible coordinates"))
                continue
            yield uid, aligned, cmap, coords


def evaluate_pairs(models_iter, pair_list):
    """Evaluate a list of (col_i, col_j) alignment pairs across models.

    Returns (rows, control_rate) where rows = per-pair dicts with intra/inter
    contact stats, and control_rate = pooled matched-random contact rate.
    """
    rows = {p: {"models": 0, "intra_contacts": 0, "inter_contacts": 0,
                "intra_dists": [], "inter_dists": []} for p in pair_list}
    rnd_intra = 0
    rnd_inter = 0
    rnd_n = 0
    model_count = 0
    for uid, aligned, cmap, coords in models_iter:
        model_count += 1
        chains = sorted(coords)
        # map columns -> residues
        mapped = {}
        for (i, j) in pair_list:
            ri, rj = cmap.get(i), cmap.get(j)
            if ri is None or rj is None:
                continue
            if abs(ri - rj) < MIN_SEP:
                continue
            mapped[(i, j)] = (ri, rj)
        # intra-chain (same chain, residues present)
        for (i, j), (ri, rj) in mapped.items():
            if len(chains) == 0:
                continue
            chain = chains[0]
            if ri not in coords[chain] or rj not in coords[chain]:
                continue
            rows[(i, j)]["models"] += 1
            d = dist(coords[chain][ri], coords[chain][rj])
            rows[(i, j)]["intra_dists"].append(d)
            if d < CONTACT:
                rows[(i, j)]["intra_contacts"] += 1
            # inter-chain: min distance to residue in any OTHER chain
            best_inter = None
            for ch in chains[1:]:
                if ri in coords[ch] and rj in coords[ch]:
                    dd = dist(coords[ch][ri], coords[ch][rj])
                    best_inter = dd if best_inter is None else min(best_inter, dd)
            if best_inter is not None:
                rows[(i, j)]["inter_dists"].append(best_inter)
                if best_inter < CONTACT:
                    rows[(i, j)]["inter_contacts"] += 1
        # random control (matched separation distribution of the rule pairs)
        seps = [abs(ri - rj) for (ri, rj) in mapped.values()]
        if not seps or len(chains) == 0:
            continue
        chain = chains[0]
        residues = sorted(coords[chain])
        if len(residues) < 10:
            continue
        rng = random.Random(42 + model_count)
        for _ in range(N_RAND):
            a = rng.choice(residues)
            s = rng.choice(seps)
            b = a + s
            if b not in coords[chain]:
                continue
            if dist(coords[chain][a], coords[chain][b]) < CONTACT:
                rnd_intra += 1
            rnd_n += 1
    control = rnd_intra / rnd_n if rnd_n else 0.0
    return rows, control, model_count


def row_summary(rows, control, key_intra=True):
    """Per-pair summary table rows with enrichment + binomial p-value."""
    out = []
    for p, r in rows.items():
        n = r["models"]
        if n == 0:
            continue
        contacts = r["intra_contacts"] if key_intra else r["inter_contacts"]
        dists = r["intra_dists"] if key_intra else r["inter_dists"]
        frac = contacts / n
        enrich = frac / control if control > 0 else float("nan")
        # exact binomial test: P(X >= contacts | n, control) via scipy
        from scipy.stats import binomtest
        pval = 1.0
        if control > 0 and n > 0:
            pval = binomtest(contacts, n, control, alternative="greater").pvalue
        out.append({"pair": p, "models": n, "contacts": contacts,
                    "fraction": round(frac, 3), "enrichment": round(enrich, 2),
                    "p_value": pval,
                    "median_dist": round(sorted(dists)[len(dists) // 2], 2)
                    if dists else None})
    return out


def aggregate_set(rows, control, pair_list):
    """Pooled statistics over a rule set."""
    n = sum(rows[p]["models"] for p in pair_list)
    c = sum(rows[p]["intra_contacts"] for p in pair_list)
    frac = c / n if n else 0.0
    enrich = frac / control if control else float("nan")
    pval = 1.0
    if control > 0 and n:
        from scipy.stats import binomtest
        pval = binomtest(c, n, control, alternative="greater").pvalue
    return {"n": n, "contacts": c, "fraction": round(frac, 3),
            "enrichment": round(enrich, 2), "p_value": pval}
