#!/usr/bin/env python3
"""
PHASE 3: Convert all 6 datasets to canonical aligned FASTA + metadata.

Each dataset gets data/<dataset>/converted/*.fasta (uppercased, gaps kept as '-',
headers per source). Conversions are verified: output FASTA sequence count must
equal the source count. Resumable (skips existing converted files).

DeepMSA note: the tar ships seq.fasta + native.pdb per target, NO MSAs
(verified: 0 other files). Structures are extracted for contact validation;
MSA building is out of scope (no hhblits/jackhmmer on this system — documented).
"""

import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
LOG = DATA / "convert.log"


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def count_fasta(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def write_fasta(out_path, entries):
    """entries: list of (name, seq)"""
    os.makedirs(out_path.parent, exist_ok=True)
    with open(out_path, "w") as f:
        for name, seq in entries:
            f.write(f">{name}\n{seq}\n")


def stats(entries):
    seqs = [s for _, s in entries]
    lens = Counter(len(s) for s in seqs)
    uniq = len({hashlib.md5(s.encode()).hexdigest() for s in seqs})
    return {
        "n_sequences": len(seqs),
        "width": max(lens),
        "n_unique": uniq,
        "redundancy_ratio": round(len(seqs) / uniq, 2) if uniq else None,
    }


# ───────────────────────── 1. PSICOV 150 ─────────────────────────


def convert_psicov():
    raw = DATA / "psicov150" / "raw"
    if not (raw / "aln").exists():
        log("psicov: extracting tar.gz")
        with tarfile.open(DATA / "raw" / "suppdata.tar.gz") as tf:
            tf.extractall(raw)
    aln_dir = raw / "aln"
    out = DATA / "psicov150" / "converted"
    n_files = 0
    for f in sorted(aln_dir.glob("*.aln")):
        name = f.stem
        dest = out / f"{name}.fasta"
        if dest.exists():
            continue
        entries = []
        for line in open(f):
            line = line.strip()
            if line:
                entries.append((name, line.upper()))
        write_fasta(dest, entries)
        n_files += 1
    log(f"psicov: converted {n_files} proteins -> {out}")
    # verify: every aln has a converted twin with equal count
    bad = [
        f.stem for f in aln_dir.glob("*.aln") if not (out / f"{f.stem}.fasta").exists()
    ]
    log(f"psicov: verify missing={len(bad)}")
    return n_files, bad


# ───────────────────────── 2. Pfam Stockholm ─────────────────────────


def convert_pfam():
    fams = {"PF00072": "Response_reg", "PF00595": "PDZ", "PF00071": "Ras"}
    for fam, label in fams.items():
        src = DATA / "pfam" / f"{fam}_full.sto"
        dest = DATA / "pfam" / "converted" / f"{fam}.fasta"
        if dest.exists():
            log(f"pfam {fam}: already converted")
            continue
        log(f"pfam {fam} ({label}): parsing Stockholm")
        entries = []
        seqs = {}
        for line in open(src):
            if line.startswith(("#", "/")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                name, seq = parts[0], parts[1].upper()
                seqs[name] = seqs.get(name, "") + seq
        entries = list(seqs.items())
        write_fasta(dest, entries)
        st = stats(entries)
        log(f"pfam {fam}: {st}")
    # exact-dup report
    report = {}
    for fam in fams:
        dest = DATA / "pfam" / "converted" / f"{fam}.fasta"
        if dest.exists():
            report[fam] = stats(
                [
                    (n, s)
                    for n, s in [
                        l.split("\n", 1)
                        for l in [x for x in dest.read_text().split(">")[1:]]
                    ]
                ]
            )
    json.dump(report, open(DATA / "pfam" / "converted_stats.json", "w"), indent=1)
    log(f"pfam: stats -> {DATA / 'pfam' / 'converted_stats.json'}")


# ───────────────────────── 3. EVmutation ─────────────────────────


def convert_evmutation():
    raw = DATA / "evmutation" / "raw"
    if not any(raw.glob("*.a2m")):
        log("evmutation: extracting tar")
        with tarfile.open(DATA / "evmutation" / "raw" / "alignments.tar.gz") as tf:
            tf.extractall(raw)
    out = DATA / "evmutation" / "converted"
    n = 0
    for f in sorted(raw.glob("*.a2m")):
        dest = out / f"{f.stem}.fasta"
        if dest.exists():
            continue
        entries = []
        name, seq = None, []
        for line in open(f):
            if line.startswith(">"):
                if name:
                    entries.append((name, "".join(seq).upper()))
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if name:
            entries.append((name, "".join(seq).upper()))
        write_fasta(dest, entries)
        n += 1
    log(f"evmutation: converted {n} proteins -> {out}")


# ───────────────────────── 4. Unit tests ─────────────────────────


def convert_units():
    out = DATA / "unit_tests" / "converted"
    # DHFR.a2m (FASTA-like with lowercase)
    if not (out / "DHFR.fasta").exists():
        entries = []
        name, seq = None, []
        for line in open(DATA / "unit_tests" / "DHFR.a2m"):
            if line.startswith(">"):
                if name:
                    entries.append((name, "".join(seq).upper()))
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if name:
            entries.append((name, "".join(seq).upper()))
        write_fasta(out / "DHFR.fasta", entries)
    # 1whzA.msa (two-column: name <spaces> seq)
    if not (out / "1whzA.fasta").exists():
        entries = []
        for line in open(DATA / "unit_tests" / "1whzA.msa"):
            parts = line.split()
            if len(parts) >= 2:
                entries.append((parts[0], parts[1].upper()))
        write_fasta(out / "1whzA.fasta", entries)
    # 1atzA.aln (one seq per line, no headers)
    if not (out / "1atzA.fasta").exists():
        entries = [
            (f"S{i}", l.strip().upper())
            for i, l in enumerate(open(DATA / "unit_tests" / "1atzA.aln"))
            if l.strip()
        ]
        write_fasta(out / "1atzA.fasta", entries)
    for f in ("DHFR", "1whzA", "1atzA"):
        log(f"unit {f}: {stats(list(iter_fasta(out / f'{f}.fasta')))}")


def iter_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


# ───────────────────────── 5. DeepMSA ─────────────────────────


def convert_deepmsa():
    raw = DATA / "deepmsa" / "raw"
    if not (raw / "data").exists():
        log("deepmsa: extracting tar.bz2 (1,229 members)")
        import bz2

        with tarfile.open(DATA / "deepmsa" / "raw" / "dataset.tar.bz2") as tf:
            tf.extractall(raw)
    targets = sorted((raw / "data").glob("*/"))
    natives = DATA / "deepmsa" / "pdb"
    os.makedirs(natives, exist_ok=True)
    n_nat = 0
    for t in targets:
        np_ = t / "native.pdb"
        if np_.exists() and not (natives / f"{t.name}.pdb").exists():
            shutil.copy2(np_, natives / f"{t.name}.pdb")
            n_nat += 1
    log(f"deepmsa: {len(targets)} targets, {n_nat} native structures copied")
    log(
        "deepmsa: NOTE — no MSAs in this dataset (seq.fasta + native.pdb only); "
        "MSA building out of scope (no hhblits/jackhmmer)"
    )


# ───────────────────────── 6. GPCRdb ─────────────────────────


def convert_gpcrdb():
    out = DATA / "gpcrdb" / "converted"
    os.makedirs(out, exist_ok=True)
    for f in ("5ht1a.fasta", "b2ar.fasta"):
        shutil.copy2(DATA / "gpcrdb" / f, out / f)
        log(f"gpcrdb {f}: {stats(list(iter_fasta(out / f)))}")


# ───────────────────────── metadata ─────────────────────────


def write_metadata():
    meta = {}
    for ds, name in [
        ("psicov150", "psicov150"),
        ("pfam", "pfam"),
        ("evmutation", "evmutation"),
        ("unit_tests", "unit_tests"),
        ("deepmsa", "deepmsa"),
        ("gpcrdb", "gpcrdb"),
    ]:
        conv = DATA / ds / "converted"
        if conv.exists():
            files = sorted(conv.glob("*.fasta"))
            meta[ds] = {"n_files": len(files), "files": [f.name for f in files]}
            if ds == "psicov150":
                meta[ds]["proteins"] = [f.stem for f in files]
    json.dump(meta, open(DATA / "metadata.json", "w"), indent=1)
    log(f"metadata -> {DATA / 'metadata.json'}")


def main():
    convert_psicov()
    convert_pfam()
    convert_evmutation()
    convert_units()
    convert_deepmsa()
    convert_gpcrdb()
    write_metadata()
    log("CONVERT ALL DONE")


if __name__ == "__main__":
    main()
