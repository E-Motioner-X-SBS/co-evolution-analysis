#!/usr/bin/env python3
"""
COVERAGE AUDIT: prove every script used ALL sequences and FULL alignment length.

For each completed run (results/<ds>/<target>/<script>/):
  - input FASTA stats: n_seqs, max length (full width)
  - reported stats from the run log tail (stdout) and summary JSONs:
      "N sequences" / "Total: N sequences" / "Positions: L" / "Variable positions: X / L"
  - FLAGS any run where reported sequences < input n_seqs, or reported
    positions < input width (i.e., truncation).

Usage: python scripts/verify_coverage.py [--dataset X]
Output: results/coverage_report.json + summary table
"""

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
RESULTS = REPO / "results"
DATA = REPO / "data"


def fasta_stats(fasta_path):
    n = 0
    width = 0
    cur = []
    for line in open(fasta_path):
        if line.startswith(">"):
            if cur:
                n += 1
                width = max(width, len("".join(cur)))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        n += 1
        width = max(width, len("".join(cur)))
    return n, width


SEQ_RE = re.compile(r"(\d+)[ ,]*sequences?", re.I)
POS_RE = re.compile(r"([\d,]+) positions?|Positions: (\d+)|(\d+) / (\d+)")
TOTAL_RE = re.compile(r"Total: (\d+) sequences")


def parse_tail(tail):
    """Extract reported n_seqs and n_positions from stdout tail (totals only).

    Handles: 'Total: N sequences', 'Sequences: N', 'over X of Y sequences',
    'Train: A sequences, Test: B sequences' (sums to total), and position
    reports. Sub-counts (cluster sizes) are ignored by taking the max.
    """
    n_seq = None
    n_pos = None
    m = TOTAL_RE.search(tail)
    if m:
        n_seq = int(m.group(1))
    m = re.search(r"over (\d+) of (\d+) sequences", tail)
    if m:
        n_seq = int(m.group(2))
    m = re.search(r"Train: (\d+) sequences, Test: (\d+) sequences", tail)
    if m:
        n_seq = int(m.group(1)) + int(m.group(2))
    for m in re.finditer(r"Sequences: (\d+)", tail):
        n_seq = int(m.group(1))
    # generic fallback: max of all 'N sequences' mentions (totals dominate)
    if n_seq is None:
        vals = [int(x) for x in re.findall(r"(\d+) sequences?", tail)]
        if vals:
            n_seq = max(vals)
    for m in re.finditer(r"Positions: (\d+)", tail):
        n_pos = int(m.group(1))
    for m in re.finditer(r"positions 0-(\d+)", tail):
        n_pos = int(m.group(1)) + 1
    for m in re.finditer(r"(\d+) positions", tail):
        n_pos = int(m.group(1))
    return n_seq, n_pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset")
    args = ap.parse_args()

    rows = []
    for run_log in sorted(RESULTS.glob("*/*/*/run_log.json")):
        ds, tgt, script_dir = run_log.parts[-4], run_log.parts[-3], run_log.parts[-2]
        if args.dataset and ds != args.dataset:
            continue
        d = json.load(open(run_log))
        if d["status"] != "ok":
            rows.append(
                {
                    "dataset": ds,
                    "target": tgt,
                    "script": script_dir,
                    "status": d["status"],
                    "note": "not ok — skipped",
                }
            )
            continue
        fasta = DATA / ds / "converted" / f"{tgt}.fasta"
        if not fasta.exists():
            continue
        n_in, w_in = fasta_stats(str(fasta))
        n_rep, p_rep = parse_tail(d.get("tail", ""))
        flags = []
        if n_rep is not None and n_rep < n_in:
            flags.append(f"SEQS {n_rep}<{n_in}")
        if p_rep is not None and w_in > 0 and p_rep < w_in:
            flags.append(f"LEN {p_rep}<{w_in}")
        rows.append(
            {
                "dataset": ds,
                "target": tgt,
                "script": script_dir,
                "status": d["status"],
                "n_input": n_in,
                "width_input": w_in,
                "n_reported": n_rep,
                "len_reported": p_rep,
                "flags": flags,
            }
        )

    out = RESULTS / "coverage_report.json"
    json.dump(rows, open(out, "w"), indent=1)
    flagged = [r for r in rows if r.get("flags")]
    n_checked = len(rows)
    print(f"coverage audit: {n_checked} runs, {len(flagged)} flagged")
    for r in flagged:
        print(f"  FLAG {r['dataset']}/{r['target']}/{r['script']}: {r['flags']}")
    if not flagged:
        print("ALL RUNS: no sequence or length truncation detected.")
    return 0 if not flagged else 1


if __name__ == "__main__":
    sys.exit(main())
