#!/usr/bin/env python3
"""
PHASE 2: Parameterize the 23 analysis scripts for multi-dataset runs.

Adds two environment variables with UNCHANGED defaults:
  COEVO_FASTA   -> FASTA/MSA input  (default: datasets/co-evolution/Spike_protein.aln-fasta)
  COEVO_RESULTS -> output dir       (default: each script's original results dir)

Every replacement is counted and asserted; the script fails loudly if a pattern
did not match. Reversible via `git checkout -- <files>` (repo is git-tracked).
"""

import re
import subprocess
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
PY = "from __future__ import annotations\n"  # not needed; placeholder
OS_ENV = '__import__("os").environ.get'


def patch_file(rel, rules, label):
    p = REPO / rel
    src = p.read_text()
    if "COEVO_FASTA" in src or "COEVO_RESULTS" in src:
        print(f"  [skip] {rel}: already parameterized")
        return
    total = 0
    for pattern, repl, count in rules:
        new_src, n = re.subn(pattern, repl, src, count=count, flags=re.M)
        if n == 0:
            print(f"  !! {rel}: NO MATCH for {pattern[:70]}")
            sys.exit(1)
        src = new_src
        total += n
    p.write_text(src)
    print(f"  [OK] {rel}: {total} replacements ({label})")


print("=== A. shared loader ===")
patch_file(
    "coevolution_shared.py",
    [
        (
            r'fasta_path = Path\(__file__\)\.resolve\(\)\.parent / "Spike_protein\.aln-fasta"',
            'fasta_path = os.environ.get("COEVO_FASTA") or (Path(__file__).resolve().parent / "Spike_protein.aln-fasta")',
            1,
        ),
    ],
    "loader",
)

print("=== B. uniform 20-script patch ===")
UNIFORM = [
    "advanced_co-evolution_analysis.py",
    "allseq_constraint_function.py",
    "boolean_co-evolution.py",
    "dca_boolean_coevolution.py",
    "flipped_boolean_coevolution.py",
    "full_length_analysis.py",
    "kmap_boolean_coevolution.py",
    "master_boolean.py",
    "nary_kmap_co-evolution.py",
    "perplexity_coevolution.py",
    "position_kmap_coevolution.py",
    "predictive_constraint_function.py",
    "run_allseq_analysis.py",
    "run_kmap_analysis.py",
    "variable_position_coevolution.py",
]
for f in UNIFORM:
    patch_file(
        f,
        [
            (
                r'^(\s*)fasta_file = base_dir / "Spike_protein\.aln-fasta"',
                r"\1fasta_file = Path("
                + OS_ENV
                + r'("COEVO_FASTA") or (base_dir / "Spike_protein.aln-fasta"))',
                1,
            ),
            (
                r'^(\s*)results_dir = base_dir / "([^"]+)"',
                r"\1results_dir = Path("
                + OS_ENV
                + r'("COEVO_RESULTS") or (base_dir / "\2"))',
                1,
            ),
            (
                r"^(\s*)results_dir\.mkdir\(exist_ok=True\)",
                r"\1results_dir.mkdir(parents=True, exist_ok=True)",
                1,
            ),
        ],
        "uniform",
    )

print("=== C. create_mi_heatmap ===")
patch_file(
    "create_mi_heatmap.py",
    [
        (
            r'^(\s*)fasta_file = base_dir / "Spike_protein\.aln-fasta"',
            r"\1fasta_file = Path("
            + OS_ENV
            + r'("COEVO_FASTA") or (base_dir / "Spike_protein.aln-fasta"))',
            1,
        ),
        (
            r'output_dir = base_dir / "mi_heatmap"',
            "output_dir = Path("
            + OS_ENV
            + r'("COEVO_RESULTS") or (base_dir / "mi_heatmap"))',
            1,
        ),
        (
            r"output_dir\.mkdir\(exist_ok=True\)",
            "output_dir.mkdir(parents=True, exist_ok=True)",
            1,
        ),
    ],
    "heatmap",
)

print("=== D. gpu_full_analysis ===")
patch_file(
    "gpu_full_analysis.py",
    [
        (
            r'FASTA = BASE / "Spike_protein\.aln-fasta"',
            "FASTA = Path("
            + OS_ENV
            + r'("COEVO_FASTA") or (BASE / "Spike_protein.aln-fasta"))',
            1,
        ),
        (
            r'OUT = BASE / "full_gpu_results"',
            "OUT = Path("
            + OS_ENV
            + r'("COEVO_RESULTS") or (BASE / "full_gpu_results"))',
            1,
        ),
        (r"OUT\.mkdir\(exist_ok=True\)", "OUT.mkdir(parents=True, exist_ok=True)", 1),
    ],
    "gpu_full",
)

print("=== E. dca_mf_analysis ===")
patch_file(
    "dca_mf_analysis.py",
    [
        (
            r'FASTA = BASE / "Spike_protein\.aln-fasta"',
            "FASTA = Path("
            + OS_ENV
            + r'("COEVO_FASTA") or (BASE / "Spike_protein.aln-fasta"))',
            1,
        ),
        (
            r'OUT = BASE / "dca_results"',
            "OUT = Path(" + OS_ENV + r'("COEVO_RESULTS") or (BASE / "dca_results"))',
            1,
        ),
        (r"OUT\.mkdir\(exist_ok=True\)", "OUT.mkdir(parents=True, exist_ok=True)", 1),
    ],
    "dca_mf",
)

print("=== F. generate_co-evolution_md ===")
patch_file(
    "generate_co-evolution_md.py",
    [
        (
            r'sequences = parse_fasta\(base_dir / "Spike_protein\.aln-fasta"\)',
            "sequences = parse_fasta("
            + OS_ENV
            + r'("COEVO_FASTA") or (base_dir / "Spike_protein.aln-fasta"))',
            1,
        ),
        (
            r'open\(base_dir / "master_boolean" / "master_boolean_summary\.json"\)',
            "open("
            + OS_ENV
            + r'("COEVO_RESULTS") or (base_dir / "master_boolean") / "master_boolean_summary.json")',
            1,
        ),
        (
            r'output_path = base_dir / "kmap_boolean_coevolution" / "COEVOLUTION_KMAP_BOOLEAN\.md"',
            "output_path = ("
            + OS_ENV
            + r'("COEVO_RESULTS") or (base_dir / "kmap_boolean_coevolution")) / "COEVOLUTION_KMAP_BOOLEAN.md"',
            1,
        ),
    ],
    "coev_md",
)

print("=== G. generate_full_analysis_md ===")
patch_file(
    "generate_full_analysis_md.py",
    [
        (
            r'OUT = BASE / "FULL_COEVOLUTION_ANALYSIS\.md"',
            "RESULTS = Path("
            + OS_ENV
            + r'("COEVO_RESULTS") or BASE)\nOUT = Path('
            + OS_ENV
            + r'("COEVO_OUT") or (BASE / "FULL_COEVOLUTION_ANALYSIS.md"))',
            1,
        ),
        (r"load_json\(BASE / ", "load_json(RESULTS / ", 0),  # all reads
    ],
    "full_analysis_md",
)

print("=== H. generate_full_pipeline_doc ===")
patch_file(
    "generate_full_pipeline_doc.py",
    [
        (
            r'OUT = BASE / "FULL_PIPELINE_ANALYSIS\.md"',
            "RESULTS = Path("
            + OS_ENV
            + r'("COEVO_RESULTS") or BASE)\nOUT = Path('
            + OS_ENV
            + r'("COEVO_OUT") or (BASE / "FULL_PIPELINE_ANALYSIS.md"))',
            1,
        ),
        (r"open\(BASE / ", "open(RESULTS / ", 0),
    ],
    "pipeline_doc",
)

# final sanity: compile all
print("=== compile check ===")
files = (
    ["coevolution_shared.py"]
    + UNIFORM
    + [
        "create_mi_heatmap.py",
        "gpu_full_analysis.py",
        "dca_mf_analysis.py",
        "generate_co-evolution_md.py",
        "generate_full_analysis_md.py",
        "generate_full_pipeline_doc.py",
    ]
)
ok = True
for f in files:
    r = subprocess.run(
        [sys.executable, "-m", "py_compile", str(REPO / f)], capture_output=True
    )
    if r.returncode != 0:
        ok = False
        print(f"  COMPILE FAIL {f}: {r.stderr.decode()[:300]}")
print("ALL COMPILE OK" if ok else "COMPILE ERRORS")
sys.exit(0 if ok else 1)
