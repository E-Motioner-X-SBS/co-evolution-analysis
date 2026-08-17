#!/usr/bin/env python3
"""
PHASE 4: Run the 23 analysis scripts on every dataset.

Layout: results/<dataset>/<target>/<script>/  (script outputs + run_log.json)

Usage:
  python scripts/run_dataset.py --dataset unit_tests --target DHFR
  python scripts/run_dataset.py --all --parallel 4
  python scripts/run_dataset.py --psicov-mode core   # 150 proteins, core scripts
  python scripts/run_dataset.py --psicov-mode full   # 150 proteins, all scripts (slow)

Every run: env COEVO_FASTA + COEVO_RESULTS, cwd = results dir, timeout, exit code,
stdout/stderr tails captured into run_log.json. Resumable (skips done runs).
"""

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
PY = "/store/shuvam/.venv/bin/python"
RESULTS = REPO / "results"
TIMEOUT_S = None  # user directive: NO timeouts — let every job run to completion


# scripts that need >30s sometimes (mfDCA, full boolean, LOO-CV)
SLOW = {
    "dca_mf_analysis.py",
    "boolean_co-evolution.py",
    "allseq_constraint_function.py",
}

CORE_SCRIPTS = [
    "full_length_analysis.py",
    "master_boolean.py",
    "kmap_boolean_coevolution.py",
    "flipped_boolean_coevolution.py",
    "nary_kmap_co-evolution.py",
    "perplexity_coevolution.py",
    "run_kmap_analysis.py",
    "variable_position_coevolution.py",
    "position_kmap_coevolution.py",
    "run_allseq_analysis.py",
    "advanced_co-evolution_analysis.py",
]
FULL_SCRIPTS = CORE_SCRIPTS + [
    "boolean_co-evolution.py",
    "create_mi_heatmap.py",
    "gpu_full_analysis.py",
    "dca_boolean_coevolution.py",
    "dca_mf_analysis.py",
    "allseq_constraint_function.py",
    "predictive_constraint_function.py",
    "generate_co-evolution_md.py",
    "generate_full_analysis_md.py",
    "generate_full_pipeline_doc.py",
]

# Canonical results-dir names (report generators read these fixed names)
SCRIPT_DIRNAME = {
    "master_boolean.py": "master_boolean",
    "boolean_co-evolution.py": "boolean_results",
    "nary_kmap_co-evolution.py": "nary_kmap_results",
    "position_kmap_coevolution.py": "position_kmap_results",
    "run_allseq_analysis.py": "full_position_results",
    "run_kmap_analysis.py": "kmap_results",
    "flipped_boolean_coevolution.py": "flipped_boolean_results",
    "kmap_boolean_coevolution.py": "kmap_boolean_coevolution",
    "variable_position_coevolution.py": "variable_position_results",
    "predictive_constraint_function.py": "constraint_function_results",
    "allseq_constraint_function.py": "allseq_constraint_results",
    "dca_boolean_coevolution.py": "dca_boolean_results",
    "dca_mf_analysis.py": "dca_results",
    "perplexity_coevolution.py": "perplexity_results",
    "advanced_co-evolution_analysis.py": "advanced_analysis_results",
    "full_length_analysis.py": "full_length_results",
    "gpu_full_analysis.py": "full_gpu_results",
    "create_mi_heatmap.py": "mi_heatmap",
}


# dataset -> list of (target, fasta)
def dataset_targets():
    D = {}
    for f in ["DHFR", "1whzA", "1atzA"]:
        p = REPO / "data" / "unit_tests" / "converted" / f"{f}.fasta"
        if p.exists():
            D.setdefault("unit_tests", []).append((f, p))
    for f in ["5ht1a", "b2ar"]:
        p = REPO / "data" / "gpcrdb" / "converted" / f"{f}.fasta"
        if p.exists():
            D.setdefault("gpcrdb", []).append((f, p))
    ev = REPO / "data" / "evmutation" / "converted"
    if ev.exists():
        D["evmutation"] = [(p.stem, p) for p in sorted(ev.glob("*.fasta"))]
    pf = REPO / "data" / "pfam" / "converted"
    if pf.exists():
        D["pfam"] = [(p.stem, p) for p in sorted(pf.glob("*.fasta"))]
    ps = REPO / "data" / "psicov150" / "converted"
    if ps.exists():
        D["psicov150"] = [(p.stem, p) for p in sorted(ps.glob("*.fasta"))]
    return D


REPORT_GENERATORS = {
    "generate_co-evolution_md.py",
    "generate_full_analysis_md.py",
    "generate_full_pipeline_doc.py",
}

# GPU-heavy scripts: constrained to 2 concurrent (one runner at a time).
# Light scripts may touch CUDA transiently (memory-fraction fail-fast protects).
GPU_SCRIPTS = {
    "master_boolean.py", "gpu_full_analysis.py", "create_mi_heatmap.py",
    "dca_mf_analysis.py", "boolean_co-evolution.py", "full_length_analysis.py",
    "run_allseq_analysis.py", "dca_boolean_coevolution.py",
}
_GPU_SEM = __import__("threading").BoundedSemaphore(8)  # adaptive chunks bound memory; 8 is safe


def run_one(job):
    dataset, target, script, fasta, results_dir = job
    if script in GPU_SCRIPTS:
        _GPU_SEM.acquire()
    try:
        return _run_one(job)
    finally:
        if script in GPU_SCRIPTS:
            _GPU_SEM.release()


def _run_one(job):
    dataset, target, script, fasta, results_dir = job
    dirname = SCRIPT_DIRNAME.get(script, script.replace(".py", ""))
    outdir = results_dir / dataset / target / dirname
    outdir.mkdir(parents=True, exist_ok=True)
    logf = outdir / "run_log.json"
    if logf.exists():
        return (dataset, target, script, "skipped")
    env = dict(os.environ)
    env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    # CRITICAL: torch's OpenMP pins the process to CPU 0 unless OMP_PROC_BIND=FALSE
    # (verified: import torch -> affinity {0}; with OMP_PROC_BIND=FALSE -> all 24).
    env["OMP_PROC_BIND"] = "FALSE"
    # LOO-CV scripts are CPU-bound algorithms; keep them OFF CUDA entirely
    # (their initial GPU MI causes allocator spins on shared GPUs).
    if script in ("allseq_constraint_function.py", "predictive_constraint_function.py"):
        env["CUDA_VISIBLE_DEVICES"] = ""
    elif os.environ.get("COEVO_FORCE_CPU") and script in GPU_SCRIPTS:
        env["CUDA_VISIBLE_DEVICES"] = ""  # CPU fallback for memory-heavy scripts
    env["COEVO_FASTA"] = str(fasta)
    if script in REPORT_GENERATORS:
        # generators read sibling result dirs at the TARGET level
        env["COEVO_RESULTS"] = str(results_dir / dataset / target)
        env["COEVO_OUT"] = str(outdir / script.replace(".py", ".md"))
    else:
        env["COEVO_RESULTS"] = str(outdir)
    t0 = time.time()
    try:
        r = subprocess.run(
            [PY, str(REPO / script)],
            cwd=str(outdir),
            env=env,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_S,
        )
        status = "ok" if r.returncode == 0 else f"exit{r.returncode}"
        tail = (r.stdout + "\n" + r.stderr)[-2000:]
    except subprocess.TimeoutExpired:
        status = "timeout"
        tail = f"timeout >{TIMEOUT_S}s"
    except Exception as e:
        status = "exc"
        tail = str(e)[:2000]
    json.dump(
        {
            "dataset": dataset,
            "target": target,
            "script": script,
            "status": status,
            "seconds": round(time.time() - t0, 1),
            "tail": tail,
            "outputs": [f.name for f in outdir.iterdir() if f.name != "run_log.json"],
        },
        open(logf, "w"),
        indent=1,
    )
    return (dataset, target, script, status)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", help="limit to one dataset")
    ap.add_argument("--target", help="limit to one target")
    ap.add_argument("--psicov-mode", choices=["core", "full", "none"], default="core")
    ap.add_argument("--parallel", type=int, default=4)
    ap.add_argument("--max-targets", type=int, default=0, help="limit targets (test)")
    ap.add_argument(
        "--retry-failed",
        action="store_true",
        help="quarantine failed-run dirs so they re-run",
    )
    ap.add_argument(
        "--retry-failed-cpu",
        action="store_true",
        help="retry failed GPU scripts forcing CPU",
    )
    args = ap.parse_args()

    scripts = FULL_SCRIPTS if args.psicov_mode == "full" else CORE_SCRIPTS
    targets = dataset_targets()
    if args.retry_failed_cpu:
        os.environ["COEVO_FORCE_CPU"] = "1"
    if args.retry_failed or args.retry_failed_cpu:
        import shutil as _sh

        trash = REPO.parent / ".rigor-trash" / "retried"
        for lf in RESULTS.glob("*/*/*/run_log.json"):
            if args.dataset and lf.parts[-4] != args.dataset:
                continue
            d = json.load(open(lf))
            if d.get("status") not in ("ok",):
                dest = trash / lf.parent.name
                dest.mkdir(parents=True, exist_ok=True)
                _sh.move(
                    str(lf.parent),
                    str(dest / f"{lf.parent.parent.name}_{lf.parent.name}_{int(time.time())}"),
                )
                print(f"  retry: {lf.parent}", flush=True)
    jobs = []
    for ds, tlist in targets.items():
        if args.dataset and ds != args.dataset:
            continue
        for tgt, fasta in tlist:
            if args.target and tgt != args.target:
                continue
            for s in scripts:
                jobs.append((ds, tgt, s, fasta, RESULTS))
    if args.max_targets:
        seen = set()
        jobs = [j for j in jobs if not (seen.add(j[1]) or len(seen) > args.max_targets)]
    # Two waves: analysis scripts first, report generators after (they read
    # the analysis outputs, so they must run last).
    gen_jobs = [j for j in jobs if j[2] in REPORT_GENERATORS]
    analysis_jobs = [j for j in jobs if j[2] not in REPORT_GENERATORS]
    print(
        f"jobs: {len(jobs)} (analysis={len(analysis_jobs)}, "
        f"generators={len(gen_jobs)}, parallel={args.parallel})",
        flush=True,
    )
    t0 = time.time()
    done = {"ok": 0, "skipped": 0, "fail": 0}
    for wave in (analysis_jobs, gen_jobs):
        if not wave:
            continue
        with ThreadPoolExecutor(max_workers=args.parallel) as ex:
            for ds, tgt, s, st in ex.map(run_one, wave):
                done[st if st in done else "fail"] += 1
                if st != "ok" and st != "skipped":
                    print(f"  [{st}] {ds}/{tgt}/{s}", flush=True)
    print(f"summary: {done} in {round(time.time() - t0, 1)}s", flush=True)


if __name__ == "__main__":
    main()
