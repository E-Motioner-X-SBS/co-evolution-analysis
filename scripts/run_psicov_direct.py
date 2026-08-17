#!/usr/bin/env python3
"""Direct driver for remaining psicov150 jobs (bypasses runner executor)."""
import json, os, subprocess, sys, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
PY = "/store/shuvam/.venv/bin/python"
RESULTS = REPO / "results"

DM = {"master_boolean.py":"master_boolean","boolean_co-evolution.py":"boolean_results",
"nary_kmap_co-evolution.py":"nary_kmap_results","position_kmap_coevolution.py":"position_kmap_results",
"run_allseq_analysis.py":"full_position_results","run_kmap_analysis.py":"kmap_results",
"flipped_boolean_coevolution.py":"flipped_boolean_results","kmap_boolean_coevolution.py":"kmap_boolean_coevolution",
"variable_position_coevolution.py":"variable_position_results","predictive_constraint_function.py":"constraint_function_results",
"allseq_constraint_function.py":"allseq_constraint_results","dca_boolean_coevolution.py":"dca_boolean_results",
"dca_mf_analysis.py":"dca_results","perplexity_coevolution.py":"perplexity_results",
"advanced_co-evolution_analysis.py":"advanced_analysis_results","full_length_analysis.py":"full_length_results",
"gpu_full_analysis.py":"full_gpu_results","create_mi_heatmap.py":"mi_heatmap"}
GEN = ["generate_co-evolution_md.py","generate_full_analysis_md.py","generate_full_pipeline_doc.py"]
FULL = list(DM) + GEN
GEN_DIRS = {"generate_co-evolution_md.py":"generate_co-evolution_md","generate_full_analysis_md.py":"generate_full_analysis_md","generate_full_pipeline_doc.py":"generate_full_pipeline_doc"}

jobs = []
for t in sorted((RESULTS/"psicov150").iterdir()):
    if not t.is_dir(): continue
    for s in FULL:
        dn = DM.get(s, GEN_DIRS.get(s))
        logf = t / dn / "run_log.json"
        if not logf.exists():
            jobs.append((t.name, s))

print(f"missing jobs: {len(jobs)}", flush=True)

def run_one(job):
    tgt, script = job
    dn = DM.get(script, GEN_DIRS.get(script))
    outdir = RESULTS / "psicov150" / tgt / dn
    outdir.mkdir(parents=True, exist_ok=True)
    logf = outdir / "run_log.json"
    if logf.exists():
        return "skipped"
    env = dict(os.environ)
    env["OMP_PROC_BIND"] = "FALSE"
    env["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
    if script in ("allseq_constraint_function.py","predictive_constraint_function.py"):
        env["CUDA_VISIBLE_DEVICES"] = ""
    env["COEVO_FASTA"] = str(REPO/"data"/"psicov150"/"converted"/f"{tgt}.fasta")
    if script in GEN:
        env["COEVO_RESULTS"] = str(RESULTS/"psicov150"/tgt)
        env["COEVO_OUT"] = str(outdir / script.replace(".py",".md"))
    else:
        env["COEVO_RESULTS"] = str(outdir)
    t0 = time.time()
    try:
        r = subprocess.run([PY, str(REPO/script)], cwd=str(outdir), env=env,
                           capture_output=True, text=True, timeout=None)
        status = "ok" if r.returncode == 0 else f"exit{r.returncode}"
        tail = (r.stdout + "\n" + r.stderr)[-2000:]
    except Exception as e:
        status, tail = "exc", str(e)[:2000]
    json.dump({"dataset":"psicov150","target":tgt,"script":script,"status":status,
               "seconds":round(time.time()-t0,1),"tail":tail,
               "outputs":[f.name for f in outdir.iterdir() if f.name!="run_log.json"]},
              open(logf,"w"), indent=1)
    return status

done = {"ok":0,"fail":0,"skip":0}
with ThreadPoolExecutor(max_workers=12) as ex:
    for st in ex.map(run_one, jobs):
        done[st if st in done else "fail"] += 1
print(f"summary: {done}", flush=True)
