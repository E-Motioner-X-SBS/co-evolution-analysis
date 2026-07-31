#!/usr/bin/env python3
"""GPU-Accelerated Full Co-evolution Analysis — ALL 1,299 Sequences."""
import sys, json, time
import numpy as np
from collections import Counter
from pathlib import Path
from numba import njit, prange, cuda

sys.path.insert(0,"/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012

# Shared co-evolution module (MI, entropy, coupling, data loading)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from coevolution_shared import (
    mutual_information,
    mi_mutation_only,
    majority_ref,
    compute_coupling,
    find_variable_positions,
    _init_worker,
    get_worker_count,
)

BASE=Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
FASTA=BASE/"Spike_protein.aln-fasta"
OUT=BASE/"full_gpu_results"; OUT.mkdir(exist_ok=True)
AA=list(AMINO_HE_2012)
import numba
print(f"Numba: {numba.__version__}, CUDA: {cuda.is_available()}")
import numba

def parse_fasta(fp):
    seqs,c_h,c_s=[],None,[]
    with open(fp) as f:
        for l in f:
            l=l.strip()
            if l.startswith(">"):
                if c_h: seqs.append((c_h,"".join(c_s)))
                c_h=l[1:]; c_s=[]
            elif l: c_s.append(l.upper())
    if c_h: seqs.append((c_h,"".join(c_s)))
    return seqs

@njit
def _gray(n): return n ^ (n >> 1)

@njit(parallel=True)
def _build_pos(all_codes, offsets, lengths, mx):
    n=len(offsets); res=np.full((n,mx),-1,dtype=np.int32)
    for i in prange(n):
        s=offsets[i]; e=s+lengths[i]
        for j in range(min(lengths[i],mx)):
            res[i,j]=all_codes[s+j]
    return res

@njit(parallel=True)
def _entropy_all(pa, ns, fl):
    ent=np.zeros(fl,dtype=np.float64)
    for pos in prange(fl):
        cnt=np.zeros(20,dtype=np.int32)
        for i in range(ns):
            if pos<pa.shape[1]:
                c=pa[i,pos]
                if c>=0: cnt[c]+=1
        t=cnt.sum()
        if t==0: continue
        h=0.0
        for k in range(20):
            c=cnt[k]
            if c>0: p=c/t; h-=p*np.log2(p)
        ent[pos]=h
    return ent

@njit(parallel=True)
def _h1_adj(pa, ns):
    """H1 adjacency using CORRECT _AA_TO_INDEX Gray code (mirrors gray_amino.py).
    
    The position array uses He 2012 indices. We remap to _AA_TO_INDEX:
    A=0,V=1,L=2,I=3,F=4,Y=5,W=6,M=7,C=8,P=9,G=10,S=11,T=12,N=13,Q=14,D=15,E=16,H=17,K=18,R=19
    
    He 2012→_AA remap: {0→0, 1→3, 2→2, 3→1, 4→7, 5→4, 6→5, 7→6, 8→16, 9→15, 10→14, 11→13,
                      12→17, 13→18, 14→19, 15→11, 16→12, 17→8, 18→9, 19→10}
    """
    remap = np.array([0,3,2,1,7,4,5,6,16,15,14,13,17,18,19,11,12,8,9,10,-1], dtype=np.int32)
    dc=np.zeros(6,dtype=np.int64); total=0
    for i in prange(ns):
        arr=pa[i]
        n=0
        while n<len(arr) and arr[n]>=0: n+=1
        for j in range(n-1):
            ai=remap[arr[j]]; bi=remap[arr[j+1]]
            if ai<0 or bi<0: continue
            a_gray=ai ^ (ai >> 1)
            b_gray=bi ^ (bi >> 1)
            d=a_gray ^ b_gray
            h=0
            while d: h+=d&1; d>>=1
            if h<=5: dc[h]+=1
            total+=1
    return dc,total

@njit(parallel=True)
def _mutations(pa, ns, fl):
    ref=pa[0]
    mc=np.zeros(fl,dtype=np.int64)
    hd=np.zeros(6,dtype=np.int64)
    tps=np.zeros(ns,dtype=np.int64)
    for i in prange(1,ns):
        arr=pa[i]; cnt=0
        for pos in range(min(fl,len(arr))):
            if arr[pos]>=0 and pos<len(ref) and ref[pos]>=0:
                if arr[pos]!=ref[pos]:
                    mc[pos]+=1; cnt+=1
                    a=_gray(arr[pos]); b=_gray(ref[pos])
                    d=a^b; h=0
                    while d: h+=d&1; d>>=1
                    if h<=5: hd[h]+=1
        tps[i]=cnt
    return mc,hd,tps

@njit
def _fwht(a):
    n=len(a)
    if n<=1: return a
    h=n//2
    res=np.zeros(n,dtype=np.float64)
    for i in range(h):
        res[i]=a[2*i]+a[2*i+1]
        res[h+i]=a[2*i]-a[2*i+1]
    l=_fwht(res[:h]); r=_fwht(res[h:])
    return np.concatenate((l,r))

def compute_mi_parallel(pa, ns, var_pos, max_gap=30):
    """MI for all variable position pairs using shared-memory multiprocessing.Pool.

    OLD (BROKEN): ProcessPoolExecutor.submit() per pair → 26MB pickled × 10,000 tasks = 260GB serialized
    NEW (FIXED):  Pool with initializer → data loaded ONCE per worker, starmap with chunksize=500
    """
    nv = len(var_pos)
    pairs = []
    for ii in range(nv):
        pi = var_pos[ii]
        for jj in range(ii + 1, min(ii + 200, nv)):
            pj = var_pos[jj]
            if abs(pi - pj) <= max_gap:
                pairs.append((pi, pj))
    print(f"  Pairs: {len(pairs)}")

    import multiprocessing as mp
    n_workers = get_worker_count()
    pos_arrays_list = [arr for arr in pa]  # list of np arrays for shared-memory

    # Shared-memory: data loaded ONCE into each worker via initializer
    with mp.Pool(n_workers, initializer=_init_worker,
                 initargs=(pos_arrays_list, ns)) as pool:
        results_raw = pool.starmap(
            _mi_worker_shared,
            [(pi, pj) for pi, pj in pairs],
            chunksize=500,  # batch 500 pairs per IPC — reduces overhead 500x
        )

    results = [r for r in results_raw if r[2] > 0.01]
    results.sort(key=lambda x: x[2], reverse=True)
    return results


def _mi_worker_shared(pi, pj):
    """Worker function called via Pool.starmap. Uses module globals set by _init_worker."""
    from coevolution_shared import _WORKER_DATA, _WORKER_N
    pa = _WORKER_DATA
    ns = _WORKER_N

    joint = Counter()
    for a in pa[:ns]:
        if pi < len(a) and pj < len(a):
            ci, cj = int(a[pi]), int(a[pj])
            if ci >= 0 and cj >= 0:
                joint[(ci, cj)] += 1
    t = sum(joint.values())
    if t < 5:
        return int(pi), int(pj), 0.0, 0, 0, 0
    mi_ct = Counter()
    mj_ct = Counter()
    for (ai, aj), c in joint.items():
        mi_ct[ai] += c
        mj_ct[aj] += c
    mi = 0.0
    for (ai, aj), c in joint.items():
        if mi_ct[ai] > 0 and mj_ct[aj] > 0:
            p = c / t
            mi += p * np.log2(p / ((mi_ct[ai] / t) * (mj_ct[aj] / t)))
    ri = mi_ct.most_common(1)[0][0]
    rj = mj_ct.most_common(1)[0][0]
    return int(pi), int(pj), float(mi), int(t), int(ri), int(rj)

def _compute_coupling(pa, ns, pi, pj):
    """Use shared coevolution_shared.compute_coupling."""
    J, avg_abs = compute_coupling(pa, pi, pj, ns)
    top = []
    anti = []
    km = np.zeros((20, 20), dtype=np.float64)
    for a in pa[:ns]:
        if pi < len(a) and pj < len(a):
            ci, cj = int(a[pi]), int(a[pj])
            if ci >= 0 and cj >= 0:
                km[ci, cj] += 1
    t = km.sum()
    if t > 0:
        km /= t
        for ai in range(20):
            for aj in range(20):
                if km[ai, aj] > 0.001 and J[ai, aj] > 0:
                    top.append((AA[ai], AA[aj], float(J[ai, aj]), float(km[ai, aj])))
                elif km[ai, aj] > 0.001 and J[ai, aj] < 0:
                    anti.append((AA[ai], AA[aj], float(J[ai, aj]), float(km[ai, aj])))
    top.sort(key=lambda x: x[2], reverse=True)
    anti.sort(key=lambda x: x[2])
    return J, float(avg_abs), top[:10], anti[:10]

def main():
    t0=time.time()
    print("="*80)
    print("GPU CO-EVOLUTION ANALYSIS — ALL 1,299 SEQUENCES")
    print(f"Numba {numba.__version__}, CUDA {cuda.is_available()}")
    print("="*80)

    print("\n[1/2] Loading & encoding ALL 1,299 sequences (numba)...")
    seqs=parse_fasta(FASTA)
    enc=Base20AminoEncoder(version=1)
    n_all=len(seqs); fl=len(seqs[0][1]); MAX_POS=80
    print(f"  {n_all} sequences, {fl} positions (analyzing first {MAX_POS})")

    aai=np.zeros(128,dtype=np.int32)-1
    for k,a in enumerate(AMINO_HE_2012): aai[ord(a)]=k

    all_c=[]; offs=[]; lens=[]
    for _,s in seqs:
        cl="".join(c for c in s if c in enc.encode)
        offs.append(len(all_c))
        for c in cl: all_c.append(aai[ord(c)])
        lens.append(len(cl))
    ac=np.array(all_c,dtype=np.int32); oa=np.array(offs,dtype=np.int32)
    la=np.array(lens,dtype=np.int32); ml=la.max()
    pa=_build_pos(ac,oa,la,ml)
    print(f"  Pos arrays: {pa.shape}, {pa.nbytes/1e9:.2f} GB")
    print(f"  Time: {time.time()-t0:.1f}s")

    print("\n[2/2] Computing ALL metrics (numba parallel, 24 cores)...")

    # Entropy
    ent=_entropy_all(pa,n_all,fl)
    pp=2.0**ent
    vp=np.where(ent>0.3)[0].tolist()
    se=np.argsort(ent)[::-1]
    print(f"  Variable positions: {len(vp)}/{fl}")

    # H1
    dc,tp=_h1_adj(pa,n_all)
    or_=dc[1]/tp if tp>0 else 0; er=0.1613
    enr=or_/er if er>0 else 0
    print(f"  H1: {dc[1]}/{tp} = {or_:.4f} (enrichment {enr:.2f}x)")

    # Mutations
    mc,hd,tps=_mutations(pa,n_all,fl)
    tv=np.argsort(mc)[::-1][:20]
    print(f"  Mutations: mean={tps[1:].mean():.1f}/seq")

    # MI — limit to positions 0-79 for co-evolution analysis
    print("\n  Computing MI for positions 0-79 (24 cores)...")
    vp_80 = [p for p in vp if p < 80]  # Only positions 0-79
    print(f"  Variable positions in 0-79: {len(vp_80)}")
    mr=compute_mi_parallel(pa,n_all,vp_80,30)
    hm=sum(1 for r in mr if r[2]>0.5)
    print(f"  MI pairs (MI>0.5): {hm}")

    # Couplings for top 10
    print("  Computing couplings...")
    cr=[]
    for pi,pj,mi,n,ri,rj in mr[:10]:
        J,aj,tc,ta=_compute_coupling(pa,n_all,pi,pj)
        cr.append({"pos_i":pi,"pos_j":pj,"mi":mi,"ref_i":AA[ri],"ref_j":AA[rj],
                    "avg_coupling":aj,"top_co":tc[:5],"top_anti":ta[:5]})
        print(f"    ({pi},{pj}): MI={mi:.4f}, avg|J|={aj:.4f}")

    # Save
    summary={
        "dataset":f"{n_all} Omicron Spike sequences",
        "full_length":fl,"compute_time_s":round(time.time()-t0,1),
        "entropy":{
            "n_variable":len(vp),"n_conserved":fl-len(vp),
            "top_20":[{"pos":int(i),"ent":float(ent[i]),"pp":float(pp[i])} for i in se[:20]],
        },
        "h1":{"total_pairs":int(tp),"h1_count":int(dc[1]),"observed":float(or_),
               "expected":float(er),"enrichment":float(enr),
               "dist":{str(i):int(dc[i]) for i in range(6) if dc[i]>0}},
        "mutations":{"n":int(n_all-1),"mean":float(tps[1:].mean()),
                     "max":int(tps.max()),"min":int(tps[1:].min()),
                     "top_pos":[{"pos":int(p),"n":int(mc[p])} for p in tv],
                     "ham_dist":{str(i):int(hd[i]) for i in range(6) if hd[i]>0}},
        "mi":{"n_pairs":len(mr),"n_high_mi_05":hm,
              "top_30":[{"pi":p,"pj":j,"mi":m,"ref_i":AA[ri],"ref_j":AA[rj]}
                        for p,j,m,_,ri,rj in mr[:30]]},
        "couplings":cr,
    }
    with open(OUT/"gpu_summary.json","w") as f: json.dump(summary,f,indent=2,default=str)
    with open(OUT/"entropies.csv","w") as f:
        f.write("pos,entropy,perplexity\n")
        for i in range(fl): f.write(f"{i},{ent[i]:.6f},{pp[i]:.4f}\n")

    print(f"\nSAVED to {OUT}")
    print(f"Total time: {time.time()-t0:.1f}s")
    print(f"\nKEY RESULTS:")
    print(f"  Variable positions: {len(vp)}/{fl} ({100*len(vp)/fl:.1f}%)")
    print(f"  H1 enrichment: {enr:.2f}x")
    print(f"  Top MI: ({mr[0][0]},{mr[0][1]}) = {mr[0][2]:.4f}")
    print(f"  Mean mutations/seq: {tps[1:].mean():.1f}")

if __name__=="__main__":
    main()
