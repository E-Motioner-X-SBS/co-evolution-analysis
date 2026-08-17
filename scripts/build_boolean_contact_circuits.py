#!/usr/bin/env python3
"""
build_boolean_contact_circuits.py
=================================
Build Boolean circuits from FASTA sequence data that predict 3D contacts.

Phase 1: Per-pair contact circuits (residue identity → contact)
  For each variable position pair, build a 32×32 contact K-map from 299 models.
  QM minimize → Boolean expression. 5-fold CV.

Phase 2: Cross-pair circuit (physicochemical group + separation → contact)
  Pool all variable pairs. 8-bit K-map (3+3+2). Train/test by position pair.

Phase 3: Co-evolution transfer circuit (MI + perplexity + sep → contact)
  6-bit K-map (2+2+2). LOO-CV. Tests co-evolution → structure transfer.

Output: results/contacts/boolean_contact_circuits.json + printed report.
"""
import sys, json, math, random
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")

from rules3d_common import OmicronModels, CANON, dist, CONTACT, MIN_SEP
from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap
from nkmap.encoding.bio_sequences import AMINO_HE_2012
from coevolution_shared import (load_position_arrays, compute_entropy_vectorized,
                                 mutual_information, perplexity_ratio,
                                 dense_from_arrays)

AA_LIST = list(AMINO_HE_2012)
N_AA = 20

# Physicochemical groups (3-bit, 8 groups)
AA_GROUPS = {'A':0,'I':0,'L':0,'V':0,'M':0, 'F':1,'Y':1,'W':1,
             'E':2,'D':2, 'Q':3,'N':3, 'H':4,'K':4,'R':4,
             'S':5,'T':5, 'C':6, 'P':7,'G':7}
GROUP_NAMES = ['hydrophobic','aromatic','acidic','amide','basic','hydroxyl','sulfur','special']
SEP_NAMES = ['3-5','6-10','11-20','21+']

def gray_code(n): return n ^ (n >> 1)

def pi_to_residues(pi, offset, n_bits=5):
    vals = list(pi["values"]); mask = list(pi["mask"])
    while len(vals) < offset + n_bits: vals.append(0); mask.append(False)
    matching = []
    for code in range(20):
        g = gray_code(code)
        bits = [(g >> (n_bits-1-k)) & 1 for k in range(n_bits)]
        if all(bits[k] == vals[offset+k] for k in range(n_bits) if not mask[offset+k]):
            matching.append(AA_LIST[code])
    return matching

def format_pi_residues(pi):
    ri = pi_to_residues(pi, 0)
    rj = pi_to_residues(pi, 5)
    return f"aa_i∈{{{','.join(ri)}}} ∧ aa_j∈{{{','.join(rj)}}}"

def sep_bin(s): return 0 if s<=5 else (1 if s<=10 else (2 if s<=20 else 3))

def metrics(preds, acts):
    tp=sum(1 for p,a in zip(preds,acts) if p==1 and a==1)
    tn=sum(1 for p,a in zip(preds,acts) if p==0 and a==0)
    fp=sum(1 for p,a in zip(preds,acts) if p==1 and a==0)
    fn=sum(1 for p,a in zip(preds,acts) if p==0 and a==1)
    n=len(preds); acc=(tp+tn)/n if n else 0
    prec=tp/(tp+fp) if (tp+fp) else 0
    rec=tp/(tp+fn) if (tp+fn) else 0
    f1=2*prec*rec/(prec+rec) if (prec+rec) else 0
    s=(tp+fn)*(tp+fp)*(tn+fn)*(tn+fp)
    mcc=(tp*tn-fp*fn)/math.sqrt(s) if s>0 else 0
    return {"accuracy":round(acc,4),"precision":round(prec,4),"recall":round(rec,4),
            "f1":round(f1,4),"mcc":round(mcc,4),"tp":tp,"tn":tn,"fp":fp,"fn":fn,"n":n}

def build_32x32(votes, min_obs=2, thresh=0.5):
    km = np.full((32,32), -1, dtype=np.int32); km[:20,:20] = 0
    for (ci,cj),(c,t) in votes.items():
        if t >= min_obs: km[ci,cj] = 1 if c/t > thresh else 0
    return km

def eval_32x32(km_flat, test_data):
    preds=[]; acts=[]
    for aa_i, aa_j, ic in test_data:
        ci = AA_LIST.index(aa_i) if aa_i in AA_LIST else -1
        cj = AA_LIST.index(aa_j) if aa_j in AA_LIST else -1
        if ci<0 or cj<0: continue
        v = int(km_flat[ci*32+cj]); preds.append(0 if v==-1 else v); acts.append(ic)
    return preds, acts

def decode_p2_cell(idx):
    gi=(idx>>5)&7; gj=(idx>>2)&7; sb=idx&3
    return f"{GROUP_NAMES[gi]}-{GROUP_NAMES[gj]} sep={SEP_NAMES[sb]}"

def decode_p3_cell(idx):
    mb=(idx>>4)&3; rb=(idx>>2)&3; sb=idx&3
    mi_n=['MI<0.1','MI.1-.3','MI.3-.6','MI>.6']
    ra_n=['r<1.0','r1-1.5','r1.5-2','r>2']
    return f"{mi_n[mb]} {ra_n[rb]} sep={SEP_NAMES[sb]}"

def main():
    print("="*80); print("BOOLEAN CONTACT CIRCUITS FROM FASTA SEQUENCE DATA"); print("="*80)
    # 1. Load models + variable positions
    print("\n[1] Loading 299 models + variable positions...")
    models = OmicronModels()
    pa, na, fl = load_position_arrays(aligned=True)
    ent = compute_entropy_vectorized(pa, na, fl)
    vp = [p for p in range(fl) if ent[p] > 0.3]
    pairs = [(i,j) for a,i in enumerate(vp) for j in vp[a+1:] if j-i <= 30]
    print(f"  Variable positions: {len(vp)} | Pairs: {len(pairs)}")

    # 2. Collect per-pair per-model data
    print("\n[2] Collecting per-model contact data...")
    pdata = {p: [] for p in pairs}
    mc = 0
    for uid, aligned, cmap, coords in models.iter_models():
        mc += 1; ch = sorted(coords)
        if not ch: continue
        ch = ch[0]
        for (i,j) in pairs:
            ri, rj = cmap.get(i), cmap.get(j)
            if ri is None or rj is None or abs(ri-rj) < MIN_SEP: continue
            if ri not in coords[ch] or rj not in coords[ch]: continue
            ai = aligned[i] if i < len(aligned) else '?'
            aj = aligned[j] if j < len(aligned) else '?'
            if ai not in CANON or aj not in CANON: continue
            d = dist(coords[ch][ri], coords[ch][rj])
            pdata[(i,j)].append((ai, aj, 1 if d < CONTACT else 0, d))
    print(f"  Models: {mc}")

    # 3. Pair-level features
    print("\n[3] Computing MI + perplexity ratio...")
    dense = dense_from_arrays(pa, na)
    pfeat = {}
    for (i,j) in pairs:
        mi = mutual_information(pa, i, j, na)
        rat = perplexity_ratio(pa, i, j, na, ent)
        pfeat[(i,j)] = {"mi": mi, "ratio": rat if rat else 0, "sep": j-i}

    # ===== PHASE 1 =====
    print("\n" + "="*80); print("PHASE 1: PER-PAIR CONTACT CIRCUITS"); print("="*80)
    p1 = {}
    for pair in pairs:
        data = pdata[pair]
        if len(data) < 10: continue
        rp = Counter((ai,aj) for ai,aj,_,_ in data)
        contacts = sum(c for _,_,c,_ in data)
        cr = contacts / len(data)
        if cr == 0 or cr == 1:
            p1[str(pair)] = {"trivial": True, "contact_rate": round(cr,3),
                             "n_models": len(data), "n_residue_pairs": len(rp),
                             "circuit": "CONSTANT " + str(int(cr))}
            continue
        votes = defaultdict(lambda: [0,0])
        for ai,aj,ic,_ in data:
            ci=AA_LIST.index(ai); cj=AA_LIST.index(aj)
            votes[(ci,cj)][1] += 1
            if ic: votes[(ci,cj)][0] += 1
        km = build_32x32(votes); kmf = km.flatten().astype(int)
        try:
            res = boolean_minimize_kmap(kmf)
            pis = res.get("prime_implicants",[])
            ess = res.get("essential_prime_implicants",[])
        except Exception as e:
            p1[str(pair)] = {"error": str(e)}; continue
        exprs = [format_pi_residues(p) for p in ess]
        td = [(ai,aj,c) for ai,aj,c,_ in data]
        pr, ar = eval_32x32(kmf, td)
        tm = metrics(pr, ar)
        # 5-fold CV
        random.seed(42); idx = list(range(len(data))); random.shuffle(idx)
        fs = len(idx)//5; cp=[]; ca=[]
        for f in range(5):
            ti = set(idx[f*fs:(f+1)*fs])
            tr = [data[k] for k in range(len(data)) if k not in ti]
            te = [data[k] for k in ti]
            cv = defaultdict(lambda: [0,0])
            for ai,aj,ic,_ in tr:
                ci=AA_LIST.index(ai); cj=AA_LIST.index(aj)
                cv[(ci,cj)][1]+=1
                if ic: cv[(ci,cj)][0]+=1
            ckm = build_32x32(cv); ckf = ckm.flatten().astype(int)
            ctd = [(ai,aj,c) for ai,aj,c,_ in te]
            p,a = eval_32x32(ckf, ctd); cp.extend(p); ca.extend(a)
        cvm = metrics(cp, ca)
        baseline = max(cr, 1-cr)
        p1[str(pair)] = {"trivial": False, "n_models": len(data),
                         "n_residue_pairs": len(rp), "contact_rate": round(cr,3),
                         "n_pis": len(pis), "n_essential": len(ess),
                         "essential_expressions": exprs,
                         "all_expressions": [format_pi_residues(p) for p in pis[:10]],
                         "train_metrics": tm, "cv_metrics": cvm,
                         "baseline_accuracy": round(baseline,3)}
        print(f"\n  Pair {pair}: {len(data)} models, {len(rp)} residue pairs, CR={cr:.3f}")
        print(f"    PIs: {len(pis)} ({len(ess)} essential)")
        for e in exprs[:5]: print(f"    → {e}")
        print(f"    Train: acc={tm['accuracy']:.3f} prec={tm['precision']:.3f} rec={tm['recall']:.3f} F1={tm['f1']:.3f} MCC={tm['mcc']:.3f}")
        print(f"    5-fold CV: acc={cvm['accuracy']:.3f} prec={cvm['precision']:.3f} rec={cvm['recall']:.3f} F1={cvm['f1']:.3f} MCC={cvm['mcc']:.3f}")
        print(f"    Baseline: {baseline:.3f}")

    # ===== PHASE 2 =====
    print("\n" + "="*80); print("PHASE 2: CROSS-PAIR CIRCUIT (group+sep → contact)"); print("="*80)
    all_dp = []
    for pair in pairs:
        data = pdata[pair]
        sb = sep_bin(pair[1]-pair[0])
        for ai,aj,ic,_ in data:
            all_dp.append((AA_GROUPS.get(ai,7), AA_GROUPS.get(aj,7), sb, ic, pair))
    print(f"  Data points: {len(all_dp)}")
    votes = defaultdict(lambda: [0,0])
    for gi,gj,sb,ic,_ in all_dp:
        idx = (gi<<5)|(gj<<2)|sb; votes[idx][1]+=1
        if ic: votes[idx][0]+=1
    truth = np.full(256, -1, dtype=np.int32)
    for idx,(c,t) in votes.items():
        if t >= 3: truth[idx] = 1 if c/t > 0.5 else 0
    try:
        res = boolean_minimize_kmap(truth)
        gpis = res.get("prime_implicants",[]); gess = res.get("essential_prime_implicants",[])
        print(f"  Global circuit: {len(gpis)} PIs ({len(gess)} essential)")
    except Exception as e:
        print(f"  QM failed: {e}"); gpis=[]; gess=[]
    # Report on-set cells
    on_cells = [(idx, votes[idx]) for idx in range(256) if truth[idx]==1]
    print(f"  On-set cells (contact rules): {len(on_cells)}")
    for idx,(c,t) in sorted(on_cells, key=lambda x:-x[1][0])[:15]:
        print(f"    {decode_p2_cell(idx)}: {c}/{t} contact ({c/t:.1%})")
    # Train/test by pair
    random.seed(42); up = list(set(p for *_,p in all_dp)); random.shuffle(up)
    ntr = int(0.8*len(up)); trp=set(up[:ntr]); tep=set(up[ntr:])
    tvotes = defaultdict(lambda: [0,0])
    for gi,gj,sb,ic,pair in all_dp:
        if pair not in trp: continue
        idx=(gi<<5)|(gj<<2)|sb; tvotes[idx][1]+=1
        if ic: tvotes[idx][0]+=1
    ttruth = np.full(256,-1,dtype=np.int32)
    for idx,(c,t) in tvotes.items():
        if t>=3: ttruth[idx] = 1 if c/t>0.5 else 0
    tpreds=[]; tacts=[]
    for gi,gj,sb,ic,pair in all_dp:
        if pair not in tep: continue
        idx=(gi<<5)|(gj<<2)|sb; v=int(ttruth[idx]); tpreds.append(0 if v==-1 else v); tacts.append(ic)
    tm2 = metrics(tpreds, tacts)
    tcr = sum(tacts)/len(tacts) if tacts else 0
    print(f"\n  Train pairs: {len(trp)} | Test pairs: {len(tep)}")
    print(f"  Test CR: {tcr:.3f} | Baseline: {1-tcr:.3f}")
    print(f"  Test: acc={tm2['accuracy']:.3f} prec={tm2['precision']:.3f} rec={tm2['recall']:.3f} F1={tm2['f1']:.3f} MCC={tm2['mcc']:.3f}")

    # ===== PHASE 3 =====
    print("\n" + "="*80); print("PHASE 3: CO-EVOLUTION TRANSFER (MI+ratio+sep → contact)"); print("="*80)
    pl = []
    for pair in pairs:
        data = pdata[pair]
        if len(data) < 10: continue
        contacts = sum(c for _,_,c,_ in data)
        cf = contacts/len(data)
        f = pfeat[pair]
        mb = 0 if f["mi"]<0.1 else (1 if f["mi"]<0.3 else (2 if f["mi"]<0.6 else 3))
        rb = 0 if f["ratio"]<1.0 else (1 if f["ratio"]<1.5 else (2 if f["ratio"]<2.0 else 3))
        sb = sep_bin(f["sep"])
        label = 1 if cf > 0.3 else 0
        pl.append({"pair":pair,"mi":f["mi"],"ratio":f["ratio"],"sep":f["sep"],
                    "mi_bin":mb,"ratio_bin":rb,"sep_bin":sb,"contact_frac":cf,"label":label})
    print(f"  Pairs: {len(pl)} | Contact (frac>0.3): {sum(1 for p in pl if p['label']==1)}")
    cvotes = defaultdict(lambda: [0,0])
    for p in pl:
        idx=(p["mi_bin"]<<4)|(p["ratio_bin"]<<2)|p["sep_bin"]
        cvotes[idx][1]+=1
        if p["label"]==1: cvotes[idx][0]+=1
    ctruth = np.full(64,-1,dtype=np.int32)
    for idx,(c,t) in cvotes.items():
        if t>=1: ctruth[idx] = 1 if c/t>0.5 else 0
    try:
        cres = boolean_minimize_kmap(ctruth)
        cpis = cres.get("prime_implicants",[]); cess = cres.get("essential_prime_implicants",[])
        print(f"  Circuit: {len(cpis)} PIs ({len(cess)} essential)")
    except Exception as e:
        print(f"  QM failed: {e}"); cpis=[]; cess=[]
    on_c = [(idx, cvotes[idx]) for idx in range(64) if ctruth[idx]==1]
    print(f"  On-set cells (contact rules): {len(on_c)}")
    for idx,(c,t) in sorted(on_c, key=lambda x:-x[1][0]):
        print(f"    {decode_p3_cell(idx)}: {c}/{t} contact")
    # LOO-CV
    lp=[]; la=[]
    for k in range(len(pl)):
        tr=[pl[i] for i in range(len(pl)) if i!=k]; te=pl[k]
        cv=defaultdict(lambda: [0,0])
        for p in tr:
            idx=(p["mi_bin"]<<4)|(p["ratio_bin"]<<2)|p["sep_bin"]
            cv[idx][1]+=1
            if p["label"]==1: cv[idx][0]+=1
        ct=np.full(64,-1,dtype=np.int32)
        for idx,(c,t) in cv.items():
            if t>=1: ct[idx]=1 if c/t>0.5 else 0
        idx=(te["mi_bin"]<<4)|(te["ratio_bin"]<<2)|te["sep_bin"]
        v=int(ct[idx]); lp.append(0 if v==-1 else v); la.append(te["label"])
    lm = metrics(lp, la)
    lcr = sum(la)/len(la) if la else 0
    print(f"\n  LOO-CV: acc={lm['accuracy']:.3f} prec={lm['precision']:.3f} rec={lm['recall']:.3f} F1={lm['f1']:.3f} MCC={lm['mcc']:.3f}")
    print(f"  Baseline (no contact): {1-lcr:.3f}")

    # Save
    OUT = REPO/"results"/"contacts"/"boolean_contact_circuits.json"
    json.dump({"phase1": p1,
               "phase2": {"n_data":len(all_dp),"n_train_pairs":len(trp),"n_test_pairs":len(tep),
                           "n_pis":len(gpis),"n_essential":len(gess),
                           "on_set_cells":[{"cell":decode_p2_cell(idx),"contacts":c,"total":t,"rate":round(c/t,3)} for idx,(c,t) in on_cells],
                           "test_metrics": tm2, "test_cr": round(tcr,3), "baseline": round(1-tcr,3)},
               "phase3": {"n_pairs":len(pl),"n_contact":sum(1 for p in pl if p['label']==1),
                           "n_pis":len(cpis),"n_essential":len(cess),
                           "on_set_cells":[{"cell":decode_p3_cell(idx),"contacts":c,"total":t} for idx,(c,t) in on_c],
                           "loo_metrics": lm, "contact_rate": round(lcr,3), "baseline": round(1-lcr,3)}},
              open(OUT,"w"), indent=1, default=str)
    print(f"\nSaved -> {OUT}")

if __name__ == "__main__":
    main()
