#!/usr/bin/env python3
"""Aggregate the long run over all joinable 2_2->3_2 rules (results/rules_joinable*; started by scripts/run_all_joinable.sh).
Writes paper/generated/tab_joinable.tex and joinable_stats.tex.  Safe to run on partial data (reports how many rules/runs are complete)."""
import csv, os
from collections import defaultdict, Counter
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
res = os.path.join(root, "results"); gen = os.path.join(root, "paper", "generated")
items = [l.strip().split("|") for l in open(os.path.join(res, "rules_joinable.txt")) if l.strip()]
cls = {r["rule"]: r for r in csv.DictReader(open(os.path.join(res, "classification", "classification.csv")))}
nec = {r["rule"]: r for r in csv.DictReader(open(os.path.join(res, "necessity", "necessity_summary.csv")))}

def rows(d, rid):
    f = os.path.join(res, d, rid, "evolution.csv")
    if not os.path.exists(f): return []
    return list(csv.DictReader(open(f)))

def runs(rs):
    g = defaultdict(list)
    for r in rs: g[(r["background"], r["experiment"], r["defect"], r["order"])].append(r)
    return g

per = {}
cor_tot = cor_ok = 0
for rid, rule in items:
    d = int(cls[rule]["delta_R"]); info = dict(rule=rule, bucket=cls[rule]["bucket"], sm=nec[rule]["simple_minor_preserving"])
    cone = runs(rows("rules_joinable_runs", rid)); model = runs(rows("rules_joinable_runs_model", rid))
    cone3 = runs(rows("rules_joinable_runs3d", rid)); model3 = runs(rows("rules_joinable_runs3d_model", rid))
    for g in (cone, model, cone3, model3):     # each directory separately (3D n=4 and n=7 runs must not be merged)
        for k, rs in g.items():
            b0 = int(rs[0]["beta1_multi_global"]); ev = 0; ok = True
            for r in rs[1:]:
                ev += int(r["events"]); ok &= int(r["beta1_multi_global"]) == b0 + d * ev
            cor_tot += 1; cor_ok += ok
    s2 = [rs for k, rs in cone.items() if k[0] == "grid2d" and k[1] == "single"]
    c2 = [rs for k, rs in cone.items() if k[0] == "grid2d" and k[1] == "collision"]
    m2 = [rs for k, rs in model.items() if k[0] == "grid2d" and k[1] == "collision"]
    info["n_single2d"] = len(s2)
    info["D_nonplanar_all"] = all(all(r["D1_planar"] == "0" for r in rs) for rs in s2) if len(s2) == 4 else None
    info["cone_sigma_all"] = all(all(r["Q_sum_conserved"] == "1" for r in rs) for rs in c2) if len(c2) == 4 else None
    info["model_sigma_all"] = all(all(r["Q_sum_conserved"] == "1" for r in rs) for rs in m2) if len(m2) == 4 else None
    info["model_iota0"] = all(all(int(r["iota_union_minus_sum"]) == 0 for r in rs) for rs in m2) if len(m2) == 4 else None
    per[rid] = info

def cnt(sel, key):
    have = [i for i in sel if i[key] is not None]
    return "%d/%d" % (sum(1 for i in have if i[key]), len(have)) if have else "--"
order = ["protected_frozen", "protected_cycle_growing", "dissolving"]
with open(os.path.join(gen, "tab_joinable.tex"), "w") as f:
    f.write("\\begin{tabular}{lrcccc}\n\\toprule\nclass & rules & $D$ non-planar $\\forall t$ (2D, cone) & $\\Sigma$ conserved (cone) & $\\Sigma$ conserved (model) & $\\iota_t\\equiv0$ (model)\\\\\n\\midrule\n")
    for b in order:
        sel = [i for i in per.values() if i["bucket"] == b]
        f.write("%s & %d & %s & %s & %s & %s\\\\\n" % (b.replace("_", "\\_"), len(sel), cnt(sel, "D_nonplanar_all"), cnt(sel, "cone_sigma_all"), cnt(sel, "model_sigma_all"), cnt(sel, "model_iota0")))
    f.write("\\bottomrule\n\\end{tabular}\n")
open(os.path.join(gen, "joinable_stats.tex"), "w").write("%d of %d evolutions" % (cor_ok, cor_tot))
nrules = Counter(i["bucket"] for i in per.values())
print(dict(nrules)); print("Cor1", cor_ok, cor_tot)
print(open(os.path.join(gen, "tab_joinable.tex")).read())
# consistency flags
bad = [r for r, i in per.items() if i["bucket"] != "dissolving" and i["D_nonplanar_all"] is False]
print("minor-preserving rules whose cone-labelled D became planar:", bad)
