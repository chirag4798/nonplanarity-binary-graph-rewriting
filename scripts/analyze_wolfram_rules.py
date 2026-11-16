#!/usr/bin/env python3
"""Tables for the section 'Relation to the Wolfram model' from results/wolfram_rules (lattice runs + self-loop evolutions)."""
import csv, os, sys
from collections import defaultdict
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(root, "code"))
import oligon_experiment as oe
W = os.path.join(root, "results", "wolfram_rules"); gen = os.path.join(root, "paper", "generated")
rules = [l.strip().split("|") for l in open(os.path.join(W, "rules.txt")) if l.strip()]
src = {"W1": "Wolfram Cloud working notebooks; visual summary", "W2": "Gorard (2020a), Fig.~14", "W4": "used in several Gorard papers (web record)"}
disp = {"W1": r"$\{\{x,y\},\{x,z\}\}\to\{\{x,z\},\{x,w\},\{y,w\},\{z,w\}\}$",
        "W2": r"$\{\{x,y\},\{x,z\}\}\to\{\{x,y\},\{x,w\},\{y,w\},\{z,w\}\}$",
        "W4": r"$\{\{x,y\},\{y,z\}\}\to\{\{w,y\},\{y,z\},\{z,w\},\{x,w\}\}$"}

def runs(mode, rid):
    g = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(W, "%s-%s" % (rid, mode), "evolution.csv"))):
        g[(r["background"], r["experiment"], r["defect"], r["order"])].append(r)
    return g

with open(os.path.join(gen, "tab_wolfram_rules.tex"), "w") as f:
    f.write("\\begin{tabular}{lccccc}\n\\toprule\nrule & $\\delta_R$ & $t_{\\max}$ (2D / 3D) & $D$ non-planar $\\forall t$ (cone; 2D, 3D) & $\\Sigma$ conserved (cone; 2D, 3D) & $\\Sigma$ conserved, $\\iota\\equiv0$ (model; 2D, 3D)\\\\\n\\midrule\n")
    law_tot = law_ok = 0
    for rid, s in rules:
        r = oe.parse_rule(s); d = oe.delta_R(r)
        cone, model = runs("cone", rid), runs("model", rid)
        tm = {bg: [int(rs[-1]["t"]) for k, rs in cone.items() if k[0] == bg] for bg in ("grid2d", "grid3d")}
        def cnt(g, pred, bg, exp):
            sel = [rs for k, rs in g.items() if k[0] == bg and k[1] == exp]
            return sum(1 for rs in sel if pred(rs)), len(sel)
        s2 = cnt(cone, lambda rs: all(x["D1_planar"] == "0" for x in rs), "grid2d", "single")
        s3 = cnt(cone, lambda rs: all(x["D1_planar"] == "0" for x in rs), "grid3d", "single")
        c2 = cnt(cone, lambda rs: all(x["Q_sum_conserved"] == "1" for x in rs), "grid2d", "collision")
        c3 = cnt(cone, lambda rs: all(x["Q_sum_conserved"] == "1" for x in rs), "grid3d", "collision")
        m2 = cnt(model, lambda rs: all(x["Q_sum_conserved"] == "1" and int(x["iota_union_minus_sum"]) == 0 for x in rs), "grid2d", "collision")
        m3 = cnt(model, lambda rs: all(x["Q_sum_conserved"] == "1" and int(x["iota_union_minus_sum"]) == 0 for x in rs), "grid3d", "collision")
        for g in (cone, model):
            for k, rs in g.items():
                b0 = int(rs[0]["beta1_multi_global"]); ev = 0; ok = True
                for x in rs[1:]:
                    ev += int(x["events"]); ok &= int(x["beta1_multi_global"]) == b0 + d * ev
                law_tot += 1; law_ok += ok
        f.write("%s & %d & %d--%d / %d--%d & %d/%d, %d/%d & %d/%d, %d/%d & %d/%d, %d/%d\\\\\n" % (
            disp[rid], d, min(tm["grid2d"]), max(tm["grid2d"]), min(tm["grid3d"]), max(tm["grid3d"]),
            s2[0], s2[1], s3[0], s3[1], c2[0], c2[1], c3[0], c3[1], m2[0], m2[1], m3[0], m3[1]))
    f.write("\\bottomrule\n\\end{tabular}\n")
open(os.path.join(gen, "wolfram_law.tex"), "w").write("%d of %d" % (law_ok, law_tot))
open(os.path.join(gen, "macros.tex"), "w").write("\\newcommand{\\wolframlaw}{%d of %d}\n" % (law_ok, law_tot))

rows = list(csv.DictReader(open(os.path.join(W, "selfloop_evolution.csv"))))
g = defaultdict(list)
for r in rows: g[(r["rule"], r["order"])].append(r)
with open(os.path.join(gen, "tab_selfloop.tex"), "w") as f:
    f.write("\\begin{tabular}{llrrrrc}\n\\toprule\nrule & order & generations & $|E|$ at end & first non-planar generation & non-injective events & cycle-rank law\\\\\n\\midrule\n")
    for (rid, order), v in g.items():
        first = [int(x["t"]) for x in v if x["planar_simple"] == "0"]
        f.write("%s & %s & %s & %s & %s & %d & %s\\\\\n" % (rid, order, v[-1]["t"], v[-1]["E"], first[0] if first else "none",
                sum(int(x["noninj_events"]) for x in v), "holds" if all(x["law_ok"] == "1" for x in v) else "VIOLATED"))
    f.write("\\bottomrule\n\\end{tabular}\n")
open(os.path.join(gen, "noninjective_check.txt"), "w").write(open(os.path.join(W, "noninjective_check.txt")).read())
print(open(os.path.join(gen, "tab_wolfram_rules.tex")).read()); print(open(os.path.join(gen, "wolfram_law.tex")).read())
