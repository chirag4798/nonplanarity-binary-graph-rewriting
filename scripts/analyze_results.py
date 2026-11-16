#!/usr/bin/env python3
"""Aggregate results/ into LaTeX tables (paper/generated/*.tex) and figures (paper/figs/*.pdf).
Every number in the generated tables is computed from the CSV files; nothing is typed by hand."""
import csv, os, math
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
res = os.path.join(root, "results")
gen = os.path.join(root, "paper", "generated"); figs = os.path.join(root, "paper", "figs")
os.makedirs(gen, exist_ok=True); os.makedirs(figs, exist_ok=True)

sel = [l.strip().split("|") for l in open(os.path.join(res, "rules_selected.txt")) if l.strip()]
ids = [i for i, _ in sel]; rule_of = dict(sel)
cls = {r["rule"]: r for r in csv.DictReader(open(os.path.join(res, "classification", "classification.csv")))}


def readcsv(*p):
    f = os.path.join(res, *p)
    return {r["rule"]: r for r in csv.DictReader(open(f))} if os.path.exists(f) else {}


cen = readcsv("census", "rule_census.csv")
fl = lambda x: float(x) if x not in ("", None) else float("nan")


def esc(s):
    return s.replace("{", r"\{").replace("}", r"\}").replace("->", r"$\to$").replace("_", r"\_")


def tt(rule):
    return r"\texttt{" + esc(rule).replace(" ", "") + "}"


def bucket_short(rule):
    return {"protected_frozen": "frozen", "protected_cycle_growing": "growing", "dissolving": "dissolving"}.get(cls[rule]["bucket"], cls[rule]["bucket"])


def build(d2, d3):
    R = defaultdict(list)
    for dname, tag3 in ((d2, "grid3d_n4"), (d3, "grid3d_n7")):
        for rid in ids:
            f = os.path.join(res, dname, rid, "evolution.csv")
            if not os.path.exists(f):
                continue
            for r in csv.DictReader(open(f)):
                if dname == d3 and r["background"] != "grid3d":
                    continue
                tag = "grid2d_n8" if r["background"] == "grid2d" else tag3
                R[(rid, tag, r["experiment"], r["defect"], r["order"])].append(r)
    return dict(R)


RC = build("runs", "runs3d")                # causal-cone labelling
RM = build("runs_model", "runs3d_model")    # minor-model label transport
delta = lambda rid: int(cls[rule_of[rid]]["delta_R"])


def first_t(rows, col, val):
    for r in rows:
        if r[col] != "" and int(r[col]) == val:
            return int(r["t"])
    return None


# ------------------------------------------------------------ rules table
with open(os.path.join(gen, "tab_rules.tex"), "w") as f:
    f.write("\\begin{tabular}{llcccl}\n\\toprule\nID & rule & $\\delta_R$ & path-pres. & minor-pres. & class / joinability\\\\\n\\midrule\n")
    for i in ids:
        r = rule_of[i]; c = cls[r]; ce = cen.get(r, {})
        f.write("%s & %s & %s & %s & %s & %s / %s\\\\\n" % (i, tt(r), c["delta_R"], "yes" if c["path_preserving"] == "1" else "no",
                "yes" if c["minor_preserving"] == "1" else "no", bucket_short(r), esc(ce.get("conf_status", "n/a"))))
    f.write("\\bottomrule\n\\end{tabular}\n")

# ------------------------------------------------------------ Corollary 1 check (all runs, both labellings)
tot = ok_ = 0; bad = []
for R in (RC, RM):
    for k, rows in R.items():
        b0 = int(rows[0]["beta1_multi_global"]); ev = 0; good = True
        for r in rows[1:]:
            ev += int(r["events"])
            good &= int(r["beta1_multi_global"]) == b0 + delta(k[0]) * ev
        tot += 1; ok_ += good
        if not good:
            bad.append(k)
open(os.path.join(gen, "cor1_check.tex"), "w").write("%d of %d evolutions" % (ok_, tot) + ("" if not bad else "; violations: " + esc(str(bad[:3]))))
print("Cor1", ok_, tot)


# ------------------------------------------------------------ single-defect table
def single(R, rids, fname):
    def summ(tags):
        o = {}
        for rid in rids:
            n = Dnp = gn = gnp = 0; fp = []
            for tag in tags:
                for d in ("K5", "K33"):
                    for ordr in ("random", "canonical"):
                        rows = R.get((rid, tag, "single", d, ordr))
                        if not rows:
                            continue
                        n += 1
                        t1 = first_t(rows, "D1_planar", 1)
                        Dnp += t1 is None
                        if t1 is not None:
                            fp.append(t1)
                        if rows[0]["planar_simple_global"] != "":
                            gn += 1; gnp += first_t(rows, "planar_simple_global", 1) is None
            o[rid] = (n, Dnp, fp, gn, gnp)
        return o
    s2 = summ(["grid2d_n8"]); s3 = summ(["grid3d_n4", "grid3d_n7"])
    with open(os.path.join(gen, fname), "w") as f:
        f.write("\\begin{tabular}{lcccccc}\n\\toprule\n & & \\multicolumn{3}{c}{2D grid ($n=8$)} & \\multicolumn{2}{c}{3D grids ($n=4,7$)}\\\\\n"
                "\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}\nID & $\\delta_R$ & $D$ non-planar $\\forall t$ & first planar $t$ & global non-planar $\\forall t$ & $D$ non-planar $\\forall t$ & first planar $t$\\\\\n\\midrule\n")
        for rid in rids:
            n, Dnp, fp, gn, gnp = s2[rid]; n3, Dnp3, fp3, _, _ = s3[rid]
            f.write("%s & %d & %d/%d & %s & %d/%d & %d/%d & %s\\\\\n" % (rid, delta(rid), Dnp, n, ("%d" % min(fp) if fp else "--"), gnp, gn, Dnp3, n3, ("%d" % min(fp3) if fp3 else "--")))
        f.write("\\bottomrule\n\\end{tabular}\n")


def collision(R, rids, fname):
    def coll(rid, tag):
        n = ok = 0; iota = []; ds = []; d1 = []
        for d in ("K5", "K33"):
            for ordr in ("random", "canonical"):
                rows = R.get((rid, tag, "collision", d, ordr))
                if not rows:
                    continue
                n += 1
                ok += all(r["Q_sum_conserved"] == "1" for r in rows)
                iota += [int(r["iota_union_minus_sum"]) for r in rows]
                ds.append(int(rows[-1]["beta1_sum"]) - int(rows[0]["beta1_sum"]))
        return n, ok, iota, ds
    with open(os.path.join(gen, fname), "w") as f:
        f.write("\\begin{tabular}{lcccccc}\n\\toprule\n & & \\multicolumn{3}{c}{2D grid ($n=8$)} & \\multicolumn{2}{c}{3D grid ($n=7$)}\\\\\n"
                "\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}\nID & $\\delta_R$ & $\\Sigma$ conserved $\\forall t$ & $\\iota_t$ range & $\\Delta\\Sigma$ at end & $\\Sigma$ conserved $\\forall t$ & $\\Delta\\Sigma$ at end\\\\\n\\midrule\n")
        for rid in rids:
            n, ok, io, ds = coll(rid, "grid2d_n8"); n3, ok3, io3, ds3 = coll(rid, "grid3d_n7")
            rg = lambda x: ("%d..%d" % (min(x), max(x))) if x else "--"
            f.write("%s & %d & %d/%d & %s & %s & %d/%d & %s\\\\\n" % (rid, delta(rid), ok, n, rg(io), rg(ds), ok3, n3, rg(ds3)))
        f.write("\\bottomrule\n\\end{tabular}\n")


def curv(R, rids, fname):
    with open(os.path.join(gen, fname), "w") as f:
        f.write("\\begin{tabular}{lccccc}\n\\toprule\nID & $|E|$ at end & growth exp.\\ (control) & growth exp.\\ (with $K_5$) & $\\kappa$ defect, $t{=}0$ & $\\kappa$ defect, end\\\\\n\\midrule\n")
        for rid in rids:
            c = R.get((rid, "grid2d_n8", "control", "none", "random")); s = R.get((rid, "grid2d_n8", "single", "K5", "random"))
            if not c or not s:
                continue
            e = fl(s[-1]["orc_defect_mean"])
            f.write("%s & %d & %.2f & %.2f & %.3f & %s\\\\\n" % (rid, int(s[-1]["E"]), fl(c[-1]["dim_est"]), fl(s[-1]["dim_est"]), fl(s[0]["orc_defect_mean"]), "%.3f" % e if not math.isnan(e) else "--"))
        f.write("\\bottomrule\n\\end{tabular}\n")
    c0 = R.get((rids[0], "grid2d_n8", "control", "none", "random"))
    return fl(c0[0]["orc_global_mean"]) if c0 else float("nan")


mrids = [i for i in ids if i[0] in "FG"]
single(RC, ids, "tab_single.tex"); single(RM, mrids, "tab_single_model.tex")
collision(RC, ids, "tab_collision.tex"); collision(RM, mrids, "tab_collision_model.tex")
g0 = curv(RC, ids, "tab_curv.tex")
open(os.path.join(gen, "orc_global0.tex"), "w").write("%.3f" % g0)

# ------------------------------------------------------------ classification & census summaries
counts = list(csv.DictReader(open(os.path.join(res, "classification", "classification_counts.csv"))))
with open(os.path.join(gen, "tab_class.tex"), "w") as f:
    f.write("\\begin{tabular}{llrr}\n\\toprule\nBucket & $\\delta_R$ & $2_2\\to3_2$ & $2_2\\to4_2$\\\\\n\\midrule\n")
    keys = sorted({(r["bucket"], r["delta_R"]) for r in counts}, key=lambda x: (x[0] == "open_rule", x[0], x[1]))
    for b, d in keys:
        v = {r["signature"]: r["count"] for r in counts if r["bucket"] == b and r["delta_R"] == d}
        f.write("%s & %s & %s & %s\\\\\n" % (esc(b), d if d != "n/a" else "--", v.get("2_3", "0"), v.get("2_4", "0")))
    f.write("\\bottomrule\n\\end{tabular}\n")


def census_table():
    cl = {}
    for sig, path in (("2_3", "census"), ("2_4", "census_2_4")):
        rows = list(csv.DictReader(open(os.path.join(res, path, "rule_census.csv")))) if os.path.exists(os.path.join(res, path, "rule_census.csv")) else []
        c = Counter()
        for r in rows:
            c[(cls[r["rule"]]["bucket"] if r["rule"] in cls else "?", r["conf_status"])] += 1
        cl[sig] = (len(rows), c)
    stats = ["joinable_all_tested", "unjoined_within_J", "inconclusive"]
    with open(os.path.join(gen, "tab_census.tex"), "w") as f:
        f.write("\\begin{tabular}{llrrr}\n\\toprule\nsignature & class & joinable (all tested) & unjoined within $J$ & inconclusive\\\\\n\\midrule\n")
        for sig, (n, c) in cl.items():
            for b in ["protected_frozen", "protected_cycle_growing", "dissolving", "inapplicable_to_simple_H0"]:
                if sum(c[(b, s)] for s in stats) == 0:
                    continue
                f.write("%s & %s & %d & %d & %d\\\\\n" % ("$2_2\\to %s_2$" % sig[-1], esc(b), *[c[(b, s)] for s in stats]))
        f.write("\\bottomrule\n\\end{tabular}\n")
census_table()

# ------------------------------------------------------------ figures
fig_ids = [i for i in ("F1", "F3", "G1", "D1", "D4", "RC") if i in ids]
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
for rid in fig_ids:
    rows = RC.get((rid, "grid2d_n8", "single", "K5", "canonical"))
    if not rows:
        continue
    t = [int(r["t"]) for r in rows]
    ax[0].plot(t, [int(r["beta1_multi_global"]) for r in rows], marker="o", ms=3, label=rid)
    ax[1].plot(t, [int(r["beta1_D1"]) for r in rows], marker="o", ms=3, label=rid)
    ax[2].plot(t, [1 - int(r["D1_planar"]) for r in rows], marker="o", ms=3, label=rid)
ax[0].set_title("global $\\beta_1(H_t)$ (log)"); ax[0].set_yscale("log")
ax[1].set_title("cone-labelled $\\beta_1(D_1(t))$"); ax[2].set_title("$D_1(t)$ non-planar (1/0)")
for a in ax: a.set_xlabel("generation $t$")
ax[0].legend(fontsize=7, ncol=2); plt.tight_layout(); plt.savefig(os.path.join(figs, "fig_beta1.pdf")); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.4))
for rid in fig_ids:
    rows = RC.get((rid, "grid2d_n8", "single", "K5", "canonical"))
    if not rows:
        continue
    t = [int(r["t"]) for r in rows]
    ax[0].plot(t, [fl(r["orc_global_mean"]) for r in rows], marker="o", ms=3, label=rid)
    ax[1].plot(t, [fl(r["orc_defect_mean"]) for r in rows], marker="o", ms=3, label=rid)
ax[0].set_title("mean $\\kappa_{1/2}$, sampled edges"); ax[1].set_title("mean $\\kappa_{1/2}$, defect-labelled edges")
for a in ax: a.set_xlabel("generation $t$")
ax[0].legend(fontsize=7, ncol=2); plt.tight_layout(); plt.savefig(os.path.join(figs, "fig_curv.pdf")); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.4), sharey=True)
for rid in [i for i in ("F1", "F3", "G1", "G3") if i in ids]:
    for a, R, ttl in ((ax[0], RC, "causal-cone labelling"), (ax[1], RM, "minor-model transport")):
        rows = R.get((rid, "grid2d_n8", "collision", "K5", "canonical"))
        if rows:
            a.plot([int(r["t"]) for r in rows], [int(r["beta1_sum"]) for r in rows], marker="o", ms=3, label=rid)
        a.set_title("$\\beta_1(D_1)+\\beta_1(D_2)$, " + ttl, fontsize=9); a.set_xlabel("generation $t$")
ax[0].set_yscale("symlog"); ax[0].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(figs, "fig_collision.pdf")); plt.close()

# ------------------------------------------------------------ necessity certificates
nec = list(csv.DictReader(open(os.path.join(res, "necessity", "necessity_summary.csv"))))
vlog = open(os.path.join(res, "necessity", "verify_log.txt")).read()
with open(os.path.join(gen, "tab_necessity.tex"), "w") as f:
    f.write("\\begin{tabular}{lrrrr}\n\\toprule\nsignature & closed rules & simple-minor-preserving & minor-preserving (multigraph) & not simple-minor, certificate on simple $K_5$\\\\\n\\midrule\n")
    for sig in ("2_3", "2_4"):
        rr = [r for r in nec if r["signature"] == sig]
        f.write("$2_2\\to %s_2$ & %d & %d & %d & %d\\\\\n" % (sig[-1], len(rr), sum(r["simple_minor_preserving"] == "1" for r in rr),
                sum(r["minor_preserving"] == "1" for r in rr), sum(r["simple_minor_preserving"] == "0" and r["certificate"] == "A5" for r in rr)))
    f.write("\\bottomrule\n\\end{tabular}\n")
open(os.path.join(gen, "nec_stats.tex"), "w").write(
    "%d certificates, %d verified" % tuple(int(x) for x in __import__("re").search(r"certificates: (\d+) verified: (\d+)", vlog).groups()))
open(os.path.join(gen, "nec_gaps.tex"), "w").write("%d" % len([l for l in open(os.path.join(res, "necessity", "gaps.txt")).read().split("\n") if l]))
nviol = sum(int(r["bad"]) for r in nec if r["simple_minor_preserving"] == "1")
open(os.path.join(gen, "nec_viol.tex"), "w").write("%d" % nviol)

# ------------------------------------------------------------ control runs: when does the planar 2D grid become non-planar by itself?
fn = {}
for rid in ids:
    rows = RC.get((rid, "grid2d_n8", "control", "none", "canonical"))
    t = [int(r["t"]) for r in rows if r["planar_simple_global"] == "0"] if rows else []
    fn[rid] = min(t) if t else None
open(os.path.join(gen, "control_nonplanar.tex"), "w").write(", ".join("%s: %s" % (k, v if v is not None else "never") for k, v in fn.items()))

# ------------------------------------------------------------ general-signature search
import sys as _sys
_sys.path.insert(0, os.path.join(root, "code"))
import oligon_experiment as _oe
import networkx as _nx
gen_rows = []
for d in ("general_1", "general_3", "general_2k"):
    f = os.path.join(res, d, "general_signature_summary.csv")
    if os.path.exists(f):
        gen_rows += list(csv.DictReader(open(f)))

def explain(rule_s):
    r = _oe.parse_rule(rule_s)
    if not _oe.is_planar(_oe.simple_graph(list(r.K))):
        return "kuratowski_gadget"
    Lp = {frozenset(e) for e in r.L}
    G = _nx.Graph(); G.add_edges_from(r.K)
    star = G.number_of_nodes() == 4 and sorted(d for _, d in G.degree()) == [1, 1, 1, 3]
    if len(Lp) == 3 and r.n_old == 3 and star:
        return "delta_to_Y"
    return "unexplained"

order = []
for r in gen_rows:
    if r["sig"] not in order: order.append(r["sig"])
with open(os.path.join(gen, "tab_general.tex"), "w") as f:
    f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\nsignature $m_2\\to k_2$ ($n_{\\rm new}\\le$) & closed rules & simple-minor-pres. & not s.m., certificate & not s.m., no certificate: $\\Delta\\!\\to\\!Y$ & Kuratowski gadget & unexplained\\\\\n\\midrule\n")
    for sg in order:
        rr = [r for r in gen_rows if r["sig"] == sg]; m, k, mx = sg.split(":")
        gaps = [r for r in rr if r["sm"] == "0" and not r["cert"]]
        ex = Counter(explain(r["rule"]) for r in gaps)
        f.write("$%s_2\\to%s_2$ (%s) & %d & %d & %d & %d & %d & %d\\\\\n" % (m, k, mx, len(rr), sum(r["sm"] == "1" for r in rr), sum(1 for r in rr if r["sm"] == "0" and r["cert"]),
                ex["delta_to_Y"], ex["kuratowski_gadget"], ex["unexplained"]))
    f.write("\\bottomrule\n\\end{tabular}\n")
open(os.path.join(gen, "general_bad.tex"), "w").write("%d" % sum(int(r["bad"]) for r in gen_rows))

# ------------------------------------------------------------ audit appendix: enumeration table and certificate example
import json as _json
cnt = defaultdict(lambda: [0, 0, 0])
for r in nec:
    rr = _oe.parse_rule(r["rule"])
    c = cnt[(r["signature"], rr.n_old, rr.n_new)]
    c[0] += 1; c[1] += r["simple_minor_preserving"] == "1"; c[2] += r["simple_minor_preserving"] == "0"
with open(os.path.join(gen, "tab_enum.tex"), "w") as f:
    f.write("\\begin{tabular}{lcrrr}\n\\toprule\nsignature & $(n_{\\rm old},n_{\\rm new})$ & closed rules & simple-minor-preserving & not (certificate)\\\\\n\\midrule\n")
    for (sig, no, nn), (a, b, c) in sorted(cnt.items()):
        f.write("$2_2\\to%s_2$ & (%d,%d) & %d & %d & %d\\\\\n" % (sig[-1], no, nn, a, b, c))
    f.write("\\bottomrule\n\\end{tabular}\n")
first = _json.loads(open(os.path.join(res, "necessity", "certificates.jsonl")).readline())
import textwrap as _tw
open(os.path.join(gen, "cert_example.txt"), "w").write("\n".join(_tw.wrap(_json.dumps(first), 105)) + "\n")
print("done")
