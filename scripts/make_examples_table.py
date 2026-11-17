#!/usr/bin/env python3
"""Writes paper/generated/tab_examples.tex: five running-example rules and which of the notions of the paper each satisfies.
Every entry is computed by the code of the archive (oligon_experiment.py, planarity_preserving.py); nothing is typed by hand."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "code")); sys.path.insert(0, HERE)
import oligon_experiment as oe
import planarity_preserving as pp

EX = [("E1", "{{1,2},{2,3}}->{{1,4},{4,2},{2,3}}", "subdivides an edge"),
      ("E2", "{{1,2},{1,3}}->{{1,4},{2,4},{3,4}}", "vertex expansion (star rule)"),
      ("E3", "{{1,2},{1,2}}->{{1,2},{1,3},{1,3}}", "parallel pair, one edge kept"),
      ("E4", "{{1,2},{2,3}}->{{1,2},{1,3},{2,3}}", "closes a triangle"),
      ("E5", "{{1,2},{1,3}}->{{1,2},{2,3},{2,4}}", r"the rule $R_c$ of Proposition~\ref{prop:counter}")]
yes, no = r"\checkmark", "--"
out = [r"\begin{tabular}{llllcccccc}", r"\toprule",
       r"id & rule & description & $\delta_R$ & path & minor & simple-minor & planarity-pres. & guaranteed in $H_t$\\", r"\midrule"]
for name, s, d in EX:
    r = oe.parse_rule(s)
    path, mn, sm, P = oe.is_path_preserving(r), oe.is_minor_preserving(r), oe.is_simple_minor_preserving(r), pp.satisfies_P(r)
    g = ("a subdivision of $H_0$" if path else "$H_0$ as a minor" if mn else "simple graph of $H_0$ as a minor" if sm else "no guarantee")
    rs = "$" + s.replace("->", r"\to").replace("{", r"\{").replace("}", r"\}") + "$"
    out.append("%s & %s & %s & %d & %s & %s & %s & %s & %s\\\\" % (name, rs, d, oe.delta_R(r), *(yes if x else no for x in (path, mn, sm, P)), g))
out += [r"\bottomrule", r"\end{tabular}"]
open(os.path.join(HERE, "..", "paper", "generated", "tab_examples.tex"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
