#!/usr/bin/env python3
"""Run the evolution/collision experiments for each rule in results/rules_selected.txt (ID|rule), in parallel.
usage: run_selected.py JOBS MODE      MODE in {main, 3d, model, 3dmodel}
  main     cone labelling, 2D n=8 and 3D n=4                        -> results/runs
  3d       cone labelling, 3D n=7 only (room for collisions)         -> results/runs3d
  model    minor-model label transport (minor-preserving rules), 2D n=8 and 3D n=4 -> results/runs_model
  3dmodel  minor-model label transport, 3D n=7 only                  -> results/runs3d_model"""
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
code = os.path.join(root, "code")
jobs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
mode = sys.argv[2] if len(sys.argv) > 2 else "main"
rules_file = sys.argv[3] if len(sys.argv) > 3 else "rules_selected.txt"   # e.g. rules_joinable.txt (outputs go to results/<name>/ instead of results/runs*)
sub = "" if rules_file == "rules_selected.txt" else rules_file.replace(".txt", "") + "_"
cfg = {"main": ("runs", [], "4", "cone"), "3d": ("runs3d", ["--backgrounds", "grid3d"], "7", "cone"),
       "model": ("runs_model", [], "4", "model"), "3dmodel": ("runs3d_model", ["--backgrounds", "grid3d"], "7", "model")}[mode]
items = [l.strip().split("|") for l in open(os.path.join(root, "results", rules_file)) if l.strip()]
if cfg[3] == "model":
    sys.path.insert(0, code)
    import oligon_experiment as _oe
    items = [it for it in items if _oe.is_minor_preserving(_oe.parse_rule(it[1]))]   # minor-preserving rules only
def go(it):
    rid, rule = it
    out = os.path.join(root, "results", sub + cfg[0], rid)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out + ".log", "w") as log:
        return subprocess.call([sys.executable, "oligon_experiment.py", "run", "--out", out, "--rule", rule,
                                "--tmax", "10", "--n2d", "8", "--n3d", cfg[2], "--max-edges", "15000", "--orc-sample", "100",
                                "--tracking", cfg[3]] + cfg[1], cwd=code, stdout=log, stderr=subprocess.STDOUT)
with ThreadPoolExecutor(jobs) as ex:
    print(list(ex.map(go, items)))
