#!/bin/sh
# Long job: all joinable (bounded filter) 2_2->3_2 rules that apply to simple H0. Estimated several hours.
cd "$(dirname "$0")/.."
for m in main 3d model 3dmodel; do python3 scripts/run_selected.py 5 $m rules_joinable.txt; done > results/run_joinable.log 2>&1
echo FINISHED >> results/run_joinable.log
