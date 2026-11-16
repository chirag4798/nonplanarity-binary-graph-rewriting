#!/bin/sh
# Lattice experiments for the three Wolfram-model binary rules of Section "Relation to the Wolfram model".
cd "$(dirname "$0")/../code"
for spec in "W1|{{1,2},{1,3}}->{{1,3},{1,4},{2,4},{3,4}}" "W2|{{1,2},{1,3}}->{{1,2},{1,4},{2,4},{3,4}}" "W4|{{1,2},{2,3}}->{{4,2},{2,3},{3,4},{1,4}}"; do
  id=${spec%%|*}; rule=${spec#*|}
  for mode in cone model; do
    python3 oligon_experiment.py run --out ../results/wolfram_rules/$id-$mode --rule "$rule" --tracking $mode --tmax 10 --n2d 8 --n3d 7 --max-edges 15000 --orc-sample 60 --no-dim > ../results/wolfram_rules/$id-$mode.log 2>&1 &
  done
done
wait
echo FINISHED > ../results/wolfram_rules/DONE
