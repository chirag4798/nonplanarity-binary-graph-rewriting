# Oligon paper: manuscript, code, results

    paper/main.tex, references.bib   manuscript; paper/main.pdf is the compiled draft (tectonic)
    paper/generated, paper/figs      tables and figures generated from results/ by scripts/analyze_results.py
    code/oligon_experiment.py        experiments + selftest (networkx, numpy, scipy)
    code/oligon_wolframmodel.wl      UNTESTED WolframModel cross-check (no Wolfram Language available)
    scripts/classify_rules.py        exhaustive rule classification on isolated K5/K33
    scripts/run_selected.py          lattice evolution/collision runs (modes: main, 3d, model, 3dmodel)
    scripts/necessity_certificates.py   computer-assisted necessity check (certificates in results/necessity)
    scripts/verify_certificates.py      independent checker for those certificates
    scripts/run_all_joinable.sh         LONG job (hours): all 126 joinable 2_2->3_2 rules -> results/rules_joinable_runs*
    scripts/analyze_joinable.py         aggregates that long run (safe on partial data) -> paper/generated/tab_joinable.tex
    scripts/general_signature_search.py necessity search in other signatures (results/general_*)
    scripts/verify_independent.py       Burnside count + second decision procedure + semantic contraction test
    scripts/verify_counterexamples.py   checks the Delta->Y and Kuratowski-gadget counterexamples
    code/general_rules.py               m-edge left-side enumeration/matching (reproduces the 2-edge enumeration exactly)
    scripts/wolfram_bridge_study.py, run_wolfram_rules.sh, analyze_wolfram_rules.py, check_incidence_formula.py
                                        Section 'Relation to the Wolfram model' (code/wolfram_bridge.py: non-injective matching)
    scripts/analyze_results.py       CSV -> LaTeX tables and figures
    results/                         all raw CSV/log outputs

## Reproduce (seed 20240601)
    pip install -r code/requirements.txt
    cd code && python oligon_experiment.py selftest > ../results/selftest_log.txt
    python oligon_experiment.py census --out ../results/census --sigs 2_3 --jobs 6
    python oligon_experiment.py census --out ../results/census_2_4 --sigs 2_4 --jobs 4
    cd .. && python scripts/classify_rules.py --out results/classification --jobs 8
    for m in main 3d model 3dmodel; do python scripts/run_selected.py 4 $m; done
    python scripts/necessity_certificates.py --jobs 3 && python scripts/verify_certificates.py
    python scripts/analyze_results.py
    cd paper && tectonic main.tex     # set TECTONIC_CACHE_DIR to a writable dir if needed

The joinability filter is bounded, not a proof of confluence; see the paper's Limitations.

## Quick check of the headline claim (Theorem 8)
    python3 scripts/verify_certificates.py     # standard library only; reads results/necessity/*
See supplement/DEPOSIT.md for depositing the archive.
