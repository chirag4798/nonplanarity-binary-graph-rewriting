# Persistence of non-planarity under minimal binary graph rewriting: code, data and certificates

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24310939.svg)](https://doi.org/10.5281/zenodo.24310939) Archived release v1.0.0 on Zenodo (all versions: https://doi.org/10.5281/zenodo.24310938).

Supplementary material for the paper *Persistence of Non-Planarity under Minimal Binary Graph Rewriting: Cycle Rank, Minor Preservation,
and an Exhaustive Classification* (Chirag Shetty, 2026). The paper source and PDF are in `paper/`.

**Licences.** Code (`code/`, `scripts/`): MIT (`LICENSE`). Data, results, figures and text (`results/`, `paper/`): CC BY 4.0 (`LICENSE-DATA`).

## Check the central claim in one minute (Python 3 only, no packages)
    python3 scripts/verify_certificates.py
This verifies the 2388 certificates behind Theorem 8 of the paper (host, event, planar rotation system, and that the certified rules are
exactly the rules that are not simple-minor-preserving among all 5225), reading `results/necessity/certificates.jsonl` and `necessity_summary.csv`.

## Where everything is
| Item | In the paper | Files | Scripts (`scripts/` unless a path is given) |
|---|---|---|---|
| Rule enumeration and counts | App. A, Table 3 | `results/necessity/necessity_summary.csv` | `necessity_certificates.py`, `code/general_rules.py` |
| Certificates (Theorem 8) | Sec. 6, App. A | `results/necessity/certificates.jsonl`, `verify_log.txt` | `necessity_certificates.py`; checker `verify_certificates.py` |
| Independent checks, counterexamples | App. A | `results/selftest_log.txt`, `verify_independent_log.txt`, `counterexamples_log.txt` | `code/oligon_experiment.py selftest`, `verify_independent.py`, `verify_counterexamples.py` |
| Classification on isolated K5, K3,3 | App. B | `results/classification/` | `classify_rules.py` |
| Joinability census | App. B | `results/census/`, `results/census_2_4/` | `code/oligon_experiment.py census` |
| Search beyond the signatures of Theorem 8 | Sec. 6.1 | `results/general_1`, `general_2k`, `general_3` | `general_signature_search.py` |
| Lattice runs, all 126 joinable rules | Sec. 8 | `results/rules_joinable.txt`, `results/rules_joinable_runs*/` | `run_all_joinable.sh`, `run_selected.py`, `analyze_joinable.py` |
| Lattice runs, selected rules | App. B | `results/rules_selected.txt`, `results/runs*/` | `run_selected.py`, `analyze_results.py` |
| Wolfram-model rules, self-loop evolutions | Sec. 9 | `results/wolfram_rules/` | `run_wolfram_rules.sh`, `wolfram_bridge_study.py`, `analyze_wolfram_rules.py`, `code/wolfram_bridge.py` |
| Planarity-preserving rules (Thm 14, Table 2) | Sec. 6.4 | `results/planarity/` | `planarity_preserving.py` |
| Non-injective matches (Prop. 13 check) | Sec. 4 | `results/noninjective_minor_check.txt` | `check_noninjective_minor.py` |
| Literature search log | App. C | `LITERATURE_LOG.md` | – |
| Incidence-graph check | Sec. 9 | `results/incidence_check.txt` | `check_incidence_formula.py` |
| Tables and figures | all | `paper/generated/`, `paper/figs/` | `analyze_results.py`, `analyze_joinable.py`, `analyze_wolfram_rules.py` |

## Reproduce everything (seed 20240601; Python 3.13.5, networkx 3.6.1, numpy 2.3.3, scipy 1.16.3, matplotlib 3.11.2)
    pip install -r code/requirements.txt
    cd code && python oligon_experiment.py selftest > ../results/selftest_log.txt
    python oligon_experiment.py census --out ../results/census --sigs 2_3 --jobs 6
    python oligon_experiment.py census --out ../results/census_2_4 --sigs 2_4 --jobs 4
    cd .. && python scripts/classify_rules.py --out results/classification --jobs 8
    python scripts/necessity_certificates.py --jobs 3 && python scripts/verify_certificates.py
    python scripts/verify_independent.py && python scripts/verify_counterexamples.py
    python scripts/general_signature_search.py --sigs 1:2:1 1:3:2 1:4:2 1:5:3 --out results/general_1
    python scripts/general_signature_search.py --sigs 3:3:1 3:4:1 --out results/general_3
    python scripts/general_signature_search.py --sigs 2:5:3 2:6:2 2:7:1 --out results/general_2k
    for m in main 3d model 3dmodel; do python scripts/run_selected.py 4 $m; done      # selected rules
    scripts/run_all_joinable.sh                                                       # several hours
    scripts/run_wolfram_rules.sh && python scripts/wolfram_bridge_study.py
    python scripts/planarity_preserving.py          # ~2 min on 8 cores
    python scripts/check_noninjective_minor.py
    python scripts/analyze_results.py && python scripts/analyze_joinable.py && python scripts/analyze_wolfram_rules.py
The long jobs (all 126 rules, the 2->6 search) take hours on 8 cores. Absolute local paths were removed from run manifests and log files are not included.

## Archive
Each release is archived on Zenodo (see the DOI badge/link in the repository description once minted). `python3 scripts/make_supplement.py`
rebuilds the deterministic archive and `supplement/SHA256SUMS`.

## Citation
See `CITATION.cff`.
