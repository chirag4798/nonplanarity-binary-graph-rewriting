# Literature search log (16 November 2026)

This is the query log behind Appendix C of the paper (protocol and inclusion criteria are described there). A negative result means "not found by these queries". MathSciNet, Scopus and Web of Science were not searched (institutional access was not available); zbMATH Open, OpenAlex, Crossref and arXiv were searched through their public APIs, and a forward citation search was made (see the sections below).

| # | Query | Outcome |
|---|---|---|
| 1 | web: oligons, particles, topological defects, Wolfram model | Wolfram essays and technical introduction; Winter School 2021 $K_5$-in-tiling project (repository abstract); Gorard's tangle remarks. Cited. |
| 2 | Gorard (2020b), full text, tangles | Read via HTML rendering: §2.4 conditional heuristic on planarity-preserving rules, annihilation, charge analogy; Noether analogue deferred; "oligon" not used in that paper. |
| 3 | web: causal invariance, confluence, critical pairs | Gorard 2020 (arXiv:2011.12174); Grossberndt notebook (experimental prover). Cited. |
| 4 | arXiv API: "Wolfram model" with planar/defect/particle/topological | 2 records, all by Gorard on relativity/numerical GR; none on non-planar structures, defects or invariants. |
| 5 | arXiv API: "hypergraph rewriting" with topological/homology/planarity/invariant | 3 records (Gorard); none relevant to (a)–(c) of our criteria. |
| 6 | INSPIRE-HEP: "hypergraph rewriting" or "Wolfram model" in title/abstract | 8 records; closest: Arsiwalla–Gorard homotopy types (published Int. J. Theor. Phys. 63(4), 2024); none on defects, non-planarity or Betti numbers by title. Confirmed DOI of Gorard (2020b). |
| 7 | INSPIRE-HEP: "non-planar" and "graph rewriting" in abstract | 0 records. |
| 8 | zbMATH Open API, Semantic Scholar | Failed (404/403/429); no records. |
| 9 | Lin–Lu–Yau; Bourne et al.; Wagner (1937); Ehrig et al.; Plump | Verified by web records (Ehrig book authors are Ehrig, Ehrig, Prange, Taentzer; Plump 2005 start page 280, end page unconfirmed; a second Plump title could not be verified and is not cited). |
| 10 | Bilson-Thompson–Markopoulou–Smolin; Chester–Arsiwalla–Kauffman | Verified (Class. Quantum Grav. 24, 3975; arXiv:2501.03260). |
| 11 | web: Betti numbers as conserved quantity in hypergraph rewriting | No source found. |
| 12 | Patarroyo post; Summer School 2020 stable-defect project; technical-introduction page | Post exists, content not read; 2020 project not found; page listed by search but not fetchable (DNS failure) and not read. |

## zbMATH Open queries (17 November 2026)

The public zbMATH Open API (`https://api.zbmath.org`) was queried with `scripts/zbmath_search.py`; the full output, with the first eight hits per query, is in `results/literature/zbmath_queries.txt`. Queries: "graph transformation planarity preserving"; "graph grammar planar graphs"; "double pushout minor" (request failed); "hyperedge replacement planar"; "graph rewriting non-planar persistence" (no result); "graph rewriting" planar; "graph transformation" planar; "double pushout"; "graph grammar" planarity preserving (no result); "graph rewriting" minor; "hypergraph rewriting" Wolfram (no result); "graph rewriting" "Betti" (no result); "graph transformation" "graph minor".

Items cited from these records (bibliographic record, and for Habel-Muller-Plump the zbMATH review, read; full texts not read): Habel, Muller and Plump (2001), Joshi and Konig (2008), Montanari (1970), Pascu and Paun (1979). No record concerned persistence of non-planarity or Betti numbers under rewriting. MathSciNet and Scopus were not searched.


## Multi-source search and citation chain (18 November 2026)

Script: `scripts/literature_search.py` (modes `search`, `arxiv`, `chain`, `targeted`). Outputs in `results/literature/`:

| File | Content |
|---|---|
| `multisource_queries.txt` | 20 queries on graph rewriting/transformation with planarity, minors, Betti numbers, Wolfram model; first 8 hits from OpenAlex and Crossref. The Semantic Scholar column is empty (rate-limited, HTTP 429); the arXiv column of this file is unreliable (empty responses) and was redone in `arxiv_queries.txt`. |
| `arxiv_queries.txt` | the same 20 queries as conjunctive keyword queries on arXiv (0 to 7 records each; none relevant). |
| `citation_chain.txt` | forward citation search in OpenAlex from 8 seed works, citing titles screened by keywords (at most the first 500 citing works per seed). |
| `citation_chain_targeted.txt` | OpenAlex `cites:` filter combined with a title-and-abstract search on ten topic terms, for four seeds. |

Outcome: no record concerns persistence of non-planarity or Betti numbers under graph rewriting. Items added as related work, after confirming title, authors and venue in Crossref or the arXiv API: Llabres and Rossello (1997), Pascual et al. (2022), Arrighi and Dowek (2012). Records seen and not cited: Montanari/Courcelle-style graph-grammar papers beyond Montanari (1970), Pascual and Le Gall's earlier work, Gorard/Arsiwalla/Leuenberger papers not about defects or invariants.

Full texts of Habel-Muller-Plump (2001), Joshi-Konig (2008), Montanari (1970) and Pascu-Paun (1979) could not be retrieved (login walls, a 404 on the author copy, a bot check on the Duisburg-Essen repository). Google Scholar returned a bot check and was not used.
