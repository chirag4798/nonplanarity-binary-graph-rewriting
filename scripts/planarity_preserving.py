#!/usr/bin/env python3
"""
planarity_preserving.py -- which rules of S preserve planarity of the underlying simple graph?

Deterministic and exhaustive over S (5225 rules), bounded in the host family.

Condition P (purely combinatorial, on the rule only).  Let G = simple(K) on var(K), O = old variables,
N = new variables, pairs(L) = distinct unordered pairs of simple(L).  R satisfies P iff
  (P1) every edge of G between two old variables is a pair of pairs(L);
  (P2) every connected component C of G[N] has attachment set A(C) (old variables adjacent to C) with
       |A(C)| <= 2, and if |A(C)| = 2 then A(C) is a pair of pairs(L);
  (P3) for |A(C)| = 2, A(C) = {a,b}: G[C u {a,b}] + ab is planar (always true here, but checked).
(Proposition in the paper: P implies that planarity is preserved by every event, including non-injective matches.)

Certificate search (the converse, bounded).  For a rule violating P we look for a planar host and a match whose
output is non-planar.  Host: a planar simple graph Hs (a base graph from a fixed list, with the pairs of L added
if absent and the result planar), every edge present with 4 copies (both orientations twice) so that every L pattern
matches and removing two copies never removes an adjacency.  The output's simple graph is then Hs + image(K); since
dropping the L-edges can only delete edges, this is the most favourable case for non-planarity, and a non-planar
output is a certificate that the rule does not preserve planarity.  Rules violating P without certificate would be
counterexamples to the equivalence; the run reports them.

Also reported: P implies simple-minor-preserving (checked), and how P overlaps with the dichotomy of Theorem 8.
Usage: python3 planarity_preserving.py [--out ../results/planarity]
"""
import argparse, itertools, os, sys, json, csv
from multiprocessing import Pool
from collections import Counter
import networkx as nx

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe


def simple_of(edges):
    G = nx.Graph()
    G.add_edges_from((a, b) for a, b in edges if a != b)
    return G


def pairs_L(r):
    return {frozenset(e) for e in r.L if e[0] != e[1]}


def satisfies_P(r):
    G = simple_of(r.K)
    old = set(range(r.n_old)); new = set(range(r.n_old, r.n_old + r.n_new))
    PL = pairs_L(r)
    for a, b in G.edges():
        if a in old and b in old and frozenset((a, b)) not in PL:
            return False
    comps = list(nx.connected_components(G.subgraph(new))) if new else []
    for C in comps:
        A = {o for v in C for o in G[v] if o in old}
        if len(A) > 2:
            return False
        if len(A) == 2 and frozenset(A) not in PL:
            return False
        if len(A) == 2:
            B = G.subgraph(set(C) | A).copy()
            B.add_edge(*A)
            if not nx.check_planarity(B)[0]:
                return False
    return True


BASES = {
    "octahedron": nx.octahedral_graph(),
    "cube": nx.cubical_graph(),
    "icosahedron": nx.icosahedral_graph(),
    "K4": nx.complete_graph(4),
    "wheel6": nx.wheel_graph(6),
    "grid3x3": nx.convert_node_labels_to_integers(nx.grid_2d_graph(3, 3)),
}


def certificate(r):
    """First (base, binding) whose output simple graph is non-planar, or None."""
    PL = [tuple(sorted(p)) for p in pairs_L(r)]
    nid0 = 1000
    for name, G in BASES.items():
        V = sorted(G.nodes())
        if len(V) < r.n_old:
            continue
        for sigma in itertools.permutations(V, r.n_old):
            Hs = G.copy()
            Hs.add_edges_from((sigma[a], sigma[b]) for a, b in PL)
            if not nx.check_planarity(Hs)[0]:
                continue
            m = {i: sigma[i] for i in range(r.n_old)}
            for k in range(r.n_new):
                m[r.n_old + k] = nid0 + k
            Hp = Hs.copy()
            Hp.add_edges_from((m[u], m[v]) for u, v in r.K if u != v)
            if not nx.check_planarity(Hp)[0]:
                return name, sigma
    return None


def analyse(r):
    sm = oe.is_simple_minor_preserving(r)
    P = satisfies_P(r)
    cert = certificate(r)          # searched for every rule; for P rules it must find nothing
    sig = "2_2->%d_2" % len(r.K)
    return (dict(rule=oe.rule_str(r), sig=sig, n_new=r.n_new, delta_R=oe.delta_R(r), simple_minor=int(sm), P=int(P),
                 cert_base=(cert[0] if cert else ""), cert_binding=(",".join(map(str, cert[1])) if cert else "")),)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "planarity"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    S = oe.enumerate_rules(3, 1) + oe.enumerate_rules(4, 2)
    assert len(S) == 5225, len(S)
    with Pool() as pool:
        res = pool.map(analyse, S, chunksize=20)
    rows = [x[0] for x in res]
    P_viol = [x[0]["rule"] for x in res if x[0]["P"] and x[0]["cert_base"]]   # P rules with a non-planar output: must be none
    with open(os.path.join(a.out, "planarity_classification.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    n = len(rows)
    nP = sum(x["P"] for x in rows)
    nSM = sum(x["simple_minor"] for x in rows)
    P_not_sm = [x for x in rows if x["P"] and not x["simple_minor"]]
    nonP_nocert = [x for x in rows if not x["P"] and not x["cert_base"]]
    by_base = Counter(x["cert_base"] for x in rows if x["cert_base"])
    by_sig = {}
    for sig in sorted({x["sig"] for x in rows}):
        sub = [x for x in rows if x["sig"] == sig]
        by_sig[sig] = dict(rules=len(sub), simple_minor=sum(x["simple_minor"] for x in sub), P=sum(x["P"] for x in sub),
                           P_delta0=sum(1 for x in sub if x["P"] and x["delta_R"] == 0),
                           P_with_cycle_growth=sum(1 for x in sub if x["P"] and x["delta_R"] > 0))
    summ = dict(rules=n, simple_minor_preserving=nSM, satisfy_P=nP, P_not_simple_minor=len(P_not_sm),
                nonP_without_certificate=len(nonP_nocert), P_with_nonplanar_output_found=len(P_viol),
                certificates_by_base=dict(by_base), by_signature=by_sig,
                nonP_without_certificate_rules=[x["rule"] for x in nonP_nocert][:50],
                P_violations=P_viol[:50])
    json.dump(summ, open(os.path.join(a.out, "planarity_summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))
    gen = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper", "generated")
    if os.path.isdir(gen):
        tot = lambda k: sum(v[k] for v in by_sig.values())
        with open(os.path.join(gen, "tab_planarity.tex"), "w") as f:
            f.write("\\begin{tabular}{lrrrrrr}\n\\toprule\n")
            f.write("signature & rules & simple-minor-pres. & planarity-pres. ($P$) & $P$, $\\delta_R=0$ & $P$, $\\delta_R>0$ & s.m.-pres., not planarity-pres.\\\\\n\\midrule\n")
            for sig, v in by_sig.items():
                nm = sig.replace("2_2->", "$2_2\\to").replace("_2", "_2$", 1) if False else "$2_2\\to%s_2$" % sig.split("->")[1].split("_")[0]
                f.write("%s & %d & %d & %d & %d & %d & %d\\\\\n" % (nm, v["rules"], v["simple_minor"], v["P"], v["P_delta0"], v["P_with_cycle_growth"], v["simple_minor"] - v["P"]))
            f.write("\\midrule\nall of $\\mathcal S$ & %d & %d & %d & %d & %d & %d\\\\\n\\bottomrule\n\\end{tabular}\n" % (n, nSM, nP, tot("P_delta0"), tot("P_with_cycle_growth"), nSM - nP))
        with open(os.path.join(gen, "planarity_nums.tex"), "w") as f:
            for k, v in (("plnAll", n), ("plnSM", nSM), ("plnP", nP), ("plnNonP", n - nP), ("plnSMnotP", nSM - nP), ("plnNotSM", n - nSM),
                         ("plnPzero", tot("P_delta0")), ("plnPgrow", tot("P_with_cycle_growth"))):
                f.write("\\newcommand{\\%s}{%d}\n" % (k, v))
    ok = (not P_not_sm) and (not nonP_nocert) and (not P_viol)
    print("PLANARITY CLASSIFICATION CONSISTENT" if ok else "INCONSISTENT: see summary")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
