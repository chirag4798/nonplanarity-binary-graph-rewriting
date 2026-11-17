#!/usr/bin/env python3
"""
check_noninjective_minor.py -- numerical check of Proposition (non-injective matches) for Theorem 4.

For every rule of S that is minor-preserving (multigraph version) or simple-minor-preserving, and for matches in which
distinct pattern variables may bind the SAME vertex (hosts: K5 and K33 with tight multiplicities - each edge once, in both orientations once, or twice in one
orientation - and a loop at every vertex, plus the standard double self-loop), apply one event, perform the contraction used in the
proof (new vertices of each branch set are contracted into the image of their old variable, discarded vertices are
deleted) and check that the result contains the host
  - as a multigraph (loops included) for minor-preserving rules,
  - as a simple graph for simple-minor-preserving rules.
Matches are sampled with a fixed seed (up to --cap per rule and host); the number of non-injective matches is reported.
No sampling bias is claimed: this checks the proof step on many cases, it is not a substitute for the proof.
"""
import argparse, os, random, sys
from collections import Counter
import networkx as nx

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe
import wolfram_bridge as wb


def hosts():
    """Hosts with TIGHT multiplicities, so that removing the matched edges really removes adjacencies:
    once oriented, both orientations once, same orientation twice; each with one loop per vertex."""
    out = {}
    for name, G in (("K5", nx.complete_graph(5)), ("K33", nx.complete_bipartite_graph(3, 3))):
        base = sorted((min(u, v), max(u, v)) for u, v in G.edges())
        loops = [(v, v) for v in sorted(G.nodes())]
        out[name + "_once"] = base + loops
        out[name + "_both"] = base + [(b, a) for a, b in base] + loops
        out[name + "_twice"] = base + base + loops
    out["selfloop"] = [(0, 0), (0, 0)]
    return out


def contract(edges_after, nid0, rule, bind, assign):
    mp = {}
    for k, a in enumerate(assign):
        mp[nid0 + k] = None if a == -1 else bind[a]
    res = []
    for u, v in edges_after:
        u2 = mp.get(u, u); v2 = mp.get(v, v)
        if u2 is None or v2 is None:
            continue
        res.append((u2, v2))
    return res


def conn_ok(r, assign):
    br = {a: a for a in range(r.n_old)}; br.update({r.n_old + k: x for k, x in enumerate(assign)})
    for a in range(r.n_old):
        S_ = {v for v, b in br.items() if b == a}
        g = nx.Graph(); g.add_nodes_from(S_); g.add_edges_from((u, v) for u, v in r.K if u in S_ and v in S_ and u != v)
        if not nx.is_connected(g):
            return False
    return True


def negative_control(S, H, limit=600):
    """Rules that are NOT simple-minor-preserving must fail the contraction test for every valid (connected-branch-set)
    assignment on some event: shows that the test above is not vacuous."""
    import itertools
    tot = passing = 0
    for r in S:
        if oe.is_simple_minor_preserving(r):
            continue
        tot += 1; found = False
        for assign in itertools.product(list(range(r.n_old)) + [-1], repeat=r.n_new):
            if not conn_ok(r, assign):
                continue
            ok = True
            for hn, E in H.items():
                nid0 = 1 + max(max(e) for e in E)
                for (i, j, b) in wb.find_matches_general(E, r):
                    C = contract(wb.apply_events(E, [(i, j, b)], r), nid0, r, b, assign)
                    Gs = nx.Graph(); Gs.add_edges_from(e for e in E if e[0] != e[1])
                    Gc = nx.Graph(); Gc.add_edges_from(e for e in C if e[0] != e[1])
                    if not all(Gc.has_edge(*e) for e in Gs.edges()):
                        ok = False; break
                if not ok:
                    break
            if ok:
                found = True; break
        passing += found
        if tot >= limit:
            break
    return tot, passing


def msig(edges):
    return Counter(tuple(sorted(e)) for e in edges)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20240601)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    S = oe.enumerate_rules(3, 1) + oe.enumerate_rules(4, 2)
    H = hosts()
    stats = Counter(); viol = []
    for r in S:
        mm = oe.minor_models(r)
        sm = mm or oe.minor_models(r, simple=True)
        if not sm:
            continue
        multi = bool(mm)
        model = sm[0]["assign"]
        for hname, edges in H.items():
            ms = wb.find_matches_general(edges, r)
            if len(ms) > a.cap:
                ms = rng.sample(ms, a.cap)
            for (i, j, bind) in ms:
                noninj = not wb.injective(bind)
                after = wb.apply_events(edges, [(i, j, bind)], r)
                nid0 = 1 + max(max(e) for e in edges)
                C = contract(after, nid0, r, bind, model)
                stats["events"] += 1; stats["noninjective_events"] += noninj
                if multi:
                    need = msig(edges); have = msig(C)
                    ok = all(have[k] >= v for k, v in need.items())
                    stats["multigraph_checked"] += 1
                else:
                    Gs = nx.Graph(); Gs.add_edges_from(e for e in edges if e[0] != e[1])
                    Gc = nx.Graph(); Gc.add_edges_from(e for e in C if e[0] != e[1])
                    ok = all(Gc.has_edge(*e) for e in Gs.edges())
                    stats["simple_checked"] += 1
                if not ok:
                    viol.append((oe.rule_str(r), hname, bind)); stats["violations"] += 1
    print(dict(stats))
    for v in viol[:10]:
        print("VIOLATION", v)
    tot, passing = negative_control(S, H)
    print("negative control: %d rules that are not simple-minor-preserving, %d pass the test for some valid assignment (must be 0)" % (tot, passing))
    gen = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper", "generated")
    if os.path.isdir(gen):
        with open(os.path.join(gen, "noninj_minor_nums.tex"), "w") as f:
            for k, v in (("nimEvents", stats["events"]), ("nimNonInj", stats["noninjective_events"]), ("nimCtrl", tot)):
                f.write("\\newcommand{\\%s}{%d}\n" % (k, v))
    good = (not viol) and passing == 0
    print("NON-INJECTIVE MINOR CHECK PASSED" if good else "NON-INJECTIVE MINOR CHECK FAILED")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
