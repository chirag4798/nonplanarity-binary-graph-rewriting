#!/usr/bin/env python3
"""Search for counterexamples to the necessity direction (iii)=>(i) of the dichotomy in signatures beyond 2_2->{3,4}_2.
For each closed rule m_2 -> k_2 (given n_new bound): decide simple-minor preservation (oe.minor_models); for rules that are NOT
simple-minor-preserving look for a one-event planarisation of a non-planar host (simple/doubled K5, K33); report rules with no
certificate (candidate counterexamples) together with whether K itself is non-planar (then the rule is a PROVEN counterexample:
every output contains K).  For simple-minor-preserving rules a sample of matches is checked for consistency with Theorem 4."""
import argparse, csv, os, random, sys, time
from collections import Counter
from multiprocessing import Pool
import networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe
import general_rules as gr


def host(variant):
    G = nx.complete_graph(5) if variant.endswith("5") else nx.complete_bipartite_graph(3, 3)
    es = []
    for a, b in sorted((min(u, v), max(u, v)) for u, v in G.edges()):
        es.append((a, b, 0))
        if variant.startswith("B"):
            es += [(a, b, 0), (b, a, 0)]
    return es


HOSTS = {v: host(v) for v in ("A5", "A33", "B5", "B33")}


def simple_planar(edges):
    return oe.is_planar(oe.simple_graph([(a, b) for a, b, *_ in edges]))


def planted_hosts(r, rng, max_hosts):
    """K5 / K33 with L planted: for an injective placement of the old variables on host vertices, every host edge between
    two placed vertices carries exactly the edges (multiplicity, orientation) of L for that pair (if any); all other host
    edges are single oriented edges.  Guarantees a match with exactly L's multiplicities."""
    pairs_L = {}
    for (a, b) in r.L:
        pairs_L.setdefault(frozenset((a, b)), []).append((a, b))
    out = []
    for gname, G in (("K5", nx.complete_graph(5)), ("K33", nx.complete_bipartite_graph(3, 3))):
        places = list(__import__("itertools").permutations(range(G.number_of_nodes()), r.n_old))
        rng.shuffle(places)
        for pl in places:
            inv = {v: i for i, v in enumerate(pl)}      # host vertex -> old variable
            es = []
            for u, v in sorted((min(x, y), max(x, y)) for x, y in G.edges()):
                key = frozenset((inv.get(u, -1 - u), inv.get(v, -1 - v)))
                if u in inv and v in inv and key in pairs_L:
                    for (a, b) in pairs_L[key]:
                        es.append((pl[a], pl[b], 0))
                else:
                    es.append((u, v, 0))
            out.append((gname, es))
            if len(out) >= max_hosts:
                return out
    return out


def analyse(task):
    r, seed, sample = task
    rng = random.Random(seed)
    sm = oe.is_simple_minor_preserving(r)
    Kplanar = oe.is_planar(oe.simple_graph(list(r.K)))
    out = dict(rule=oe.rule_str(r), m=len(r.L), k=len(r.K), n_new=r.n_new, sm=int(sm), K_planar=int(Kplanar), cert="", bad=0)
    hosts = planted_hosts(r, rng, 4 if sm else 400) + [(v, h) for v, h in HOSTS.items()]
    for v, base in hosts:
        ms = gr.find_matches(base, r)
        if sm and len(ms) > sample:
            ms = rng.sample(ms, sample)
        for mt in ms:
            new = gr.apply_event(base, mt, r, gr.next_vid(base))
            if simple_planar(new):
                if sm:
                    out["bad"] += 1
                else:
                    out["cert"] = v
                    return out
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigs", nargs="+", required=True, help="m:k:max_new, e.g. 1:3:2 2:5:3 3:4:2")
    ap.add_argument("--out", default="results/general")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--sample", type=int, default=20)
    ap.add_argument("--seed", type=int, default=oe.DEFAULT_SEED)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for sg in a.sigs:
        m, k, mx = map(int, sg.split(":"))
        t0 = time.time()
        rules = gr.enumerate_rules(m, k, mx)
        t1 = time.time()
        with Pool(a.jobs) as p:
            res = p.map(analyse, [(r, a.seed + i, a.sample) for i, r in enumerate(rules)], chunksize=32)
        c = Counter()
        gaps = [o for o in res if (not o["sm"]) and not o["cert"]]
        for o in res:
            o["sig"] = sg
        rows += res
        print("signature %d_2->%d_2 (n_new<=%d): closed rules %d | simple-minor-preserving %d | not sm: certificate %d, NO certificate %d "
              "(of which K non-planar: %d) | planarisations among sm rules (sampled): %d | enum %.0fs total %.0fs"
              % (m, k, mx, len(res), sum(o["sm"] for o in res), sum(1 for o in res if not o["sm"] and o["cert"]), len(gaps),
                 sum(1 for o in gaps if not o["K_planar"]), sum(o["bad"] for o in res), t1 - t0, time.time() - t0), flush=True)
        for o in gaps[:5]:
            print("   no-certificate example:", o["rule"], "K planar" if o["K_planar"] else "K NON-planar")
    with open(os.path.join(a.out, "general_signature_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
