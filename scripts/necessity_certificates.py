#!/usr/bin/env python3
"""
Computer-assisted necessity check for the bounded signatures.

Claim checked (finite, exhaustive over the rule class):
  for every closed rule R in {2_2->3_2 (n_new<=1), 2_2->4_2 (n_new<=2)} that is NOT simple-minor-preserving there is a
  non-planar H0 and ONE injective match whose single event gives a planar underlying simple graph.
H0 candidates (all non-planar, all with underlying simple graph K5 or K33):
  A5, A33 : simple, edges oriented low->high
  B5, B33 : every undirected edge {a<b} present as (a,b),(a,b),(b,a)   (so parallel / antiparallel L match)
For each such rule a CERTIFICATE is written (H0 variant, match, resulting edge list, planar embedding as cyclic
neighbour orders).  Certificates are verified by an independent routine that does not call the planarity tester:
faces of the rotation system are traced and Euler's formula n - m + f = 1 + c is checked (this characterises a valid
planar embedding of the simple graph).
For minor-preserving rules (planarity preserved, Theorem 4 / simple version) a random sample of matches is checked as a consistency test.
Rules with no certificate are listed explicitly as open cases.
"""
import argparse, csv, json, os, random, sys
from multiprocessing import Pool
import networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe


def h0(variant):
    G = nx.complete_graph(5) if variant.endswith("5") else nx.complete_bipartite_graph(3, 3)
    es = []
    for a, b in sorted((min(u, v), max(u, v)) for u, v in G.edges()):
        es.append((a, b, 0))
        if variant.startswith("B"):
            es += [(a, b, 0), (b, a, 0)]
    return es


VARIANTS = ["A5", "A33", "B5", "B33"]
BASES = {v: h0(v) for v in VARIANTS}


def simple(edges):
    return oe.simple_graph([(a, b) for a, b, _ in edges])


def verify_embedding(G, rot):
    """rot: vertex -> cyclic list of neighbours. True iff it is a valid planar rotation system of G (Euler check)."""
    for v in G:
        if sorted(rot.get(v, [])) != sorted(G[v]):
            return False
    half = {(u, v) for u in G for v in G[u]}
    pos = {u: {w: i for i, w in enumerate(rot[u])} for u in rot}
    seen, faces = set(), 0
    for (u, v) in half:
        if (u, v) in seen:
            continue
        faces += 1
        a, b = u, v
        while (a, b) not in seen:
            seen.add((a, b))
            lst = rot[b]
            c = lst[(pos[b][a] + 1) % len(lst)]   # next edge in rotation around b after the incoming edge
            a, b = b, c
    n, m, c = G.number_of_nodes(), G.number_of_edges(), nx.number_connected_components(G)
    return n - m + faces == 1 + c


def rotation(G):
    ok, emb = nx.check_planarity(G)
    if not ok:
        return None
    return {v: list(emb.neighbors_cw_order(v)) for v in G}


def analyse(task):
    r, nsample, seed = task
    rule = oe.rule_str(r)
    mp = oe.is_simple_minor_preserving(r)
    mp_multi = oe.is_minor_preserving(r)
    rng = random.Random(seed)
    out = dict(rule=rule, minor_preserving=int(mp_multi), simple_minor_preserving=int(mp), n_new=r.n_old and r.n_new, delta_R=oe.delta_R(r), certificate="", matches_checked=0, bad=0)
    if not mp:
        for v in VARIANTS:
            base = BASES[v]
            for m in oe.find_matches(base, r):
                new, _ = oe.apply_matches(base, [m], r, oe.next_vertex_id(base))
                out["matches_checked"] += 1
                G = simple(new)
                rot = rotation(G)
                if rot is not None:
                    assert verify_embedding(G, rot)
                    cert = dict(rule=rule, H0=v, match={str(k): val for k, val in m[2].items()}, edges_after=[(a, b) for a, b, _ in new], rotation={str(k): val for k, val in rot.items()})
                    out["certificate"] = v
                    out["_cert"] = cert
                    return out
    else:
        for v in VARIANTS:
            base = BASES[v]
            ms = oe.find_matches(base, r)
            if len(ms) > nsample:
                ms = rng.sample(ms, nsample)
            for m in ms:
                new, _ = oe.apply_matches(base, [m], r, oe.next_vertex_id(base))
                out["matches_checked"] += 1
                if oe.is_planar(simple(new)):
                    out["bad"] += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/necessity")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--sample", type=int, default=60, help="matches sampled per H0 variant for minor-preserving rules")
    ap.add_argument("--seed", type=int, default=oe.DEFAULT_SEED)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rules = [("2_3", r) for r in oe.enumerate_rules(3, 1)] + [("2_4", r) for r in oe.enumerate_rules(4, 2)]
    tasks = [(r, a.sample, a.seed) for _, r in rules]
    with Pool(a.jobs) as p:
        res = p.map(analyse, tasks, chunksize=16)
    with open(os.path.join(a.out, "certificates.jsonl"), "w") as f:
        for o in res:
            if "_cert" in o:
                f.write(json.dumps(o["_cert"]) + "\n")
    rows = []
    for (sig, _), o in zip(rules, res):
        o = {k: v for k, v in o.items() if k != "_cert"}
        o["signature"] = sig
        rows.append(o)
    with open(os.path.join(a.out, "necessity_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    from collections import Counter
    c = Counter()
    for o in rows:
        if o["simple_minor_preserving"]:
            c[(o["signature"], "simple_minor_preserving", "planarised_in_sample" if o["bad"] else "never_planarised_in_sample")] += 1
        else:
            c[(o["signature"], "not_simple_minor_preserving", o["certificate"] or "NO_CERTIFICATE")] += 1
    for k in sorted(c):
        print(k, c[k])
    gaps = [o["rule"] for o in rows if (not o["simple_minor_preserving"]) and not o["certificate"]]
    open(os.path.join(a.out, "gaps.txt"), "w").write("\n".join(gaps))
    print("rules without certificate:", len(gaps), "| minor-preserving violations:", sum(o["bad"] for o in rows if o["simple_minor_preserving"]))


if __name__ == "__main__":
    main()
