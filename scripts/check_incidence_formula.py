#!/usr/bin/env python3
"""Numerical check of the incidence-graph version of the cycle-rank accounting for hyperedges of arbitrary arity:
   Delta beta1(incidence graph) = sum_{e in K}(ar(e)-1) - sum_{e in L}(ar(e)-1) - n_new   for closed rules (L connected, var(L) in var(K), K connected).
Random closed rules with hyperedges of arity 2 and 3, random connected hosts, random matches (distinct pattern variables may bind equal vertices)."""
import random, itertools, networkx as nx
rng = random.Random(20240601)

def beta1_inc(hedges):
    G = nx.MultiGraph()
    for i, e in enumerate(hedges):
        G.add_node(("e", i))
        for v in e:
            G.add_edge(("e", i), ("v", v))
    return G.number_of_edges() - G.number_of_nodes() + nx.number_connected_components(G)

def connected(hedges):
    G = nx.Graph()
    for e in hedges:
        for a in e: G.add_node(a)
        for a, b in zip(e, e[1:]): G.add_edge(a, b)
    return G.number_of_nodes() > 0 and nx.is_connected(G)

def random_rule():
    while True:
        nl = rng.randint(1, 2); n_old = rng.randint(2, 4); n_new = rng.randint(0, 2)
        L = [tuple(rng.randrange(n_old) for _ in range(rng.choice((2, 3)))) for _ in range(nl)]
        if {v for e in L for v in e} != set(range(n_old)) or not connected(L): continue
        nk = rng.randint(1, 4); nv = n_old + n_new
        K = [tuple(rng.randrange(nv) for _ in range(rng.choice((2, 3)))) for _ in range(nk)]
        if {v for e in K for v in e} != set(range(nv)) or not connected(K): continue
        return L, K, n_old, n_new

def matches(host, L):
    out = []
    def rec(i, used, bind):
        if i == len(L): out.append((tuple(used), dict(bind))); return
        for j, e in enumerate(host):
            if j in used or len(e) != len(L[i]): continue
            b = dict(bind); ok = True
            for var, val in zip(L[i], e):
                if var in b and b[var] != val: ok = False; break
                b[var] = val
            if ok: rec(i + 1, used + [j], b)
    rec(0, [], {}); return out

checked = bad = 0
for _ in range(3000):
    L, K, n_old, n_new = random_rule()
    host = [tuple(rng.randrange(6) for _ in range(rng.choice((2, 3)))) for _ in range(rng.randint(4, 8))]
    if not connected(host): continue
    ms = matches(host, L)
    if not ms: continue
    used, bind = rng.choice(ms)
    nid = 1 + max(v for e in host for v in e); m = dict(bind)
    for k in range(n_old, n_old + n_new): m[k] = nid; nid += 1
    new = [e for j, e in enumerate(host) if j not in used] + [tuple(m[v] for v in e) for e in K]
    d = sum(len(e) - 1 for e in K) - sum(len(e) - 1 for e in L) - n_new
    checked += 1
    if beta1_inc(new) - beta1_inc(host) != d or not connected(new): bad += 1
print("incidence-graph cycle-rank accounting: %d random events (hyperedge arity 2 or 3, closed rules), violations %d" % (checked, bad))
