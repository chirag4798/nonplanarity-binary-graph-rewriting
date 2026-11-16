"""Matching that follows the Wolfram-model convention more closely than oligon_experiment: distinct pattern variables MAY bind the same vertex,
states may contain loops and parallel edges.  (Update order and the exact event ordering of WolframModel are NOT reproduced.)"""
import random
from collections import defaultdict
import oligon_experiment as oe


def find_matches_general(edges, rule):
    (p, q), (r, s) = rule.L
    inc = defaultdict(list)
    for i, (a, b) in enumerate(edges):
        inc[a].append(i); inc[b].append(i)
    out = []
    for i, (a, b) in enumerate(edges):
        bind = {p: a, q: b}
        cand = set()
        for v in (r, s):
            if v in bind:
                cand.update(inc[bind[v]])
        for j in sorted(cand):
            if j == i:
                continue
            c, d = edges[j]
            b2 = dict(bind); ok = True
            for var, val in ((r, c), (s, d)):
                if var in b2:
                    if b2[var] != val: ok = False; break
                else:
                    b2[var] = val
            if ok:
                out.append((i, j, b2))
    return out


def injective(bind):
    return len(set(bind.values())) == len(bind)


def apply_events(edges, chosen, rule):
    nid = 1 + max((max(a, b) for a, b in edges), default=-1)
    consumed, produced = set(), []
    for (i, j, bind) in chosen:
        m = dict(bind)
        for k in range(rule.n_old, rule.n_old + rule.n_new):
            m[k] = nid; nid += 1
        produced += [(m[u], m[v]) for (u, v) in rule.K]
        consumed.update((i, j))
    return [e for k, e in enumerate(edges) if k not in consumed] + produced


def generation(edges, rule, rng, order="canonical"):
    ms = find_matches_general(edges, rule)
    if order == "random":
        rng.shuffle(ms)
    used, chosen = set(), []
    for m in ms:
        if m[0] in used or m[1] in used: continue
        used.update((m[0], m[1])); chosen.append(m)
    return apply_events(edges, chosen, rule), chosen


def observe(edges):
    b1, V, E, c = oe.beta1_of(edges)
    loops = sum(1 for a, b in edges if a == b)
    return dict(V=V, E=E, loops=loops, beta1=b1, planar_simple=int(oe.is_planar(oe.simple_graph(edges))))
