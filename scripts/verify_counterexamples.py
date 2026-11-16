#!/usr/bin/env python3
"""Checks of the two counterexamples to the necessity direction in general signatures:
 (1) Delta->Y rule  {{1,2},{1,3},{2,3}} -> {{1,4},{2,4},{3,4}}  (3_2 -> 3_2): not simple-minor-preserving, yet an event never turns a
     non-planar host planar (proof: Y->Delta preserves planarity and Y->Delta(H') = H).  Checked here on random hosts, all matches.
 (2) Kuratowski-gadget rule (2_2 -> 12_2): K contains a K_{3,3}, so every output is non-planar, whatever the host; L's pair {1,3} cannot
     be covered by any branch model because variable 1 only touches variable 2 in K."""
import os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import networkx as nx
import oligon_experiment as oe
import general_rules as gr

rng = random.Random(oe.DEFAULT_SEED)

# (1) Delta -> Y
dy = oe.Rule(((0, 1), (0, 2), (1, 2)), ((0, 3), (1, 3), (2, 3)), 3, 1)
print("Delta->Y: closed", oe.is_closed(dy), "| simple-minor-preserving", oe.is_simple_minor_preserving(dy), "| minor-preserving", oe.is_minor_preserving(dy))
viol = events = nonplanar_hosts = made_nonplanar = 0
for _ in range(4000):
    n = rng.randint(5, 9)
    G = nx.gnm_random_graph(n, rng.randint(n + 2, min(n * (n - 1) // 2, 3 * n)), seed=rng.randint(0, 10**9))
    host = [((a, b) if rng.random() < .5 else (b, a)) + (0,) for a, b in G.edges()]
    hp = oe.is_planar(oe.simple_graph([(a, b) for a, b, _ in host]))
    for mt in gr.find_matches(host, dy):
        new = gr.apply_event(host, mt, dy, gr.next_vid(host))
        npn = oe.is_planar(oe.simple_graph([(a, b) for a, b, _ in new]))
        events += 1
        if not hp:
            nonplanar_hosts += 1
            if npn: viol += 1
        elif not npn:
            made_nonplanar += 1
print("  events %d | on non-planar hosts %d: events giving a planar graph = %d | planar hosts made non-planar (Delta->Y does NOT preserve planarity): %d" % (events, nonplanar_hosts, viol, made_nonplanar))

# (2) Kuratowski gadget
K = [(0, 1), (1, 3), (2, 3)] + [(a, b) for a in (3, 4, 5) for b in (6, 7, 8)]
gad = oe.Rule(((0, 1), (0, 2)), tuple(K), 3, 6)
print("Gadget rule 2_2->%d_2: closed %s | K planar %s | simple-minor-preserving %s | minor-preserving %s | delta_R %d" % (
    len(K), oe.is_closed(gad), oe.is_planar(oe.simple_graph(K)), oe.is_simple_minor_preserving(gad), oe.is_minor_preserving(gad), oe.delta_R(gad)))
planar_out = events = 0
for _ in range(500):
    n = rng.randint(4, 8)
    G = nx.gnm_random_graph(n, rng.randint(n, 2 * n), seed=rng.randint(0, 10**9))
    if not nx.is_connected(G): continue
    host = [((a, b) if rng.random() < .5 else (b, a)) + (0,) for a, b in G.edges()]
    for mt in gr.find_matches(host, gad)[:5]:
        new = gr.apply_event(host, mt, gad, gr.next_vid(host)); events += 1
        planar_out += oe.is_planar(oe.simple_graph([(a, b) for a, b, _ in new]))
print("  events %d | outputs that are planar: %d (hosts include planar ones)" % (events, planar_out))
