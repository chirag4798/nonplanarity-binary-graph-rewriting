#!/usr/bin/env python3
"""(1) Evolve representative Wolfram-model binary rules from the standard double self-loop initial condition with Wolfram-style
(non-injective) matching; check the cycle-rank law and record planarity of the simple projection.
(2) Check Theorem 1 (cycle-rank law) and persistence of non-planarity for NON-INJECTIVE matches on random hosts with loops and
parallel edges, for the Wolfram rules and for samples of minor-preserving rules.  Outputs: results/wolfram_rules/*.csv, *.txt."""
import csv, os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import oligon_experiment as oe
import wolfram_bridge as wb

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "wolfram_rules")
rules = [l.strip().split("|") for l in open(os.path.join(root, "rules.txt")) if l.strip()]
rows = []
for rid, s in rules:
    r = oe.parse_rule(s); d = oe.delta_R(r)
    for order in ("canonical", "random"):
        rng = random.Random(oe.DEFAULT_SEED)
        edges = [(0, 0), (0, 0)]
        o0 = wb.observe(edges); events = 0; noninj_total = 0
        rows.append(dict(rule=rid, order=order, t=0, events=0, noninj_events=0, cum_events=0, law_ok=1, **o0))
        for t in range(1, 21):
            edges, chosen = wb.generation(edges, r, rng, order)
            events += len(chosen); ni = sum(1 for m in chosen if not wb.injective(m[2])); noninj_total += ni
            o = wb.observe(edges)
            rows.append(dict(rule=rid, order=order, t=t, events=len(chosen), noninj_events=ni, cum_events=events,
                             law_ok=int(o["beta1"] == o0["beta1"] + d * events), **o))
            if o["E"] > 6000 or not chosen:
                break
with open(os.path.join(root, "selfloop_evolution.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# (2) non-injective check
rng = random.Random(oe.DEFAULT_SEED)
minor_rules = [r for r in oe.enumerate_rules(3, 1) + oe.enumerate_rules(4, 2) if oe.is_minor_preserving(r)]
sample = rng.sample(minor_rules, 150)
named = [oe.parse_rule(s) for _, s in rules]
tot = law_bad = planar_bad = noninj = 0
for r in named + sample:
    d = oe.delta_R(r)
    for _ in range(12):
        base = [(a, b) for a in range(5) for b in range(a + 1, 5)]               # K5
        extra = [(rng.randrange(5), rng.randrange(5)) for _ in range(rng.randint(0, 6))]   # loops / parallel edges allowed
        edges = [((a, b) if rng.random() < .5 else (b, a)) for a, b in base] + extra
        b0 = oe.beta1_of(edges)[0]
        for step in range(6):
            ms = wb.find_matches_general(edges, r)
            if not ms: break
            m = rng.choice(ms)
            edges = wb.apply_events(edges, [m], r)
            tot += 1; noninj += not wb.injective(m[2])
            if oe.beta1_of(edges)[0] != b0 + d * (step + 1): law_bad += 1
            if oe.is_planar(oe.simple_graph(edges)): planar_bad += 1
msg = ("non-injective-match check: %d events on random non-planar hosts with loops and parallel edges (%d of them with a non-injective match), "
       "%d Wolfram rules + 150 sampled minor-preserving rules: cycle-rank law violations %d, events producing a planar graph %d\n" % (tot, noninj, len(named), law_bad, planar_bad))
open(os.path.join(root, "noninjective_check.txt"), "w").write(msg); print(msg)
