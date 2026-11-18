#!/usr/bin/env python3
"""Multi-source literature search for Appendix C (queries are fixed below; results change over time).

Sources (public APIs, no login): OpenAlex, arXiv, Crossref, Semantic Scholar.
  python3 literature_search.py search  > results/literature/multisource_queries.txt
  python3 literature_search.py chain   > results/literature/citation_chain.txt
'search' prints, per source and query, the number of hits and the first hits (year, authors, title, venue/DOI).
'chain' resolves the seed works in OpenAlex, then lists the works that cite each seed (forward chain) and, for the seed
list, the works they cite (backward chain), filtered by a keyword screen on the title (printed unfiltered counts too)."""
import json, subprocess, sys, time, urllib.parse

UA = "research-script (literature check)"


def get(url, tries=6):
    for k in range(tries):
        out = subprocess.run(["curl", "-s", "-m", "60", "-A", UA, url], capture_output=True, text=True).stdout
        try:
            d = json.loads(out)
        except Exception:
            time.sleep(5 * (k + 1)); continue
        if isinstance(d, dict) and (d.get("error") == "Rate limit exceeded" or d.get("code") == "429"):
            time.sleep(min(30, d.get("retryAfter", 10) + 3)); continue
        return d
    return None


def raw(url, tries=3):
    for k in range(tries):
        out = subprocess.run(["curl", "-s", "-m", "60", "-A", UA, url], capture_output=True, text=True).stdout
        if out.strip():
            return out
        time.sleep(5)
    return ""


QUERIES = [
    "graph rewriting planarity preserved",
    "graph transformation planar graphs preserved by rules",
    "graph transformation minor closed class",
    "graph minor rewriting rule preserves",
    "double pushout rewriting cycle rank",
    "graph rewriting first Betti number",
    "hypergraph rewriting topology invariant",
    "Wolfram model hypergraph rewriting non-planar",
    "Wolfram model particles topological defects",
    "hypergraph rewriting homology",
    "preservation of graph properties under graph transformation",
    "graph grammar planarity",
    "graph grammars Kuratowski",
    "local graph rewriting minor theory",
    "edge replacement grammar planar graphs",
    "graph rewriting K5 K3,3 subdivision",
    "graph rewriting dynamic graphs planarity maintenance",
    "monadic second order logic graph transformation minors",
    "homeomorphic embedding graph rewriting well-quasi-order",
    "rewriting of graphs and cyclomatic number",
]


def openalex_search(q):
    d = get("https://api.openalex.org/works?search=%s&per-page=8&select=id,title,publication_year,doi,cited_by_count,authorships,primary_location" % urllib.parse.quote(q))
    if not d:
        print("  [openalex] no response"); return
    print("  [openalex] total %s" % d.get("meta", {}).get("count"))
    for r in d.get("results", []):
        au = "; ".join(a["author"]["display_name"] for a in r.get("authorships", [])[:3])
        ven = ((r.get("primary_location") or {}).get("source") or {}).get("display_name", "")
        print("   - %s | %s | %s | %s | %s | cited %s" % (r.get("publication_year"), au, (r.get("title") or "")[:100], (ven or "")[:50], r.get("doi") or "", r.get("cited_by_count")))


def arxiv_search(q):
    import re
    url = "https://export.arxiv.org/api/query?search_query=%s&max_results=8" % "+AND+".join("all:" + urllib.parse.quote(t) for t in q.replace(",", "").split())
    x = ""
    for k in range(5):                     # arXiv asks for one request per 3 s; retry on rate-limit answers
        x = raw(url)
        if "totalResults" in x:
            break
        time.sleep(10 * (k + 1))
    tot = re.search(r"<opensearch:totalResults[^>]*>(\d+)<", x)
    print("  [arxiv] total %s" % (tot.group(1) if tot else "no response"))
    for e in re.findall(r"<entry>(.*?)</entry>", x, re.S)[:8]:
        t = re.sub(r"\s+", " ", re.search(r"<title>(.*?)</title>", e, re.S).group(1))
        i = re.search(r"<id>(.*?)</id>", e).group(1)
        y = re.search(r"<published>(\d{4})", e).group(1)
        au = "; ".join(re.findall(r"<name>(.*?)</name>", e)[:3])
        print("   - %s | %s | %s | %s" % (y, au, t[:100], i))
    time.sleep(4)


def crossref_search(q):
    d = get("https://api.crossref.org/works?query=%s&rows=8&select=title,DOI,issued,author,container-title" % urllib.parse.quote(q))
    if not d:
        print("  [crossref] no response"); return
    print("  [crossref] total %s" % d["message"].get("total-results"))
    for r in d["message"]["items"]:
        au = "; ".join(a.get("family", "") for a in r.get("author", [])[:3])
        print("   - %s | %s | %s | %s | %s" % ((r.get("issued", {}).get("date-parts") or [[None]])[0][0], au, (r.get("title") or [""])[0][:100], (r.get("container-title") or [""])[0][:50], r["DOI"]))


def s2_search(q):
    d = get("https://api.semanticscholar.org/graph/v1/paper/search?query=%s&limit=8&fields=title,year,externalIds,authors,venue" % urllib.parse.quote(q))
    if not d or "data" not in d:
        print("  [semanticscholar] no usable response"); return
    print("  [semanticscholar] total %s" % d.get("total"))
    for r in d["data"]:
        au = "; ".join(a["name"] for a in r.get("authors", [])[:3])
        print("   - %s | %s | %s | %s | %s" % (r.get("year"), au, (r.get("title") or "")[:100], (r.get("venue") or "")[:40], (r.get("externalIds") or {}).get("DOI", "")))


SEEDS = [  # title fragments resolved in OpenAlex
    "Double-pushout graph transformation revisited",
    "Graph minors. XX. Wagner's conjecture",
    "Graph minors. I. Excluding a forest",
    "A Class of Models with the Potential to Represent Fundamental Physics",
    "Confluence of graph transformation revisited",
    "Fundamentals of Algebraic Graph Transformation",
    "Algebraic approaches to graph transformation",
    "Some Physical Properties of Wolfram Models",
    "Some Relativistic and Gravitational Properties of the Wolfram Model",
]
KEY = ("planar", "minor", "betti", "cycle", "homolog", "topolog", "invariant", "preserv", "wolfram", "hypergraph", "kuratowski", "genus", "embedd", "rewrit", "transformation", "grammar")


def chain():
    for s in SEEDS:
        d = get("https://api.openalex.org/works?search=%s&per-page=3&select=id,title,publication_year,cited_by_count,referenced_works" % urllib.parse.quote(s))
        if not d or not d.get("results"):
            print("== seed %r: not resolved" % s); continue
        w = d["results"][0]
        wid = w["id"].split("/")[-1]
        print("== seed %r -> %s (%s) cited by %s; references %d" % (s, (w["title"] or "")[:90], w["publication_year"], w["cited_by_count"], len(w.get("referenced_works", []))))
        hits = total = 0
        page = 1
        while page <= 5:
            c = get("https://api.openalex.org/works?filter=cites:%s&per-page=100&page=%d&select=title,publication_year,doi,authorships" % (wid, page))
            if not c or not c.get("results"):
                break
            for r in c["results"]:
                total += 1
                t = (r.get("title") or "")
                if any(k in t.lower() for k in KEY):
                    hits += 1
                    au = "; ".join(a["author"]["display_name"] for a in r.get("authorships", [])[:2])
                    print("   * %s | %s | %s | %s" % (r.get("publication_year"), au, t[:110], r.get("doi") or ""))
            if len(c["results"]) < 100:
                break
            page += 1
            time.sleep(1)
        print("   forward chain: %d citing works screened, %d title-keyword hits" % (total, hits))
        time.sleep(1)


def targeted():
    """Forward chain restricted by OpenAlex title/abstract search, for seeds with many citing works (the plain chain screens only the 500 most relevant)."""
    seeds = {"Double-pushout graph transformation revisited": None, "Fundamentals of Algebraic Graph Transformation": None,
             "ALGEBRAIC APPROACHES TO GRAPH TRANSFORMATION PART I BASIC CONCEPTS AND DOUBLE PUSHOUT APPROACH": None,
             "Graph Minors. XX. Wagner's conjecture": None}
    topics = ["planar", "planarity", "Betti", "cycle rank", "cyclomatic", "homology", "Kuratowski", "graph minor", "topological invariant", "genus"]
    for s in seeds:
        d = get("https://api.openalex.org/works?search=%s&per-page=1&select=id,title" % urllib.parse.quote(s))
        if not d or not d.get("results"):
            print("== seed %r not resolved" % s); continue
        wid = d["results"][0]["id"].split("/")[-1]
        print("== seed %s (%s)" % (d["results"][0]["title"][:80], wid))
        seen = set()
        for t in topics:
            c = get("https://api.openalex.org/works?filter=cites:%s,title_and_abstract.search:%s&per-page=25&select=title,publication_year,doi,authorships" % (wid, urllib.parse.quote(t)))
            n = c.get("meta", {}).get("count") if c else None
            print("  topic %r: %s citing works" % (t, n))
            for r in (c or {}).get("results", []):
                if r["doi"] in seen: continue
                seen.add(r["doi"])
                au = "; ".join(a["author"]["display_name"] for a in r.get("authorships", [])[:2])
                print("   * %s | %s | %s | %s" % (r.get("publication_year"), au, (r.get("title") or "")[:110], r.get("doi") or ""))
            time.sleep(1)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "chain":
        chain(); return
    if len(sys.argv) > 1 and sys.argv[1] == "targeted":
        targeted(); return
    if len(sys.argv) > 1 and sys.argv[1] == "arxiv":
        for q in QUERIES:
            print("== %s" % q); arxiv_search(q); time.sleep(4)
        return
    for q in QUERIES:
        print("== %s" % q)
        for f in (openalex_search, arxiv_search, crossref_search, s2_search):
            try:
                f(q)
            except Exception as e:
                print("  [%s] error %s" % (f.__name__, e))
            time.sleep(2)


if __name__ == "__main__":
    main()
