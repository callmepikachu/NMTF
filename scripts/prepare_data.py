"""Turn the raw files into the gene network and GO annotations used by the analysis.

Reads   data/raw/string/9606.protein.physical.links.v12.0.txt.gz   STRING physical links
        data/raw/string/9606.protein.info.v12.0.txt.gz             STRING protein -> gene symbol
        data/raw/go/goa_human.gaf.gz                               GO annotations of human genes
        config/config.yaml                                         thresholds and network size
Writes  data/processed/edges.tsv        gene_a, gene_b, score      (undirected, a < b)
        data/processed/genes.tsv        gene, degree               (the genes of the network)
        data/processed/annotations.tsv  gene, go_id                (GO terms of those genes)
        data/processed/summary.json     counts for the report

Selection: keep links with score >= string_min_score; keep the annotated genes; take the
n_genes with the most links (ties broken by name); keep the largest connected component of
the network among them.

    python scripts/prepare_data.py [--config config/config.yaml] [--raw data/raw] [--out data/processed]
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent.parent

GAF_COLUMNS = ["db", "db_id", "symbol", "qualifier", "go_id", "reference", "evidence", "with",
               "aspect", "name", "synonym", "type", "taxon", "date", "assigned_by", "extension",
               "gene_product_form"]


def read_string_links(path, min_score):
    """STRING links file (space separated, every link listed in both directions)."""
    links = pd.read_csv(path, sep=" ", dtype={"protein1": str, "protein2": str})
    n_rows = len(links)
    links = links[links["combined_score"] >= min_score]
    return links, n_rows


def read_string_symbols(path):
    """STRING protein id -> preferred gene name."""
    info = pd.read_csv(path, sep="\t", dtype=str, usecols=[0, 1])
    info.columns = ["protein", "symbol"]
    return dict(zip(info["protein"], info["symbol"]))


def read_go_annotations(path, aspect, exclude_evidence):
    """GAF file -> unique (gene symbol, GO id) pairs of one aspect, without NOT annotations."""
    gaf = pd.read_csv(path, sep="\t", comment="!", header=None, names=GAF_COLUMNS, dtype=str,
                      usecols=["symbol", "qualifier", "go_id", "evidence", "aspect"])
    n_rows = len(gaf)
    gaf = gaf[(gaf["aspect"] == aspect)
              & ~gaf["qualifier"].fillna("").str.contains("NOT")
              & ~gaf["evidence"].isin(exclude_evidence)]
    return gaf[["symbol", "go_id"]].drop_duplicates(), n_rows


def build_network(links, symbols, annotated_genes, n_genes):
    """Undirected gene network, restricted to the n best-connected annotated genes."""
    edges = pd.DataFrame({"a": links["protein1"].map(symbols), "b": links["protein2"].map(symbols),
                          "score": links["combined_score"]}).dropna()
    edges = edges[edges["a"] != edges["b"]]
    edges = edges[edges["a"].isin(annotated_genes) & edges["b"].isin(annotated_genes)]
    lo = edges[["a", "b"]].min(axis=1)
    hi = edges[["a", "b"]].max(axis=1)
    edges = pd.DataFrame({"a": lo, "b": hi, "score": edges["score"]})
    edges = edges.groupby(["a", "b"], as_index=False)["score"].max()
    n_edges_all, n_genes_all = len(edges), pd.unique(edges[["a", "b"]].values.ravel()).size

    degree = pd.concat([edges["a"], edges["b"]]).value_counts()
    top = degree.sort_index().sort_values(ascending=False, kind="stable").head(n_genes).index
    edges = edges[edges["a"].isin(top) & edges["b"].isin(top)]

    genes = sorted(pd.unique(edges[["a", "b"]].values.ravel()))
    index = {g: i for i, g in enumerate(genes)}
    r = edges["a"].map(index).to_numpy()
    c = edges["b"].map(index).to_numpy()
    A = sp.coo_matrix((np.ones(len(edges)), (r, c)), shape=(len(genes), len(genes)))
    _, comp = connected_components(A, directed=False)
    biggest = np.flatnonzero(comp == np.bincount(comp).argmax())
    keep = {genes[i] for i in biggest}
    edges = edges[edges["a"].isin(keep) & edges["b"].isin(keep)].sort_values(["a", "b"])
    return edges.reset_index(drop=True), n_edges_all, n_genes_all


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default=ROOT / "config" / "config.yaml")
    parser.add_argument("--raw", default=ROOT / "data" / "raw")
    parser.add_argument("--out", default=ROOT / "data" / "processed")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))["data"]
    raw, out = Path(args.raw), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    links, n_link_rows = read_string_links(
        raw / "string" / "9606.protein.physical.links.v12.0.txt.gz", cfg["string_min_score"])
    symbols = read_string_symbols(raw / "string" / "9606.protein.info.v12.0.txt.gz")
    annot, n_gaf_rows = read_go_annotations(raw / "go" / "goa_human.gaf.gz", cfg["go_aspect"],
                                            cfg["go_exclude_evidence"])
    n_go_pairs = len(annot)
    edges, n_edges_all, n_genes_all = build_network(links, symbols, set(annot["symbol"]),
                                                    cfg["n_genes"])
    genes = sorted(pd.unique(edges[["a", "b"]].values.ravel()))
    degree = pd.concat([edges["a"], edges["b"]]).value_counts()
    annot = annot[annot["symbol"].isin(genes)].sort_values(["symbol", "go_id"])

    edges.rename(columns={"a": "gene_a", "b": "gene_b"}).to_csv(out / "edges.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": genes, "degree": degree.reindex(genes).to_numpy()}).to_csv(
        out / "genes.tsv", sep="\t", index=False)
    annot.rename(columns={"symbol": "gene"}).to_csv(out / "annotations.tsv", sep="\t", index=False)

    per_gene = annot.groupby("symbol").size()
    summary = {
        "string_links_in_file": int(n_link_rows),
        "string_links_above_threshold": int(len(links)),
        "string_min_score": int(cfg["string_min_score"]),
        "annotated_genes_in_network_before_selection": int(n_genes_all),
        "edges_before_selection": int(n_edges_all),
        "gaf_rows": int(n_gaf_rows),
        "go_pairs_after_filters": int(n_go_pairs),
        "n_genes": len(genes),
        "n_edges": int(len(edges)),
        "median_degree": float(degree.median()),
        "genes_with_annotation": int(per_gene.size),
        "n_terms": int(annot["go_id"].nunique()),
        "median_terms_per_gene": float(per_gene.median()),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
