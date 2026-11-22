"""Embed the gene network with orthogonal NMTF and test the clusters for GO enrichment.

Question: does the random-walk PPMI of the network give functionally more coherent clusters
than the plain adjacency matrix? For each input and each random seed:
    X (PPMI or adjacency)  ->  X ~ P S G^T  ->  clusters (argmax of P)  ->  GO enrichment

Reads   data/processed/{edges,genes,annotations}.tsv, config/config.yaml
Writes  results/metrics.csv          one row per input and seed
        results/convergence.csv      error after every iteration
        results/stability.csv        adjusted Rand index between the clusterings of two seeds
        results/top_terms.csv        best GO term of every enriched cluster (seed 0)

    python scripts/run_embeddings.py [--config config/config.yaml] [--data data/processed] [--out results]
"""

import argparse
import itertools
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

import nmtf_lib as lib

ROOT = Path(__file__).resolve().parent.parent


def load_inputs(data_dir):
    """Adjacency matrix (dense 0/1), gene list, and the sparse gene x GO-term matrix."""
    genes = pd.read_csv(data_dir / "genes.tsv", sep="\t")["gene"].tolist()
    index = {g: i for i, g in enumerate(genes)}
    edges = pd.read_csv(data_dir / "edges.tsv", sep="\t")
    r, c = edges["gene_a"].map(index).to_numpy(), edges["gene_b"].map(index).to_numpy()
    A = sp.coo_matrix((np.ones(len(edges)), (r, c)), shape=(len(genes),) * 2).toarray()
    A = np.maximum(A, A.T)
    annotations = pd.read_csv(data_dir / "annotations.tsv", sep="\t")
    terms = sorted(annotations["go_id"].unique())
    tindex = {t: j for j, t in enumerate(terms)}
    rows = annotations["gene"].map(index).to_numpy()
    cols = annotations["go_id"].map(tindex).to_numpy()
    annot = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(genes), len(terms))).tocsr()
    return A, genes, annot, terms


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--config", default=ROOT / "config" / "config.yaml")
    parser.add_argument("--data", default=ROOT / "data" / "processed")
    parser.add_argument("--out", default=ROOT / "results")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    A, genes, annot, terms = load_inputs(Path(args.data))
    net, nm, en = cfg["network"], cfg["nmtf"], cfg["enrichment"]
    inputs = {"PPMI": lib.ppmi(A, net["ppmi_window"], net["ppmi_neg"]), "Adjacency": A}
    print(f"{len(genes)} genes, {int(A.sum() / 2)} edges, {len(terms)} GO terms", flush=True)

    metrics, convergence, labels, top_terms = [], [], {}, []
    for name, X in inputs.items():
        for seed in nm["seeds"]:
            start = time.time()
            P, S, G, history = lib.nmtf_orth(X, nm["k"], nm["n_iter"], seed)
            labels[(name, seed)] = lib.clusters_from_factor(P)
            axes = lib.clusters_from_factor(G)
            enr = lib.enrichment(labels[(name, seed)], annot, terms, en["min_term_size"],
                                 en["max_term_size"], en["fdr"])
            enr_axes = lib.enrichment(axes, annot, terms, en["min_term_size"],
                                      en["max_term_size"], en["fdr"])
            row = {"input": name, "seed": seed,
                   "rel_error": history["rel_error"].iloc[-1],
                   "orth_error": history["orth_error"].iloc[-1],
                   "seconds": round(time.time() - start, 1)}
            row.update(lib.enrichment_summary(enr))
            row.update({f"axes_{k}": v for k, v in lib.enrichment_summary(enr_axes).items()})
            metrics.append(row)
            convergence.append(history.assign(input=name, seed=seed))
            if seed == nm["seeds"][0]:
                sig = enr[enr["significant"]].groupby("cluster", as_index=False).first()
                top_terms.append(sig.assign(input=name))
            print(f"  {name:9s} seed {seed}: error {row['rel_error']:.3f}, "
                  f"{row['pct_genes_in_enriched_clusters']:.1f}% of genes in enriched clusters "
                  f"({row['seconds']} s)", flush=True)

    stability = []
    for name in inputs:
        for s1, s2 in itertools.combinations(nm["seeds"], 2):
            stability.append({"input": name, "seed_a": s1, "seed_b": s2,
                              "ari": lib.adjusted_rand_index(labels[(name, s1)], labels[(name, s2)])})

    pd.DataFrame(metrics).to_csv(out / "metrics.csv", index=False)
    pd.concat(convergence).to_csv(out / "convergence.csv", index=False)
    pd.DataFrame(stability).to_csv(out / "stability.csv", index=False)
    pd.concat(top_terms).to_csv(out / "top_terms.csv", index=False)
    print(f"wrote {out}/metrics.csv, convergence.csv, stability.csv, top_terms.csv")


if __name__ == "__main__":
    main()
