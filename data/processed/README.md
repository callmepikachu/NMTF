# data/processed

Made from `data/raw/` by `scripts/prepare_data.py` (`pixi run prepare`). Not committed.

- `edges.tsv`: `gene_a`, `gene_b`, `score`. One row per undirected link, STRING combined score >= 700, among the
  analysed genes.
- `genes.tsv`: `gene`, `degree`. The genes of the network (the best-connected annotated genes, largest
  connected component).
- `annotations.tsv`: `gene`, `go_id`. GO biological-process annotations of those genes, without IEA and NOT.
- `summary.json`: the counts that the report quotes.
