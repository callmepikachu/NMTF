# results

Made by `scripts/run_embeddings.py` (`pixi run embed`). Not committed.

- `metrics.csv`: one row per input (PPMI, Adjacency) and seed: reconstruction and orthogonality error,
  GO enrichment summaries of the clusters (`pct_*`) and of the axes (`axes_pct_*`), run time.
- `convergence.csv`: reconstruction and orthogonality error after every iteration.
- `stability.csv`: adjusted Rand index between the clusterings of two seeds.
- `top_terms.csv`: the best GO term of every enriched cluster (seed 0, both inputs).
