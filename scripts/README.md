# scripts

- `download_data.py`: download the raw files listed in `metadata/sources.tsv` and verify their SHA-256 (`pixi run download`).
- `prepare_data.py`: raw files -> gene network, GO annotations, `summary.json` in `data/processed/` (`pixi run prepare`).
- `nmtf_lib.py`: PPMI, orthogonal NMTF, clustering, adjusted Rand index and GO enrichment. Pure functions, tested in `tests/`.
- `run_embeddings.py`: for each input (PPMI, adjacency) and seed, factorise, cluster and test the clusters; writes `results/` (`pixi run embed`).
