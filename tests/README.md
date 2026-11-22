# tests

`pixi run test` runs them in a few seconds.

- `test_nmtf_lib.py`: PPMI, NMTF, ARI, BH and enrichment on small networks with known structure.
- `test_pipeline.py`: runs `prepare_data.py` and `run_embeddings.py` on tiny fake raw files in the formats of the real ones,
  and checks the filters (score, GO aspect, evidence, NOT).
