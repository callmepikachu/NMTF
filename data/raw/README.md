# data/raw

Raw data exactly as downloaded.

- **Never edit these files.** If something needs fixing, fix it in code and write the result to
  `data/processed/`. The download script makes them read-only.
- Get them with `pixi run download`, not by hand, so anyone gets the same files.
- Where each file comes from, its license and its SHA-256 are in `metadata/sources.tsv`.
- The files are not committed to Git (see `.gitignore`); the download script and the checksums are.

Layout after the download:

- `string/9606.protein.physical.links.v12.0.txt.gz`: STRING v12.0 physical links, human (9 MB)
- `string/9606.protein.info.v12.0.txt.gz`: STRING protein identifiers and gene names, human (2 MB)
- `go/goa_human.gaf.gz`: UniProt-GOA annotations of human proteins, GO release 2026-05-19 (15 MB)
