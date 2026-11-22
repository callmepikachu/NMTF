# metadata

Small text files that describe the data. They are committed to Git.

- `sources.tsv`: where each raw file is downloaded from.
  - `name`: short label of the download
  - `url`: download address (versioned or dated, so it keeps pointing to the same file)
  - `sha256`: SHA-256 of the file (computed 2026-11-22)
  - `dest`: where the file goes in the project
  - `license`: license of the data
- `checksums.sha256`: the same SHA-256 values in the format of `sha256sum -c metadata/checksums.sha256`
  (macOS before 14: `shasum -a 256 -c metadata/checksums.sha256`).
