# Decision log

- 2026-11-22: choose this project for the capstone because it is small, public and CPU only. It is a cut-down version
  of an NMTF project on a server (PPI network, GO annotations, PPMI input, orthogonal NMTF, enrichment of the
  clusters); the full version needs a GPU, 19,000 genes and data that need accounts or are very large.
- 2026-11-22: data are STRING v12 physical links (combined score >= 700) and UniProt-GOA annotations, all downloadable
  without an account. BioGRID, the source of the original project, is large; MSigDB needs a registration.
- 2026-11-22: GO release pinned to the dated archive 2026-05-19, because `current` changes every month. The newer
  releases 2026-06-19 and 2026-08-05 did not have `annotations/goa_human.gaf.gz` at the time of writing.
- 2026-11-22: use the 2,000 best-connected annotated genes, so one run takes about 15 s and the whole analysis about
  3 minutes on 2 CPU cores.
- 2026-11-22: GO terms are used as annotated (no propagation along the ontology), so the 31 MB `go-basic.obo` is not
  needed. Term names are not shown, only GO identifiers.
- 2026-11-22: main comparison is the share of GO terms enriched and the agreement between seeds (ARI). The share of
  enriched clusters and genes is about 100% for both inputs, so it cannot separate them (see the report).
- 2026-11-22: platforms are linux-64, osx-arm64 and osx-64 only; Windows is not tested.
