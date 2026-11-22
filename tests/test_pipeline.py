"""End-to-end test of prepare_data.py and run_embeddings.py on tiny fake raw files.

The fake files copy the formats of the real ones (STRING links and protein info, GOA GAF),
so this also tests that the filters (score, GO aspect, evidence, NOT) do what they say.
"""

import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
N_MODULES, PER_MODULE = 3, 40


def write_fake_raw(raw):
    rng = np.random.default_rng(0)
    (raw / "string").mkdir(parents=True)
    (raw / "go").mkdir(parents=True)
    n = N_MODULES * PER_MODULE
    module = np.repeat(np.arange(N_MODULES), PER_MODULE)
    protein = [f"9606.ENSP{i:05d}" for i in range(n)]
    symbol = [f"GENE{i}" for i in range(n)]

    lines = ["protein1 protein2 combined_score"]
    for i in range(n):
        for j in range(i + 1, n):
            same = module[i] == module[j]
            if rng.random() < (0.4 if same else 0.01):
                score = 900 if same else 950
                lines += [f"{protein[i]} {protein[j]} {score}", f"{protein[j]} {protein[i]} {score}"]
    lines += [f"{protein[0]} {protein[1]} 300", f"{protein[1]} {protein[0]} 300"]  # below threshold
    with gzip.open(raw / "string" / "9606.protein.physical.links.v12.0.txt.gz", "wt") as fh:
        fh.write("\n".join(lines) + "\n")

    info = ["#string_protein_id\tpreferred_name\tprotein_size\tannotation"]
    info += [f"{p}\t{s}\t100\tfake protein" for p, s in zip(protein, symbol)]
    with gzip.open(raw / "string" / "9606.protein.info.v12.0.txt.gz", "wt") as fh:
        fh.write("\n".join(info) + "\n")

    gaf = ["!gaf-version: 2.2", "!another comment line"]

    def row(sym, go, qualifier="enables", evidence="IDA", aspect="P"):
        return "\t".join(["UniProtKB", "P0", sym, qualifier, go, "PMID:1", evidence, "", aspect,
                          "name", "", "protein", "taxon:9606", "20260101", "UniProt", "", ""])

    for i in range(n):
        gaf.append(row(symbol[i], f"GO:000{module[i]}"))            # the module's own term
        gaf.append(row(symbol[i], "GO:9999", evidence="IEA"))        # removed: electronic
        gaf.append(row(symbol[i], "GO:8888", aspect="F"))            # removed: not biological process
        gaf.append(row(symbol[i], "GO:7777", qualifier="NOT|involved_in"))  # removed: NOT
    with gzip.open(raw / "go" / "goa_human.gaf.gz", "wt") as fh:
        fh.write("\n".join(gaf) + "\n")


def run(script, *args):
    return subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
                          check=True, capture_output=True, text=True)


def test_prepare_and_embed_on_fake_data(tmp_path):
    raw, processed, results = tmp_path / "raw", tmp_path / "processed", tmp_path / "results"
    write_fake_raw(raw)
    cfg = {"data": {"string_min_score": 700, "go_aspect": "P", "go_exclude_evidence": ["IEA"],
                    "n_genes": 100},
           "network": {"ppmi_window": 5, "ppmi_neg": 1},
           "nmtf": {"k": 3, "n_iter": 60, "seeds": [0, 1]},
           "enrichment": {"min_term_size": 5, "max_term_size": 200, "fdr": 0.05}}
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))

    run("prepare_data.py", "--config", config, "--raw", raw, "--out", processed)
    summary = json.loads((processed / "summary.json").read_text())
    assert summary["n_genes"] <= 100                      # n_genes limit
    assert summary["string_links_above_threshold"] < summary["string_links_in_file"]
    annotations = pd.read_csv(processed / "annotations.tsv", sep="\t")
    assert set(annotations["go_id"]) <= {"GO:0000", "GO:0001", "GO:0002"}   # only IDA, aspect P, no NOT
    edges = pd.read_csv(processed / "edges.tsv", sep="\t")
    assert (edges["gene_a"] < edges["gene_b"]).all()      # undirected, one row per pair
    assert not edges.duplicated(["gene_a", "gene_b"]).any()

    run("run_embeddings.py", "--config", config, "--data", processed, "--out", results)
    metrics = pd.read_csv(results / "metrics.csv")
    assert set(metrics["input"]) == {"PPMI", "Adjacency"} and len(metrics) == 4
    # the fake modules carry their own GO term, so clusters must be enriched
    assert metrics.loc[metrics["input"] == "PPMI", "pct_genes_in_enriched_clusters"].min() > 80
    stability = pd.read_csv(results / "stability.csv")
    assert stability["ari"].between(-1, 1).all()
    assert len(pd.read_csv(results / "convergence.csv")) == 2 * 2 * 61
    assert len(pd.read_csv(results / "top_terms.csv")) > 0
