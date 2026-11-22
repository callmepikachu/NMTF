"""Tests of the NMTF building blocks on small synthetic networks with known structure."""

import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import nmtf_lib as lib  # noqa: E402


def planted_network(n_per=25, n_blocks=4, p_in=0.5, p_out=0.02, seed=0):
    """Stochastic block model: dense inside the blocks, sparse between them."""
    rng = np.random.default_rng(seed)
    n = n_per * n_blocks
    truth = np.repeat(np.arange(n_blocks), n_per)
    prob = np.where(truth[:, None] == truth[None, :], p_in, p_out)
    A = (rng.random((n, n)) < prob).astype(float)
    A = np.triu(A, 1)
    return A + A.T, truth


def test_ppmi_is_symmetric_and_non_negative():
    A, _ = planted_network()
    X = lib.ppmi(A, window=5)
    assert X.shape == A.shape
    assert X.min() >= 0
    assert np.allclose(X, X.T)


def test_ppmi_isolated_node_gives_zero_row():
    A, _ = planted_network()
    A[0, :] = 0
    A[:, 0] = 0
    X = lib.ppmi(A, window=5)
    assert np.all(X[0] == 0) and np.isfinite(X).all()


def test_nmtf_lowers_the_error_and_keeps_factors_non_negative():
    A, _ = planted_network()
    X = lib.ppmi(A, window=5)
    P, S, G, history = lib.nmtf_orth(X, k=4, n_iter=100, seed=0)
    assert P.min() >= 0 and S.min() >= 0 and G.min() >= 0
    assert history["rel_error"].iloc[-1] < history["rel_error"].iloc[0]
    assert history["orth_error"].iloc[-1] < history["orth_error"].iloc[0]


def test_nmtf_is_deterministic_for_a_seed_and_differs_between_seeds():
    A, _ = planted_network()
    X = lib.ppmi(A, window=5)
    P1, *_ = lib.nmtf_orth(X, k=4, n_iter=20, seed=3)
    P2, *_ = lib.nmtf_orth(X, k=4, n_iter=20, seed=3)
    P3, *_ = lib.nmtf_orth(X, k=4, n_iter=20, seed=4)
    assert np.array_equal(P1, P2)
    assert not np.allclose(P1, P3)


def test_nmtf_recovers_planted_modules():
    A, truth = planted_network()
    X = lib.ppmi(A, window=5)
    P, S, G, _ = lib.nmtf_orth(X, k=4, n_iter=200, seed=0)
    ari = lib.adjusted_rand_index(truth, lib.clusters_from_factor(P))
    assert ari > 0.8


def test_adjusted_rand_index_extremes():
    a = np.array([0, 0, 1, 1, 2, 2])
    assert lib.adjusted_rand_index(a, a) == 1.0
    assert lib.adjusted_rand_index(a, np.array([5, 5, 7, 7, 9, 9])) == 1.0  # names do not matter
    rng = np.random.default_rng(0)
    big = rng.integers(0, 5, 2000)
    assert abs(lib.adjusted_rand_index(big, rng.integers(0, 5, 2000))) < 0.05


def test_bh_is_monotone_and_bounded():
    p = np.array([0.001, 0.04, 0.03, 0.5, 0.9])
    q = lib.bh(p)
    assert np.all(q >= p) and np.all(q <= 1)
    order = np.argsort(p)
    assert np.all(np.diff(q[order]) >= 0)
    # textbook value: smallest of 5 p-values scaled by 5/1
    assert np.isclose(q[0], 0.005)


def test_enrichment_finds_the_term_that_matches_a_cluster():
    n_per, n_blocks = 30, 4
    labels = np.repeat(np.arange(n_blocks), n_per)
    n = len(labels)
    # one term per block, containing 25 of that block's 30 genes, plus 2 stray genes elsewhere
    cols = []
    for b in range(n_blocks):
        col = np.zeros(n)
        col[b * n_per: b * n_per + 25] = 1
        col[(b * n_per + 40) % n] = 1
        col[(b * n_per + 70) % n] = 1
        cols.append(col)
    annot = sp.csr_matrix(np.array(cols).T)
    terms = [f"GO:{b}" for b in range(n_blocks)]
    df = lib.enrichment(labels, annot, terms, min_size=5, max_size=100)
    top = df[df["significant"]]
    assert set(top["cluster"]) == set(range(n_blocks))
    best = top.groupby("cluster")["term"].first()  # rows are sorted by q, so this is the best term
    assert all(term == f"GO:{cluster}" for cluster, term in best.items())
    summary = lib.enrichment_summary(df)
    assert summary["pct_clusters_enriched"] == 100
    assert summary["pct_genes_in_enriched_clusters"] == 100


def test_enrichment_random_annotation_finds_little():
    rng = np.random.default_rng(1)
    labels = rng.integers(0, 6, 300)
    annot = sp.csr_matrix((rng.random((300, 40)) < 0.1).astype(float))
    df = lib.enrichment(labels, annot, [f"T{i}" for i in range(40)], min_size=5, max_size=200)
    assert lib.enrichment_summary(df)["pct_clusters_enriched"] <= 34
