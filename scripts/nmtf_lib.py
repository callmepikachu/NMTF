"""The numerical building blocks of the analysis: PPMI, orthogonal NMTF, GO enrichment.

Everything here is a pure function (no files, no global state), so it is easy to test
(tests/test_nmtf_lib.py) and to reuse.

Method (Xenos et al. 2021, Bioinformatics; Doria-Belenguer et al. 2023, Bioinformatics):
    1. network A  ->  PPMI matrix X (random-walk co-occurrence, closed form of NetMF)
    2. X ~ P S G^T  with P, S, G >= 0 and G^T G = I   (orthogonal NMTF, Ding et al. 2006)
    3. every gene gets a cluster: the column where its row of P is largest
    4. every cluster is tested for enrichment in GO terms (hypergeometric test, BH correction)
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import hypergeom


def ppmi(A, window=10, neg=1.0):
    """Random-walk PPMI of a symmetric adjacency matrix (NetMF closed form).

    vol  = sum of all entries of A,  D = diag(row sums)
    S    = P + P^2 + ... + P^window,  P = D^-1 A
    PPMI = max(0, log( vol / (neg * window) * S * D^-1 ))   (entries below 1 become 0)
    The result is symmetrised to remove floating-point asymmetry. Isolated nodes get 0.
    """
    A = np.asarray(A, dtype=np.float64)
    vol = A.sum()
    d = A.sum(axis=1)
    d_inv = np.divide(1.0, d, out=np.zeros_like(d), where=d > 0)
    P = d_inv[:, None] * A
    S = np.zeros_like(A)
    Pr = P.copy()
    for r in range(window):
        S += Pr
        if r < window - 1:
            Pr = Pr @ P
    M = (vol / (neg * window)) * S * d_inv[None, :]
    M = np.log(np.maximum(M, 1.0))
    return (M + M.T) / 2


def _normalize_columns(P, S):
    """Give P unit-length columns and move the scale into S (P S G^T does not change)."""
    norms = np.linalg.norm(P, axis=0) + 1e-12
    P /= norms
    S *= norms[:, None]


def nmtf_orth(X, k, n_iter=300, seed=0, eps=1e-10):
    """Orthogonal NMTF by multiplicative updates:  X ~ P S G^T,  P, S, G >= 0,  G^T G = I.

    Updates (Ding et al. 2006; square roots keep the steps stable):
        G <- G * sqrt( X^T P S      / (G G^T X^T P S) )         (keeps G^T G close to I)
        P <- P * sqrt( X G S^T      / (P S G^T G S^T) )
        S <- S * sqrt( P^T X G      / (P^T P S G^T G) )
    P and G start random (this is where the seed matters); S starts as the least-squares
    solution. Columns of P are rescaled to unit length after each round.

    Returns P (n x k), S (k x k), G (n x k) and a DataFrame with one row per iteration:
    relative reconstruction error ||X - P S G^T||_F / ||X||_F and orthogonality error
    ||G^T G - I||_F.
    """
    X = np.asarray(X, dtype=np.float64)
    if X.min() < 0:
        raise ValueError("nmtf_orth: X must be non-negative")
    n = X.shape[0]
    if X.shape != (n, n):
        raise ValueError("nmtf_orth: X must be square")
    rng = np.random.default_rng(seed)
    P = rng.random((n, k)) + 1e-3
    G = rng.random((n, k)) + 1e-3
    G /= np.linalg.norm(G, axis=0) + 1e-12
    S = np.linalg.pinv(P.T @ P) @ (P.T @ X @ G) @ np.linalg.pinv(G.T @ G)
    S = np.maximum(S, 1e-6 * max(S.max(), 1e-12))
    _normalize_columns(P, S)
    x_norm = np.linalg.norm(X)
    eye = np.eye(k)

    def errors(it):
        rel = np.linalg.norm(X - P @ S @ G.T) / x_norm
        orth = np.linalg.norm(G.T @ G - eye)
        return {"iter": it, "rel_error": float(rel), "orth_error": float(orth)}

    history = [errors(0)]
    for it in range(1, n_iter + 1):
        A = X.T @ (P @ S)
        G *= np.sqrt(A / (G @ (G.T @ A) + eps))
        GtG = G.T @ G
        P *= np.sqrt((X @ (G @ S.T)) / (P @ (S @ GtG @ S.T) + eps))
        S *= np.sqrt((P.T @ (X @ G)) / ((P.T @ P) @ S @ GtG + eps))
        _normalize_columns(P, S)
        history.append(errors(it))
    return P, S, G, pd.DataFrame(history)


def clusters_from_factor(F):
    """Hard clustering: every row goes to the column where it is largest (-1 if the row is 0)."""
    F = np.asarray(F)
    labels = F.argmax(axis=1)
    labels[F.max(axis=1) <= 0] = -1
    return labels


def adjusted_rand_index(a, b):
    """Adjusted Rand index between two clusterings (1 = identical, about 0 = random)."""
    a, b = np.asarray(a), np.asarray(b)
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    table = np.zeros((ai.max() + 1, bi.max() + 1))
    np.add.at(table, (ai, bi), 1)
    comb2 = lambda x: x * (x - 1) / 2  # noqa: E731
    sum_cells = comb2(table).sum()
    sum_rows, sum_cols = comb2(table.sum(1)).sum(), comb2(table.sum(0)).sum()
    expected = sum_rows * sum_cols / comb2(len(a))
    maximum = (sum_rows + sum_cols) / 2
    if maximum == expected:
        return 1.0
    return float((sum_cells - expected) / (maximum - expected))


def bh(p):
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def enrichment(labels, annot, terms, min_size=5, max_size=200, fdr=0.05, min_overlap=2):
    """Hypergeometric enrichment of every cluster for every annotation term.

    labels: (n,) cluster number per gene, -1 = not clustered
    annot:  (n x T) sparse 0/1 matrix, genes x terms (rows aligned with labels)
    terms:  (T,) term names
    Only terms with min_size..max_size genes in the universe are tested; the universe is
    the set of clustered genes with at least one annotation. Cluster-term pairs that
    share fewer than min_overlap genes are not reported, but they count as tests for the
    BH correction. Returns a DataFrame (one row per reported pair) with its attrs.
    """
    labels = np.asarray(labels)
    A = sp.csr_matrix(annot).astype(bool).astype(np.float64)
    in_universe = (labels >= 0) & (np.asarray(A.sum(1)).ravel() > 0)
    A = A[in_universe]
    lab = labels[in_universe]
    N = int(in_universe.sum())
    term_size = np.asarray(A.sum(0)).ravel()
    keep = np.flatnonzero((term_size >= min_size) & (term_size <= max_size))
    A = A[:, keep]
    clusters = np.unique(lab)
    M = sp.csr_matrix((np.ones(N), (np.searchsorted(clusters, lab), np.arange(N))),
                      shape=(len(clusters), N))
    overlap = (M @ A).toarray()
    csize = np.asarray(M.sum(1)).ravel()
    tsize = term_size[keep]
    n_tests = overlap.size
    ci, ti = np.nonzero(overlap >= min_overlap)
    k = overlap[ci, ti]
    p = hypergeom.sf(k - 1, N, tsize[ti], csize[ci])
    # pairs that are not reported count as p = 1 in the BH correction
    q = bh(np.r_[p, np.ones(n_tests - len(p))])[: len(p)]
    names = np.asarray(terms)[keep]
    df = pd.DataFrame({"cluster": clusters[ci], "term": names[ti], "overlap": k.astype(int),
                       "cluster_size": csize[ci].astype(int), "term_size": tsize[ti].astype(int),
                       "universe": N, "p": p, "q": q, "significant": q <= fdr})
    df = df.sort_values("q").reset_index(drop=True)
    df.attrs.update(n_tests=n_tests, universe=N, tested_terms=len(keep),
                    clusters=clusters.tolist(),
                    cluster_sizes=dict(zip(clusters.tolist(), csize.astype(int).tolist())))
    return df


def enrichment_summary(df):
    """Three summary numbers per clustering (the ones Pržulj-group papers report):

    pct_clusters_enriched          clusters enriched in at least one term, %
    pct_genes_in_enriched_clusters genes that sit in those clusters, %
    pct_terms_enriched             tested terms enriched in at least one cluster, %
    """
    sig = df[df["significant"]]
    clusters, sizes = df.attrs["clusters"], df.attrs["cluster_sizes"]
    enriched = set(sig["cluster"])
    return {"n_clusters": len(clusters),
            "pct_clusters_enriched": 100 * len(enriched) / max(len(clusters), 1),
            "pct_genes_in_enriched_clusters":
                100 * sum(sizes[c] for c in enriched) / max(sum(sizes.values()), 1),
            "pct_terms_enriched": 100 * sig["term"].nunique() / max(df.attrs["tested_terms"], 1),
            "n_significant_pairs": int(len(sig)),
            "tested_terms": df.attrs["tested_terms"]}
