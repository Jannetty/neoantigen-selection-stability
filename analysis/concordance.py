"""Concordance between neoantigen candidate sets from different predictors.

The selection step lives here rather than in the workflow on purpose: selection
logic is the object of study, so it should be explicit, versioned and tested
rather than buried in a tool's defaults.

Candidate identity is (chromosome, position, ref, alt, hla_allele, peptide).
Two predictors "agree" on a candidate when both place it in their top N.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.stats import kendalltau

#: Candidates encoded by a real individualized neoantigen product.
DEFAULT_TOP_N = 34

#: Columns that jointly identify a candidate epitope.
KEY_COLUMNS = ["chromosome", "start", "reference", "variant", "hla_allele", "mt_epitope_seq"]


@dataclass(frozen=True)
class OverlapResult:
    """Top-N agreement between two predictors for one sample."""

    sample: str
    predictor_a: str
    predictor_b: str
    top_n: int
    n_shared: int

    @property
    def jaccard(self) -> float:
        union = 2 * self.top_n - self.n_shared
        return self.n_shared / union if union else float("nan")

    @property
    def fraction_shared(self) -> float:
        return self.n_shared / self.top_n if self.top_n else float("nan")


def load_pvacseq_tsv(path: str | Path, predictor: str, sample: str) -> pd.DataFrame:
    """Read one pVACseq ``*.filtered.tsv`` into a normalized frame.

    pVACseq column names vary across versions, so normalize to lowercase with
    underscores and fail loudly on a missing key column rather than silently
    producing an empty intersection later.
    """
    frame = pd.read_csv(path, sep="\t")
    frame.columns = [c.strip().lower().replace(" ", "_") for c in frame.columns]

    missing = [c for c in KEY_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(
            f"{path}: missing expected column(s) {missing}. "
            f"Found: {sorted(frame.columns)[:20]}..."
        )

    frame["predictor"] = predictor
    frame["sample"] = sample
    return frame


def rank_candidates(frame: pd.DataFrame, score_column: str = "median_mt_score") -> pd.DataFrame:
    """Add a 1-based rank by predicted affinity, lower score ranking better."""
    if score_column not in frame.columns:
        raise ValueError(f"no score column {score_column!r} in {sorted(frame.columns)[:20]}")
    ranked = frame.sort_values(score_column, ascending=True).copy()
    ranked["rank"] = range(1, len(ranked) + 1)
    return ranked


def candidate_keys(frame: pd.DataFrame, top_n: int | None = None) -> set[tuple]:
    """Identity tuples for the top ``top_n`` candidates (all, if ``None``)."""
    subset = frame if top_n is None else frame.nsmallest(top_n, "rank")
    return set(map(tuple, subset[KEY_COLUMNS].astype(str).itertuples(index=False, name=None)))


def pairwise_overlap(
    ranked: dict[str, pd.DataFrame], sample: str, top_n: int = DEFAULT_TOP_N
) -> list[OverlapResult]:
    """Top-N overlap for every predictor pair. Answers question 1."""
    results = []
    for a, b in combinations(sorted(ranked), 2):
        shared = candidate_keys(ranked[a], top_n) & candidate_keys(ranked[b], top_n)
        results.append(OverlapResult(sample, a, b, top_n, len(shared)))
    return results


def rank_concordance(ranked: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Kendall tau over candidates both predictors scored. Answers question 2.

    Restricted to the shared candidate set: a predictor that never emitted a
    peptide has no opinion on its rank, which is not the same as ranking it last.
    """
    rows = []
    for a, b in combinations(sorted(ranked), 2):
        left = ranked[a].set_index(
            pd.MultiIndex.from_frame(ranked[a][KEY_COLUMNS].astype(str))
        )["rank"]
        right = ranked[b].set_index(
            pd.MultiIndex.from_frame(ranked[b][KEY_COLUMNS].astype(str))
        )["rank"]
        shared = left.index.intersection(right.index)
        if len(shared) < 3:
            rows.append({"predictor_a": a, "predictor_b": b, "n": len(shared),
                         "kendall_tau": float("nan"), "p_value": float("nan")})
            continue
        tau, p = kendalltau(left.loc[shared], right.loc[shared])
        rows.append({"predictor_a": a, "predictor_b": b, "n": len(shared),
                     "kendall_tau": tau, "p_value": p})
    return pd.DataFrame(rows)


def overlap_by_allele(
    ranked: dict[str, pd.DataFrame], top_n: int = DEFAULT_TOP_N
) -> pd.DataFrame:
    """Top-N overlap split by HLA allele. Answers question 3.

    Expect lower concordance for alleles with less training data behind them.
    """
    rows = []
    alleles = sorted({a for f in ranked.values() for a in f["hla_allele"].unique()})
    for allele in alleles:
        per_allele = {
            name: rank_candidates(frame[frame["hla_allele"] == allele].drop(columns="rank"))
            for name, frame in ranked.items()
            if (frame["hla_allele"] == allele).any()
        }
        for a, b in combinations(sorted(per_allele), 2):
            shared = candidate_keys(per_allele[a], top_n) & candidate_keys(per_allele[b], top_n)
            rows.append({"hla_allele": allele, "predictor_a": a, "predictor_b": b,
                         "n_shared": len(shared), "top_n": top_n})
    return pd.DataFrame(rows)


def filter_sensitivity(
    unfiltered: pd.DataFrame, filtered: pd.DataFrame, top_n: int = DEFAULT_TOP_N
) -> int:
    """Candidates that change when the expression filter is applied.

    Compare against the predictor-swap numbers from :func:`pairwise_overlap` to
    answer question 4: which choice moves the design more?
    """
    before = candidate_keys(rank_candidates(unfiltered), top_n)
    after = candidate_keys(rank_candidates(filtered), top_n)
    return len(before - after)


def release_impact(result: OverlapResult) -> str:
    """One sentence framing an overlap number as a change-control statement."""
    changed = result.top_n - result.n_shared
    return (
        f"Switching {result.predictor_a} to {result.predictor_b} for {result.sample} "
        f"changes {changed} of {result.top_n} selected candidates "
        f"({changed / result.top_n:.0%} of the design)."
    )
