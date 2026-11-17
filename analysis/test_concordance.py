"""Tests for the selection and concordance logic.

The selection step is the object of study, so it gets tests. Synthetic frames
only: no pVACseq run required, so CI stays fast and does not need the 9.4 GB
image.
"""

from __future__ import annotations

import pandas as pd
import pytest

from concordance import (
    KEY_COLUMNS,
    OverlapResult,
    candidate_keys,
    filter_sensitivity,
    overlap_by_allele,
    pairwise_overlap,
    rank_candidates,
    rank_concordance,
    release_impact,
)


def make_frame(peptides, scores, allele="HLA-A*02:01"):
    """Minimal frame with the columns the module requires."""
    return pd.DataFrame(
        {
            "chromosome": ["1"] * len(peptides),
            "start": range(100, 100 + len(peptides)),
            "reference": ["A"] * len(peptides),
            "variant": ["T"] * len(peptides),
            "hla_allele": [allele] * len(peptides),
            "mt_epitope_seq": peptides,
            "median_mt_score": scores,
        }
    )


class TestRankCandidates:
    def test_ranks_lower_score_first(self):
        frame = make_frame(["AAA", "BBB", "CCC"], [500.0, 10.0, 100.0])
        ranked = rank_candidates(frame)
        assert list(ranked["mt_epitope_seq"]) == ["BBB", "CCC", "AAA"]
        assert list(ranked["rank"]) == [1, 2, 3]

    def test_missing_score_column_raises(self):
        frame = make_frame(["AAA"], [1.0]).drop(columns="median_mt_score")
        with pytest.raises(ValueError, match="no score column"):
            rank_candidates(frame)


class TestCandidateKeys:
    def test_respects_top_n(self):
        ranked = rank_candidates(make_frame(["AAA", "BBB", "CCC"], [1.0, 2.0, 3.0]))
        assert len(candidate_keys(ranked, top_n=2)) == 2

    def test_none_returns_all(self):
        ranked = rank_candidates(make_frame(["AAA", "BBB"], [1.0, 2.0]))
        assert len(candidate_keys(ranked, top_n=None)) == 2

    def test_key_is_all_identity_columns(self):
        ranked = rank_candidates(make_frame(["AAA"], [1.0]))
        (key,) = candidate_keys(ranked)
        assert len(key) == len(KEY_COLUMNS)


class TestPairwiseOverlap:
    def test_identical_predictors_agree_completely(self):
        frame = make_frame(["AAA", "BBB", "CCC"], [1.0, 2.0, 3.0])
        ranked = {"p1": rank_candidates(frame), "p2": rank_candidates(frame.copy())}
        (result,) = pairwise_overlap(ranked, "S1", top_n=3)
        assert result.n_shared == 3
        assert result.fraction_shared == 1.0
        assert result.jaccard == 1.0

    def test_reversed_ranking_shares_only_the_middle(self):
        # Same three candidates, opposite order. Top-1 sets are disjoint.
        frame_a = make_frame(["AAA", "BBB", "CCC"], [1.0, 2.0, 3.0])
        frame_b = make_frame(["AAA", "BBB", "CCC"], [3.0, 2.0, 1.0])
        ranked = {"p1": rank_candidates(frame_a), "p2": rank_candidates(frame_b)}
        (result,) = pairwise_overlap(ranked, "S1", top_n=1)
        assert result.n_shared == 0

    def test_one_result_per_pair(self):
        frame = make_frame(["AAA", "BBB"], [1.0, 2.0])
        ranked = {name: rank_candidates(frame.copy()) for name in ("p1", "p2", "p3")}
        assert len(pairwise_overlap(ranked, "S1", top_n=2)) == 3


class TestOverlapResult:
    def test_jaccard_on_half_agreement(self):
        # 5 shared out of 10 each: union is 15, so Jaccard is 1/3.
        result = OverlapResult("S1", "a", "b", top_n=10, n_shared=5)
        assert result.jaccard == pytest.approx(1 / 3)
        assert result.fraction_shared == 0.5

    def test_empty_selection_is_nan(self):
        result = OverlapResult("S1", "a", "b", top_n=0, n_shared=0)
        assert pd.isna(result.jaccard)


class TestRankConcordance:
    def test_identical_ranking_gives_tau_one(self):
        frame = make_frame(["AAA", "BBB", "CCC", "DDD"], [1.0, 2.0, 3.0, 4.0])
        ranked = {"p1": rank_candidates(frame), "p2": rank_candidates(frame.copy())}
        result = rank_concordance(ranked)
        assert result["kendall_tau"].iloc[0] == pytest.approx(1.0)

    def test_reversed_ranking_gives_tau_minus_one(self):
        frame_a = make_frame(["AAA", "BBB", "CCC", "DDD"], [1.0, 2.0, 3.0, 4.0])
        frame_b = make_frame(["AAA", "BBB", "CCC", "DDD"], [4.0, 3.0, 2.0, 1.0])
        ranked = {"p1": rank_candidates(frame_a), "p2": rank_candidates(frame_b)}
        result = rank_concordance(ranked)
        assert result["kendall_tau"].iloc[0] == pytest.approx(-1.0)

    def test_too_few_shared_candidates_is_nan_not_error(self):
        ranked = {
            "p1": rank_candidates(make_frame(["AAA", "BBB"], [1.0, 2.0])),
            "p2": rank_candidates(make_frame(["CCC", "DDD"], [1.0, 2.0])),
        }
        result = rank_concordance(ranked)
        assert pd.isna(result["kendall_tau"].iloc[0])
        assert result["n"].iloc[0] == 0


class TestOverlapByAllele:
    def test_splits_rows_by_allele(self):
        frame_a = pd.concat([
            make_frame(["AAA", "BBB"], [1.0, 2.0], allele="HLA-A*02:01"),
            make_frame(["CCC", "DDD"], [1.0, 2.0], allele="HLA-B*07:02"),
        ])
        ranked = {"p1": rank_candidates(frame_a), "p2": rank_candidates(frame_a.copy())}
        result = overlap_by_allele(ranked, top_n=2)
        assert set(result["hla_allele"]) == {"HLA-A*02:01", "HLA-B*07:02"}


class TestFilterSensitivity:
    def test_counts_candidates_dropped_by_the_filter(self):
        unfiltered = make_frame(["AAA", "BBB", "CCC"], [1.0, 2.0, 3.0])
        filtered = unfiltered[unfiltered["mt_epitope_seq"] != "BBB"].copy()
        assert filter_sensitivity(unfiltered, filtered, top_n=3) == 1

    def test_no_change_when_filter_removes_nothing(self):
        frame = make_frame(["AAA", "BBB"], [1.0, 2.0])
        assert filter_sensitivity(frame, frame.copy(), top_n=2) == 0


class TestReleaseImpact:
    def test_states_change_as_a_fraction_of_the_design(self):
        sentence = release_impact(OverlapResult("S1", "NetMHCpan", "MHCflurry", 34, 20))
        assert "14 of 34" in sentence
        assert "41%" in sentence
