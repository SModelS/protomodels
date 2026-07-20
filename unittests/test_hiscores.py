"""Tests for walker.hiscores — hiscore list management."""

import pytest
from walker.hiscores import Hiscores


class TestSimilarDicts:
    """Tests for Hiscores.similarDicts."""

    def _make_hiscores(self):
        h = object.__new__(Hiscores)
        h.hiscores = [None, None, None]
        return h

    def test_identical_dicts(self):
        h = self._make_hiscores()
        a = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        assert h.similarDicts(a, a) is True

    def test_different_k(self):
        h = self._make_hiscores()
        a = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        b = {"K": 6.0, "TL": 1.0, "masses": {1000021: 500.0}}
        assert h.similarDicts(a, b) is False

    def test_different_tl(self):
        h = self._make_hiscores()
        a = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        b = {"K": 5.0, "TL": 2.0, "masses": {1000021: 500.0}}
        assert h.similarDicts(a, b) is False

    def test_different_masses(self):
        h = self._make_hiscores()
        a = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        b = {"K": 5.0, "TL": 1.0, "masses": {1000021: 600.0}}
        assert h.similarDicts(a, b) is False

    def test_different_mass_keys(self):
        h = self._make_hiscores()
        a = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        b = {"K": 5.0, "TL": 1.0, "masses": {1000022: 500.0}}
        assert h.similarDicts(a, b) is False

    def test_one_k_none(self):
        h = self._make_hiscores()
        a = {"K": None, "TL": 1.0, "masses": {1000021: 500.0}}
        b = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        assert h.similarDicts(a, b) is False

    def test_both_k_none(self):
        h = self._make_hiscores()
        a = {"K": None, "TL": 1.0, "masses": {1000021: 500.0}}
        b = {"K": None, "TL": 1.0, "masses": {1000021: 500.0}}
        assert h.similarDicts(a, b) is True


class TestInsertHiscore:
    """Tests for Hiscores.insertHiscore."""

    def _make_hiscores(self):
        h = object.__new__(Hiscores)
        h.hiscores = [None, None, None]
        h.logdir = "/tmp"
        h.walkerid = 0
        h.module = "test"
        h.printLogMessages = False
        return h

    def test_insert_into_empty(self):
        h = self._make_hiscores()
        hi = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        result, added = h.insertHiscore([], hi)
        assert added is True
        assert len(result) == 1
        assert result[0] == hi

    def test_insert_higher_k(self):
        h = self._make_hiscores()
        existing = [{"K": 3.0, "TL": 0.5, "masses": {1000021: 500.0}}]
        hi = {"K": 5.0, "TL": 1.0, "masses": {1000021: 600.0}}
        result, added = h.insertHiscore(existing, hi)
        assert added is True
        assert result[0] == hi
        assert result[1]["K"] == 3.0

    def test_insert_lower_k(self):
        h = self._make_hiscores()
        existing = [{"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}]
        hi = {"K": 3.0, "TL": 0.5, "masses": {1000021: 600.0}}
        result, added = h.insertHiscore(existing, hi)
        assert added is True
        assert result[0]["K"] == 5.0
        assert result[1]["K"] == 3.0

    def test_duplicate_not_inserted(self):
        h = self._make_hiscores()
        existing = [{"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}]
        hi = {"K": 5.0, "TL": 1.0, "masses": {1000021: 500.0}}
        result, added = h.insertHiscore(existing, hi)
        assert added is False
        assert len(result) == 1

    def test_truncated_to_ten(self):
        h = self._make_hiscores()
        existing = [{"K": float(i), "TL": 0.0, "masses": {1000021: 500.0 + i}} for i in range(10)]
        hi = {"K": 100.0, "TL": 5.0, "masses": {1000021: 1500.0}}
        result, added = h.insertHiscore(existing, hi)
        assert added is True
        assert len(result) <= 10
