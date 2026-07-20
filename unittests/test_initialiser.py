"""Tests for the initialiser.TopNDict — bounded dictionary."""

import pytest
from walker.initialiser import TopNDict


class TestTopNDict:
    """Tests for TopNDict (keeps only top N entries by key)."""

    def test_basic_insert(self):
        d = TopNDict(nmax=5)
        d[1] = "a"
        d[2] = "b"
        assert len(d) == 2
        assert d[1] == "a"
        assert d[2] == "b"

    def test_evicts_smallest_when_full(self):
        d = TopNDict(nmax=3)
        d[10] = "a"
        d[20] = "b"
        d[30] = "c"
        d[5] = "d"  # 5 is smallest, should be evicted
        assert len(d) == 3
        assert 5 not in d
        assert 10 in d
        assert 20 in d
        assert 30 in d

    def test_preserves_larger_keys(self):
        d = TopNDict(nmax=2)
        d[100] = "a"
        d[200] = "b"
        d[50] = "c"  # evicts 50? No, 100 is smallest
        # nmax=2, we have 100,200, then add 50 -> 3 items -> evict smallest (50)
        # Wait, 50 < 100, so 50 should be evicted
        assert len(d) == 2
        assert 50 not in d
        assert 100 in d

    def test_update_method(self):
        d = TopNDict(nmax=5)
        d.update({1: "a", 2: "b", 3: "c"})
        assert len(d) == 3
        assert d[1] == "a"
        assert d[2] == "b"
        assert d[3] == "c"

    def test_default_nmax(self):
        d = TopNDict()
        assert d.nmax == 20

    def test_overwrite_existing_key(self):
        d = TopNDict(nmax=2)
        d[1] = "a"
        d[1] = "b"
        assert len(d) == 1
        assert d[1] == "b"

    def test_many_insertions(self):
        d = TopNDict(nmax=5)
        for i in range(20):
            d[i] = f"val{i}"
        assert len(d) == 5
        # Only the 5 largest keys should remain
        assert set(d.keys()) == {15, 16, 17, 18, 19}
