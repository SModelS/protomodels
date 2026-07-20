"""Tests for builder.manipulator — manipulator utility methods."""

from unittest.mock import MagicMock
import pytest
from types import SimpleNamespace


class TestGetClosestPair:
    """Tests for Manipulator.getClosestPair (O(n log n) implementation)."""

    def _make_manipulator(self, masses: dict):
        from builder.manipulator import Manipulator
        manip = object.__new__(Manipulator)
        manip.M = SimpleNamespace(masses=masses)
        return manip

    def test_two_pids(self):
        manip = self._make_manipulator({1: 100.0, 2: 200.0})
        pair, dmin = manip.getClosestPair([1, 2])
        assert pair == (1, 2)
        assert dmin == 100.0

    def test_three_pids_closest_middle(self):
        manip = self._make_manipulator({1: 100.0, 2: 150.0, 3: 300.0})
        pair, dmin = manip.getClosestPair([1, 2, 3])
        assert pair == (1, 2)
        assert dmin == 50.0

    def test_returns_none_for_single_pid(self):
        manip = self._make_manipulator({1: 100.0})
        assert manip.getClosestPair([1]) is None

    def test_returns_none_for_empty(self):
        manip = self._make_manipulator({})
        assert manip.getClosestPair([]) is None

    def test_unordered_input(self):
        """PIDs can come in any order; result should still find the closest pair."""
        manip = self._make_manipulator({10: 500.0, 20: 100.0, 30: 150.0})
        pair, dmin = manip.getClosestPair([10, 20, 30])
        assert pair == (20, 30)
        assert dmin == 50.0

    def test_equal_masses(self):
        manip = self._make_manipulator({1: 100.0, 2: 100.0})
        pair, dmin = manip.getClosestPair([1, 2])
        assert dmin == 0.0

    def test_many_pids(self):
        masses = {i: float(i * 10) for i in range(1, 11)}
        manip = self._make_manipulator(masses)
        pair, dmin = manip.getClosestPair(list(range(1, 11)))
        assert pair == (1, 2)
        assert dmin == 10.0
