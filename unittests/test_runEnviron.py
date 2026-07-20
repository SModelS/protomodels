"""Tests for base.runEnviron — RunEnvironment configuration."""

import os
import pytest
import tempfile
import shutil

from base.runEnviron import RunEnviron, dict_diff


class TestDictDiff:
    """Tests for the dict_diff utility function."""

    def test_identical_dicts(self):
        d = {"a": 1, "b": 2}
        assert dict_diff(d, d) == {}

    def test_added_key(self):
        d1 = {"a": 1}
        d2 = {"a": 1, "b": 2}
        diff = dict_diff(d1, d2)
        assert "b" in diff
        assert diff["b"] == ("<missing>", 2)

    def test_removed_key(self):
        d1 = {"a": 1, "b": 2}
        d2 = {"a": 1}
        diff = dict_diff(d1, d2)
        assert "b" in diff
        assert diff["b"] == (2, "<missing>")

    def test_changed_value(self):
        d1 = {"a": 1}
        d2 = {"a": 2}
        diff = dict_diff(d1, d2)
        assert "a" in diff
        assert diff["a"] == (1, 2)

    def test_empty_dicts(self):
        assert dict_diff({}, {}) == {}


class TestRunEnvironDefaults:
    """Tests for RunEnviron.defaults()."""

    def test_returns_dict(self):
        defaults = RunEnviron.defaults()
        assert isinstance(defaults, dict)

    def test_has_expected_keys(self):
        defaults = RunEnviron.defaults()
        expected_keys = {"dbpath", "select", "do_srcombine", "forbiddenparticles",
                         "templateSLHA", "allowN1N1Prod", "susy_mode", "rundir",
                         "strategy", "use_initialiser", "dbversion",
                         "extrapolation_acceptance"}
        assert expected_keys.issubset(set(defaults.keys()))

    def test_default_values(self):
        defaults = RunEnviron.defaults()
        assert defaults["dbpath"] == "official"
        assert defaults["select"] == "all"
        assert defaults["do_srcombine"] is True
        assert defaults["allowN1N1Prod"] is False
        assert defaults["susy_mode"] is False
        assert defaults["forbiddenparticles"] == []
