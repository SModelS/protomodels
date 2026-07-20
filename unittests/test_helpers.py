"""Tests for ptools.helpers — utility functions."""

import math
import os
import pytest
import numpy as np

from ptools.helpers import (
    formatObject,
    mkdir,
    nround,
    py_dump,
    py_dumps,
    prettyPrint,
    computeZFromP,
    roughZValue,
    computePAnalytically,
    computeP,
    simplifyUnixPath,
    lrEquiv,
    reorder_list_by_dict,
    repr_double_quotes,
    stripUnits,
)


class TestFormatObject:
    """Tests for formatObject."""

    def test_none_returns_string(self):
        assert formatObject(None) == "None"

    def test_float_format(self):
        assert formatObject(3.14159, ".2f") == "3.14"

    def test_int_format_string(self):
        assert formatObject(42, "1") == "42"

    def test_int_format_int(self):
        assert formatObject(42, 1) == "42.0"

    def test_large_number(self):
        assert formatObject(1234567.89, ".0f") == "1234568"


class TestMkdir:
    """Tests for mkdir."""

    def test_creates_directory(self, tmp_dir):
        new_dir = os.path.join(tmp_dir, "newdir")
        result = mkdir(new_dir)
        assert result is True
        assert os.path.isdir(new_dir)

    def test_existing_directory(self, tmp_dir):
        result = mkdir(tmp_dir)
        assert result is False

    def test_empty_string(self):
        result = mkdir("")
        assert result is False

    def test_nested_creation_fails(self, tmp_dir):
        """mkdir does not use os.makedirs, so nested creation raises."""
        nested = os.path.join(tmp_dir, "a", "b")
        with pytest.raises(FileNotFoundError):
            mkdir(nested)


class TestNround:
    """Tests for nround."""

    def test_rounds_float(self):
        assert nround(3.14159, 2) == 3.14

    def test_none_passthrough(self):
        assert nround(None, 2) is None

    def test_rounds_to_zero_decimals(self):
        assert nround(3.7, 0) == 4.0


class TestPyDumps:
    """Tests for py_dumps (pretty-print nested structures)."""

    def test_empty_dict(self):
        assert py_dumps({}) == "{}"

    def test_empty_list(self):
        assert py_dumps([]) == "[]"

    def test_empty_tuple(self):
        assert py_dumps(()) == "()"

    def test_simple_dict(self):
        result = py_dumps({"a": 1})
        assert "'a': 1" in result or '"a": 1' in result

    def test_nested_structure(self):
        result = py_dumps({"a": [1, 2], "b": {"c": 3}})
        assert "'a'" in result or '"a"' in result
        assert "'b'" in result or '"b"' in result

    def test_tuple_key(self):
        result = py_dumps({(1, 2): "value"})
        assert "value" in result

    def test_stop_at_level(self):
        result = py_dumps({"a": {"b": {"c": 1}}}, stop_at_level=1)
        # Should not have deeply nested indentation
        assert result.count("\n") < 5

    def test_nan_handling(self):
        result = py_dumps(float("nan"))
        assert "nan" in result.lower()

    def test_inf_handling(self):
        result = py_dumps(float("inf"))
        assert "inf" in result.lower()


class TestPyDump:
    """Tests for py_dump (write to file)."""

    def test_write_to_filename(self, tmp_file):
        py_dump({"key": "value"}, tmp_file)
        with open(tmp_file) as f:
            content = f.read()
        assert "key" in content

    def test_write_to_handle(self, tmp_file):
        with open(tmp_file, "w") as f:
            py_dump({"hello": 42}, f)
        with open(tmp_file) as f:
            content = f.read()
        assert "42" in content


class TestPrettyPrint:
    """Tests for prettyPrint."""

    def test_none(self):
        assert prettyPrint(None) == "None"

    def test_float(self):
        result = prettyPrint(3.14159, ndecimals=2)
        assert result == "3.14"

    def test_integer(self):
        assert prettyPrint(42) == "42"

    def test_list(self):
        result = prettyPrint([1.0, 2.0, 3.0])
        assert "1.00" in result

    def test_list_with_maxrows(self):
        result = prettyPrint([1.0, 2.0, 3.0], maxrows=2)
        assert "3.00" not in result

    def test_dict(self):
        result = prettyPrint({"a": 1.0})
        assert "1.00" in result

    def test_string(self):
        assert prettyPrint("hello") == "hello"


class TestComputeZFromP:
    """Tests for computeZFromP."""

    def test_known_value(self):
        # p=0.05 -> Z ≈ 1.645
        z = computeZFromP(0.05)
        assert abs(z - 1.645) < 0.01

    def test_5sigma(self):
        # p ≈ 2.87e-7 for 5 sigma
        z = computeZFromP(2.87e-7)
        assert abs(z - 5.0) < 0.01

    def test_one_sigma(self):
        z = computeZFromP(0.1587)
        assert abs(z - 1.0) < 0.01


class TestRoughZValue:
    """Tests for roughZValue."""

    def test_no_excess(self):
        z = roughZValue(10, 10, 1)
        assert abs(z) < 0.1

    def test_positive_excess(self):
        z = roughZValue(20, 10, 1)
        assert z > 0

    def test_negative_excess(self):
        z = roughZValue(5, 10, 1)
        assert z < 0


class TestComputePAnalytically:
    """Tests for the analytical p-value computation."""

    def test_symmetric_case(self):
        p = computePAnalytically(10, 10, 1)
        assert 0.3 < p < 0.7

    def test_high_observation(self):
        p = computePAnalytically(20, 10, 1)
        assert p < 0.05

    def test_low_observation(self):
        p = computePAnalytically(2, 10, 1)
        assert p > 0.95

    def test_zero_observation(self):
        p = computePAnalytically(0, 5, 1)
        assert p > 0.5

    def test_returned_value_is_float(self):
        p = computePAnalytically(5, 5, 1)
        assert isinstance(p, float)


class TestComputeP:
    """Tests for computeP (mixed analytical/numerical)."""

    def test_returns_tuple(self):
        result = computeP(5, 5, 1)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_pvalue_and_method(self):
        p, method = computeP(5, 5, 1)
        assert isinstance(p, float)
        assert method in ("analytical", "numerical")

    def test_analytical_forced(self):
        p, method = computeP(5, 5, 1, force="analytical")
        assert method == "analytical"

    def test_extreme_values(self):
        p, _ = computeP(0, 0.007, 0.002)
        assert isinstance(p, float)


class TestSimplifyUnixPath:
    """Tests for simplifyUnixPath."""

    def test_replaces_cwd(self, monkeypatch):
        monkeypatch.chdir("/tmp")
        result = simplifyUnixPath("/tmp/test.py")
        assert result == "./test.py"

    def test_replaces_home(self, monkeypatch):
        monkeypatch.setenv("HOME", "/home/user")
        result = simplifyUnixPath("/home/user/file.py")
        assert result == "~/file.py"

    def test_removes_double_slashes(self):
        result = simplifyUnixPath("/tmp//test.py")
        assert "//" not in result


class TestLrEquiv:
    """Tests for lrEquiv."""

    def test_same_string(self):
        assert lrEquiv("abc", "abc") is True

    def test_prefix_pm(self):
        assert lrEquiv("+-abc", "+-abc") is True

    def test_different_strings(self):
        assert lrEquiv("abc", "xyz") is False

    def test_non_string_input(self):
        assert lrEquiv(123, "abc") is False
        assert lrEquiv("abc", 123) is False

    def test_single_char_equivalent(self):
        assert lrEquiv("+-a", "+-a") is True


class TestReorderByDict:
    """Tests for reorder_list_by_dict."""

    def test_basic_ordering(self):
        result = reorder_list_by_dict(["b", "a"], {"a": "b"})
        assert result.index("a") < result.index("b")

    def test_no_mapping(self):
        lst = [3, 1, 2]
        result = reorder_list_by_dict(lst, {})
        # Should maintain relative order for unmapped items
        assert set(result) == set(lst)

    def test_preserves_all_elements(self):
        lst = ["x", "y", "z"]
        result = reorder_list_by_dict(lst, {"x": "y"})
        assert set(result) == set(lst)

    def test_nonexistent_elements_ignored(self):
        lst = ["a", "b"]
        result = reorder_list_by_dict(lst, {"a": "c"})
        assert set(result) == set(lst)


class TestReprDoubleQuotes:
    """Tests for repr_double_quotes."""

    def test_string(self):
        result = repr_double_quotes("hello")
        assert result == '"hello"'

    def test_int(self):
        assert repr_double_quotes(42) == "42"

    def test_float(self):
        assert repr_double_quotes(3.14) == "3.14"

    def test_nan(self):
        result = repr_double_quotes(float("nan"))
        assert "nan" in result

    def test_inf(self):
        result = repr_double_quotes(float("inf"))
        assert "inf" in result

    def test_list(self):
        result = repr_double_quotes([1, 2])
        assert result.startswith("[")
        assert result.endswith("]")

    def test_dict(self):
        result = repr_double_quotes({"a": 1})
        assert result.startswith("{")

    def test_tuple(self):
        result = repr_double_quotes((1, 2))
        assert result.startswith("(")


class TestStripUnits:
    """Tests for stripUnits."""

    def test_none_input(self):
        with pytest.raises(TypeError):
            stripUnits(None)

    def test_numeric_values(self):
        result = stripUnits([[1.0, 2.0], [3.0]])
        assert result == [[1.0, 2.0], [3.0]]

    def test_mixed_types(self):
        result = stripUnits([[1, "a", 2.0]])
        assert result == [[1, 2.0]]
