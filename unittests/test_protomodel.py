"""Tests for builder.protomodel — ProtoModel data class.

These tests cover the pure-Python aspects of ProtoModel that do not
require external physics libraries (smodels, pyslha, etc.).
"""

import pytest


class TestProtoModelModuleStructure:
    """Verify module-level attributes and imports are sane."""

    def test_import_class(self):
        from builder.protomodel import ProtoModel
        assert ProtoModel is not None

    def test_has_expected_class_attributes(self):
        from builder.protomodel import ProtoModel
        assert ProtoModel.LSP == 1000022
        assert isinstance(ProtoModel.SLHATEMPDIR, str)

    def test_has_key_methods(self):
        from builder.protomodel import ProtoModel
        assert callable(getattr(ProtoModel, "unFrozenParticles", None))
        assert callable(getattr(ProtoModel, "frozenParticles", None))
        assert callable(getattr(ProtoModel, "hasAntiParticle", None))
        assert callable(getattr(ProtoModel, "toTuple", None))
        assert callable(getattr(ProtoModel, "dict", None))
        assert callable(getattr(ProtoModel, "copy", None))
        assert callable(getattr(ProtoModel, "almostSameAs", None))
        assert callable(getattr(ProtoModel, "relevantSSMultipliers", None))
        assert callable(getattr(ProtoModel, "describe", None))


class TestProtoModelToTuple:
    """Tests for the static toTuple method."""

    def test_sorting(self):
        from builder.protomodel import ProtoModel
        result = ProtoModel.toTuple(None, 1000024, 1000022)
        # toTuple is an instance method but we can test via a mock
        # Actually it's `self, pid1, pid2` - let's test differently

    def test_sorted_output(self):
        from builder.protomodel import ProtoModel
        # Create a minimal mock to test the logic
        class MockModel:
            pass
        m = MockModel()
        m.toTuple = ProtoModel.toTuple.__get__(m)
        assert m.toTuple(1000024, 1000022) == (1000022, 1000024)
        assert m.toTuple(1000022, 1000022) == (1000022, 1000022)
        assert m.toTuple(1000001, 1000021) == (1000001, 1000021)


class TestProtoModelHasAntiParticle:
    """Tests for hasAntiParticle."""

    def test_self_conjugate_returns_false(self):
        from builder.protomodel import ProtoModel
        for pid in [1000021, 1000022, 1000023, 1000025, 1000035,
                    1000012, 1000014, 1000016,
                    2000012, 2000014, 2000016, 2000021]:
            class MockModel:
                pass
            m = MockModel()
            m.hasAntiParticle = ProtoModel.hasAntiParticle.__get__(m)
            assert m.hasAntiParticle(pid) is False, f"PID {pid} should be self-conjugate"

    def test_charged_returns_true(self):
        from builder.protomodel import ProtoModel
        for pid in [1000024, 1000037, 1000001, 1000006]:
            class MockModel:
                pass
            m = MockModel()
            m.hasAntiParticle = ProtoModel.hasAntiParticle.__get__(m)
            assert m.hasAntiParticle(pid) is True, f"PID {pid} should have antiparticle"
