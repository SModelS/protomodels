"""Tests for ptools.sparticleNames — particle name resolution."""

import pytest
from ptools.sparticleNames import SParticleNames


@pytest.fixture
def xids_namer():
    """SParticleNames in default X-notation."""
    return SParticleNames(susy=False)


@pytest.fixture
def susy_namer():
    """SParticleNames in SUSY notation."""
    return SParticleNames(susy=True)


class TestInit:
    """Tests for __init__."""

    def test_default_is_xids(self):
        namer = SParticleNames()
        assert namer.susy is False

    def test_susy_mode(self, susy_namer):
        assert susy_namer.susy is True


class TestIsSM:
    """Tests for isSM."""

    def test_sm_particles(self, xids_namer):
        for pid in [1, 2, 3, 4, 5, 6, 11, 13, 15, 21, 22, 23, 24, 25]:
            assert xids_namer.isSM(pid) is True, f"PID {pid} should be SM"

    def test_bsm_particles(self, xids_namer):
        for pid in [1000021, 1000022, 1000023, 1000024, 1000001]:
            assert xids_namer.isSM(pid) is False, f"PID {pid} should be BSM"

    def test_negative_pids(self, xids_namer):
        assert xids_namer.isSM(-6) is True
        assert xids_namer.isSM(-1000022) is False


class TestName:
    """Tests for name."""

    def test_known_sm_pid(self, xids_namer):
        assert xids_namer.name(24) == "W"

    def test_known_bsm_pid(self, xids_namer):
        name = xids_namer.name(1000022)
        assert "X" in name or "chi" in name

    def test_tuple_of_pids(self, xids_namer):
        result = xids_namer.name((1000022, 1000024))
        assert isinstance(result, str)
        assert "," in result or " " in result

    def test_string_input(self, xids_namer):
        assert xids_namer.name("W") == "W"

    def test_unknown_pid_returns_str(self, xids_namer):
        result = xids_namer.name(9999999)
        assert result == "9999999"

    def test_none_returns_question_mark(self, xids_namer):
        assert xids_namer.name(None) == "?"


class TestAsciiName:
    """Tests for asciiName."""

    def test_gluino(self, xids_namer):
        name = xids_namer.asciiName(1000021)
        assert "g" in name.lower() or "Xg" in name

    def test_lsp(self, xids_namer):
        name = xids_namer.asciiName(1000022)
        assert "X" in name

    def test_list_input(self, xids_namer):
        result = xids_namer.asciiName([1000022, 1000024])
        assert isinstance(result, str)

    def test_no_braces(self, xids_namer):
        """asciiName should strip LaTeX formatting."""
        name = xids_namer.asciiName(1000022)
        assert "{" not in name
        assert "}" not in name
        assert "\\" not in name


class TestTexName:
    """Tests for texName."""

    def test_basic(self, xids_namer):
        name = xids_namer.texName(1000022, addDollars=True)
        assert name.startswith("$")
        assert name.endswith("$")

    def test_no_dollars(self, xids_namer):
        name = xids_namer.texName(1000022, addDollars=False)
        assert not name.startswith("$")

    def test_with_sign(self, xids_namer):
        name = xids_namer.texName(1000024, addSign=True)
        assert "^{+}" in name or "chi" in name or "X" in name

    def test_tuple_input(self, xids_namer):
        result = xids_namer.texName((1000022, 1000024), addDollars=True)
        assert isinstance(result, str)


class TestPid:
    """Tests for pid lookup by name."""

    def test_known_name(self, xids_namer):
        pid = xids_namer.pid("W")
        assert pid == 24

    def test_integer_passthrough(self, xids_namer):
        assert xids_namer.pid(42) == 42

    def test_unknown_name(self, xids_namer):
        result = xids_namer.pid("nonexistent_particle_xyz")
        assert result is None

    def test_unsigned(self, xids_namer):
        pid = xids_namer.pid("t-", signed=False)
        assert pid == 6

    def test_comma_separated(self, xids_namer):
        result = xids_namer.pid("W,Z")
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestHas:
    """Tests for has."""

    def test_existing_name(self, xids_namer):
        assert xids_namer.has("W") is True

    def test_existing_pid(self, xids_namer):
        assert xids_namer.has(24) is True

    def test_nonexistent(self, xids_namer):
        assert xids_namer.has("nonexistent") is False


class TestParticleType:
    """Tests for particleType."""

    def test_squark(self, xids_namer):
        assert xids_namer.particleType(1000001) == "q"

    def test_sbottom(self, xids_namer):
        assert xids_namer.particleType(1000005) == "b"

    def test_stop(self, xids_namer):
        assert xids_namer.particleType(1000006) == "t"

    def test_gluino(self, xids_namer):
        assert xids_namer.particleType(1000021) == "g"

    def test_neutralino(self, xids_namer):
        assert xids_namer.particleType(1000022) == "n"

    def test_chargino(self, xids_namer):
        assert xids_namer.particleType(1000024) == "n"

    def test_slepton(self, xids_namer):
        assert xids_namer.particleType(1000011) == "l"

    def test_sneutrino(self, xids_namer):
        assert xids_namer.particleType(1000012) == "l"


class TestSusyNames:
    """Tests for SUSY notation mode."""

    def test_gluino_susy(self, susy_namer):
        name = susy_namer.asciiName(1000021)
        assert "~" in name or "tilde" in name.lower() or "g" in name

    def test_neutralino_susy(self, susy_namer):
        name = susy_namer.asciiName(1000022)
        assert isinstance(name, str)
        assert len(name) > 0
