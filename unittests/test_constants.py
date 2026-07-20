"""Tests for base.constants — Standard Model particle data."""

import pytest
from base.constants import particleNames, smMasses, smWidths


class TestParticleNames:
    """Tests for the particleNames mapping."""

    def test_contains_expected_particles(self):
        expected = {"W", "Z", "t", "h", "tau", "c", "b"}
        assert set(particleNames.keys()) == expected

    def test_pdg_ids_are_integers(self):
        for name, pid in particleNames.items():
            assert isinstance(pid, int), f"PID for {name} should be int"

    def test_known_pdg_values(self):
        assert particleNames["W"] == 24
        assert particleNames["Z"] == 23
        assert particleNames["t"] == 6
        assert particleNames["h"] == 25
        assert particleNames["tau"] == 15


class TestSmMasses:
    """Tests for the smMasses dictionary."""

    def test_string_keys_present(self):
        for name in particleNames:
            assert name in smMasses, f"Missing mass for {name}"

    def test_pid_keys_populated(self):
        """After module init, pid-based keys should be added."""
        for name, pid in particleNames.items():
            if name in smMasses:
                assert pid in smMasses, f"Missing PID {pid} for {name}"

    def test_mass_values_positive(self):
        for key, mass in smMasses.items():
            assert mass > 0, f"Mass for {key} should be positive, got {mass}"

    def test_top_mass_reasonable(self):
        assert 160 < smMasses["t"] < 190

    def test_w_mass_reasonable(self):
        assert 75 < smMasses["W"] < 90

    def test_z_mass_reasonable(self):
        assert 85 < smMasses["Z"] < 100

    def test_higgs_mass_reasonable(self):
        assert 120 < smMasses["h"] < 130


class TestSmWidths:
    """Tests for the smWidths dictionary."""

    def test_w_width(self):
        assert smWidths["W"] == pytest.approx(2.14)

    def test_z_width(self):
        assert smWidths["Z"] == pytest.approx(2.5)

    def test_pid_based_widths(self):
        assert 24 in smWidths, "W width by PID should exist"
        assert 23 in smWidths, "Z width by PID should exist"

    def test_width_values_positive(self):
        for key, width in smWidths.items():
            assert width > 0, f"Width for {key} should be positive"
