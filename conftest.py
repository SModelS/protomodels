"""Shared fixtures and configuration for the protomodels test suite."""

import os
import sys
import tempfile
import pytest

# Ensure the project root is on sys.path for imports
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


@pytest.fixture
def tmp_dir():
    """Provide a temporary directory that is cleaned up after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def tmp_file(tmp_dir):
    """Provide a path to a temporary file inside a temporary directory."""
    return os.path.join(tmp_dir, "testfile.txt")
