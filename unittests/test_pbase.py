"""Tests for base.pbase — utility functions."""

import os
import pytest
from base.pbase import prettyFileName, openWithRetry


class TestPrettyFileName:
    """Tests for the prettyFileName path-shortening function."""

    def test_replaces_home(self, monkeypatch, tmp_path):
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.chdir(tmp_path)
        result = prettyFileName(str(tmp_path / "foo.py"))
        assert result == "./foo.py"

    def test_replaces_cwd(self, monkeypatch):
        monkeypatch.setenv("HOME", "/home/testuser")
        monkeypatch.chdir("/tmp")
        result = prettyFileName("/tmp/test.py")
        assert result == "./test.py"

    def test_no_replacement_needed(self, monkeypatch):
        monkeypatch.setenv("HOME", "/home/testuser")
        monkeypatch.chdir("/tmp")
        result = prettyFileName("/var/log/syslog")
        assert result == "/var/log/syslog"


class TestOpenWithRetry:
    """Tests for the openWithRetry file-opening function."""

    def testOpensExistingFile(self, tmp_file):
        with open(tmp_file, "w") as f:
            f.write("hello")
        with openWithRetry(tmp_file) as f:
            assert f.read() == "hello"

    def testRaisesOnMissingFile(self):
        with pytest.raises(FileNotFoundError):
            openWithRetry("/nonexistent/path/file.txt")

    def testWriteMode(self, tmp_file):
        with openWithRetry(tmp_file, "w") as f:
            f.write("world")
        assert open(tmp_file).read() == "world"

    def testRetriesEventuallyFail(self):
        """With retries=0, should fail immediately."""
        with pytest.raises(Exception):
            openWithRetry("/nonexistent/file", retries=0)
