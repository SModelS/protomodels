"""Tests for base.locker — file locking mechanism."""

import os
import pytest
import tempfile

from base.locker import lock, unlock, lockfile


class TestLockfile:
    """Tests for the lockfile naming function."""

    def test_returns_hidden_lockfile(self):
        assert lockfile("data.txt") == ".data.txt.lock"

    def test_handles_pathlike(self):
        result = lockfile("/tmp/test.dat")
        assert result == ".test.dat.lock"

    def test_strips_directory(self):
        result = lockfile("/some/deep/path/model.pkl")
        assert result == ".model.pkl.lock"


class TestLockUnlock:
    """Tests for the lock/unlock pair."""

    def test_lock_unlock_cycle(self, tmp_file):
        """Locking and unlocking should succeed."""
        with open(tmp_file, "w") as f:
            f.write("data")
        result = lock(tmp_file)
        assert result is True
        unlock(tmp_file)
        # Lock file should be removed
        lf = lockfile(tmp_file)
        assert not os.path.exists(lf)

    def test_lock_nonexistent_file_returns_false(self):
        result = lock("/nonexistent/file.txt")
        assert result is False

    def test_unlock_nonexistent_lock_returns_false(self):
        result = unlock("/nonexistent/file.txt")
        assert result is False

    def test_double_lock(self, tmp_file):
        """Locking twice should not crash."""
        with open(tmp_file, "w") as f:
            f.write("data")
        lock(tmp_file)
        lock(tmp_file)  # Should wait and then force-unlock
        unlock(tmp_file)

    def test_lock_file_created(self, tmp_file):
        with open(tmp_file, "w") as f:
            f.write("data")
        lock(tmp_file)
        lf = lockfile(tmp_file)
        assert os.path.exists(lf)
        unlock(tmp_file)

    def test_lock_file_contains_host_info(self, tmp_file):
        with open(tmp_file, "w") as f:
            f.write("data")
        lock(tmp_file)
        lf = lockfile(tmp_file)
        with open(lf) as f:
            content = f.read()
        assert "'host'" in content
        assert "'time'" in content
        unlock(tmp_file)


class TestLockIntegration:
    """Integration tests for locking mechanism."""

    def test_lock_multiple_files(self, tmp_dir):
        """Locking multiple files should work independently."""
        files = []
        for i in range(3):
            path = os.path.join(tmp_dir, f"file{i}.txt")
            with open(path, "w") as f:
                f.write(f"content{i}")
            files.append(path)

        for path in files:
            assert lock(path) is True

        # All lock files should exist
        for path in files:
            assert os.path.exists(lockfile(path))

        for path in files:
            unlock(path)

        for path in files:
            assert not os.path.exists(lockfile(path))
