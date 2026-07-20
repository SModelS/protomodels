#!/usr/bin/env python3

"""Logging mixin for all protomodels classes.

Provides a ``LoggerBase`` class that can be inherited to gain structured
logging to both files and the terminal, with rate-limited message
quenching for repetitive output.
"""

import sys
import time
import os
from typing import Dict, Union

from ptools import helpers

try:
    from smodels_utils.helper.terminalcolors import GREEN, RED, YELLOW, RESET
except ImportError:
    GREEN = RED = YELLOW = RESET = ""

__all__ = ["LoggerBase"]


class LoggerBase:
    """Mixin that provides logging to file and terminal.

    Every subclass automatically gets a per-walker log file under the
    ``logs/`` directory.  Messages can be rate-limited so that
    repetitive lines are printed at most three times.

    :param walkerid: Identifier for the walker producing the log.
    :param verbosity: Minimum severity to print (``"debug"``,
        ``"info"``, ``"warning"``/``"warn"``, ``"error"`` or an int).
    """

    __slots__ = ["walkerid", "module", "logdir"]

    _VERBOSITY_MAP: Dict[str, int] = {
        "error": 40,
        "err": 40,
        "warning": 30,
        "warn": 30,
        "info": 20,
        "debug": 10,
    }

    def __init__(
        self,
        walkerid: Union[str, int] = 0,
        verbosity: Union[int, str] = "info",
    ) -> None:
        self.verbose = self.getVerbosity(verbosity)
        self.walkerid = walkerid
        self.countLogs: Dict[str, int] = {}
        self.printLogMessages: bool = False
        self.logdir: str = "logs/"
        # Derive short module name from fully-qualified class name
        module = str(type(self)).replace("<class '", "").replace("'>", "")
        p1 = module.find(".")
        p2 = module.rfind(".")
        if module.count(".") == 2:
            self.module = module[p1 + 1 : p2]
        else:
            self.module = module[p2 + 1 :]
        helpers.mkdir(self.logdir)

    def getVerbosity(self, verbosity: Union[int, str]) -> int:
        """Translate a verbosity label or integer into a numeric level.

        :param verbosity: e.g. ``"info"`` (→ 20) or ``30``.
        :returns: Numeric verbosity level.
        """
        if isinstance(verbosity, int):
            return verbosity
        key = verbosity.lower()
        if key not in self._VERBOSITY_MAP:
            self.error(f"verbosity {verbosity} unknown")
            sys.exit()
        return self._VERBOSITY_MAP[key]

    def logThrice(self, *args) -> None:
        """Log a message at most three times, then quench."""
        txt = " ".join(map(str, args))
        if txt not in self.countLogs:
            self.countLogs[txt] = 0
        if self.countLogs[txt] < 3:
            self.log(*args)
        if self.countLogs[txt] == 3:
            self.log("(quenching repeating log messages)")
        self.countLogs[txt] += 1

    def error(self, *args) -> None:
        """Log an error message (highlighted in red)."""
        self.highlight("error", *args)

    def warn(self, *args) -> None:
        """Log a warning message (highlighted in yellow)."""
        self.highlight("warn", *args)

    def warning(self, *args) -> None:
        """Alias for :meth:`warn`."""
        self.highlight("warn", *args)

    def info(self, *args) -> None:
        """Log to file and print to screen if verbosity ≥ 20."""
        self.log(*args)
        if self.verbose > 19:
            print(f"[logger] {' '.join(map(str, args))}")

    def highlight(self, msgType: str = "info", *args) -> None:
        """Log a coloured message to screen and to file."""
        col = GREEN
        if msgType.lower() in ("error", "red"):
            col = RED
        elif msgType.lower() in ("warn", "warning", "yellow"):
            col = YELLOW
        elif msgType.lower() in ("green", "info"):
            col = GREEN
        else:
            self.highlight("red", "called highlight without msg type")
        print(
            f"{col}[{self.module}:{time.strftime('%H:%M:%S')}] "
            f"{' '.join(map(str, args))}{RESET}"
        )
        self.log(*args)

    def debug(self, *args) -> None:
        """Log a debug-level message (written to file only)."""
        tmp = list(args)
        for i, arg in enumerate(tmp):
            if isinstance(arg, str):
                tmp[i] = "DEBUG: " + arg
                break
        args = tuple(tmp)
        self.log(*args)

    def pprint(self, *args) -> None:
        """Pretty-print log with quenching after 3 identical messages."""
        line = " ".join(map(str, args))
        if line not in self.countLogs:
            self.countLogs[line] = 0
        self.countLogs[line] += 1
        if self.countLogs[line] == 3:
            print(f"[{self.module}:{self.walkerid}] skipping repeating messages")
            return
        if self.countLogs[line] > 3:
            return
        print(f"[{self.module}:{self.walkerid}] {line}")
        self.prevMessage = line
        self.log(*args)

    def cprint(self, color: str, *args) -> None:
        """Coloured pretty-print with quenching."""
        line = " ".join(map(str, args))
        from smodels_utils.helper.terminalcolors import colordict, RESET as RST

        if line not in self.countLogs:
            self.countLogs[line] = 0
        self.countLogs[line] += 1
        if self.countLogs[line] == 3:
            print(f"[{self.module}:{self.walkerid}] skipping repeating messages")
            return
        if self.countLogs[line] > 3:
            return
        print(
            f"[{self.module}:{self.walkerid}] {colordict[color]}{line}{RST}"
        )
        self.prevMessage = line
        self.log(*args)

    def log(self, *args) -> None:
        """Append a timestamped message to the walker's log file.

        Retries up to 10 times on ``OSError`` to handle transient
        network filesystem failures.
        """
        helpers.mkdir(self.logdir)
        ctr = 0
        while True:
            try:
                with open(
                    f"{self.logdir}/walker_{self.walkerid}.log", "a"
                ) as f:
                    f.write(
                        f"[{self.module}-{time.strftime('%H:%M:%S')}] "
                        f"{' '.join(map(str, args))}\n"
                    )
                if self.printLogMessages:
                    print(f"[{self.module}-log] {' '.join(map(str, args))}")
                return
            except OSError as e:
                ctr += 1
                time.sleep(ctr**2)
                if ctr > 10:
                    raise e


if __name__ == "__main__":
    logger = LoggerBase(0)
    logger.pprint("what now!")
