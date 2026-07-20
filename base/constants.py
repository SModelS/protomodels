"""Physics constants used throughout protomodels.

Provides Standard Model particle masses, widths, and PDG ID mappings
that are referenced by the builder, manipulator, and initialiser modules.
"""

from typing import Dict, Final

## vector boson masses and widths

particleNames: Final[Dict[str, int]] = {
    "W": 24, "Z": 23, "t": 6, "h": 25, "tau": 15, "c": 4, "b": 5,
}

smMasses: Final[Dict] = {
    "t": 173.0, "W": 80.377, "Z": 91.1876, "h": 125.0,
    "tau": 1.77, "c": 1.2, "b": 5.0,
}

smWidths: Final[Dict[str, float]] = {"W": 2.14, "Z": 2.5}

for _name, _pid in particleNames.items():
    if _name in smMasses:
        smMasses[_pid] = smMasses[_name]
    if _name in smWidths:
        smWidths[_pid] = smWidths[_name]
