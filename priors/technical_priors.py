""" Priors for technical constraints like convex hulls
"""

from pbase.constants import LSP, X2Z, X1W, X3Z, X2W
from ptools.sparticleNames import namer

def outsideCMS18004Penalty ( protomodel ) -> tuple:
    large_penalty = 20.
    if protomodel.masses[LSP] < 70:
        return large_penalty, f"{large_penalty} for LSP mass at {protomodel.masses[LSP]} (CMS-SUS-18-004)"
    if not X2Z in protomodel.masses:
        return large_penalty, f"{large_penalty} for no X2Z in {protomodel.masses} (CMS-SUS-18-004)"
    if protomodel.masses[X2Z]<120:
        return large_penalty, f"{large_penalty} for X2Z mass at {protomodel.masses[X2Z]} (CMS-SUS-18-004)"
    if not X1W in protomodel.masses:
        return large_penalty, f"{large_penalty} for no X1W particle in {protomodel.masses} (CMS-SUS-18-004)"
    if protomodel.masses[X1W]<120:
        return large_penalty, f"{large_penalty} for X1W mass at {protomodel.masses[X1W]} (CMS-SUS-18-004)"
    return 0, f"0 for staying within CMS-SUS-18-004 ranges"

