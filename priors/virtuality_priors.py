#!/usr/bin/env python3
    
from pbase.constants import LSP, X2Z, X1W, X3Z, X2W
from ptools.sparticleNames import namer

""" penalties based on the on/off-shellness of the gauge
bosons emitted in the decays
"""

def penaltyForOnshell ( protomodel ) -> tuple:
    """ compute the penalty for X1W, X2Z having onshell bosons radiating off
    and the reason for it

    :param protomodel: the model the prior is computed for
    :returns: (penalty, reason) tuple
    """
    ## the penalty (in nll) that we give to the models we do not like
    large_penalty = 20.

    mass_gaps = { X1W: 80., X2Z: 90. }

    for ewkino in mass_gaps.keys():
        if not ewkino in protomodel.masses:
            return large_penalty, f"no ewkino {ewkino} in the model"
    dms = {}
    for ewkino, gap in mass_gaps.items():
        dm = protomodel.masses[ewkino]-protomodel.masses[LSP]
        if dm > gap:
            return large_penalty, f"{namer.asciiName(ewkino)} with mass gap {dm:.2f}"
    return 0, f"0 for N2,C1 offshell"
