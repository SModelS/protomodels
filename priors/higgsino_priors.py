#!/usr/bin/env python3
    
import numpy as np
from pbase.constants import LSP, X2Z, X1W, X3Z, X2W
from ptools.sparticleNames import namer

""" A collection of possibly penalty terms for higgsino type scenarios
"""

def mass_splitting(width):
    """Estimate chargino-neutralino mass splitting (GeV).

    Valid below the pion threshold, for a pure wino approximation.
    width: chargino width in GeV.
    """
    return (width / 5.8e-13)**0.2  # GeV

def higgsinoXSecPenalty ( protomodel ) -> tuple:
    """ a penalty for ewkino cross sections away from higgsino hypotheses 
    """
    prod_modes = [ (LSP,LSP), (LSP,X2Z),
                   (X2Z,X2Z), (-X1W,X1W),
                   (LSP,X1W), ( X2Z, X1W ) ]
    penalty = 0
    ssms = protomodel.ssmultipliers
    hino = .25
    wino = 1.0
    reasons = []
    for prod_mode in prod_modes:
        if not prod_mode in ssms:
            penalty += 20
            reasons.append ( f"{20} for {prod_mode} not in ssms" )
            continue
        deviation = ssms[prod_mode]-wino
        if deviation>0.:
            p = float ( abs ( np.log ( 1 + deviation ) ) )
            penalty += p
            reasons.append ( f"{p:.3g} for {namer.asciiName(prod_mode)} deviating from higgsino xsec" )
        else:
            p = float ( abs ( np.log ( (1 - deviation)/4. ) ) )
            penalty += p
            reasons.append ( f"{p:.3g} for {namer.asciiName(prod_mode)} deviating from higgsino xsec" )
    return penalty, reasons

def degenerateEWKinos ( protomodel ) -> tuple:
    """ punish if N2 and C1 arent mass degenerate """
    if not X2Z in protomodel.masses:
        return large_penalty, f"{large_penalty} for no X2Z in {protomodel.masses}"
    if not X1W in protomodel.masses:
        return large_penalty, f"{large_penalty} for no X2Z in {protomodel.masses}"
    dm = abs(protomodel.masses[X1W]-protomodel.masses[X2Z])/\
            (protomodel.masses[X2Z]+protomodel.masses[X1W])
    if dm > 1e-5:
        return large_penalty, f"{large_penalty} for dm(X1W,X2Z)={dm:.3g}"
    return 0, f"0 for dm(X1W,X2Z) zero"

def degenerateNeutralinos ( protomodel ) -> tuple:
    """ punish if N2 and N1 arent mass degenerate """
    if not X2Z in protomodel.masses:
        return large_penalty, f"{large_penalty} for no X2Z in {protomodel.masses}"
    if not LSP in protomodel.masses:
        return large_penalty, f"{large_penalty} for no X2Z in {protomodel.masses}"
    dm = abs(protomodel.masses[LSP]-protomodel.masses[X2Z])
    p = 10*( dm - .0004 )**2
    return p, f"{p:.3g} for dm(LSP,X2Z)"

def charginoMassWidthPenalty ( protomodel ) -> tuple:
    """ penalize if the chargino width doesnt really
    match with the mass gap to the LSP """
    if not X1W in protomodel.masses:
        return large_penalty, f"{large_penalty} for no X1W in {protomodel.masses}"
    if not X1W in protomodel.widths:
        return large_penalty, f"{large_penalty} for no X1W in {protomodel.widths}"
    wanted_splitting = mass_splitting ( protomodel.widths[X1W] )
    actual_splitting = protomodel.masses[X1W] - protomodel.masses[LSP]
    discrepancy = abs ( actual_splitting - wanted_splitting )
    p = 2*discrepancy
    return p, f"{p:.3g} for mass-to-width C1 discrepancy {discrepancy:.3g} GeV"

