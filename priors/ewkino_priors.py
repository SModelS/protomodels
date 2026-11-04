#!/usr/bin/env python3
    
import numpy as np
from pbase.constants import LSP, X2Z, X1W, X3Z, X2W

""" A collection of possibly penalty terms for ewkino type scenarios
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
            reasons.append ( f"{p} for {prod_mode} deviating from higgsino" )
        else:
            p = float ( abs ( np.log ( (1 - deviation)/4. ) ) )
            penalty += p
            reasons.append ( f"{p} for {prod_mode} deviating from higgsino" )
    return penalty, ", ".join ( reasons )

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

def ewkinoPenalty ( protomodel ) -> tuple:
    """ compute the penalty and the reason for it

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
            return large_penalty, f"{ewkino} with mass gap {dm:.2f}"
    return 0, f"0 for N2,C1 offshell"

def tpredsBasedPrior ( protomodel ) -> tuple[float,str]:
    """ a prior that penalizes based on the theorypredictions """
    from tester.predictor import Predictor
    pr = Predictor ( "myprior", protomodel.environ )
    pr.obtainPredictions ( protomodel, keep_predictions = True )
    n_predictions = len ( pr.predictions )
    return 0, "not yet in use"
