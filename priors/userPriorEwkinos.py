#!/usr/bin/env python3
    
import numpy as np
from pbase.constants import LSP, X2Z, X1W, X3Z, X2W

""" Example of a user-provided prior: we only like models with an
*off-shell* ewkino.

The walkers can be told to use a user-provided prior, i.e. a python file
that contains a function 'myprior'. That function takes the protomodel and
the nll flag and returns a float:

- for nll=True, the returned value is *added* to the (negative log likelihood
  of the) prior
- for nll=False, the returned value is *multiplied* to the prior

So, for nll=True, the user function should return a penalty (>=0 means the
model is disfavoured), while for nll=False it should return the
corresponding multiplicative factor (i.e. exp(-penalty) for the example below).

The rules implemented here are:

- if the model contains X2Z (ewkino N2) or X1W (ewkino C1) that
  decays *off-shell*, there is no penalty at all
- if the model contains neither X2Z nor X1W, there is a large penalty
- if the model contains them, but they decay via *on-shell* W, Z or h,
  there is a large penalty as well
- (if one of the two ewkinos decays off-shell, but the other one decays via
  an on-shell boson, the model is penalised, i.e. the on-shell decays win)
- (if the ewkino is there, but does not decay at all, the model is
  penalised as well)

Note how the off-shell and on-shell decays are told apart: in this framework
the off-shell decays are the 3-body decays (the intermediate W, Z or h is
virtual, e.g. ewkino -> LSP f f' via a virtual W), while the on-shell
decays are the 2-body ones (e.g. ewkino -> LSP Z with a real Z). The same
can be seen from the masses: the ewkino is in the off-shell regime if the
mass difference to the LSP is smaller than the mass of the intermediate
boson (plus its width).

Use it like this:

    slurm_walk.py --userprior /path/to/userPriorewkinos.py ...

The path is then put into the run.dict (run environment), from where the
walkers pick it up. If you run a walker yourself, you can either use the
run environment:

    environ = RunEnviron.create ( userprior = "/path/to/userPriorewkinos.py" )
    walker = RandomWalker ( environ, rvars )

or set it directly:

    environ = RunEnviron()
    environ.userprior = "/path/to/userPriorewkinos.py"
"""
    

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

def myprior ( protomodel, nll : bool = False ) -> float:
    """ the user prior

    :param protomodel: the model the prior is computed for
    :param nll: if True, return a penalty, else return a factor
    :returns: a float
    """
    reasons = []
    penalty, reason = ewkinoPenalty ( protomodel )
    reasons.append ( reason )
    if False:
        p1, r1 = tpredsBasedPrior ( protomodel )
        penalty += p1
        reasons.append ( reason )
    if True:
        p2, r2 = higgsinoXSecPenalty ( protomodel )
        penalty += p2
        reasons.append ( reason )
    if True:
        p3, r3 = outsideCMS18004Penalty ( protomodel )
        penalty += p3
        reasons.append ( reason )
    if False:
        p4, r4 = degenerateEWKinos ( protomodel )
        penalty += p4
        reasons.append ( reason )
    sreasons = " - " + "\n - ".join ( reasons )
    print ( f"[userPriorEwkinos] penalty={penalty} reason: '{sreasons}'" )
    if nll:
        return penalty ## added to the nll
    return np.exp ( - penalty ) ## multiplied to the prior
