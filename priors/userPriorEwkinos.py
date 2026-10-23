#!/usr/bin/env python3
    
import numpy as np

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

- if the model contains 1000023 (ewkino N2) or 1000024 (ewkino C1) that
  decays *off-shell*, there is no penalty at all
- if the model contains neither 1000023 nor 1000024, there is a large penalty
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
    
LSP = 1000022

def higgsinoXSecPenalty ( protomodel ) -> tuple:
    """ a penalty for ewkino cross sections away from higgsino hypotheses 
    """
    prod_modes = [ (1000022,1000022), (1000023,1000023),
                   (1000023,1000023), (-1000024,1000024),
                   (1000022,1000024), ( 1000023, 1000024 ) ]
    penalty = 0
    ssms = protomodel.ssmultipliers
    hino = .25
    wino = 1.0
    for prod_mode in prod_modes:
        if not prod_mode in ssms:
            penalty += 20
            continue
        deviation = ssms[prod_mode]-wino
        if deviation>0.:
            p = float ( np.log ( 1 + deviation ) )
            penalty += p
        else:
            p = float ( np.log ( (1 - deviation)/4. ) )
    return penalty, f"{penalty} from deviation form higgsino xsec"

def degenerateEWKinos ( protomodel ) -> tuple:
    """ punish if N2 and C1 arent mass degenerate """
    if not 1000023 in protomodel.masses:
        return large_penalty, f"no X2Z in {protomodel.masses}"
    if not 1000024 in protomodel.masses:
        return large_penalty, f"no X2Z in {protomodel.masses}"
    dm = abs(protomodel.masses[1000024]-protomodel.masses[1000023])/\
            (protomodel.masses[1000023]+protomodel.masses[1000024])
    if dm > 1e-5:
        return large_penalty, f"dm(X1W,X2Z) is {dm:.3g}"
    return 0, f"dm(X1W,X2Z) zero"

def outsideCMS18004Penalty ( protomodel ) -> tuple:
    large_penalty = 20.
    if protomodel.masses[LSP] < 70:
        return large_penalty, f"LSP mass at {protomodel.masses[LSP]}"
    if not 1000023 in protomodel.masses:
        return large_penalty, f"no X2Z in {protomodel.masses}"
    if protomodel.masses[1000023]<120:
        return large_penalty, f"X2Z mass at {protomodel.masses[1000023]}"
    if not 1000024 in protomodel.masses:
        return large_penalty, f"no X1W particle in {protomodel.masses}"
    if protomodel.masses[1000024]<120:
        return large_penalty, f"X1W mass at {protomodel.masses[1000024]}"
    return 0, f"within CMS-SUS-18-004"

def ewkinoPenalty ( protomodel ) -> tuple:
    """ compute the penalty and the reason for it

    :param protomodel: the model the prior is computed for
    :returns: (penalty, reason) tuple
    """
    ## the penalty (in nll) that we give to the models we do not like
    large_penalty = 20.

    mass_gaps = { 1000024: 80., 1000023: 90. }

    for ewkino in mass_gaps.keys():
        if not ewkino in protomodel.masses:
            return large_penalty, f"no ewkino {ekwino} in the model"
    dms = {}
    for ewkino, gap in mass_gaps.items():
        dm = protomodel.masses[ewkino]-protomodel.masses[LSP]
        if dm > gap:
            return large_penalty, f"{ewkino} with mass gap {dm:.2f}"
    return 0, f"N2,C1 offshell"

def myprior ( protomodel, nll : bool = False ) -> float:
    """ the user prior

    :param protomodel: the model the prior is computed for
    :param nll: if True, return a penalty, else return a factor
    :returns: a float
    """
    penalty, reason = ewkinoPenalty ( protomodel )
    print ( f"[userPriorEwkinos] penalty {penalty} reason {reason}" )
    if True:
        p2, r2 = higgsinoXSecPenalty ( protomodel )
        penalty += p2
        reason = f"{reason}, {r2}"
    if True:
        p3, r3 = outsideCMS18004Penalty ( protomodel )
        penalty += p3
        reason = f"{reason}, {r3}"
    if False:
        p4, r4 = degenerateEWKinos ( protomodel )
        penalty += p4
        reason = f"{reason}, {r4}"
    if nll:
        return penalty ## added to the nll
    return np.exp ( - penalty ) ## multiplied to the prior
