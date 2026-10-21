#!/usr/bin/env python3

""" Example of a user-provided prior: we only like models with an
*off-shell* chargino.

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

- if the model contains 1000023 (chargino N2) or 1000024 (chargino C1) that
  decays *off-shell*, there is no penalty at all
- if the model contains neither 1000023 nor 1000024, there is a large penalty
- if the model contains them, but they decay via *on-shell* W, Z or h,
  there is a large penalty as well
- (if one of the two charginos decays off-shell, but the other one decays via
  an on-shell boson, the model is penalised, i.e. the on-shell decays win)
- (if the chargino is there, but does not decay at all, the model is
  penalised as well)

Note how the off-shell and on-shell decays are told apart: in this framework
the off-shell decays are the 3-body decays (the intermediate W, Z or h is
virtual, e.g. chargino -> LSP f f' via a virtual W), while the on-shell
decays are the 2-body ones (e.g. chargino -> LSP Z with a real Z). The same
can be seen from the masses: the chargino is in the off-shell regime if the
mass difference to the LSP is smaller than the mass of the intermediate
boson (plus its width).

Use it like this:

    slurm_walk.py --userprior /path/to/userPriorCharginos.py ...

The path is then put into the run.dict (run environment), from where the
walkers pick it up. If you run a walker yourself, you can either use the
run environment:

    environ = RunEnviron.create ( userprior = "/path/to/userPriorCharginos.py" )
    walker = RandomWalker ( environ, rvars )

or set it directly:

    environ = RunEnviron()
    environ.userprior = "/path/to/userPriorCharginos.py"
"""

import numpy

## the particles we are interested in: chargino N2 and chargino C1
CHARGINOS = [ 1000023, 1000024 ]
## the SM bosons we do not want to see on-shell
BOSONS = { 24: "W", 23: "Z", 25: "h" }
## the intermediate boson of the off-shell decays of each chargino
INTERMEDIATE = { 1000023: "Z", 1000024: "W" }
## the penalty (in nll) that we give to the models we do not like
LARGE_PENALTY = 10.

try:
    from pbase.constants import smMasses, smWidths
except ImportError: ## if the protomodels package is not importable, we
    smMasses, smWidths = {}, {} ## just cannot do the mass-based cross-check


def isOffshellDecay ( decay ) -> bool:
    """ check if a decay is an off-shell decay

    In this framework the off-shell decays are the 3-body decays, where the
    intermediate boson (W, Z or h) is virtual (e.g. chargino -> LSP f f'),
    while the on-shell decays are the 2-body ones (e.g. chargino -> LSP Z).

    :param decay: the decay (a tuple of particle ids)
    :returns: True if the decay is an off-shell decay
    """
    return len ( decay ) == 3

def onShellBoson ( decay ):
    """ check if a decay goes through an on-shell W, Z or h

    :param decay: the decay (a tuple of particle ids)
    :returns: the name of the boson ('W', 'Z' or 'h') or None
    """
    for dpid in decay:
        if abs ( dpid ) in BOSONS:
            return BOSONS [ abs ( dpid ) ]
    return None

def offshellRegime ( protomodel, pid ) -> bool:
    """ check if the chargino pid is in the off-shell regime, i.e. if it is
    lighter than the intermediate boson (plus its width) above the LSP

    :param protomodel: the model
    :param pid: the pid of the chargino
    :returns: True if the chargino is in the off-shell regime
    """
    lsp = getattr ( protomodel, "LSP", None )
    if lsp is None or lsp not in protomodel.masses or pid not in protomodel.masses:
        return False
    boson = INTERMEDIATE.get ( pid, None )
    if boson is None or boson not in smMasses:
        return False
    dm = protomodel.masses[pid] - protomodel.masses[lsp]
    return dm < ( smMasses[boson] + smWidths.get ( boson, 0. ) )

def charginoPenalty ( protomodel ):
    """ compute the penalty and the reason for it

    :param protomodel: the model the prior is computed for
    :returns: (penalty, reason) tuple
    """
    ## the charginos that are actually in the model (frozen particles have
    ## been removed from the masses)
    charginos = [ pid for pid in CHARGINOS if pid in protomodel.masses and
                  protomodel.masses[pid] > 0 ]

    if len ( charginos ) == 0:
        return LARGE_PENALTY, "no chargino (1000023 or 1000024) in the model"

    offshell = [] ## decays via a virtual W/Z/h
    onshell = []  ## decays via a real W/Z/h
    for pid in charginos:
        for decay, br in protomodel.decays.get ( pid, {} ).items():
            if br is None or br <= 0.:
                continue ## ignore closed channels
            if isOffshellDecay ( decay ):
                offshell.append ( ( pid, decay ) )
            else:
                boson = onShellBoson ( decay )
                if boson != None:
                    onshell.append ( ( pid, decay, boson ) )

    if len ( onshell ) > 0:
        ## the chargino decays via a real W, Z or h -> we do not like it
        decays = ", ".join ( f"{pid} -> {' '.join(map(str,decay))} ({boson})" for
                            pid, decay, boson in onshell )
        return LARGE_PENALTY, f"on-shell boson decays: {decays}"

    if len ( offshell ) == 0:
        ## the chargino is there, but it does not decay off-shell
        regime = ", ".join ( f"{pid} ({'off-shell' if offshellRegime(protomodel,pid) else 'on-shell'})" for
                            pid in charginos )
        return LARGE_PENALTY, f"the charginos do not decay off-shell ({regime})"

    decays = ", ".join ( f"{pid} -> {' '.join(map(str,decay))}" for pid, decay in offshell )
    return 0., f"off-shell chargino decays: {decays}"

def myprior ( protomodel, nll : bool = False ) -> float:
    """ the user prior

    :param protomodel: the model the prior is computed for
    :param nll: if True, return a penalty, else return a factor
    :returns: a float
    """
    penalty, reason = charginoPenalty ( protomodel )
    print ( f"@@userPriorCharginos penalty {penalty} reason {reason}" )
    if nll:
        return penalty ## added to the nll
    return numpy.exp ( - penalty ) ## multiplied to the prior
