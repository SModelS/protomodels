#!/usr/bin/env python3

""" Example of a user-configurable prior built from the registered priors.

Run it from the repository root with

    python -m priors.example_configurable_prior

The example discovers all *_priors.py files in this directory, prints the
registered priors and their descriptions, then builds a single user prior
from a configuration dictionary.
"""

from pbase.constants import LSP
from priors.prior_tools import PRIOR_REGISTRY, discoverPriors, \
        makeUserPrior, registerPrior

def lspMassPenalty ( protomodel, min_mass : float = 70. ) -> tuple:
    """ penalise models with an LSP lighter than min_mass (in GeV)

    This shows how a prior with configurable keyword arguments can be added
    to the registry by hand, e.g. when it does not live in a *_priors.py file.

    :param protomodel: the model the prior is computed for
    :param min_mass: the minimum allowed LSP mass in GeV
    :returns: (penalty, reason) tuple
    """
    lsp_mass = protomodel.masses [ LSP ]
    if lsp_mass < min_mass:
        penalty = 20.
        return penalty, f"{penalty} for LSP mass {lsp_mass:.1f} < {min_mass:.1f} GeV"
    return 0., f"0 for LSP mass {lsp_mass:.1f} >= {min_mass:.1f} GeV"

def main ():
    """ discover the priors and build an example user prior """
    discoverPriors ()

    print ( "registered priors:" )
    for name, entry in sorted ( PRIOR_REGISTRY.items () ):
        description = entry [ "description" ].splitlines ()
        first_line = description [ 0 ] if description else "(no description)"
        print ( f"  {name}: {first_line}" )

    ## priors from other files can also be registered by hand
    registerPrior ( lspMassPenalty )

    priorConfig = {
        ## values that are dicts are passed as keyword arguments
        "lspMassPenalty": { "min_mass": 100. },
        ## True means: use this prior without extra arguments
        "higgsinoXSecPenalty": True,
        "outsideCMS18004Penalty": {},
        "degenerateNeutralinos": {},
        "charginoMassWidthPenalty": {},
        ## disabled priors can simply stay in the configuration
        "tpredsBasedPrior": False,
        "degenerateEWKinos": None,
    }

    myprior = makeUserPrior ( priorConfig, verbose = True )

    ## The result can be used exactly like the 'myprior' functions in the
    ## other examples, e.g. by a walker or the Combiner:
    ##
    ##     penalty = myprior ( protomodel, nll = True )
    ##     factor = myprior ( protomodel, nll = False )
    ##
    print ( myprior.__doc__ )

if __name__ == "__main__":
    main ()
