#!/usr/bin/env python3
    
""" an example of a user based prior, compased of pre-defined
    penalty terms"""

import ewkinoPenalty, higgsinoXSecPenalty, outsideCMS18004Penalty, \
       degenerate, tpredsBasedPrior

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
