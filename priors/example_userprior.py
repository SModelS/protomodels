#!/usr/bin/env python3
    
from priors.ewkino_priors import ewkinoPenalty, tpredsBasedPrior, \
         higgsinoXSecPenalty, outsideCMS18004Penalty, \
         degenerateEWKinos, degenerateNeutralinos, \
         charginoMassWidthPenalty, mergeReasons

def myprior ( protomodel, nll : bool = False ) -> float:
    """ the user prior

    :param protomodel: the model the prior is computed for
    :param nll: if True, return a penalty, else return a factor
    :returns: a float
    """
    reasons = []
    penalty, reason = ewkinoPenalty ( protomodel )
    mergeReasons ( reasons, reason )
    if False:
        p1, r1 = tpredsBasedPrior ( protomodel )
        penalty += p1
        mergeReasons ( reasons, r1 )
    if True:
        p2, r2 = higgsinoXSecPenalty ( protomodel )
        penalty += p2
        mergeReasons ( reasons, r2 )
    if True:
        p3, r3 = outsideCMS18004Penalty ( protomodel )
        penalty += p3
        mergeReasons ( reasons, r3 )
    if False:
        p4, r4 = degenerateEWKinos ( protomodel )
        penalty += p4
        mergeReasons ( reasons, r4 )
    if True:
        p5, r5 = degenerateNeutralinos ( protomodel )
        penalty += p5
        mergeReasons ( reasons, r5 )
    if True:
        p6, r6 = charginoMassWidthPenalty ( protomodel )
        penalty += p6
        mergeReasons ( reasons, r6 )
    sreasons = " - " + "\n - ".join ( reasons )
    print ( f"[userPriorEwkinos] penalty={penalty:.3g} reason:" )
    print ( f"{sreasons}" )
    if nll:
        return penalty ## added to the nll
    import numpy as np
    return np.exp ( - penalty ) ## multiplied to the prior
