""" tools, convenience functions that come in handy when dealing with user
defined priors
"""

def mergeReasons ( reasons : list, more_reasons : str|list ):
    """ very simple convenience method to merge reasons """
    if type(more_reasons) == str:
        reasons.append ( more_reasons )
        return
    for r in more_reasons:
        reasons.append ( r )
    
def return_penalty ( p : float, is_nll : bool ) -> float:
    """ very simple convenience to exponentiate penalty if needed """
    import numpy as np
    if is_nll:
        return p
    return np.exp ( -p )
