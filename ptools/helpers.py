#!/usr/bin/env python3

""" various helper functions that do not fit in any of the more
specific modules. This module is meant to contain the helper functions
that do not depend on protomodels functionality.
The "moreHelpers" module is meant to contain those that do depend on protomodels
code.
"""

import copy, math, time, random, subprocess, os, unum, numpy
from smodels.experiment.datasetObj import DataSet
from smodels.experiment.expResultObj import ExpResult
from smodels.experiment.infoObj import Info
from smodels.base.physicsUnits import GeV
from smodels.matching.theoryPrediction import TheoryPrediction
import scipy.stats
from os import PathLike
from typing import Union, Set, List
import numpy as np

def repr_double_quotes(obj):
    import json, math
    if isinstance(obj, str):
        # Use json.dumps to get a double-quoted string with escapes handled
        return json.dumps(obj)
    elif isinstance(obj, (list, tuple, set)):
        # Preserve type and recursively format elements
        open_bracket, close_bracket = {
            list: ("[", "]"),
            tuple: ("(", ")"),
            set: ("{", "}")
        }[type(obj)]
        inner = ", ".join(repr_double_quotes(x) for x in obj)
        if isinstance(obj, tuple) and len(obj) == 1:  # special case for (x,)
            inner += ","
        return f"{open_bracket}{inner}{close_bracket}"
    elif isinstance(obj, dict):
        items = (f"{repr_double_quotes(k)}: {repr_double_quotes(v)}" for k, v in obj.items())
        return "{" + ", ".join(items) + "}"
    elif isinstance(obj, float):
        if math.isnan(obj):
            return 'float("nan")'
        elif math.isinf(obj):
            if obj > 0.:
                return 'float("inf")'
            else:
                return 'float("-inf")'
        else:
            return repr(obj)
    else:
        return repr(obj)

def py_dump ( obj, handle, **args ):
    """ equivalent to json.dump, convenience method
    :param handle: is either a file handle, or a file name
    """
    ds = py_dumps ( obj, **args )
    if type(handle) == str: # a pathname was given
        with open ( handle, "wt" ) as f:
            f.write ( ds + "\n" )
            f.close()
    else:
        handle.write ( ds + "\n" )

def reorder_list_by_dict( lst : list, mapping : dict ) -> list:
    """ reorder the list, but make sure the keys in mapping
    come before the values in mapping """
    from collections import defaultdict, deque
    present = set(lst)

    edges = defaultdict(set)
    indegree = defaultdict(int)

    for k, v in mapping.items():
        if k in present and v in present:
            if v not in edges[k]:
                edges[k].add(v)
                indegree[v] += 1
                indegree.setdefault(k, 0)

    queue = deque([x for x in lst if indegree[x] == 0])
    result = []

    while queue:
        node = queue.popleft()
        result.append(node)
        for nxt in edges[node]:
            indegree[nxt] -= 1
            if indegree[nxt] == 0:
                queue.append(nxt)

    if len(result) != len(present):
        raise ValueError("Cycle detected: ordering impossible")

    # preserve original duplicates & unrelated items
    order = {v: i for i, v in enumerate(result)}
    return sorted(lst, key=lambda x: order.get(x, float("inf")))

def py_dumps( obj, indent : int = 4, level : int = 0, stop_at_level : int = -1,
              a_before : dict = {}, double_quotes : bool = True ) -> str:
    """ equivalent to json.dumps (ie it pretty prints a given nested structure)
    but tuples are allowed as keys.

    :param indent: number of spaces used for an indentation
    :param level: how many indentations are we in?
    :param stop_at_level: stop indentation at that level, if positive number
    :param a_before: dictionary, if e.g. { "a": "b", "d": "c" }, then
    :param double_quotes: use double quotes, as required for json.
    FIXME maybe not even ask, its always true
    key a comes before key b

    :returns: formatted string
    """
    sp = ' ' * (level * indent)
    sp_next = ' ' * ((level + 1) * indent)
    mrepr = repr
    if double_quotes:
        mrepr = repr_double_quotes

    if isinstance(obj, dict):
        if not obj:
            return '{}'
        items = []
        keys = list ( obj.keys() )
        if len(a_before)>0:
            keys = reorder_list_by_dict ( keys, a_before )
        if stop_at_level > 0 and level >= stop_at_level:
            for k in keys:
                v = obj[k]
                value = f"{py_dumps(v, indent, level + 1, stop_at_level, a_before, double_quotes )}"
                items.append(f"{mrepr(k)}: {value}")
            return '{ ' + ', '.join(items) + ' }'
        for k in keys:
            v = obj[k]
            value = f"{py_dumps(v, indent, level + 1, stop_at_level, a_before, double_quotes )}"
            items.append(f"{sp_next}{mrepr(k)}: {value}")
        return '{\n' + ',\n'.join(items) + '\n' + sp + '}'

    elif isinstance(obj, list):
        if not obj:
            return '[]'
        items = [f"{sp_next}{py_dumps(i, indent, level + 1, stop_at_level, a_before, double_quotes )}" for i in obj]
        if stop_at_level > 0 and level >= stop_at_level:
            return '[ ' + ', '.join(items) + ' ]'
        return '[\n' + ',\n'.join(items) + '\n' + sp + ']'

    elif isinstance(obj, tuple):
        if not obj:
            return '()'
        items = [f"{sp_next}{py_dumps(i, indent, level + 1, stop_at_level, a_before, double_quotes )}" for i in obj]
        if stop_at_level > 0 and level >= stop_at_level:
            return '( ' + ', '.join(items) + ' )'
        return '(\n' + ',\n'.join(items) + '\n' + sp + ')'

    return mrepr(obj)

def mkdir ( dirname : os.PathLike ) -> bool:
    """ make a directory, gracefully """
    if dirname == "":
        return False
    dirname = os.path.expanduser ( dirname )
    if os.path.exists ( dirname ):
        return False
    try:
        os.mkdir ( dirname )
        return True
    except FileExistsError as e:
        pass
    return False

def formatObject ( obj, fmt_str : Union[int,str] = ".2f" ) -> str:
    """ format an object like a number, e.g. a test statistic

    :param fmt_str: either e.g. .2f, or '1' which gets translated to .1f
    """
    if obj == None:
        return "None"
    if type(fmt_str) == int:
        fmt_str = f".{fmt_str}f"
    return f"{obj:{fmt_str}}"

def getJsonFileName(dset: DataSet) -> str:
    "get file name of json used by the combined dataset dset"

    jsonFileDict = dset.globalInfo.jsonFiles
    dsId = [ds.getID() for ds in dset._datasets]            #get the dataset ids in the combined dataset dset

    for file, dslist in jsonFileDict.items():
        for ds in dslist:
            if ds['smodels'] in dsId:                               #check which json file has the corresponding datasets
                file = file.split(".")[0]                       #get only name of json file, not the .json part
                return file

    # if no file got matched with dataset
    print(f"JSON file present for {dset.globalInfo.id} but combined dataset does not match to any JSON file")

    return "NoJsonFound"

def experimentalId(pred : TheoryPrediction) -> str:
    """
    Return Id of tpred's expresult
        - anaId:upperLimit      if dataType = upperLimit
        - anaId:dataId          if dataType = efficiencyMap
        - anaId:combined        if dataType = combined,SLv1,v2
        - anaId:jsonFileName    if dataType = combined,pyhf
    """

    anaId = pred.analysisId()
    dtype = pred.dataType()

    if dtype == "upperLimit":
        return f"{anaId}:{dtype}"

    elif dtype == "combined":
        if pred._statsComputer.dataType == "pyhf":
            jfile = getJsonFileName(pred.dataset)
            return f"{anaId}:{jfile}"
        else:
            return f"{anaId}:{dtype}"                       #SLv1,v2
    else:
        dsId = pred.dataId()                                #for em-type results
        return f"{anaId}:{dsId}"

def getAllProdModesOfTheoryPred (pred : TheoryPrediction) -> Set:
    """ Get all Prod Modes assoaciated with a theory pred i.e (PV->PID1,PID2..)"""
    prod_modes = set()
    smses = pred.smsList
    for sms in smses:
        for momIndex,dIndex in sms.genIndexIterator():
            pdg1, pdg2 = [],[]
            mass_comp = False
            if momIndex == sms.rootIndex:  # Get Primary Vertex
                daughter = sms.indexToNode(dIndex)
                dppdg = [d.pdg for d in daughter]
                for item in dppdg:
                    if isinstance(item, list): pdg2 += item; mass_comp = True
                    else: pdg1.append(item)
                if mass_comp:                               #Ex: For cases such as: [1000023, [1000022, -1000022]]
                    import itertools
                    pprod = itertools.product(pdg1, pdg2)
                    for pids in pprod:
                        if -1000022 in pids: continue       #FIXME! Hack for now, only checking for -LSP
                        dppdg = tuple(sorted(pids))
                        prod_modes.add(dppdg)
                    continue
                dppdg = tuple(sorted(dppdg))
                prod_modes.add(dppdg)
                continue
    return prod_modes

def getAllPidsOfTheoryPred ( pred : TheoryPrediction ) -> List:
    """ get all pids that make it into a theory prediction """
    def addPDGs ( pids, pid ):
        if type(pid) == int:
            pids.add ( abs(pid) )
        if type(pid) in [ tuple, list ]:
            for p in pid:
                pids.add ( abs(p) )
    pids = set()
    smses = pred.smsList
    for sms in smses:
        for dIndex in sms.daughterIndices(sms.rootIndex):
            daughter = sms.indexToNode(dIndex)
            addPDGs ( pids, daughter.pdg )
            for nodeIndex in sms.dfsIndexIterator(dIndex):
                node = sms.indexToNode(nodeIndex)
                if node.isSM:
                    continue
                addPDGs ( pids, node.pdg )
    pids = list(pids)
    pids.sort( reverse=True )
    return pids

def prettyPrint ( value : Union[None,float,numpy.float64],
        ndecimals : int = 2, maxrows : int = 0 ) -> str:
    """ pretty print a value, but allow for it to also be None

    :param maxrows: maximum number of rows for lists and tuples. zero is all.
    """
    if type(value) in [ float, numpy.float64 ]:
        return f"{value:.{ndecimals}f}"
    if type(value) in [ list, tuple ]:
        if maxrows == 0:
            maxrows = len(value)
        ret = ', '.join ( map ( prettyPrint, value[:maxrows], [ndecimals]*len(value[:maxrows]) ) )
        return ret
    if type(value) == dict:
        ret = ', '.join(f'{k}: {prettyPrint(v,ndecimals)}' for k,v in value.items())
        return f"{{ {ret} }}"
    return str(value)

def nround ( value : Union[None,float], ndecimals : int ) -> Union[None,float]:
    if type(value) == type(None):
        return value
    return round(value,ndecimals)

def simplifyUnixPath ( path : str ) -> str:
    """ simple code to simplify file paths in printouts """
    path = path.replace( f"{os.getcwd()}/", "./" )
    if path.startswith ( os.environ["HOME" ] ):
        path = path.replace( os.environ["HOME"], "~" )
    while "//" in path:
        path = path.replace("//","/")
    return path

def computeZFromP ( pvalue : float ) -> float:
    """ compute significance Z from p-value, i.e. compute Phi^-1 ( p )

    :param pvalue: the p-value
    :returns: the corresponding significance Z
    """
    return float ( - scipy.stats.norm.ppf ( pvalue ) )

def computePForDataSet ( dataset : DataSet, obsN : Union[int,None] = None,
       nmax : int = 100000000, nmin : int = 200000 ) -> float:
    """ given a dataset, compute p for SM hypothesis
    :param obsN: if not None, compute for the observation
    :param nmax: maximum number of toys

    :returns: p-value
    """
    exp = dataset.dataInfo.expectedBG
    err = dataset.dataInfo.bgError
    if obsN is None:
        obsN = dataset.dataInfo.observedN
    thirdMoment = None
    if hasattr ( dataset.dataInfo, "thirdMoment" ):
        thirdMoment = dataset.dataInfo.thirdMoment
    if thirdMoment is None:
        p = computeP ( obsN, exp, err, nmax = nmax, nmin = nmin )
    else:
        p = computePSLv2 ( obsN, exp, err, thirdMoment, nmax = nmax, nmin = nmin )
    if p < 1e-100:
        print ( f"[helpers] {dataset.globalInfo.id}:{dataset.dataInfo.id} has p={p}" )
    return p

def computeP ( obs : float, bg : float, bgerr : float,
        lognormal : bool = False, nmax : int = 100000000,
        sigN : Union[None,float] = None, nmin : int = 50000 ) -> float:
    """ compute P value, gaussian or log-normal nuisance model, w.r.t
    SM hypothesis

    :param obs: observed number of events
    :param bg: number of expected background events
    :param bgerr: error on number of expected bg events
    :param lognormal: if true, model the enveloping nuisance parameter
    as a lognormal instead of a normal
    :param nmax: maximum number of toys
    :param nmin: minimum number of toys

    :returns: p-value
    """
    n = min ( nmin, nmax )
    ret = 0.
    while ret < .9 / nmax or ret > 1. - .9 / nmax:
        if n > nmax:
            # print ( f"[helpers] when computing p: n={n}>{nmax}. breaking off with ret={ret} obs={obs} bg={bg} bgerr={bgerr}" )
            break
        lmbda = scipy.stats.norm.rvs ( loc=[bg]*n, scale=[bgerr]*n )
        lmbda = lmbda[lmbda>0.]
        if lognormal:
            # for lognormal and signals
            central = bg
            if sigN != None:
                central = bg + sigN
            if lognormal and central > ( bgerr / 4. ):
                loc = central**2 / np.sqrt ( central**2 + bgerr**2 )
                stderr = float ( np.sqrt ( np.log ( 1 + bgerr**2 / central**2 ) ) )
                if stderr == 0.:
                    return 0.
                lmbda = scipy.stats.lognorm.rvs ( s=[stderr]*n, scale=[loc]*n )
        fakeobs = scipy.stats.poisson.rvs ( lmbda )
        ## == we count half
        ret = float ( ( sum(fakeobs>obs) + .5*sum(fakeobs==obs) ) / len(fakeobs) )
        n *= 5
    return ret

def computePSLv2 ( obs : float, bg : float, bgerr : float,
        third : float, nmax : int = 100000000,
	      nmin : int = 200000 ) -> float:
    """ compute p value, gaussian nuisance model, w.r.t SM hypothesis, for SLv2

    :param obs: observed number of events
    :param bg: number of expected background events
    :param bgerr: error on number of expected bg events
    :param third: the third moment
    :param nmax: maximum number of toys

    :returns: p-value
    """
    # return -1
    from smodels.statistics.simplifiedLikelihoods import Data
    printErr = True
    while 8*bgerr**6 - third**2 < 0.:
        if printErr:
            ## FIXME ugly hack, shrink the third momenta
            print ( f"[helpers] third moments too large (bgerr**2={bgerr**2:.3g}, third={third:.3g}). shrink them!" )
            printErr = False
        third *= 0.9
    d = Data ( obs, bg, bgerr**2, third )
    from icecream import ic
    #ic ( "FIXME needs implementation! computePSLv2" )
    n = nmin
    ret = 0.
    rhoparam = d.rho[0][0]
    # thtadbn = scipy.stats.multivariate_normal(np.zeros(self.size), rhoparam )
    while ret < .9/nmax or ret > 1. - .9/nmax:
        if n > nmax:
            print ( f"[helpers] SLv2 n={n}>{nmax}. breaking off with ret={ret} obs={obs} bg={bg} bgerr={bgerr} third={third}" )
            break
        ctr = 0
        # thtas = thtadbn.rvs( n )
        thtas = scipy.stats.norm.rvs ( loc=[0.]*n, scale=[1.]*n )
        lmbdas = d.A + d.B * thtas + d.C * thtas**2
        indices = numpy.where ( lmbdas < 0. )[0]
        while len(indices)>0:
            thta = scipy.stats.norm.rvs( loc=[0.]*len(indices), scale=[1.]*len(indices) )
            lmbdas [ indices ] = thta
            indices = numpy.where ( lmbdas < 0. )[0]
            ctr += 1
            if ctr > 20: # after trying 20 times we set to almost zero
                lmbdas [ indices ] = [0.]*len(indices)
                break
        try:
            fakeobs = scipy.stats.poisson.rvs ( lmbdas )
        except ValueError as e:
            ic ( lmbdas )
            ic ( d.A )
            ic ( d.B )
            ic ( d.C )
            ic ( d.rho )
            ic ( obs, bg, bgerr, third )
            import sys; sys.exit()
        ## == we count half
        ret = float ( ( sum(fakeobs>obs) + .5*sum(fakeobs==obs) ) / len(fakeobs) )
        n *= 5
    return ret


def stripUnits( container ):
    """ strip all units from a mass vector """
    if type(container) in [ None ]:
        return container
    ret = []
    for br in container:
        tbr = []
        for m in br:
            if type(m) in [ float, int ]:
                tbr.append ( m )
            if type(m) == type(GeV):
                tbr.append ( m.asNumber(GeV) )
            if type(m) in [ tuple, list ]:
                tbr.append ( m )
        ret.append ( tbr)
    return ret

def countSSMultipliers ():
    """ count the total number of ssmultipliers of a protomodel """
    from builder.protomodel import ProtoModel
    model = ProtoModel()
    modes = set()
    def sortMe ( p, q ):
        if p < q:
            return (p,q)
        return ( q,p )
    for p in model.particles:
        for q in model.particles:
            modes.add ( (p,q) )
            if model.hasAntiParticle(p):
                modes.add ( sortMe(-p,q ) )
            if model.hasAntiParticle(q):
                modes.add( sortMe(p,-q) )
            if model.hasAntiParticle(p) and model.hasAntiParticle(q):
                modes.add ( sortMe(-p,-q) )
    print ( f"We have {len(modes)} production modes" )
    return modes

def countDecays( templatefile = "../builder/templates/template_default.slha" ):
    """ count the number of decays in a template file """
    if not os.path.exists ( templatefile ):
        templatefile = templatefile.replace("../","./" )
        if not os.path.exists ( templatefile ):
            print ( f"Could not find template file {templatefile}" )
            return 0
    with open( templatefile ) as f:
       lines=f.readlines()
    count = []
    for line in lines:
        if "#" in line:
            p = line.find("#")
            line = line[:p]
        line = line.strip()
        if not "D" in line:
            continue
        if "DECAY" in line:
            continue
        if "BLOCK" in line:
            continue
        if line == "":
            continue
        line = line.replace("D","")
        tokens = line.split(" ")
        ids = tokens[0].split("_")
        ids = tuple ( map ( int , ids ) )
        count.append ( ids )
    print (f"I count {len(count)} decay channels" )
    return count

def seedRandomNumbers ( seed ):
    """ seed all random number generation """
    ## scipy takes random numbers from numpy.random, so
    np.random.seed ( seed )
    import scipy.stats as s
    r = s.norm.rvs()
    print(f"[helpers] seeding the random number generators with {seed}. Here is a first realization of a standard normal: {r:.3f}")

def cpPythia8 ( ):
    """ as a very ugly workaround for now, if something goes wrong with
        cross sections, cp the pythia8 install. """
    libdir = f"{os.environ['HOME']}/git/smodels/smodels/lib"
    if os.path.exists ( f"{libdir}/pythia8/pythia8226/share/Pythia8/xmldoc/Welcome.xml" ):
        return
    lockfile = f"{libdir}/lock"
    ctr = 0
    while os.path.exists ( lockfile ):
        time.sleep ( np.random.uniform ( 1, 3 ) )
        ctr += 1
        if ctr > 5:
            break
    cmd = f"touch {lockfile}"
    o = subprocess.getoutput ( cmd )
    cmd = f"chmod -R u+w {libdir}/pythia8old {libdir}/pythia8"
    o = subprocess.getoutput ( cmd )
    cmd = f"rm -rf {libdir}/pythia8old"
    o = subprocess.getoutput ( cmd )
    cmd = f"mv {libdir}/pythia8 {libdir}/pythia8old"
    o = subprocess.getoutput ( cmd )
    cmd = f"cp -r {libdir}/pythia8backup {libdir}/pythia8"
    o = subprocess.getoutput ( cmd )
    cmd = f"rm -f {lockfile} lockfile"
    o = subprocess.getoutput ( cmd )

def lrEquiv ( l, r ):
    """ check if the two strings are equivalent up to L vs R """
    if type(l) != str:
        return False
    if type(r) != str:
        return False
    if l.startswith("+-") and r.startswith("+-"):
        l=l.replace("+-","")
        r=r.replace("+-","")
    return l[1:] == r[1:]

def simplifyList ( modes ):
    """ simplify a given list of production modes """
    import itertools
    ret = copy.deepcopy ( modes )
    for combo in itertools.combinations ( modes, 2 ):
        if combo[0][0] == combo[1][0] and combo[0][1] == -combo[1][1]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( combo[0][0], f"+-{abs(combo[1][1])}" ) )
            except ValueError:
                pass
        if combo[0][0] == combo[1][1] and combo[0][1] == -combo[1][0]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( combo[0][0], f"+-{abs(combo[1][0])}" ) )
            except ValueError:
                pass
        if combo[0][1] == combo[1][1] and combo[0][0] == -combo[1][0]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( combo[0][1], f"+-{abs(combo[1][0])}" ) )
            except ValueError:
                pass
        if combo[0][0] == -combo[1][0] and combo[0][1] == -combo[1][1]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( f"+-{abs(combo[0][0])}", f"+-{abs(combo[1][0])}" ) )
            except ValueError:
                pass
        if combo[0][0] == -combo[1][1] and combo[0][1] == -combo[1][0]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( f"+-{abs(combo[0][0])}", f"+-{abs(combo[1][0])}" ) )
            except ValueError:
                pass
    modes = copy.deepcopy ( ret )
    for combo in itertools.combinations ( modes, 2 ):
        if type(combo[0][1])==str and type(combo[1][1])==str:
            if combo[0][1] == combo[1][1] and type(combo[0][0])==int and \
                            type(combo[1][0])==int and combo[0][0]==-combo[1][0]:
                try:
                    ret.remove ( combo[0] )
                    ret.remove ( combo[1] )
                    ret.append ( ( f"+-{abs(combo[0][0])}", combo[0][1] ) )
                except ValueError as e:
                    pass
        if type(combo[0][0])==int and type(combo[1][1])==str:
            c00 = abs(combo[0][0])
            if combo[0] == (-c00, -c00) and combo[1] == (c00, f'+-{c00}' ):
                try:
                    ret.remove ( combo[0] )
                    ret.remove ( combo[1] )
                    ret.append ( ( f"+-{abs(combo[0][0])}", combo[1][1] ) )
                except ValueError as e:
                    pass
            if type(combo[0][1])==str and c00 == combo[1][0] and \
                     lrEquiv ( combo[0][1], combo[1][1] ):
                ## (1000021, '+-2000006'), (1000021, '+-1000006')
                try:
                    ret.remove ( combo[0] )
                    ret.remove ( combo[1] )
                    c11 = combo[1][1].replace("+-1","+-?").replace("+-2","+-?")
                    ret.append ( ( combo[0][0], c11 ) )
                except ValueError as e:
                    pass
    # print ( "reduced to", ret )
    return ret

def lightObjCopy(obj,rmAttr=['elements','avgElement', 'computer', 'txnameList',
                          'txnames','datasets','_databaseParticles',
                          'comment','path','url','publication','contact']):

    """Tries to make a light copy of an object. The attributes in rmAttr will not be copied"""

    if obj is None:
        return obj
    elif isinstance(obj,(int,float,unum.Unum,str,numpy.float,numpy.bool_)):
        return obj
    elif isinstance(obj,list):
        return [lightObjCopy(x,rmAttr=rmAttr) for x in obj]
    elif isinstance(obj,tuple):
        return tuple([lightObjCopy(x,rmAttr=rmAttr) for x in obj])
    elif isinstance(obj,dict):
        return dict([[lightObjCopy(k,rmAttr=rmAttr),lightObjCopy(v,rmAttr=rmAttr)] for k,v in obj.items()])
    elif isinstance(obj,DataSet):
        newDS = DataSet()
        newDS.dataInfo = Info()
        for key,v in obj.dataInfo.__dict__.items():
            if key in rmAttr: continue
            setattr(newDS.dataInfo,key,v)
        return newDS
    elif isinstance(obj,ExpResult):
        newExp = ExpResult()
        newExp.globalInfo = Info()
        newExp.datasets = []
        for key,v in obj.globalInfo.__dict__.items():
            if key in rmAttr: continue
            setattr(newExp.globalInfo,key,v)
        return newExp
    else:
        newObj = obj.__class__()
        for key,val in obj.__dict__.items():
            if key in rmAttr: continue
            setattr(newObj,lightObjCopy(key,rmAttr=rmAttr),lightObjCopy(val,rmAttr=rmAttr))
        return newObj
