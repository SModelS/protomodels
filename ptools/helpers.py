#!/usr/bin/env python3

""" various helper functions that do not fit in any of the more
    specific modules """

import copy, math, time, random, subprocess, os, unum, numpy
from smodels.experiment.datasetObj import DataSet
from smodels.experiment.expResultObj import ExpResult
from smodels.experiment.infoObj import Info
from smodels.tools.physicsUnits import GeV
import scipy.stats
from os import PathLike

def readDictionaryFile ( filename : PathLike ) -> dict:
    """ read the database dictionary files, as produced by expResModifier.py -C
    :param filename: path to the dictionary file
    :returns: a dictionary with "meta" and "data" as keys
    
    .. code-block:: python3

    >>> content = readDictionaryFile ( "230.dict" )
    """
    with open( filename,"rt") as f:
        tmp=f.readlines()
    lines = []
    for line in tmp:
        if line.startswith("#"):
            continue
        lines.append ( line )
    basename = os.path.basename ( filename ).replace(".dict","")
    ret = { "meta": {}, "data": {}, "basename": basename }
    ret["meta"].update (  eval(lines[0]) )
    nan=float("nan")
    data = eval("\n".join(lines[1:]))
    ret["data"] = data
    return ret

def computeP ( obs, bg, bgerr, lognormal = False ):
    """ compute P value, gaussian nuisance model only
    :param obs: observed number of events
    :param bg: number of expected background events
    :param bgerr: error on number of expected bg events
    :returns: p-value
    """
    n = 50000
    ret = 0.
    while ret < 1e-22 or ret > 1. - 1e-22:
        lmbda = scipy.stats.norm.rvs ( loc=[bg]*n, scale=[bgerr]*n )
        lmbda = lmbda[lmbda>0.]
        if lognormal:
            # for lognormal and signals
            central = bg
            if self.signalmodel and sigN != None:
                central = bg + sigN
            if lognormal and central > ( bgerr / 4. ):
                loc = central**2 / np.sqrt ( central**2 + bgerr**2 )
                stderr = np.sqrt ( np.log ( 1 + bgerr**2 / central**2 ) )
                lmbda = scipy.stats.lognorm.rvs ( s=[stderr]*n, scale=[loc]*n )
        fakeobs = scipy.stats.poisson.rvs ( lmbda )
        ## == we count half
        ret = ( sum(fakeobs>obs) + .5*sum(fakeobs==obs) ) / len(fakeobs)
        n *= 5
        if n > 4000000:
            break
    return ret

def stripUnits( container ):
    """ strip all units from a mass vector """
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
    print ( "We have %d production modes" % len(modes) )
    return modes

def countDecays( templatefile = "../builder/templates/template1g.slha" ):
    """ count the number of decays in a template file """
    if not os.path.exists ( templatefile ):
        templatefile = templatefile.replace("../","./" )
        if not os.path.exists ( templatefile ):
            print ( "Could not find template file %s" % templatefile )
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
    print ("I count %d decay channels" % len(count) )
    return count

def seedRandomNumbers ( seed ):
    """ seed all random number generation """
    import random, numpy
    random.seed ( seed )
    ## scipy takes random numbers from numpy.random, so
    numpy.random.seed ( seed )
    import scipy.stats as s
    r = s.norm.rvs()
    print ( "[helpers] seeding the random number generators with %d. r=%.3f" % \
            ( seed, r ) )

def cpPythia8 ( ):
    """ as a very ugly workaround for now, if something goes wrong with
        cross sections, cp the pythia8 install. """
    libdir = "/users/wolfgan.waltenberger/git/smodels/smodels/lib"
    # ~/git/smodels/smodels/lib/pythia8/pythia8226/share/Pythia8/xmldoc/Welcome.xml
    if os.path.exists ( libdir + "/pythia8/pythia8226/share/Pythia8/xmldoc/Welcome.xml" ):
        return
    lockfile = libdir+"/lock"
    ctr = 0
    while os.path.exists ( lockfile ):
        time.sleep ( random.uniform ( 1, 3 ) )
        ctr += 1
        if ctr > 5:
            break
    cmd = "touch %s" % lockfile
    subprocess.getoutput ( cmd )
    cmd = "rm -rf %s/pythia8old" % libdir
    subprocess.getoutput ( cmd )
    cmd = "mv %s/pythia8 %s/pythia8old" % ( libdir, libdir )
    subprocess.getoutput ( cmd )
    cmd = "cp -r %s/pythia8backup %s/pythia8" % ( libdir, libdir )
    subprocess.getoutput ( cmd )
    cmd = "rm -f %s" % lockfile
    subprocess.getoutput ( cmd )

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
    # print ( "reducing", modes )
    import itertools
    ret = copy.deepcopy ( modes )
    for combo in itertools.combinations ( modes, 2 ):
        if combo[0][0] == combo[1][0] and combo[0][1] == -combo[1][1]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( combo[0][0], "+-%s" % abs(combo[1][1]) ) )
            except ValueError:
                pass
        if combo[0][0] == combo[1][1] and combo[0][1] == -combo[1][0]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( combo[0][0], "+-%s" % abs(combo[1][0]) ) )
            except ValueError:
                pass
        if combo[0][1] == combo[1][1] and combo[0][0] == -combo[1][0]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( combo[0][1], "+-%s" % abs(combo[1][0]) ) )
            except ValueError:
                pass
        if combo[0][0] == -combo[1][0] and combo[0][1] == -combo[1][1]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( "+-%s" % abs(combo[0][0]), "+-%s" % abs(combo[1][0]) ) )
            except ValueError:
                pass
        if combo[0][0] == -combo[1][1] and combo[0][1] == -combo[1][0]:
            try:
                ret.remove ( combo[0] )
                ret.remove ( combo[1] )
                ret.append ( ( "+-%s" % abs(combo[0][0]), "+-%s" % abs(combo[1][0]) ) )
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
                    ret.append ( ( "+-%s" % abs(combo[0][0]), combo[0][1] ) )
                except ValueError as e:
                    pass
        if type(combo[0][0])==int and type(combo[1][1])==str:
            c00 = abs(combo[0][0])
            if combo[0] == (-c00, -c00) and combo[1] == (c00, '+-%s' % c00 ):
                try:
                    ret.remove ( combo[0] )
                    ret.remove ( combo[1] )
                    ret.append ( ( "+-%s" % abs(combo[0][0]), combo[1][1] ) )
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
