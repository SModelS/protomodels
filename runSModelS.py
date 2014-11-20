#!/usr/bin/env python

"""
.. module:: Main code for running smodels
"""

import sys, os, logging
import argparse
from ConfigParser import SafeConfigParser
from smodels.tools.physicsUnits import GeV, fb, TeV
from smodels.tools import slhaChecks, ioObjects, xsecComputer, missingTopologies
from smodels.experiment import smsHelpers, smsAnalysisFactory, smsResults
from smodels.theory import slhaDecomposer, lheDecomposer
from smodels.theory.theoryPrediction import theoryPredictionFor,  _getElementsFrom
from smodels.theory import crossSection
from smodels.installation import installDirectory


def main(filename, parameterfile=None, outputfile="summary.txt"):

    log = logging.getLogger(__name__)

    inputFile = filename #get input filename

    #use ConfigParser to read input parameters
    parser = SafeConfigParser()
    #first, read defaults
    parser.read("%s/etc/parameters_default.ini" % installDirectory())
    #then, set parameters given in input file
    if parameterfile: parser.read(parameterfile)

    if os.path.exists(outputfile): #remove old output file
        log.warning("Removing old output file in "+outputfile)
        os.remove(outputfile)

    #get status, warnings from slhaChecks

    #slhastat == 0, the file is not checked
    #slhastat == 1, the check is ok
    #slhastat == -1, physical problem in scenario, e.g. charged LSP,
    #slhastat == -2  formal problems in the file, e.g. missing decay blocks
    #slhastat == -3 missing input file

    if parser.getboolean("options","doSLHAdec"):
        inputStatus = slhaChecks.SlhaStatus(inputFile, sigmacut = parser.getfloat("parameters","sigmacut"), checkXsec = not parser.getboolean("options","addMissingXsecs"), massgap = parser.getfloat("parameters","minmassgap"), maxcond = parser.getfloat("parameters","maxcond"))
    else:
        #for input lhe, check if file exists and has correct format
        inputStatus = ioObjects.LheStatus(inputFile, massgap=parser.getfloat("parameters","minmassgap"), maxcond=parser.getfloat("parameters","maxcond"))

    inStat, warnings = inputStatus.status

    #set database address
    smsHelpers.base = parser.get("path","database")

    for i in [ "$HOME", "$(HOME)", "${HOME}" ]:
        # we use a $HOME variable, for convenience
        smsHelpers.base = smsHelpers.base.replace( i,os.environ["HOME"] )
    for i in [ "$SRC", "$(SRC)", "${SRC}" ]:
        # we use a $SRC variable, for convenience
        smsHelpers.base = smsHelpers.base.replace( i,installDirectory() )

    #check if database directory exists, read database version
    try:
        databaseVersion = smsHelpers.databaseVersion()
    except:
        log.error("Database not found in %s" % os.path.realpath(smsHelpers.base))
        databaseVersion = None

    #initialize output status
    outputStatus = ioObjects.OutputStatus(inStat, warnings, databaseVersion)

    if outputStatus.status == -2 or outputStatus.status == -4:
        #do not continue decomposition for bad input files, or in case database is not found
        outputStatus.printout("file",outputfile)
        inputStatus.printout("file",outputfile)
        sys.exit()

    #FIXME can I automatize calling the printout functions and the sys.exit?

    # sqrts from parameter input file
    sqrts = parser.getfloat("options","sqrts")*TeV

    #cross section computation
    #if addMissingXsecs is set True, SModelS calculates LO cross sections
    #production processes not yet covered in the input file are added
    #else, if the input slha does not contain cross sections, SModelS will add LO cross sections
    if parser.getboolean("options","doSLHAdec") and (parser.getboolean("options","addMissingXsecs") or inputStatus.xsec[0]==-1):
        xsecs = xsecComputer.computeXSec(sqrts, 0, parser.getint("options","nevts"), inputFile)
        comment = "Nevts: " + parser.get("options","nevts")
        xsecComputer.addXSecToFile(xsecs, inputFile, comment)

    #option of adding NLO and NLL cross sections from nllFast
    #using LO cross sections written in the input slha file
    if parser.getboolean("options","doSLHAdec") and parser.getboolean("options","addnlo"):
        #read LO from input file, compute NLO
        xsecs_nlo = xsecComputer.computeXSec(sqrts, 1, parser.getint("options","nevts"), inputFile, loFromSlha=True)
        if xsecs_nlo: xsecComputer.addXSecToFile(xsecs_nlo, inputFile) #FIXME comment! (to be printed to slha file)

    if parser.getboolean("options","doSLHAdec") and parser.getboolean("options","addnll"):
        xsecs_nll = xsecComputer.computeXSec(sqrts, 2, parser.getint("options","nevts"), inputFile, loFromSlha=True)
        if xsecs_nll: xsecComputer.addXSecToFile(xsecs_nll, inputFile) #FIXME comment!

    #after cross section computation, re-evaluate slha file for displaced vertices of high cross section
    if parser.getboolean("options","doSLHAdec"):
        s, m = inputStatus.reEvaluateDisplaced()
        if s == -1:
            outputStatus.addWarning(m)
            outputStatus.updateSLHAStatus(-1)
            outputStatus.updateStatus(-2)
            outputStatus.printout("file", outputfile)
            inputStatus.printout("file", outputfile)
            sys.exit()


    #decomposition

    if parser.getboolean("options","doSLHAdec"):
        #sigmacut = minimum value of cross-section for an element to be considered eligible for decomposition.
        #Too small sigmacut leads to too large decomposition time.
        sigmacut = parser.getfloat("parameters","sigmacut")*fb

    try:
        # Decompose input SLHA file, store the output elements in smstoplist
        if parser.getboolean("options","doSLHAdec"):
            smstoplist = slhaDecomposer.decompose(inputFile, sigmacut, doCompress=parser.getboolean("options","doCompress"),
                         doInvisible=parser.getboolean("options","doInvisible"), minmassgap=parser.getfloat("parameters","minmassgap")*GeV)
        else:
            smstoplist = lheDecomposer.decompose(inputFile, doCompress=parser.getboolean("options","doCompress"),
                         doInvisible=parser.getboolean("options","doInvisible"), minmassgap=parser.getfloat("parameters","minmassgap")*GeV)   
    except:
        outputStatus.updateStatus(-1)
        outputStatus.printout("file", outputfile)
        inputStatus.printout("file", outputfile)
        sys.exit()

    # If no topologies with sigma > sigmacut are found, update status, write output file, stop running
    if not smstoplist:
        outputStatus.updateStatus(-3)
        outputStatus.printout("file", outputfile)
        inputStatus.printout("file", outputfile)
        sys.exit()

    if parser.getboolean("stdout","printGtop"):
        smstoplist.printout()

    # This is my porposed format for element tabel
    if parser.getboolean("stdout","printThEl"):
        for (i,topo) in enumerate(smstoplist):
            print '\n'
            print "A new global topoloy starts here" 
            print "====================================================================="
            for j, el in enumerate(topo.elementList):
                print "\t ........................................................................."
                el.printout()
        print "====================================================================="
        print "====================================================================="
        print "The list ends here" 
        print "====================================================================="
        print "====================================================================="

    analyses = parser.get("data","analyses")
    #in case that a list of analyses is given, retrieve list
    if "," in analyses: analyses = analyses.split(",")

    topologies = parser.get("data","topologies")
    if "," in topologies: topologies = topologies.split(",")

    # Load analyses
    listofanalyses = smsAnalysisFactory.load(analyses, topologies)

    #define result list that collects all theoryprediction objects
    #variables set to define printing options
    results = ioObjects.ResultList(bestresultonly = not parser.getboolean("file","expandedSummary"), describeTopo = parser.getboolean("file","describeTopo"))

    constrainedElements = []

    if parser.getboolean("stdout","printAnaEl"):
        for analysis in listofanalyses:
            elements = _getElementsFrom(smstoplist, analysis)
            if len(elements) == 0: continue
            # This is my porposed format for analyses elements table
            print "========================================================"
            print "Analysis Name:", analysis.label.split(":")[0]
            print "Analysis Topology:", analysis.label.split(":")[1]
            print "Analysis Sqrts:", analysis.sqrts
            print "========================================================"
            ref_el = None
            for el in elements:
                el.printout()
                print "........................................................"

    #Get theory prediction for each analysis and print basic output
    for analysis in listofanalyses:
        theorypredictions = theoryPredictionFor(analysis, smstoplist)
        if not theorypredictions:
            continue
        if parser.getboolean("stdout","printResults"):
            print "================================================================================"
            theorypredictions.printout() # again, check print function
        print "................................................................................"

        # Create a list of results, to determine the best result
        for theoryprediction in theorypredictions:
            results.addResult(theoryprediction, maxcond = parser.getfloat("parameters","maxcond"))

    # If there is no best result, this means that there are no matching experimental results for the point.
    # Decomposition status has following flags:
    #-1: "#could not run the decomposition",
    #-3: "#no cross sections above sigmacut found",
    #-2: "#bad input slha, did not run decomposition",
    #0: "#no matching experimental results",
    #1: "#decomposition was successful". 

    if results.isEmpty():
        # no experimental constraints found
        outputStatus.updateStatus(0)
    else:
        outputStatus.updateStatus(1)

    # write output file
    outputStatus.printout("file", outputfile)
    inputStatus.printout("file", outputfile)
    # add experimental constraints if found
    if outputStatus.status == 1: results.printout("file", outputfile)

    sqrtsValues = [xsec.info.sqrts for xsec in smstoplist.getTotalWeight()]
    if not sqrts in sqrtsValues: sqrts = max(sqrtsValues)

    if parser.getboolean("options","findMissingTopos"):
        #look for missing topologies, add them to the output file
        missingtopos = missingTopologies.MissingTopoList(sqrts)
        missingtopos.findMissingTopos(smstoplist, listofanalyses, parser.getfloat("parameters","minmassgap")*GeV)
        missingtopos.printout("file", outputfile)



if __name__ == "__main__":
    #get the name of input slha file (and parameter file)
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-f', '--filename', help = 'name of SLHA or LHE input file, necessary input', required = True)
    argparser.add_argument('-p', '--parameterfile', help = 'name of parameter file, optional argument, if not set, use all parameters from etc/parameters_default.ini. For directory names, "$HOME" can be used for the user home directory, "$SRC" encodes the base directory of the source code.')
    argparser.add_argument('-o', '--outputfile', help = 'name of output file, optional argument, default is: ./summary.txt', default = 'summary.txt')
    args = argparser.parse_args() # pylint: disable-msg=C0103

    # execute run function
    main(args.filename, args.parameterfile, args.outputfile)
