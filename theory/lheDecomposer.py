"""
.. module:: theory.lheDecomposer
:synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import LHEReader
import topology
import crossSection
import element
import pyslha
import branch
import ParticleNames
from tools.physicsUnits import addunit
import logging
import copy

logger = logging.getLogger(__name__)


def decompose(lhefile, inputXsecs=None, nevts=None, doCompress=False, 
              doInvisible=False, minmassgap=None):
    """
    Do LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param inputXsecs: xSectionList object with cross-sections for the mothers
    appearing in the LHE file. If None, use information from file.
    :param nevts: (maximum) number of events used in the decomposition. If
    None, all events from file are processed.
    :param doCompress: mass compression option (True/False)
    :param doInvisible: invisible compression option (True/False)
    :param minmassgap: minimum mass gap for mass compression (only used if
    doCompress=True)
    :returns: a TopologyList object 
    
    """
    reader = LHEReader.LHEReader(lhefile, nevts)
    smsTopList = topology.TopologyList ( )
    # get cross-section from file (= event weight, assuming a common weight for
    # all events)
    if not inputXsecs:
        xSectionList = crossSection.getXsecFromLHEFile(lhefile, 
                                                       addEvents=False)
    else:
        xSectionList = inputXsecs

    # Loop over events and decompose 
    for event in reader:
        momPDG = tuple(sorted(event.getMom()))  # Get mother PDGs
        eventweight = xSectionList.getXsecsFor(momPDG)        
        # Get event element        
        newElement = elementFromEvent(event, eventweight)
        allElements = [newElement] 
        # Do compression:        
        if doCompress or doInvisible:
            allElements += newElement.compressElement(doCompress, doInvisible,

        for el in allElements:
            top = topology.Topology(el)                      
            smsTopList.addList([top])                   

    return smsTopList


def elementFromEvent(event, weight=None):
    """
    Creates an element from a LHE event and the corresponding event weight.
    
    :param event: LHE event
    :param weight: event weight. Must be a XSectionList object (usually with a
    single entry) or None if not specified.
    :returns: element
    
    """
    if not event.particles:
        logger.error('Empty event!')
        return None

    brDic, massDic = getDictionariesFromEvent(event)
          
    # Creates branch list
    finalBranchList = []
    for ip,particle in enumerate(event.particles):
        # Means particle came from initial state (primary mother)
        if 1 in particle.moms:
            mombranch = branch.Branch()
            mombranch.momID = particle.pdg
            mombranch.daughterID = particle.pdg
            if weight: mombranch.maxWeight = weight.getMaxXsec()            
            # Get simple BR and Mass dictionaries for the corresponding branch
            branchBR = brDic[ip]
            branchMass = massDic[ip]
            # Generate final branches (after all R-odd particles have decayed)
            finalBranchList += branch.decayBranches([mombranch], branchBR, branchMass,
                                           sigcut=addunit(0., 'fb'))     
    
    if len(finalBranchList) != 2:
        logger.error(str(len(finalBranchList))+" branches found in event. " +
                     "R-parity violation?")
        return False
    # Finally create element from event:
    newElement = element.Element(finalBranchList)    
    if weight: newElement.weight = copy.deepcopy(weight)

    return newElement


def getDictionariesFromEvent(event):
    """
    Read an event and create simple mass and BR dictionaries 
    for each branch in the event
    
    :param event: LHE event
    :returns: BR and Mass dictionaries for the branches in the event
    
    """
    
    particles = event.particles
    
    # Identify and label individual branches:
    branchDic = {}
    for ip,particle in enumerate(particles):
        if particle.status == -1: continue
        if particles[particle.moms[0]].status == -1:
            #If a primary mother, the branch index is its own position
            initMom = ip
        else:
            #If not a primary mother, check if particle has a single parent (as it should)
            if particle.moms[0] != particle.moms[1] and min(particle.moms) != 0: 
                logger.error("More than one parent particle found!")
                return False        
            initMom = max(particle.moms)-1
            while particles[particles[initMom].moms[0]].status != -1:
                #Find primary mother (labels the branch)
                initMom = max(particles[initMom].moms)-1
        branchDic[ip] = initMom           
        
    # Get mass and BR dictionaries for all branches:
    massDic = {}
    brDic = {}
    for ibranch in branchDic.values():
        massDic[ibranch] = {}
        brDic[ibranch] = {}
    for ip,particle in enumerate(particles):        
        if particle.pdg in ParticleNames.Reven or particle.status == -1:
            # Ignore R-even particles and initial state particles
            continue
        ibranch = branchDic[ip]  #Get particle branch
        massDic[ibranch][particle.pdg] = addunit(particle.mass, 'GeV')       
    # Create empty BRs    
        brDic[ibranch][particle.pdg] = [pyslha.Decay(0., 0, [], particle.pdg)]
    
    # Get BRs from event       
    for ip,particle in enumerate(particles):        
        if particle.status == -1:
            # Ignore initial state particles
            continue
        if particles[particle.moms[0]].status == -1:
            # Ignore initial mothers
            continue        
        ibranch = branchDic[ip]
        momPdg = particles[max(particle.moms)-1].pdg
        if momPdg in ParticleNames.Reven:
            # Ignore R-even decays
            continue
        # BR = 1 always for an event        
        brDic[ibranch][momPdg][0].br = 1.
        brDic[ibranch][momPdg][0].nda += 1
        brDic[ibranch][momPdg][0].ids.append(particle.pdg)

    return brDic, massDic
    
