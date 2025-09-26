#!/usr/bin/env python3

""" Code that gives possible tpred combinations (without weights) manually without using smodels isCombinableWith """

from tester.combinationsmatrix import getMatrix
from typing import Dict,Union,List
from smodels.matching.theoryPrediction import theoryPredictionsFor, TheoryPrediction
from tester.combinationsmatrix import getYamlMatrix

class combFinder(object):

    def __init__(self, combination_matrix: Union[dict,None] = None ):
        """
        :param combination_matrix: dictionary of allowed analyses combination/strategy;
                                   if None, use default combination_matrix present
                                   in tester/combinationsmatrix.py
        """
        combinationsmatrix, status = getYamlMatrix()
        if not combinationsmatrix or status != 0:
            sys.exit("Combination matrix not loaded correctly when instantiating CombinationFinder class.")

        self.cM = combinationsmatrix
        if combination_matrix != None:
            self.cM = combination_matrix

        self._createDefaultCombinationMatrix()


    def _createDefaultCombinationMatrix(self):
        """create a True/False dict for all analyses present in the combination matrix"""

        all_ana = list(self.cM.keys())
        #print(all_ana)
        self.__eM = {ana:{combAna:False for combAna in all_ana} for ana in all_ana}

        for id1,ana in enumerate(all_ana):
            #if ana has no combinable analyses mentioned in the combination dictionary
            if self.cM.get(ana) == []: continue
            #if '2013-11' in ana: print("\n check 1 , ",ana, self.cM.get(ana))
            for combAna in all_ana[id1+1:]:
                #if combAna not present in the list of analyses combinable with ana
                if combAna not in self.cM.get(ana): continue
                #if combAna is present in the list of analyses combinable with ana
                else:
                    self.__eM[ana][combAna] = True
                    self.__eM[combAna][ana] = True

        #print("excl matrix = ", self.__eM)

    def getCombinationMatrix(self) -> Dict:
        """
        :param theoryPredictionList: List of theory predictions
        :returns: a dict of true/false combinations for the theory predictions
        """


        tp_ana = [tp.analysisId() for tp in self.listoftp if tp.analysisId() in self.__eM.keys()]
        self.listoftp = [tp for tp in self.listoftp if tp.analysisId() in tp_ana]

        #print("\n tp_ana ", tp_ana)
        #ext = ["agg","ma5","eff","adl","slv1"]
        '''
        #EXTENSIONS SHOULD BE IGNORED
        for ana in tp_ana:
            for ex in ext:
                if ex in ana: ana = ana.replace(f"-{ex}","")
        '''

        combMatrix = {ana:{combAna:self.__eM[ana][combAna] for combAna in tp_ana} for ana in tp_ana}
        '''
        for key, value in combMatrix.items():
            if '2013-11' in key:
                print("\n check = ",key, value)
        '''
        #print("combMatrix ", combMatrix)
        def checkCombinable(a1, a2):
            "Check if two analyses are combinable if not specified in combination dictionary"
            sq_s1 = a1.dataset.globalInfo.sqrts
            sq_s2 = a2.dataset.globalInfo.sqrts
            exp1 = a1.analysisId().split('-')[0]                                #split at - and get 'CMS'/'ATLAS'
            exp2 = a2.analysisId().split('-')[0]
            ana1 = a1.analysisId()
            ana2 = a2.analysisId()
            '''
            if '2017-03' in ana1:
                print(ana1, ana2)
                print("exp 1", exp1)
                print("exp 2", exp2)
                print("s1", sq_s1)
                print("s2", sq_s2)
            '''
            if exp1 != exp2: return True                                      #if diff expts return True
            elif sq_s1 != sq_s2 : return True                                   #if diff sqrts return True
            else: return False


        #check if two analyses can be combined (do we want it now?)
        for i,ana in enumerate(tp_ana):
            for j,combAna in enumerate(tp_ana[i+1:]):
                combinable = combMatrix[ana][combAna]
                #if '2013-11' in ana and '2013-12' in combAna:
                #print(ana, combAna, " combinable = ",combinable)
                if not combinable:
                    combinable = checkCombinable(self.listoftp[i],self.listoftp[i+j+1])
                    if combinable:
                        #print("combinable though not present in matrix:", ana,combAna)
                        combMatrix[ana][combAna] = True
                        combMatrix[combAna][ana] = True

        return combMatrix


    def isSubsetOf(self, newcomb_tuple, combinations) -> bool:
        """ checks if a comb already exists or is a subset of some other comb
            in the current list of combinations
            :param newcomb_tuple: tuple of combinable theory predictions
            :param combinations: set of tuples of combinable theory predictions
        """
        for comb_tuple in combinations:
            if set(newcomb_tuple).issubset(set(comb_tuple)):
                return True
        return False


    def getPossibleCombinations(self, predictions: list[TheoryPrediction]) -> List:
        """
        :param theoryPredictionList: List of theory predictions
        :returns: a list of all possible combinable theory predictions from the input theory predictions
        """
        self.listoftp = predictions
        combMatrix = self.getCombinationMatrix()

        combinables = set()
        for i,pred in enumerate(self.listoftp[:-1]):
            ana1 = pred.analysisId()
            for combPred in self.listoftp[i+1:]:
                ana2 = combPred.analysisId()
                canCombine = combMatrix[ana1][ana2]

                if canCombine:
                    #loop over existing combinations to see if it can append to any existing combination
                    update_Comb = False
                    old_comb = set()
                    new_comb = set()
                    new_combinables = sorted(list(combinables), key = lambda x: len(x), reverse=True)
                    for j,comb in enumerate(new_combinables):
                        add_pred = False
                        currentCombine = [False]
                        if pred in comb:
                            if combPred in comb:
                                update_Comb = True      #if already present in some other previous combination
                                continue
                            else:
                                currentCombine = [(True if combMatrix[tp.analysisId()][ana2] else False) for tp in comb] #check if combPred can be added to current comb

                        if False not in currentCombine:
                            old_comb.add(comb)        #add current comb to old comb
                            comb = comb + (combPred,)
                            new_comb.add(comb)        #add updated comb to new_comb
                            update_Comb = True

                        elif True in currentCombine:
                            #combPred can be combined with some analyses in current comb, store in new_comb
                            nc = [c for i,c in enumerate(comb) if currentCombine[i] == True]
                            nc.append(combPred)
                            nc = tuple(nc)
                            already_in = self.isSubsetOf(nc,combinables)        #check if nc is in combinables already
                            if not already_in:
                                already_in = self.isSubsetOf(nc,new_comb)       #check if nc is in new_comb already
                                if not already_in: new_comb.add(nc)             #add updated comb to new_comb
                            update_Comb = True

                        else: continue


                    #add new_combs to combinables and remove old_combs from combinables
                    if update_Comb:
                        combinables = combinables.difference(old_comb)
                        combinables = combinables.union(new_comb)
                    #add the two analyses in combinables as a tuple
                    else: combinables.add((pred,combPred))

        combinables = sorted(list(combinables), key = lambda x: len(x), reverse=True) #sort according to decreasing length of combinable theory predictions
        combinables = [list(comb) for comb in combinables]                            #convert tuples of combinations to list
        aids = [[c.analysisId() for c in comb] for comb in combinables]
        return combinables


        '''
        def get_weights()
        combMatrix = np.array(combMatrix)


        #get weights
        lbsm = numpy.array([preds.likelihood(1., expected = exp) for preds in self.listoftp], dtype=object)
        lsm = numpy.array([preds.likelihood(0., expected = exp) for preds in self.listoftp], dtype=object)
        weight_vector = np.log(lbsm/lsm)              #return llh ratio as discovery mode #check if sm is ever zero?

        #call pathfinder
        bam = pf.BinaryAcceptance(combMatrix, weights=np.array(weight_vector))

        #Get the allowed list of combinations with decreasing weights
        whdfs = pf.WHDFS(bam, top=self.ntop)
        whdfs.find_paths()

        #return list of theory predictions for which the combination has max weight
        top_path = whdfs.get_paths[0]  # gets indices of analyses which are best combinable
        best_comb = [self.listoftp[i] for i in top_path]
        '''

if __name__ == "__main__":
    from smodels.experiment.databaseObj import Database
    from smodels.base.model import Model
    from share.model_spec import BSMList
    from smodels.share.models.SMparticles import SMList
    from smodels.base.physicsUnits import fb, GeV
    from smodels.decomposition import decomposer

    #filename = "gluino_squarks.slha"
    filename = "ew_lwd7kc46.slha"
    #filename = "ew_3t6no481.slha"
    model = Model(BSMparticles = BSMList, SMparticles = SMList)
    model.updateParticles(inputFile = filename)

    print("Decomposing model")
    toplist = decomposer.decompose(model, 0.005*fb, doCompress=True, doInvisible=True, minmassgap=5.*GeV)
    #toplist = decomposer.decompose(model, 0.005*fb, massCompress=True, invisibleCompress=True, minmassgap=5.*GeV)

    listOfAna = ['ATLAS-SUSY-2018-32','ATLAS-SUSY-2018-41','ATLAS-SUSY-2019-02']
    listOfAna = ['all']
    comb_dict = {"ATLAS-SUSY-2018-32":['ATLAS-SUSY-2018-41'], "ATLAS-SUSY-2018-41":['ATLAS-SUSY-2018-32', 'ATLAS-SUSY-2019-02'], "ATLAS-SUSY-2019-02":['ATLAS-SUSY-2018-41']}
    db = Database ( "official" )

    print("Getting experimental Results")
    expresults = db.getExpResults(analysisIDs=listOfAna, dataTypes=['efficiencyMap','combined'])

    print("Finding theory Predictions")
    allPreds = theoryPredictionsFor(expresults, toplist, combinedResults=True)
    allPredsana = [tp.analysisId() for tp in allPreds]
    print("preds ", allPredsana)
    print("Getting combinations")
    bC = CombinationFinder()
    bestThPred = bC.getPossibleCombinations(allPreds)


    if bestThPred ==[] : print("\n Model Point: ", filename, "  , No predictions")
    else:
        for bp in bestThPred:
            ana = [b.analysisId() for b in bp]
            print("\n Model Point: ", filename, " , Combination: ", ana)
