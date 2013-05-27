#!/usr/bin/env python

import SMSglobals, sys
from prettytable import PrettyTable
from Theory import SMSDecomp, SMSmethods, SMSxsec
from Tools.PhysicsUnits import addunit, rmvunit
from Tools import SMSPrettyPrinter
from Tools.SMSPrettyPrinter import wrap
from  Experiment import SMSanalyses, SMSAnalysisFactory, SMSgetlimit


DoFactory = True

printer=SMSPrettyPrinter.SMSPrettyPrinter()

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Creat analyses list:
if DoFactory:
  ListOfAnalyses = SMSAnalysisFactory.load()
else:  
  ListOfAnalyses = SMSanalyses.load()


#Generate events and compute cross-sections:
nevts = 10000
#slhafile = "AndreSLHA/andrePT4.slha"
slhafile = "slha/DESY_stop.slha"
Wv = SMSxsec.pytinit(nevts,slhafile,rpythia = True, donlo = True)
W = Wv["Wdic"]
Xsec = Wv["Xsecdic"]
lhefile = Wv["lhefile"]
CMdic = Wv["CMdic"]
SMSglobals.CMdic = CMdic



DoCompress = True
DoInvisible = True
minmassgap = addunit(10.,'GeV')

DoSLHAdec = True
if DoSLHAdec:
  maxlum = SMSmethods.getMaxLum(ListOfAnalyses) # Maximum cross-section*BR to be included
  if rmvunit(maxlum,'fb-1'):    
    sigmacut = addunit(1./rmvunit(maxlum,'fb-1'),'fb')
  else:
    sigmacut = addunit(0.1,'fb')
  if DoCompress or DoInvisible: sigmacut = sigmacut/10.  #When compression is turned on, relax sigmacut
  SMSTopList = SMSDecomp.SLHAdecomp(slhafile,Xsec,sigmacut,DoCompress,DoInvisible,minmassgap)
else:
  SMSTopList = SMSDecomp.LHEdecomp(lhefile,W,nevts,DoCompress,DoInvisible,minmassgap)


EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])

eltot = 0
totweight = []
#Print Results:
for i in range(len(SMSTopList)):
  weightlist = [el.weight for el in SMSTopList[i].ElList]
  sumw = SMSmethods.sumweights(weightlist)
  totweight.append(sumw)
  EvTop_table.add_row([i,SMSTopList[i].vertnumb,SMSTopList[i].vertparts,len(SMSTopList[i].ElList),wrap(printer.pformat(sumw),width=30)])
  eltot += len(SMSTopList[i].ElList)

 
      
#Print element list for Topology[i]:  
  if i >= 0:       
    for j in range(len(SMSTopList[i].ElList)):
      EvElement_table.add_row([i,j,SMSTopList[i].ElList[j].B[0].particles,SMSTopList[i].ElList[j].B[1].particles,wrap(printer.pformat(SMSTopList[i].ElList[j].B[0].masses),width=25),wrap(printer.pformat(SMSTopList[i].ElList[j].B[1].masses),width=25),wrap(printer.pformat(SMSTopList[i].ElList[j].weight),width=30)])
    EvElement_table.add_row(["---","---","---","---","---","---","---"])  
      
      
print "Number of Global topologies = ",len(SMSTopList)      
print(EvTop_table)
print "Total Number of Elements = ",eltot
print "Total weight = ",SMSmethods.sumweights(totweight)
#print(EvElement_table)

print '\n \n \n'


#sys.exit()
#Add event topologies to analyses:
for Analysis in ListOfAnalyses:
  Analysis.add(SMSTopList)

#Print analyses output:
AnElement_table = PrettyTable(["Analyses","Element","Masses","Element Weight"])  

for Ana in ListOfAnalyses:
  label = Ana.label
  ifirst = True
  for i in range(len(Ana.Top.ElList)):
    ptcs = [Ana.Top.ElList[i].B[0].particles,Ana.Top.ElList[i].B[1].particles]
    for j in range(len(Ana.Top.ElList[i].B[0].masses)):
      mass = [Ana.Top.ElList[i].B[0].masses[j],Ana.Top.ElList[i].B[1].masses[j]]
      if not ifirst: label = ""
      if j != 0: ptcs = ""
      ifirst = False
      AnElement_table.add_row([label,ptcs,wrap(printer.pformat(mass),width=100),wrap(printer.pformat(Ana.Top.ElList[i].weight[j]),width=30)])
  AnElement_table.add_row(["---","---","---","---"])  
    

#print(AnElement_table)


print '\n \n \n'
#sys.exit()

    
#Compute theoretical predictions to analyses results:
for Analysis in ListOfAnalyses:
  print "---------------Analysis Label = ",Analysis.label
  Results_table = PrettyTable(["Result","Conditions","Mass","Theoretical Value","Experimental Limits"])
  
  
  if max(Analysis.Top.vertnumb) > 3:
    print "Skipping analysis",Analysis.label,"with more than 3 vertices"
    continue


  for res in Analysis.results.keys():
    theoRes = SMSmethods.EvalRes(res,Analysis)  #Theoretical values for result and conditions
    

    try:
      plot = Analysis.plots[res]
    except KeyError:  
      plot = []
        
    for imass in range(len(theoRes)):
      mass = theoRes[imass]['mass']
      tvalue = theoRes[imass]['result']
      conds = theoRes[imass]['conditions']
      sigmalimit = SMSgetlimit.GetPlotLimit(mass,plot,Analysis)
      
      if sigmalimit and len(sigmalimit) > 0:
        Results_table.add_row([wrap(printer.pformat(res),width=30),wrap(printer.pformat(conds),width=30),wrap(printer.pformat(mass),width=30),wrap(printer.pformat(tvalue),width=30),wrap(printer.pformat(sigmalimit[0]),width=30)])
        for ilim in range(1,len(sigmalimit)):
          Results_table.add_row(["","","","",wrap(printer.pformat(sigmalimit[ilim]),width=30)])
      else:
        Results_table.add_row([wrap(printer.pformat(res),width=30),wrap(printer.pformat(conds),width=30),wrap(printer.pformat(mass),width=30),wrap(printer.pformat(tvalue),width=30),wrap(printer.pformat(sigmalimit),width=30)])
      Results_table.add_row(["---------","---------","---------","---------","---------"])

  print(Results_table)


sys.exit()

