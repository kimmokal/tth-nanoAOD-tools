#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

#from  PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import *
from  tthAnalysis.NanoAODTools.postprocessing.tthModules import *
p=PostProcessor(".",inputFiles(),"",modules=[susy1lepdata(),susy_1l_Trigg2016(),susy_1l_FiltersData(),susy1lepTOPData()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
p.run()
print "DONE"
os.system("ls -lR")

