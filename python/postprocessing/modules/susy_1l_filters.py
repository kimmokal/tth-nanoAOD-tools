import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


filterName = "/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src/tthAnalysis/NanoAODTools/data/dummyFilterlist.txt" #to process 2016 data quickly...
filterList = None


def readList(fname):
    evList = set()
    with open(fname,"r") as flist:
        for line in flist.readlines():
            if ":" not in line: continue
            sline = line.split(":")
            if len(sline) != 3: continue
            evList.add((int(sline[0]),int(sline[1]),int(sline[2])))
    return evList


class susy1lepFilter(Module):
    def __init__(self, isData):
        self.branches = [
            'passFilters','passFiltersMoriond2017Tight','passCSCFilterList',
            ]
        self.isData = isData
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.inFile = inputFile.GetName()
        self.out = wrappedOutputTree

        for var in self.branches:
            self.out.branch(var, "O")
#            print "Creating branch ", var, " boolean, O"

        print self.inFile
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        
        runNr = getattr(event, "run")
        lumiNr = getattr(event, "luminosityBlock")
        eventNr = getattr(event, "event")        
        #print "Running on ", runNr, ":", lumiNr, ":", eventNr
        
        if (not self.isData):
            passCSCFilterList           = True
            passFilters                 = True
            passFiltersMoriond2017Tight = True
        else:
            global filterList
            if filterList == None: filterList = readList(filterName)

            # check MET text filter files
            if (runNr,lumiNr,eventNr) in filterList:
                #print "yes", runNr,lumiNr,eventNr
                passCSCFilterList = False
            else:
                #print "no", runNr,lumiNr,eventNr
                passCSCFilterList = True
        
            # check filters present in event
            if hasattr(event,"Flag_eeBadScFilter"):
                flag_eeBadScFilter = getattr(event, "Flag_eeBadScFilter")
                flag_HBHENoiseFilter = getattr(event, "Flag_HBHENoiseFilter")
                flag_HBHENoiseIsoFilter = getattr(event, "Flag_HBHENoiseIsoFilter")
                flag_EcalDeadCellTriggerPrimitiveFilter = getattr(event, "Flag_EcalDeadCellTriggerPrimitiveFilter")
                flag_goodVertices = getattr(event, "Flag_goodVertices")
                flag_globalTightHalo2016Filter = getattr(event, "Flag_globalTightHalo2016Filter")
                flag_badMuons = getattr(event, "Flag_badMuons")
                flag_duplicateMuons = getattr(event, "Flag_duplicateMuons")
                
                # Suggested by Isabell on 04.04.2018 meeting:
                flag_badChargedHadronSummer2016 = getattr(event, "Flag_chargedHadronTrackResolutionFilter")
                flag_badMuonSummer2016 = getattr(event, "Flag_muonBadTrackFilter")
                
                #forMoriond2017 https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
                #for Moriond 2017: use updated badChargedHadron and badPFMuon filters (Ece's Summer2016 implementation). Do NOT use Flag_badMuons and Flag_duplicateMuons (they are only to be used if new tails would appear in the metMuEGClean collection comparing to the METUncorrected collection)
                passFilters = flag_HBHENoiseFilter and flag_HBHENoiseIsoFilter and flag_EcalDeadCellTriggerPrimitiveFilter and  flag_goodVertices and flag_eeBadScFilter and flag_globalTightHalo2016Filter and flag_badChargedHadronSummer2016 and flag_badMuonSummer2016
                #also apply Flag_badMuons and Flag_duplicateMuons (they are only to be used if new tails would appear in the metMuEGClean collection comparing to the METUncorrected collection)
                passFiltersMoriond2017Tight = passFilters and not flag_badMuons and not flag_duplicateMuons
            else:
                passFilters                 = True
                passFiltersMoriond2017Tight = True
        
        self.out.fillBranch('passFilters', passFilters)
        self.out.fillBranch('passFiltersMoriond2017Tight', passFiltersMoriond2017Tight)
        self.out.fillBranch('passCSCFilterList', passCSCFilterList)
        
        return True
        

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy_1l_FiltersMC   = lambda : susy1lepFilter(False)
susy_1l_FiltersData = lambda : susy1lepFilter(True)
