import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class susy1lepTrig(Module):
    def __init__(self):
        self.branches = [
            ## Custom trigger names
            'HLT_EleOR', 'HLT_MuOR','HLT_LepOR','HLT_MetOR',
            ## Trigger efficiencies
            'TrigEff'
            ]
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.inFile = inputFile.GetName()
        self.out = wrappedOutputTree

        for var in self.branches:
            if 'HLT_' in var:
                self.out.branch(var, "O")
            else:
                self.out.branch(var, "F")
        print self.inFile
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
		
		#electron and muon collections for TrigEff:
		electrons = Collection(event, "Electron")
		muons = Collection(event, "Muon")
		
		# triggers for EleOR:
		trig_Ele105 = getattr(event, "HLT_Ele105_CaloIdVT_GsfTrkIdT")
		trig_Ele115 = getattr(event, "HLT_Ele115_CaloIdVT_GsfTrkIdT")
		trig_Ele50PFJet165 = getattr(event, "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165")
		trig_IsoEle27T = getattr(event, "HLT_Ele27_WPTight_Gsf")
		if hasattr(event, "HLT_Ele15_IsoVVVL_PFHT400"):
			trig_EleHT400 = getattr(event, "HLT_Ele15_IsoVVVL_PFHT400")
		else :trig_EleHT400 = True
		if hasattr(event, "HLT_Ele15_IsoVVVL_PFHT350"):
			trig_EleHT350 = getattr(event, "HLT_Ele15_IsoVVVL_PFHT350")
		else : trig_EleHT350 = True 
		
		# triggers for MuOR:
		trig_Mu50 = getattr(event, "HLT_Mu50")
		trig_IsoMu24 = getattr(event, "HLT_IsoMu24") or getattr(event, "HLT_IsoTkMu24")
		if hasattr(event, "HLT_Mu15_IsoVVVL_PFHT400"):
			trig_MuHT400 = getattr(event, "HLT_Mu15_IsoVVVL_PFHT400")
		else : trig_MuHT400 = True 
		if hasattr(event, "HLT_Mu15_IsoVVVL_PFHT350"): 
			trig_MuHT350 = getattr(event, "HLT_Mu15_IsoVVVL_PFHT350")
		else : trig_MuHT350 = True 
		
		# triggers for MetOR:
		if hasattr(event, "HLT_PFMET100_PFMHT100_IDTight"):
			trig_MET100MHT100 = getattr(event, "HLT_PFMET100_PFMHT100_IDTight")
		else : trig_MET100MHT100 = True
		if hasattr(event, "HLT_PFMET110_PFMHT110_IDTight"):
			trig_MET110MHT110 = getattr(event, "HLT_PFMET110_PFMHT110_IDTight")
		else :trig_MET110MHT110 = True 
		if hasattr(event, "HLT_PFMET120_PFMHT120_IDTight"):
			trig_MET120MHT120 = getattr(event, "HLT_PFMET120_PFMHT120_IDTight")
		else : trig_MET120MHT120 = True 
		
		
		if 'Ele' in self.inFile or 'Mu' in self.inFile:
			trig_Eff = 1.0
		elif electrons._len>=1 and muons._len==0: trig_Eff = 0.963 # ele efficieny (for 2016 4/fb)
		elif muons._len>=1 and electrons._len==0: trig_Eff = 0.926 # mu efficieny (for 2016 4/fb)
		else: trig_Eff = 1.0
		
		trig_EleOR = trig_Ele105 or trig_Ele115 or trig_Ele50PFJet165 or trig_IsoEle27T or trig_EleHT400 or trig_EleHT350
		trig_MuOR = trig_Mu50 or trig_IsoMu24 or trig_MuHT400 or trig_MuHT350
		trig_MetOR = trig_MET100MHT100 or trig_MET110MHT110 or trig_MET120MHT120
		
		self.out.fillBranch('TrigEff', trig_Eff)
		self.out.fillBranch('HLT_EleOR', trig_EleOR)
		self.out.fillBranch('HLT_MuOR', trig_MuOR)
		self.out.fillBranch('HLT_LepOR', trig_EleOR or trig_MuOR)
		self.out.fillBranch('HLT_MetOR', trig_MetOR)
		
		return True
        

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy_1l_Trigg = lambda : susy1lepTrig()
