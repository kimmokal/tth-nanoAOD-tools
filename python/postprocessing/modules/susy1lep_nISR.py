import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class susy1lepnISR(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		self.out.branch("nIsr","I");
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
		Jets = Collection(event, "Jet")
		jets = [j for j in Jets if j.cleanmask == True]
		genpart = Collection(event, "GenPart")
		nIsr = 0
		for jet in jets:
			if jet.pt <30.0: continue
			if abs(jet.eta )>2.4: continue
			matched = False
			for mc in genpart:
				if matched: break
				if (mc.status!=23 or abs(mc.pdgId)>5): continue
				momid = abs(mc[genpart.genPartIdxMother].pdgId)
				if not (momid==6 or momid==23 or momid==24 or momid==25 or momid>1e6): continue
					#check against daughter in case of hard initial splitting
				for idau in range(mc.numberOfDaughters()):
					dR = math.sqrt(deltaR2(jet.eta(),jet.phi(), mc.daughter(idau).p4().eta(),mc.daughter(idau).p4().phi()))
					if dR<0.3:
						matched = True
						break
			if not matched:
				nIsr+=1
		self.out.fillBranch("nIsr",nIsr)	
			
		pass
		return True
        

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy_1l_nISR  = lambda : susy1lepnISR()

