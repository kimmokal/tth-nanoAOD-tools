import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


def deltaR2( e1, p1, e2=None, p2=None):
	"""Take either 4 arguments (eta,phi, eta,phi) or two objects that have 'eta', 'phi' methods)"""
	if (e2 == None and p2 == None):
		return deltaR2(e1.eta,e1.phi, p1.eta, p1.phi)
	de = e1 - e2
	dp = deltaPhi(p1, p2)
	return de*de + dp*dp


def deltaPhi( p1, p2):
	'''Computes delta phi, handling periodic limit conditions.'''
	res = p1 - p2
	while res > math.pi:
		res -= 2*math.pi
	while res < -math.pi:
		res += 2*math.pi
	return res

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
		self.out.branch("nISRweight","F");
		self.out.branch("nISRttweightsyst_up","F");
		self.out.branch("nISRttweightsyst_down","F");
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
		Jets = Collection(event, "Jet")
		jets = [j for j in Jets if j.cleanmask == True]
		genpart = Collection(event, "GenPart")
		genParts = [l for l in genpart]
		# get the particles when they have a mother ---> getting the daughters only 
		daughters = [l for l in genpart if l.genPartIdxMother >= 0 ]
		nIsr = 0
		for jet in jets:
			if jet.pt <30.0: continue
			if abs(jet.eta )>2.4: continue
			matched = False
			
			for mc in genParts:	
				if matched: break
				if (mc.status!=23 or abs(mc.pdgId)>5): continue
				momid = abs(genParts[mc.genPartIdxMother].pdgId)
				if not (momid==6 or momid==23 or momid==24 or momid==25 or momid>1e6): continue
					#check against daughter in case of hard initial splitting
				for idau in range(len(daughters)) :#(mc.numberOfDaughters()):
					dR = math.sqrt(deltaR2(jet.eta,jet.phi, daughters[idau].eta,daughters[idau].phi))
					if dR<0.3:
						matched = True
						break
			if not matched:
				nIsr+=1
		self.out.fillBranch("nIsr",nIsr)
		nISRweight = 1
		ISRweights = { 0: 1, 1 : 0.920, 2 : 0.821, 3 : 0.715, 4 : 0.662, 5 : 0.561, 6 : 0.511}
		ISRweightssyst = { 0: 0.0, 1 : 0.040, 2 : 0.090, 3 : 0.143, 4 : 0.169, 5 : 0.219, 6 : 0.244}
		#if 'TTJets' in inputFile.GetName() or 'T1tttt' in inputFile.GetName():
		nISRforWeights = int(nIsr)
		if nIsr > 6:
			nISRforWeights = 6
		
		C_ISR = 1.090
		C_ISR_up   = 1.043
		C_ISR_down = 1.141
		nISRweight = C_ISR * ISRweights[nISRforWeights]
		nISRweightsyst_up   =  C_ISR_up   * (ISRweights[nISRforWeights] + ISRweightssyst[nISRforWeights])
		nISRweightsyst_down =  C_ISR_down * (ISRweights[nISRforWeights] - ISRweightssyst[nISRforWeights])
		
		self.out.fillBranch("nISRweight",nISRweight)
		self.out.fillBranch("nISRttweightsyst_up",nISRweightsyst_up)
		self.out.fillBranch("nISRttweightsyst_down",nISRweightsyst_down)


            # ------ Forwarded Message --------
            # Subject: Re: question for ttbar ISR reweighting
            # Date: Sat, 14 Jan 2017 20:24:14 +0100
            # From: Manuel Franco Sevilla <manuel.franco.sevilla@cern.ch>
            #The [Nom, Up, Down] values we find for the events with Nisr = 0 are:
            #[1.090, 1.043, 1.141]: TTJets_Tune
            #[1.096, 1.046, 1.151]: TTJets_SingleLeptFromT
            #[1.116, 1.055, 1.185]: TTJets_DiLept
		
		
		return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy_1l_nISR  = lambda : susy1lepnISR()

