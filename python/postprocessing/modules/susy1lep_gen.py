import ROOT
import math, os

ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class genpartsusy(Module):
	def __init__(self):#, muonSelection, electronSelection):
		#self.isMC = isMC
		#self.isSig = isSig
		#self.muSel = muonSelection
		#self.elSel = electronSelection
		pass
	def beginJob(self):
		pass
	def endJob(self):
		pass
	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
	
		self.out.branch("GenDeltaPhiLepWSum","F");
		self.out.branch("GenDeltaPhiLepWDirect","F");
		self.out.branch("GenWSumMass","F");
		self.out.branch("GenWDirectMass","F");
		self.out.branch("nidxGenWs","I");
		self.out.branch("GenmTLepNu","F");
		self.out.branch("LeptonDecayChannelFlag","I");
		self.out.branch("genTau_grandmotherId","F");
		self.out.branch("genTau_motherId","F");
		self.out.branch("genLep_grandmotherId","F");
		self.out.branch("genLep_motherId","F");
		#self.out.branch("IsDiLepEvent","O");
		#self.out.branch("IsSemiLepEvent","O");
	def mt_2(p4one, p4two):
		return sqrt(2*p4one.Pt()*p4two.Pt()*(1-cos(p4one.Phi()-p4two.Phi())))
	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass
	def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
		genpart = Collection(event, "GenPart")
		# The following variables still need to be double-checked for validity
		genLeps = [l for l in genpart if l.pdgId == (13 or -13 or 11 or -11) and l.genPartIdxMother >= 0 ]
		# for some reason TTjets_* doesnot have GenPart_statusFla; better to check it before use it 
		if hasattr(event,"GenPart_statusFlags"):
			genLepFromTau = [l for l in genLeps if l.statusFlags == 2]
		else :genLepFromTau = [l for l in genLeps]
		genTaus = [l for l in genpart if l.pdgId == (15 or -15) and l.genPartIdxMother >= 0]
		genParts = [l for l in genpart]
		#leptons from tau decay https://github.com/cms-nanoAOD/cmssw/blob/master/PhysicsTools/NanoAOD/python/genparticles_cff.py#L67
		genLepsAndLepsFromTaus = [l for l in genLeps] + [ l for l in genLepFromTau]
		
		#### for genlept mother and grandmother ID 		
		#print genLepsAndLepsFromTaus
		ngenLepFromTau = len(genLepFromTau)
		ngenLeps = len(genLeps)
		ngenTaus = len(genTaus)
		ngenParts = len (genParts)
		ngenLepsAndLepsFromTau = len(genLepsAndLepsFromTaus)
		
		GenDeltaPhiLepWSum=-999
		GenDeltaPhiLepWDirect=-999
		GenWSumMass=-999
		GenWDirectMass=-999
		GenmTLepNu=-999
		LeptonDecayChannelFlag=-999 
		idx_genWs=[]
		idx_genLeps=[]
		idx_genNus=[]
		genLep_motherId = -999
		genLep_motherIdx = -999
		genLep_grandmotherId = -999
		genTau_motherId = -999
		genTau_motherIdx = -999
		genTau_grandmotherId = -999
		# find gen-level neutrinos (status 23), calculate deltaPhi (lep, nu), and genW-masses m(lep+nu)
		# for this: start from genLeps (status 23)
		for glep in genLeps: 
			if glep.status == 23 : 
				genLep_motherId = genpart[glep.genPartIdxMother].pdgId
				genLep_motherIdx = glep.genPartIdxMother
				genLep_grandmotherIdx = genpart[genLep_motherIdx].genPartIdxMother
				if genLep_grandmotherIdx >=0 : 
					genLep_grandmotherId = genpart[genLep_grandmotherIdx].pdgId
				self.out.fillBranch("genLep_motherId",genLep_motherId)
				self.out.fillBranch("genLep_grandmotherId",genLep_grandmotherId)
		
		for gtau in genTaus:  
			genTau_motherId = genpart[gtau.genPartIdxMother].pdgId
			genTau_motherIdx = gtau.genPartIdxMother
			genTau_grandmotherIdx = genpart[genTau_motherIdx].genPartIdxMother
			if genTau_grandmotherIdx >=0:
				genTau_grandmotherId = genpart[genTau_grandmotherIdx].pdgId
			self.out.fillBranch("genTau_motherId",genLep_motherId)
			self.out.fillBranch("genTau_grandmotherId",genTau_grandmotherId)
		
		for i_lep, genLep in enumerate(genLeps):
			if genLep.status == 23 and abs(genParts[genLep.genPartIdxMother].pdgId) == 24: # genLep is outgoing and has W as mother
				W_idx = genLep.genPartIdxMother
				idx_genWs.append(W_idx)
				idx_genLeps.append(i_lep)
				for i_nu, genPart in enumerate(genParts):
					if genPart.genPartIdxMother==W_idx and genPart.status == 23: # find W as mother
						if abs(genParts[genPart.genPartIdxMother].pdgId) == 12 or abs(genParts[genPart.genPartIdxMother].pdgId) == 14 or abs(genParts[genPart.genPartIdxMother].pdgId) == 16: #check whether it is a neutrino
							idx_genNus.append(i_nu)
		
		
		if(len(idx_genLeps)>=1):
			genLepP4 = genLeps[idx_genLeps[0]].p4()
			if ngenParts >=  idx_genNus:
				genNuP4 = genParts[idx_genNus[0]].p4()
				genWSumP4 = genLepP4 + genNuP4
				genWDirectP4 = genParts[genLeps[idx_genLeps[0]].genPartIdxMother].p4()
				GenDeltaPhiLepWSum = genLepP4.DeltaPhi(genWSumP4)
				GenDeltaPhiLepWDirect = genLepP4.DeltaPhi(genWDirectP4)
				GenWSumMass = genWSumP4.M()
				GenWDirectMass = genWDirectP4.M()
				GenmTLepNu = mt_2(genLepP4,genNuP4)
		
		#print ngenLepsAndLepsFromTau, ngenLeps + ngenTaus, ngenLepFromTau+ngenLeps
		assert ngenLepsAndLepsFromTau==ngenLepFromTau+ngenLeps
		if ngenLeps + ngenTaus ==2: #looking at semileptonic events
			IsDiLepEvent = True
			IsHadTauEvent = (ngenTaus > ngenLepFromTau)
			LeptonsInAcceptance = True
			PtMax = -999
			for l in genLepsAndLepsFromTaus: 
				if l.pt>PtMax: PtMax = l.pt 
		
			if IsHadTauEvent: LeptonDecayChannelFlag = 1 # preconfigure HadTau (becaus next loop won't start for two had taus in the event)
			for l in genLepsAndLepsFromTaus:
				if PtMax>=25 and l.pt<10: LeptonsInAcceptance=False
				if PtMax<25 and l.pt<5: LeptonsInAcceptance=False
				lepEta = abs(l.eta)
				if (lepEta>2.5): LeptonsInAcceptance=False
				if (abs(l.pdgId) == 11 and lepEta >= 1.44 and lepEta < 1.57): LeptonsInAcceptance=False
				
				if IsHadTauEvent and not LeptonsInAcceptance: LeptonDecayChannelFlag = 0 # OutOfAcceptance and HadTau
				elif IsHadTauEvent: LeptonDecayChannelFlag = 1 # HadTau (only)
				elif not LeptonsInAcceptance: LeptonDecayChannelFlag = 2 # OutOfAcceptance (only)
				else: LeptonDecayChannelFlag = 3 # Rest (Id/Isolation/Resolution)

		self.out.fillBranch("GenDeltaPhiLepWSum", GenDeltaPhiLepWSum) #initialize the dictionary with a first entry
		self.out.fillBranch("GenDeltaPhiLepWDirect", GenDeltaPhiLepWDirect)
		self.out.fillBranch("GenWSumMass", GenWSumMass)
		self.out.fillBranch("GenWDirectMass", GenWDirectMass)
		self.out.fillBranch("GenmTLepNu", GenmTLepNu)
		self.out.fillBranch("nidxGenWs", len(idx_genWs))
		self.out.fillBranch("LeptonDecayChannelFlag", LeptonDecayChannelFlag)
		
		return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy_1l_gen = lambda : genpartsusy()#,
 
