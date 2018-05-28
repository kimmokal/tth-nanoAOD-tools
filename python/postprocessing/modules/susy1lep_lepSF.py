import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


#################
### Cuts and WP
#################

## Eta requirement
centralEta = 2.4
eleEta = 2.4
###########
# MUONS
###########

muID = 'medium' # 'medium'(2015) or 'ICHEPmediumMuonId' (2016)


###########
# Electrons
###########

eleID = 'CB' # 'MVA' or 'CB'

## PHYS14 IDs
## Non-triggering electron MVA id (Phys14 WP)
# Tight MVA WP
Ele_mvaPhys14_eta0p8_T = 0.73;
Ele_mvaPhys14_eta1p4_T = 0.57;
Ele_mvaPhys14_eta2p4_T = 0.05;
# Medium MVA WP  <--- UPDATE
Ele_mvaPhys14_eta0p8_M = 0.35;
Ele_mvaPhys14_eta1p4_M = 0.20;
Ele_mvaPhys14_eta2p4_M = -0.52;
# Loose MVA WP
Ele_mvaPhys14_eta0p8_L = 0.35;
Ele_mvaPhys14_eta1p4_L = 0.20;
Ele_mvaPhys14_eta2p4_L = -0.52;

## SPRING15 IDs
## Non-triggering electron MVA id (Spring15 WP):
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Electrons
# Tight MVA WP
Ele_mvaSpring15_eta0p8_T = 0.87;
Ele_mvaSpring15_eta1p4_T = 0.6;
Ele_mvaSpring15_eta2p4_T = 0.17;
# Medium MVA WP  <--- UPDATE
Ele_mvaSpring15_eta0p8_M = 0.35;
Ele_mvaSpring15_eta1p4_M = 0.20;
Ele_mvaSpring15_eta2p4_M = -0.52;
# Loose MVA WP
Ele_mvaSpring15_eta0p8_L = -0.16;
Ele_mvaSpring15_eta1p4_L = -0.65;
Ele_mvaSpring15_eta2p4_L = -0.74;
# VLoose MVA WP
Ele_mvaSpring15_eta0p8_VL = -0.11;
Ele_mvaSpring15_eta1p4_VL = -0.55;
Ele_mvaSpring15_eta2p4_VL = -0.74;

## Ele MVA check

## Isolation
ele_miniIsoCut = 0.1
muo_miniIsoCut = 0.2
Lep_miniIsoCut = 0.4
trig_miniIsoCut = 0.8

## Lepton cuts (for MVAID)
goodEl_lostHits = 0
goodEl_sip3d = 4
goodMu_sip3d = 4

def checkEleMVA(lep,WP = 'Tight', era = "Spring16" ):
		# Eta dependent MVA ID check:
		passID = False
		
		lepEta = abs(lep.eta)
		
		# eta cut
		if lepEta > eleEta:
			print "here"
			return False
		
		if era == "Spring15":
			lepMVA = lep.mvaSpring16GP
			# numbers here needed to be chaecked if we gonna use MVA ID 
			if WP == 'Tight':
				if lepEta < 0.8: passID = lepMVA > Ele_mvaSpring15_eta0p8_T
				elif lepEta < 1.44: passID = lepMVA > Ele_mvaSpring15_eta1p4_T
				elif lepEta >= 1.57: passID = lepMVA > Ele_mvaSpring15_eta2p4_T
			elif WP == 'Medium':
				if lepEta < 0.8: passID = lepMVA > Ele_mvaSpring15_eta0p8_M
				elif lepEta < 1.44: passID = lepMVA > Ele_mvaSpring15_eta1p4_M
				elif lepEta >= 1.57: passID = lepMVA > Ele_mvaSpring15_eta2p4_M
			elif WP == 'Loose':
				if lepEta < 0.8: passID = lepMVA > Ele_mvaSpring15_eta0p8_L
				elif lepEta < 1.44: passID = lepMVA > Ele_mvaSpring15_eta1p4_L
				elif lepEta >= 1.57: passID = lepMVA > Ele_mvaSpring15_eta2p4_L
			elif WP == 'VLoose':
				if lepEta < 0.8: passID = lepMVA > Ele_mvaSpring15_eta0p8_VL
				elif lepEta < 1.44: passID = lepMVA > Ele_mvaSpring15_eta1p4_VL
				elif lepEta >= 1.57: passID = lepMVA > Ele_mvaSpring15_eta2p4_VL
		
		elif era == "Phys14":
			lepMVA = lep.mvaIdPhys14
		
			if WP == 'Tight':
				if lepEta < 0.8: passID = lepMVA > Ele_mvaPhys14_eta0p8_T
				elif lepEta < 1.44: passID = lepMVA > Ele_mvaPhys14_eta1p4_T
				elif lepEta >= 1.57: passID = lepMVA > Ele_mvaPhys14_eta2p4_T
			elif WP == 'Medium':
				if lepEta < 0.8: passID = lepMVA > Ele_mvaPhys14_eta0p8_M
				elif lepEta < 1.44: passID = lepMVA > Ele_mvaPhys14_eta1p4_M
				elif lepEta >= 1.57: passID = lepMVA > Ele_mvaPhys14_eta2p4_M
			elif WP == 'Loose':
				if lepEta < 0.8: passID = lepMVA > Ele_mvaPhys14_eta0p8_L
				elif lepEta < 1.44: passID = lepMVA > Ele_mvaPhys14_eta1p4_L
				elif lepEta >= 1.57: passID = lepMVA > Ele_mvaPhys14_eta2p4_L
		
		return passID


class lepSFProducer(Module):
	def __init__(self, muonSelectionTag, electronSelectionTag):
		if muonSelectionTag=="LooseWP_2016":
			mu_f=["Mu_Trg.root","Mu_ID.root","Mu_Iso.root"]
			mu_h = ["IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio",
					"MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio",
					"LooseISO_LooseID_pt_eta/pt_abseta_ratio"]
		if electronSelectionTag=="GPMVA90_2016":
			el_f = ["EGM2D_eleGSF.root","EGM2D_eleMVA90.root"]
			el_h = ["EGamma_SF2D", "EGamma_SF2D"]
		if muonSelectionTag=="MediumWP_2016":
			mu_f=["Mu_Trg.root","Mu_ID.root","Mu_Iso.root"]
			mu_h = ["IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio",
					"MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio",
					"LooseISO_MediumID_pt_eta/pt_abseta_ratio"]
		if electronSelectionTag=="Tight_2016":
			el_f = ["EGM2D_eleGSF.root","EGM2D_eleMVA90.root"]
			el_h = ["EGamma_SF2D", "EGamma_SF2D"]
			
		mu_f = ["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in mu_f]
		el_f = ["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in el_f]
	
		self.mu_f = ROOT.std.vector(str)(len(mu_f))
		self.mu_h = ROOT.std.vector(str)(len(mu_f))
		for i in range(len(mu_f)): self.mu_f[i] = mu_f[i]; self.mu_h[i] = mu_h[i];
		self.el_f = ROOT.std.vector(str)(len(el_f))
		self.el_h = ROOT.std.vector(str)(len(el_f))
		for i in range(len(el_f)): self.el_f[i] = el_f[i]; self.el_h[i] = el_h[i];
	
		if "/LeptonEfficiencyCorrector_cc.so" not in ROOT.gSystem.GetLibraries():
			print "Load C++ Worker"
			ROOT.gROOT.ProcessLine(".L %s/src/PhysicsTools/NanoAODTools/python/postprocessing/helpers/LeptonEfficiencyCorrector.cc+" % os.environ['CMSSW_BASE'])
	def beginJob(self):
		self._worker_mu = ROOT.LeptonEfficiencyCorrector(self.mu_f,self.mu_h)
		self._worker_el = ROOT.LeptonEfficiencyCorrector(self.el_f,self.el_h)
	def endJob(self):
		pass
	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		self.out.branch("Muon_effSF", "F")
		self.out.branch("Electron_effSF", "F")
		self.out.branch("lepSF","F")

	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass
	def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
goodLep = []
		Elecs = [x for x in electrons if x.isPFcand and x.pt > 10 and abs(x.eta) < 2.4 and x.cutBased >= 1 and x.miniPFRelIso_all < 0.4]
		Mus = [x for x in muons if x.isPFcand and x.pt > 10 and abs(x.eta) < 2.4 and x.miniPFRelIso_all < 0.4  ]
		goodLep = [i for i in itertools.chain(Mus, Elecs)]
		# Clean good leptons and otherleptons 
		for Mu in goodLep :
			if abs(Mu.pdgId) == 13 :
				for El in goodLep :
					if abs(El.pdgId) == 11 : 
						if Mu.p4().DeltaR(El.p4()) < 0.05 :
							#print "overlap fonded remove Electron"
							goodLep.remove(El)
				
		leps = goodLep 
		nlep = len(leps)		
        # selected good leptons
		selectedTightLeps = []
		for idx,lep in enumerate(leps):
			# for acceptance check
			lepEta = abs(lep.eta)

            # Pt cut
			if lep.pt < 10: continue

            # Iso cut -- to be compatible with the trigger
			if lep.miniPFRelIso_all > trig_miniIsoCut: continue
			###################
			# MUONS
			###################
			if(abs(lep.pdgId) == 13):
				if lepEta > 2.4: continue
			
				#if muID=='ICHEPmediumMuonId': passID = lep.ICHEPmediumMuonId -->> not needed any more 
				passID = lep.mediumId
				passIso = lep.miniPFRelIso_all < muo_miniIsoCut
				passIP = lep.sip3d < goodMu_sip3d
			
				# selected muons
				if passID and passIso and passIP:
					selectedTightLeps.append(lep)			
			###################
			# ELECTRONS
			###################
			
			elif(abs(lep.pdgId) == 11):
			
				if lepEta > eleEta: continue
			
				# pass variables
				passIso = False
				passConv = False
			
				if eleID == 'CB':
					# ELE CutBased ID
					eidCB = lep.cutBased
			
					passTightID = (eidCB == 4 )#and lep.convVeto)
					passMediumID = (eidCB >= 3 )#and lep.convVeto)
					#passLooseID = (eidCB >= 2)
					passVetoID = (eidCB >= 1)# and lep.convVeto)
			
				elif eleID == 'MVA':
					# ELE MVA ID
					# check MVA WPs
					passTightID = checkEleMVA(lep,'Tight')
					passLooseID = checkEleMVA(lep,'VLoose')
				# selected
				if passTightID:
					# Iso check:
					if lep.miniPFRelIso_all < ele_miniIsoCut: passIso = True
					# conversion check
					if eleID == 'MVA':
						if lep.lostHits <= goodEl_lostHits and lep.convVeto and lep.sip3d < goodEl_sip3d: passConv = True
					elif eleID == 'CB':
						passConv = True # cuts already included in POG_Cuts_ID_SPRING15_25ns_v1_ConvVetoDxyDz_X

					passPostICHEPHLTHOverE = True # comment out again if (lep.hOverE < 0.04 and abs(lep.eta)>1.479) or abs(lep.eta)<=1.479 else False
			
					# fill
					if passIso and passConv and passPostICHEPHLTHOverE:
						selectedTightLeps.append(lep)

		# Fill SFs only for tight leptons 
		if len(selectedTightLeps) > 0:
			tightLeps = selectedTightLeps
			for tlep in tightLeps:
				if abs(tlep.pdgId) == 13:
					sf_mu = self._worker_mu.getSF(tlep.pdgId,tlep.pt,tlep.eta)
					self.out.fillBranch("lepSF", sf_mu)
					self.out.fillBranch("Muon_effSF", sf_mu)
				elif abs(tlep.pdgId) == 11:
					sf_el = self._worker_el.getSF(tlep.pdgId,tlep.pt,tlep.eta)
					self.out.fillBranch("lepSF", sf_el)
					self.out.fillBranch("Electron_effSF", sf_el)
				else :
					self.out.fillBranch("Muon_effSF", 1.0)
					self.out.fillBranch("Electron_effSF", 1.0)
					self.out.fillBranch("lepSF", 1.0)
							
		return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

susy_lepSF = lambda : lepSFProducer( "MediumWP_2016", "Tight_2016")

