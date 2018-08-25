import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import itertools
from deltar import bestMatch


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
	def __init__(self, muonSelectionTag, electronSelectionTag,era):
		if era == "2016":
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
				
		############################# 2017 lepSF
		elif era == "2017" : 
			if muonSelectionTag=="LooseWP_2017":
				mu_f=["Mu_Trg17.root","Mu_ID17.root","Mu_Iso17.root"]
				mu_h = ["IsoMu27_PtEtaBins/pt_abseta_ratio",
						"NUM_LooseID_DEN_genTracks_pt_abseta",
						"NUM_LooseRelIso_DEN_LooseID_pt_abseta"]
			if electronSelectionTag=="GPMVA90_2017":
				el_f = ["EGM2D_eleGSF17.root","EGM2D_eleMVA90_17.root"]
				el_h = ["EGamma_SF2D", "EGamma_SF2D"]
			if muonSelectionTag=="MediumWP_2017":
				mu_f=["Mu_Trg17.root","Mu_ID17.root","Mu_Iso17.root"]
				mu_h = ["IsoMu27_PtEtaBins/pt_abseta_ratio",
						"NUM_MediumID_DEN_genTracks_pt_abseta",
						"NUM_LooseRelIso_DEN_MediumID_pt_abseta"]
			if electronSelectionTag=="Tight_2017":
				el_f = ["EGM2D_eleGSF17.root","EGM2D_eleMVA90_17.root"]
				el_h = ["EGamma_SF2D", "EGamma_SF2D"]
			
			else : 
				raise ValueError("ERROR: Invalid era = '%s'!" % self.era)
				
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
		allelectrons = Collection(event, "Electron")
		allmuons = Collection(event, "Muon")
		# for all leptons (veto or tight)
		
		### inclusive leptons = all leptons that could be considered somewhere in the analysis, with minimal requirements (used e.g. to match to MC)
		event.inclusiveLeptons = []
		### selected leptons = subset of inclusive leptons passing some basic id definition and pt requirement
		### other    leptons = subset of inclusive leptons failing some basic id definition and pt requirement
		event.selectedLeptons = []
		event.selectedMuons = []
		event.selectedElectrons = []
		event.otherLeptons = []
		inclusiveMuons = []
		inclusiveElectrons = []
		for mu in allmuons:
			if (mu.pt>10 and abs(mu.eta)<2.4 and
					abs(mu.dxy)<0.5 and abs(mu.dz)<1.):
				inclusiveMuons.append(mu)
		for ele in allelectrons:
			if ( ele.cutBased >=1 and
					ele.pt>10 and abs(ele.eta)<2.4):# and abs(ele.dxy)<0.5 and abs(ele.dz)<1. and ele.lostHits <=1.0):
				inclusiveElectrons.append(ele)
		event.inclusiveLeptons = inclusiveMuons + inclusiveElectrons
		
		# make loose leptons (basic selection)
		for mu in inclusiveMuons :
				if (mu.pt > 10 and abs(mu.eta) < 2.4 and mu.miniPFRelIso_all < 0.4 and mu.isPFcand and abs(mu.dxy)<0.05 and abs(mu.dz)<0.5):
					event.selectedLeptons.append(mu)
					event.selectedMuons.append(mu)
				else:
					event.otherLeptons.append(mu)
		looseMuons = event.selectedLeptons[:]
		for ele in inclusiveElectrons :
			ele.looseIdOnly = ele.cutBased >=1
			if (ele.looseIdOnly and
						ele.pt>10 and abs(ele.eta)<2.4 and ele.miniPFRelIso_all < 0.4 and ele.isPFcand and ele.convVeto and # and abs(ele.dxy)<0.05 and abs(ele.dz)<0.5  and ele.lostHits <=1.0 and 
						(bestMatch(ele, looseMuons)[1] > (0.05**2))):
					event.selectedLeptons.append(ele)
					event.selectedElectrons.append(ele)
			else:
					event.otherLeptons.append(ele)		
		event.otherLeptons.sort(key = lambda l : l.pt, reverse = True)
		event.selectedLeptons.sort(key = lambda l : l.pt, reverse = True)
		event.selectedMuons.sort(key = lambda l : l.pt, reverse = True)
		event.selectedElectrons.sort(key = lambda l : l.pt, reverse = True)
		event.inclusiveLeptons.sort(key = lambda l : l.pt, reverse = True)		
		
		goodLep = [l for l in event.selectedLeptons]
		LepOther = [l for l in event.otherLeptons]				
		leps = goodLep 
		nlep = len(leps)		
		# selected good leptons
		selectedTightLeps = []
		selectedTightLepsIdx = []
		selectedVetoLeps = []

        # anti-selected leptons
		antiTightLeps = []
		antiTightLepsIdx = []
		antiVetoLeps = []
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
		otherleps = [l for l in LepOther]#(set(selectedTightLeps) or set(selectedVetoLeps) or set(antiTightLeps) or set(antiVetoLeps) )]
		nLepOther = len(otherleps)
		for idx,lep in enumerate(otherleps):
			# check acceptance
			lepEta = abs(lep.eta)
			if lepEta > 2.4: continue
			# Pt cut
			if lep.pt < 10: continue
			# Iso cut -- to be compatible with the trigger
			if lep.miniPFRelIso_all > trig_miniIsoCut: continue
		
			############
			# Muons
			if(abs(lep.pdgId) == 13):
				## Lower ID is POG_LOOSE (see cfg)
		
				# ID, IP and Iso check:
				passIso = lep.miniPFRelIso_all > muo_miniIsoCut
				#if passIso and passID and passIP:
				if passIso:
					antiTightLeps.append(lep)
					antiTightLepsIdx.append(idx)
				else:
					antiVetoLeps.append(lep)
		
			############
			# Electrons
			elif(abs(lep.pdgId) == 11):
		
				if(lepEta > eleEta): continue
		
				## Iso selection: ele should have MiniIso < 0.4 (for trigger)
				if lep.miniPFRelIso_all > Lep_miniIsoCut: continue
		
				## Set Ele IDs
				if eleID == 'CB':
					# ELE CutBased ID
					eidCB = lep.cutBased
		
					passMediumID = (eidCB >= 3)# and lep.convVeto and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1)
					passVetoID = (eidCB >= 1 )#and lep.convVeto and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1)
				else:
					passMediumID = False
					passVetoID = False
		
				# Cuts for Anti-selected electrons
				if not passMediumID:
					# should always be true for LepOther
		
					# other checks
					passOther = False
					if hasattr(lep,"hoe"):
						passOther = lep.hoe > 0.01
		
					#if not lep.conVeto:
					if passOther:
						antiTightLeps.append(lep)
						antiTightLepsIdx.append(idx);
					else:
						antiVetoLeps.append(lep)
		
				elif passVetoID: #all Medium+ eles in LepOther
					antiVetoLeps.append(lep)

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

susy_lepSF16 = lambda : lepSFProducer( 'MediumWP_2016', 'Tight_2016',"2016")
susy_lepSF17 = lambda : lepSFProducer( 'MediumWP_2017', 'Tight_2017',"2017")

