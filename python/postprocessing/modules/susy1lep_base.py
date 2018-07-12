import ROOT
import math, os
import array
import operator

ROOT.PyConfig.IgnoreCommandLineOptions = True
# for met object 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetSmearer import jetSmearer
from RunIISummer16_Moriond17 import getXsec
from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple

import itertools
from ROOT import TLorentzVector, TVector2, std
mt2obj = ROOT.heppy.Davismt2.Davismt2()
from deltar import bestMatch

#################
### Cuts and WP
#################

## Eta requirement
centralEta = 2.4
eleEta = 2.4

###########
# Jets
###########
smearJER = True

btag_LooseWP = 0.5426
btag_MediumWP = 0.8484
btag_TightWP = 0.9535

# DeepCSV (new Deep Flavour tagger)
btag_DeepLooseWP = 0.2219
btag_DeepMediumWP = 0.6324
btag_DeepTightWP = 0.8958

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


#  print cleanedValueList, min(cleanedValueList)#d, key=d.get)
class susysinglelep(Module):
	def __init__(self, isMC , isSig):#, HTFilter, LTFilter):#, muonSelection, electronSelection):
		self.isMC = isMC
		self.isSig = isSig		
		#self.HTFilt = HTFilter
		#self.LTFilt = LTFilter
		# smear jet pT to account for measured difference in JER between data and simulation.
		self.jerInputFileName = "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"
		self.jerUncertaintyInputFileName = "Spring16_25nsV10_MC_SF_AK4PFchs.txt"
		self.jetSmearer = jetSmearer("Summer16_23Sep2016V4_MC", "AK4PFchs", self.jerInputFileName, self.jerUncertaintyInputFileName)
		pass
	def beginJob(self):
		self.jetSmearer.beginJob()
		pass
	def endJob(self):
		self.jetSmearer.endJob()
		pass
	def matchLeptons(self, event):
		def plausible(rec,gen):
			if abs(rec.pdgId) == 11 and abs(gen.pdgId) != 11:   return False
			if abs(rec.pdgId) == 13 and abs(gen.pdgId) != 13:   return False
			dr = deltaR(rec.etc,rec.phi,gen.eta,gen.phi)
			if dr < 0.3: return True
			if rec.pt < 10 and abs(rec.pdgId) == 13 and gen.pdgId != rec.pdgId: return False
			if dr < 0.7: return True
			if min(rec.pt,gen.pt)/max(rec.pt,gen.pt) < 0.3: return False
			return True
		leps = event.selectedLeptons
		match = matchObjectCollectionMultiple(leps, 
									event.genleps + event.gentauleps, 
									dRmax = 1.2, presel=lambda plausible : True)
		for lep in leps:
			gen = match[lep]
			lep.mcMatchId  = (gen.sourceId if gen != None else  0)
			lep.mcMatchTau = (gen in event.gentauleps if gen else -99)
			lep.mcLep=gen
			
	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		self.out.branch("isData","I");
		self.out.branch("nLep","I");
		self.out.branch("nVeto","I");
		self.out.branch("nEl","I");
		self.out.branch("nMu","I");
		self.out.branch("nTightLeps", "I");
		self.out.branch("nTightMu", "I");
		self.out.branch("nTightEl", "I");
		self.out.branch("tightLepsIdx","I",10,"nTightLeps");
            #("tightLeps_DescFlag","I",10,"nTightLeps"),
		self.out.branch("Lep_pdgId","F");
		self.out.branch("Lep_pt","F");
		self.out.branch("Lep_eta","F");
		self.out.branch("Lep_phi","F");
		self.out.branch("Lep_Idx","I");
		self.out.branch("Lep_relIso","F");
		self.out.branch("Lep_miniIso","F");
		self.out.branch("Lep_hOverE","F");
		self.out.branch("Selected","I"); # selected (tight) or anti-selected lepton	
            # second leading lepton
		self.out.branch("Lep2_pt","F");
		self.out.branch("Selected2","I");
		self.out.branch("LT","F");
		self.out.branch("ST","F");
		self.out.branch("MT","F");
		self.out.branch("DeltaPhiLepW","F");
		self.out.branch("dPhi","F");
		self.out.branch("Lp","F");
		self.out.branch("GendPhi","F");
		self.out.branch("GenLT","F");
		self.out.branch("GenMET","F");
		 # no HF stuff
		#"METNoHF", "LTNoHF", "dPhiNoHF",
            ## jets	
		self.out.branch("ST1","F");
		self.out.branch("HT1","F");    
		self.out.branch("HT","F");
		#self.out.branch("HTphi","F");
		self.out.branch("nJets","I");
		self.out.branch("nBJet","I");
		self.out.branch("nBJetDeep","I");
		self.out.branch("nJets30","I");
		self.out.branch("Jets30Idx","I",50,"nJets30");
		self.out.branch("nBJets30","I");
		self.out.branch("nJets30Clean","I");
		self.out.branch("nJets40","I");
		self.out.branch("nBJets40","I");
		self.out.branch("htJet30j","F");
		self.out.branch("htJet30ja","F");
		self.out.branch("htJet40j","F");
		self.out.branch("Jet1_pt","F");
		self.out.branch("Jet2_pt","F");           
		 ## top tags
		self.out.branch("nHighPtTopTag","I");
		self.out.branch("nHighPtTopTagPlusTau23","I");
            ## special Vars
		self.out.branch("LSLjetptGT80","F"); # leading + subl. jet pt > 80
		self.out.branch("isSR","I"); # is it Signal or Control region
		self.out.branch("Mll","F"); #di-lepton mass
		#self.out.branch("METfilters","I"); not needed for nanoAOD 
            #Datasets

		self.out.branch("PD_JetHT","O");
		self.out.branch("PD_SingleEle","O");
		self.out.branch("PD_SingleMu","O");
		self.out.branch("PD_MET","O");
		
		self.out.branch("isDPhiSignal","I");
		self.out.branch("RA2_muJetFilter","I");
		self.out.branch("Flag_fastSimCorridorJetCleaning","I");
		self.out.branch("minMWjj","F");
		self.out.branch("minMWjjPt","F");
		self.out.branch("bestMWjj","F");
		self.out.branch("bestMWjjPt","F");
		self.out.branch("bestMTopHad","F");
		self.out.branch("bestMTopHadPt","F");   
       # self.out.branch("Jet_mhtCleaning", "b", lenVar="nJet")
       #Iso_track parameters 
		self.out.branch("iso_had","I");
		self.out.branch("iso_pt","F");
		self.out.branch("iso_MT2","F");
		self.out.branch("iso_Veto","O");
		self.out.branch("nLepGood","F");
		self.out.branch("nLepOther","F");
		self.out.branch("LepGood_Cutbased","I");
		self.out.branch("LepOther_Cutbased","I");
       # Store the Xsec 
		self.out.branch("Xsec",  "F");
		self.xs = getXsec(inputFile.GetName())
		self.filename = inputFile.GetName()
		print inputFile.GetName()
		print self.xs
		
	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass
	def met(self, met, isMC):
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_smeared which should be compared 
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
		if isMC:
			return (met.pt_smeared,met.phi_smeared)
		else:
			return (met.pt,met.phi)
	
	def analyze(self, event):
 		"""process event, return True (go to next module) or False (fail, go to next event)"""
		allelectrons = Collection(event, "Electron")
		allmuons = Collection(event, "Muon")
		Jets = Collection(event, "Jet")
		met = Object(event, "MET")
		genmet = Object(event, "GenMET")
		if self.isMC or self.isSig :
			Gen = Collection(event, "GenPart")
			event.genleps = [l for l in Gen if abs(l.pdgId) == 11 or abs(l.pdgId) == 13]
			GenTau = Collection(event, "GenVisTau")
			event.gentauleps = [l for l in GenTau ]
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
					ele.pt>10 and abs(ele.eta)<2.4 ):#and abs(ele.dxy)<0.5 and abs(ele.dz)<1. and ele.lostHits <=1.0) :
				inclusiveElectrons.append(ele)
		event.inclusiveLeptons = inclusiveMuons + inclusiveElectrons
		
		# make loose leptons (basic selection)
		for mu in inclusiveMuons :
				if (mu.pt > 10 and abs(mu.eta) < 2.4 and mu.miniPFRelIso_all < 0.4 and mu.isPFcand and abs(mu.dxy)<0.05 and abs(mu.dz)<1):
					event.selectedLeptons.append(mu)
					event.selectedMuons.append(mu)
				else:
					event.otherLeptons.append(mu)
		looseMuons = event.selectedLeptons[:]
		for ele in inclusiveElectrons :
			ele.looseIdOnly = ele.cutBased >=1
			if (ele.looseIdOnly and
						ele.pt>10 and abs(ele.eta)<2.4 and ele.miniPFRelIso_all < 0.4 and ele.isPFcand and #ele.convVeto and abs(ele.dxy)<0.05 and abs(ele.dz)<0.5  and ele.lostHits <=1.0 and 
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
		self.out.fillBranch("nLepGood",len(goodLep))
		self.out.fillBranch("nLepOther",len(LepOther))
	
		
		'''Elecs = [x for x in electrons if x.isPFcand and x.pt > 10 and abs(x.eta) < 2.4 and x.cutBased >= 1 and x.miniPFRelIso_all < 0.4]
		Mus = [x for x in muons if x.isPFcand and x.pt > 10 and abs(x.eta) < 2.4 and x.miniPFRelIso_all < 0.4  ]
		Elec_Other = [x for x in electrons if x not in set(Elecs) and abs(x.eta) < 2.4 if x.cutBased>=1]
		Mus_Other = [x for x in muons if  x not in set(Mus) and abs(x.eta) < 2.4 ]
		goodLep = [i for i in itertools.chain(Mus, Elecs)]
		LepOther = [i for i in itertools.chain(Mus_Other, Elec_Other)]
		# Clean good leptons and otherleptons 
		for Mu in goodLep :
			if abs(Mu.pdgId) == 13 :
				for El in goodLep :
					if abs(El.pdgId) == 11 : 
						if Mu.p4().DeltaR(El.p4()) < 0.05 :
							#print "overlap fonded remove Electron"
							goodLep.remove(El)
		for Mu in LepOther :
			if abs(Mu.pdgId) == 13 :
				for El in LepOther :
					if abs(El.pdgId) == 11 : 
						if Mu.p4().DeltaR(El.p4()) < 0.05 :
							#print "overlap fonded remove Electron"
							LepOther.remove(El)'''
				
		#LepOther = goodLep
		leps = goodLep 
		nlep = len(leps)
		# adding the ST HT filter 
		if nlep > 0:
			ST1 = leps[0].pt + met.pt 
			self.out.fillBranch("ST1",ST1)
		HT1 = sum([j.pt for j in Jets if (j.pt > 30 and abs(j.eta)<2.4)])
		self.out.fillBranch("HT1",HT1)
        ### LEPTONS
		Selected = False
		if self.isMC == False and self.isSig == False: self.out.fillBranch("isData",1)
		else : self.out.fillBranch("isData",0)
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
			
				## Lower ID is POG_LOOSE (see cfg)
			
				# ID, IP and Iso check:
				passID = False
				#if muID=='ICHEPmediumMuonId': passID = lep.ICHEPmediumMuonId -->> not needed any more 
				passID = lep.mediumId
				passIso = lep.miniPFRelIso_all < muo_miniIsoCut
				passIP = lep.sip3d < goodMu_sip3d
			
				# selected muons
				if passID and passIso and passIP:
					selectedTightLeps.append(lep); selectedTightLepsIdx.append(idx)
					antiVetoLeps.append(lep);
				else:
					selectedVetoLeps.append(lep)
				# anti-selected muons
				if not passIso:
					antiTightLeps.append(lep); antiTightLepsIdx.append(idx)
				else:
					antiVetoLeps.append(lep);
			
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
			
					passTightID = (eidCB == 4 )# and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1) # and lep.convVeto)
					passMediumID = (eidCB >= 3)# and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1)# and lep.convVeto)
					#passLooseID = (eidCB >= 2)
					passVetoID = (eidCB >= 1)# and abs(lep.dxy) < 0.2 and abs(lep.dz) < 0.5) # and lep.convVeto)
			
				elif eleID == 'MVA':
					# ELE MVA ID
					# check MVA WPs
					passTightID = checkEleMVA(lep,'Tight')
					passLooseID = checkEleMVA(lep,'VLoose')
				# selected
				if passTightID:
			
					# all tight leptons are veto for anti
					antiVetoLeps.append(lep)
			
					# Iso check:
					if lep.miniPFRelIso_all < ele_miniIsoCut: passIso = True
					# conversion check
					if eleID == 'MVA':
						if lep.lostHits <= goodEl_lostHits and lep.convVeto and lep.sip3d < goodEl_sip3d: passConv = True
					elif eleID == 'CB':
						passConv = True #if lep.lostHits <= goodEl_lostHits and lep.convVeto and lep.sip3d < goodEl_sip3d else False  # cuts already included in POG_Cuts_ID_SPRING15_25ns_v1_ConvVetoDxyDz_X

					passPostICHEPHLTHOverE = True #if (lep.hoe < 0.04 and abs(lep.eta)>1.479) or abs(lep.eta)<=1.479 else False
			
					# fill
					if passIso and passConv and passPostICHEPHLTHOverE:
						selectedTightLeps.append(lep); selectedTightLepsIdx.append(idx)
					else:
						selectedVetoLeps.append(lep)
			
				# anti-selected
				elif not passMediumID:#passVetoID:
			
					# all anti leptons are veto for selected
					selectedVetoLeps.append(lep)
			
					# Iso check
					passIso = lep.miniPFRelIso_all < Lep_miniIsoCut # should be true anyway
					# other checks
					passOther = False
					if hasattr(lep,"hoe"):
						passOther = lep.hoe > 0.01
			
					# fill
					if passIso and passOther:
						antiTightLeps.append(lep); antiTightLepsIdx.append(idx)
					else:
						antiVetoLeps.append(lep)
				# Veto leptons
				elif passVetoID:
					# the rest is veto for selected and anti
					selectedVetoLeps.append(lep)
					antiVetoLeps.append(lep)
        # end lepton loop

        ###################
        # EXTRA Loop for lepOther -- for anti-selected leptons
        ###################

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
            #############
            #############
		# retrive and fill branches 
            #############
        # choose common lepton collection: select selected or anti lepton
		if len(selectedTightLeps) > 0:
			tightLeps = selectedTightLeps
			tightLepsIdx = selectedTightLepsIdx
			vetoLeps = selectedVetoLeps
		
			self.out.fillBranch("nTightLeps", len(tightLeps))
			self.out.fillBranch("nTightMu",sum([ abs(lep.pdgId) == 13 for lep in tightLeps]))
			self.out.fillBranch("nTightEl", sum([ abs(lep.pdgId) == 11 for lep in tightLeps]))
			
			self.out.fillBranch("tightLepsIdx", tightLepsIdx)
		
			self.out.fillBranch("Selected", 1)
		
			# Is second leading lepton selected, too?
			if len(selectedTightLeps) > 1:
				self.out.fillBranch("Selected2", 1)
			else:
				self.out.fillBranch("Selected2", 0)
		
		elif len(antiTightLeps) > 0:
			tightLeps = antiTightLeps
			tightLepsIdx = antiTightLepsIdx
		
			vetoLeps = antiVetoLeps
		
			self.out.fillBranch("nTightLeps",0)
			self.out.fillBranch("nTightMu",0)
			self.out.fillBranch("nTightEl", 0)
			self.out.fillBranch("tightLepsIdx", [])
		
			self.out.fillBranch("Selected", -1)
		
		else:
			tightLeps = []
			tightLepsIdx = []
			vetoLeps = []
			self.out.fillBranch("nTightLeps", 0)
			self.out.fillBranch("nTightMu", 0)
			self.out.fillBranch("nTightEl", 0)
		
			self.out.fillBranch("tightLepsIdx",[])
		
			self.out.fillBranch("Selected", 0)
		
		# store Tight and Veto lepton numbers
		self.out.fillBranch("nLep",len(tightLeps))
		self.out.fillBranch("nVeto", len(vetoLeps))
		
		# get number of tight el and mu
		tightEl = [lep for lep in tightLeps if abs(lep.pdgId) == 11]
		tightMu = [lep for lep in tightLeps if abs(lep.pdgId) == 13]
		
		self.out.fillBranch("nEl", len(tightEl))
		self.out.fillBranch("nMu", len(tightMu))
		
		# save leading lepton vars
		if len(tightLeps) > 0:# leading tight lep
			self.out.fillBranch("Lep_Idx", tightLepsIdx[0])
		
			self.out.fillBranch("Lep_pt",tightLeps[0].pt)
			self.out.fillBranch("Lep_eta", tightLeps[0].eta)
			self.out.fillBranch("Lep_phi", tightLeps[0].phi)
			self.out.fillBranch("Lep_pdgId", tightLeps[0].pdgId)
		
			self.out.fillBranch("Lep_relIso", tightLeps[0].pfRelIso03_all)
			self.out.fillBranch("Lep_miniIso", tightLeps[0].miniPFRelIso_all)
			if hasattr(tightLeps[0],"hoe"):
				self.out.fillBranch("Lep_hOverE", tightLeps[0].hoe)
		
		elif len(leps) > 0: # fill it with leading lepton
			self.out.fillBranch("Lep_Idx", 0)
		
			self.out.fillBranch("Lep_pt",leps[0].pt)
			self.out.fillBranch("Lep_eta",leps[0].eta)
			self.out.fillBranch("Lep_phi", leps[0].phi)
			self.out.fillBranch("Lep_pdgId", leps[0].pdgId)
		
			self.out.fillBranch("Lep_relIso",leps[0].pfRelIso03_all)
			self.out.fillBranch("Lep_miniIso", leps[0].miniPFRelIso_all)
			if hasattr(leps[0],"hoe"):
				self.out.fillBranch("Lep_hOverE", leps[0].hoe)
		
		# save second leading lepton vars
		if len(tightLeps) > 1:# 2nd tight lep
			self.out.fillBranch("Lep2_pt", tightLeps[1].pt)
#		print " Nlep : ", len(leps) , " nTightleps " , len(tightLeps), " nVetoLEp ", len(vetoLeps)
#####################################################################
#####################################################################
#################Jets, BJets, METS and filters ######################
#####################################################################
#####################################################################
		########
		### Jets
		########
		metp4 = ROOT.TLorentzVector(0,0,0,0)
		Genmetp4 = ROOT.TLorentzVector(0,0,0,0)
		
		if self.isMC:
			Genmetp4.SetPtEtaPhiM(genmet.pt,0,genmet.phi,0)
		#self.out.fillBranch("MET", metp4.Pt())
		Jets = Collection(event, "Jet")
		jets = [j for j in Jets if j.pt > 20 and abs(j.eta) < 2.4]
		njet = len(jets)
		( met_px, met_py ) = ( met.pt*math.cos(met.phi), met.pt*math.sin(met.phi) )
		( met_px_nom, met_py_nom ) = ( met_px, met_py )
		# match reconstructed jets to generator level ones
		# (needed to evaluate JER scale factors and uncertainties)
		
		'''if self.isMC :
			rho = getattr(event,"fixedGridRhoFastjetAll")
			genJets = Collection(event, "GenJet" )
			pairs = matchObjectCollection(Jets, genJets)
			for jet in jets:
				genJet = pairs[jet]
				if smearJER==True :
					(jet_pt_jerNomVal, jet_pt_jerUpVal, jet_pt_jerDownVal) = self.jetSmearer.getSmearValsPt(jet, genJet, rho)
					jet_pt_nom = jet_pt_jerNomVal * jet.pt
					if jet.pt > 15.:
						jet_cosPhi = math.cos(jet.phi)
						jet_sinPhi = math.sin(jet.phi)
						met_px_nom = met_px_nom - (jet_pt_nom - jet.pt)*jet_cosPhi
						met_py_nom = met_py_nom - (jet_pt_nom - jet.pt)*jet_sinPhi
						met_pt_nom = math.sqrt(met_px_nom**2 + met_py_nom**2)
						met_phi_nom = math.atan2(met_py_nom, met_px_nom)
						met.pt = met_pt_nom
						met.phi = met_phi_nom
					jet.pt = jet_pt_nom'''
		metp4.SetPtEtaPhiM(met.pt,0.,met.phi,0.) # only use met vector to derive transverse quantities)	
		centralJet30 = []; centralJet30idx = []
		centralJet40 = []
		cleanJets25 = []; cleanJets25idx = [] 
		cleanJets = []; cleanJetsidx = [] 
		# fill this flage but defults to 1 and then change it after the proper selection 
		self.out.fillBranch("Flag_fastSimCorridorJetCleaning", 1)
		for i,j in enumerate(jets):
			# Cleaning up of fastsim jets (from "corridor" studies) https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsMoriond17#Cleaning_up_of_fastsim_jets_from
			if self.isSig: #only check for signals (see condition check above)
				self.out.fillBranch("isDPhiSignal",1) 
				genj = pairs[j]
				if genj is not None : 
					if j.pt>20 and abs(j.eta)<2.5 and (genj.pt == 0) and j.chHEF<0.1:	self.out.fillBranch("Flag_fastSimCorridorJetCleaning", 0  )
			if j.pt>25 :
				cleanJets25.append(j)
				cleanJets25idx.append(j)
			if j.pt > 30 and abs(j.eta)<centralEta:
				centralJet30.append(j)
				centralJet30idx.append(i)
			if j.pt>40 and abs(j.eta)<centralEta:
				centralJet40.append(j)
		
		# jets 30 (cmg cleaning only)
		nJetC = len(centralJet30)
		self.out.fillBranch("nJets", nJetC)
		self.out.fillBranch("nJets30",nJetC)
		# store indeces
		self.out.fillBranch("Jets30Idx", centralJet30idx)
		#print "nJets30:", len(centralJet30), " nIdx:", len(centralJet30idx)
		
		# jets 40
		nJet40C = len(centralJet40)
		self.out.fillBranch("nJets40",nJet40C)
		
		##############################
		## Local cleaning from leptons
		##############################
		cJet30Clean = []
		dRminCut = 0.4
		
		# Do cleaning a la CMG: clean max 1 jet for each lepton (the nearest)
		cJet30Clean = centralJet30
		cleanJets30 = centralJet30
		#clean selected leptons at First 
		'''for lep in goodLep:
			if lep.pt < 20 : continue 
			jNear, dRmin = None, 99
			# find nearest jet
			for jet in centralJet30:
				dR = jet.p4().DeltaR(lep.p4())
				if dR < dRmin:
					jNear, dRmin = jet, dR
			# remove nearest jet
			if dRmin < dRminCut:
				cJet30Clean.remove(jNear)'''
		#then clean other tight leptons 
		for lep in tightLeps:
			# don't clean LepGood, only LepOther
			#if lep not in otherleps: continue
			jNear, dRmin = None, 99
			# find nearest jet
			for jet in centralJet30:
				dR = jet.p4().DeltaR(lep.p4())
				if dR < dRmin:
					jNear, dRmin = jet, dR
			# remove nearest jet
			if dRmin < dRminCut:
				cJet30Clean.remove(jNear)
			for ijet,jet25 in enumerate(cleanJets25): 
				dR = jet25.p4().DeltaR(lep.p4())
				if dR < dRmin:
					cleanJets.append(jet25)
					cleanJetsidx.append(ijet)
		#if nJetC !=  len(cJet30Clean) and False:
		#	print "Non-clean jets: ", nJetC, "\tclean jets:", len(cJet30Clean)
		#	print jets
		#	print leps
		#	print otherleps
		
		# cleaned jets
		nJet30C = len(cJet30Clean)
		self.out.fillBranch("nJets30Clean",len(cJet30Clean))
		
		if nJet30C > 0:
			self.out.fillBranch("Jet1_pt", cJet30Clean[0].pt)
		if nJet30C > 1:
			self.out.fillBranch("Jet2_pt", cJet30Clean[1].pt)
		
		# imho, use Jet2_pt > 80 instead
		self.out.fillBranch("LSLjetptGT80", 1 if sum([j.pt>80 for j in cJet30Clean])>=2 else 0)
		
		self.out.fillBranch("htJet30j", sum([j.pt for j in cJet30Clean]))
		self.out.fillBranch("htJet30ja", sum([j.pt for j in jets if j.pt>30]))
		
		self.out.fillBranch("htJet40j", sum([j.pt for j in centralJet40]))
		
		self.out.fillBranch("HT", sum([j.pt for j in cJet30Clean]))
		
		## B tagging WPs for CSVv2 (CSV-IVF)
		## from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging#Preliminary_working_or_operating
		
		# WP defined on top
		btagWP = btag_MediumWP
		
		BJetMedium30 = []
		BJetMedium40 = []
		
		nBJetDeep = 0
		
		for i,j in enumerate(cJet30Clean):
			if j.btagCSVV2 > btagWP:
				BJetMedium30.append(j)
			if (j.btagDeepB) > btag_DeepMediumWP:
				nBJetDeep += 1
		
		for i,j in enumerate(centralJet40):
			if j.btagCSVV2 > btagWP:
				BJetMedium40.append(j)
		
		# using cleaned collection!
		self.out.fillBranch("nBJet", len(BJetMedium30))
		self.out.fillBranch("nBJets30", len(BJetMedium30))
		
		self.out.fillBranch("nBJetDeep", nBJetDeep)
		
		# using normal collection
		self.out.fillBranch("nBJets40", len(BJetMedium40))
		
		######
		# MET
		#####
		

#####################################################################
#####################################################################
#################High level variables ###############################
#####################################################################
#####################################################################
		dPhiLepW = -999 # set default value to -999 to spot "empty" entries
		GendPhiLepW = -999 # set default value to -999 to spot "empty" entries
		# LT of lepton and MET
		LT = -999
		GenLT = -999
		Lp = -99
		MT = -99
		
		if len(tightLeps) >=1:
			recoWp4 =  tightLeps[0].p4() + metp4
			GenrecoWp4 =  tightLeps[0].p4() + Genmetp4
			GendPhiLepW = tightLeps[0].p4().DeltaPhi(GenrecoWp4)
			GenLT = tightLeps[0].pt + Genmetp4.Pt()
			dPhiLepW = tightLeps[0].p4().DeltaPhi(recoWp4)
			LT = tightLeps[0].pt + metp4.Pt()
			Lp = tightLeps[0].pt / recoWp4.Pt() * math.cos(dPhiLepW)
		
			#MT = recoWp4.Mt() # doesn't work
			MT = math.sqrt(2*metp4.Pt()*tightLeps[0].pt * (1-math.cos(dPhiLepW)))
		self.out.fillBranch("DeltaPhiLepW", dPhiLepW)
		dPhi = abs(dPhiLepW) # nickname for absolute dPhiLepW
		self.out.fillBranch("dPhi", dPhi)
		self.out.fillBranch("ST", LT)
		self.out.fillBranch("LT", LT)
		self.out.fillBranch("Lp", Lp)
		self.out.fillBranch("MT", MT)
		self.out.fillBranch("GendPhi", abs(GendPhiLepW))
		self.out.fillBranch("GenLT", GenLT)
		self.out.fillBranch("GenMET", Genmetp4.Pt())
		
		#####################
		## SIGNAL REGION FLAG
		#####################
		
		## Signal region flag
		# isSR SR vs CR flag
		isSR = 0
		
		# 0-B SRs -- simplified dPhi
		if len(BJetMedium30) == 0:# check the no. of Bjets 
			if LT < 250:   isSR = 0
			elif LT > 250: isSR = dPhi > 0.75
			# BLIND data
			if (not self.isMC)  and nJet30C >= 5:
				isSR = - isSR
		# Multi-B SRs
		elif nJet30C < 99:
			if LT < 250:   isSR = 0
			elif LT < 350: isSR = dPhi > 1.0
			elif LT < 600: isSR = dPhi > 0.75
			elif LT > 600: isSR = dPhi > 0.5
		
			# BLIND data
			if (not self.isMC) and nJet30C >= 6:
				isSR = - isSR
		
		self.out.fillBranch("isSR", isSR)
		
		#############
		## Playground
		#############
		
		# di-lepton mass: opposite-sign, same flavour
		Mll = 0
		
		if len(tightLeps) > 1:
		
			lep1 = tightLeps[0]
			id1 = lep1.pdgId
		
			for lep2 in leps[1:]:
				if lep2.pdgId + lep1.pdgId == 0:
					dilepP4 = lep1.p4() + lep2.p4()
					Mll = dilepP4.M()
		
		self.out.fillBranch("Mll", Mll)
		
		# RA2 proposed filter
		self.out.fillBranch("RA2_muJetFilter", True) # normally true for now # don't know how to get the Muon energy fraction from EMEF
		#for j in cJet30Clean:
		#	if j.pt > 200 and j.chEmEF > 0.5 and abs(math.acos(math.cos(j.phi-metp4.Phi()))) > (math.pi - 0.4):
		#		self.out.fillBranch("RA2_muJetFilter", False)
		

		## MET FILTERS for data looks like the met filters are applied already for nanoAOD 
		#####################
		## Top Tagging          ------->>>>>>. to be moved to different modules keep it commented here 
		#####################		
		lightJets = [ j for j in cleanJets if not j.btagCSVV2 == btagWP ]
		bjetsLoose  = [ j for j in cleanJets if j.btagCSVV2== btag_LooseWP]
		minMWjj   = 999
		minMWjjPt = 0
		bestMWjj   = 0
		bestMWjjPt = 0
		bestMTopHad   = 0
		bestMTopHadPt = 0
		for i1,j1 in enumerate(lightJets):
			for i2 in xrange(i1+1,len(lightJets)):
				j2 = lightJets[i2]
				jjp4 = j1.p4() + j2.p4()
				mjj  = jjp4.M()
				if mjj > 30 and mjj < minMWjj:
					minMWjj = mjj
					minMWjjPt = jjp4.Pt()
					self.out.fillBranch("minMWjj",minMWjj)
					self.out.fillBranch("minMWjjPt",minMWjjPt)
				if abs(mjj-80.4) < abs(bestMWjj-80.4):
					bestMWjj = mjj
					bestMWjjPt = jjp4.Pt()
					self.out.fillBranch("bestMWjj",bestMWjj)
					self.out.fillBranch("bestMWjjPt",bestMWjjPt)
					for bj in bjetsLoose:
						if deltaR(bj.eta(),bj.phi(),j1.eta(),j1.phi()) < 0.1 or deltaR(bj.eta(),bj.phi(),j2.eta(),j2.phi()) < 0.1: continue
						tp4 = jjp4 + bj.p4()
						mtop = tp4.M()
						if abs(mtop-172) < abs(bestMTopHad - 172):
							bestMTopHad = mtop
							bestMTopHadPt = tp4.Pt()
							self.out.fillBranch("bestMTopHad",bestMTopHad)
							self.out.fillBranch("bestMTopHadPt",bestMTopHadPt)
#####################################################################
#####################################################################
#################Fill MT2 here, not point to make a separate module##
#####################################################################
#####################################################################
# isolated tracks after basic selection (((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < 2.5) && abs(dxy) < 0.2 && abs(dz) < 0.1 && ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)) and lepton veto
		# First check is the event has the IsoTrack or not 
		if hasattr(event, "nIsoTrack"):
			trks = [j for j in Collection(event,"IsoTrack","nIsoTrack")]
			tp4 = ROOT.TLorentzVector(0,0,0,0)
			
			# min dR between good lep and iso track 
			minDR = 0.1
			# MT2 cuts for hadronic and leptonic veto tracks
			hadMT2cut = 60
			lepMT2cut = 80
			if (len(tightLeps)>=1) and len(trks)>=1:
				for i,t in enumerate(trks):
					# looking for opposite charged tracks
					#if tightLeps[0].charge == t.charge: continue # not track charge is founded replace with the next copule of lines 
					if t.isHighPurityTrack == False : continue 
					#print t.miniPFRelIso_chg
					# not track mass is founded 
					tp4.SetPtEtaPhiM(t.pt,t.eta,t.phi,0.)
					dR = tp4.DeltaR(tightLeps[0].p4())
					if minDR>dR: continue
					p1=tightLeps[0].p4()
					p2=tp4
					a=array.array('d', [ p1.M(), p1.Px(), p1.Py() ])                    
					b=array.array('d', [ p2.M(), p2.Px(), p2.Py() ])                    
					c=array.array('d', [ metp4.M(), metp4.Px(), metp4.Py() ])                    
					mt2obj.set_momenta( a, b, c )
					mt2obj.set_mn(0)
					self.out.fillBranch("iso_MT2",mt2obj.get_mt2())
					self.out.fillBranch("iso_pt", p2.Pt())
					# cuts on MT2 as defined above
					if abs(t.pdgId)>10 and abs(t.pdgId)<14:
						self.out.fillBranch("iso_had", 0)  #leptonic
						cut=lepMT2cut
					else: 
						self.out.fillBranch("iso_had", 1)  #hadronic track
						cut=hadMT2cut
					if mt2obj.get_mt2()<=cut: self.out.fillBranch("iso_Veto",True)
			
		self.out.fillBranch("Xsec",self.xs)
		if 'JetHT' in self.filename:
			self.out.fillBranch("PD_JetHT",True)
		else: self.out.fillBranch("PD_JetHT",False)
		if 'SingleEle' in self.filename:
			self.out.fillBranch("PD_SingleEle",True)
		else : self.out.fillBranch("PD_SingleEle",False)
		if 'SingleMu' in self.filename:
			self.out.fillBranch("PD_SingleMu",True)
		else : self.out.fillBranch("PD_SingleMu",False)
		if 'MET' in self.filename:
			self.out.fillBranch("PD_MET",True)
		else : self.out.fillBranch("PD_MET",False)
					
		
		return True
susy1lepSIG = lambda : susysinglelep(True,True)#,
susy1lepMC = lambda : susysinglelep(True,False)#,
susy1lepdata = lambda : susysinglelep(False ,False) 
