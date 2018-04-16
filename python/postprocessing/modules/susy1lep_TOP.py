import ROOT
import math, os
import array
import operator

ROOT.PyConfig.IgnoreCommandLineOptions = True
# for met object 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetSmearer import jetSmearer

from ROOT import TLorentzVector, TVector2, std

#################
### Cuts and WP
#################

## Eta requirement
centralEta = 2.4
eleEta = 2.4

###########
# Jets
###########

corrJEC = "central" # can be "central","up","down"
JECAllowedValues = ["central","up","down"]
assert any(val==corrJEC for val in JECAllowedValues)

smearJER = "None"# can be "None","central","up","down"
JERAllowedValues = ["None","central","up","down"]
assert any(val==smearJER for val in JERAllowedValues)

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


def getPhysObjectArray(j): # https://github.com/HephySusySW/Workspace/blob/72X-master/RA4Analysis/python/mt2w.py
    px = j.pt*math.cos(j.phi )
    py = j.pt*math.sin(j.phi )
    pz = j.pt*math.sinh(j.eta )
    E = math.sqrt(px*px+py*py+pz*pz) #assuming massless particles...
    return array.array('d', [E, px, py,pz])

def mt_2(p4one, p4two):
    return math.sqrt(2*p4one.Pt()*p4two.Pt()*(1-math.cos(p4one.Phi()-p4two.Phi())))

def GetZfromM(vector1,vector2,mass):
    MT = math.sqrt(2*vector1.Pt()*vector2.Pt()*(1-math.cos(vector2.DeltaPhi(vector1))))
    if (MT<mass):
        Met2D = TVector2(vector2.Px(),vector2.Py())
        Lep2D = TVector2(vector1.Px(),vector1.Py())
        A = mass*mass/2.+Met2D*Lep2D
        Delta = vector1.E()*vector1.E()*(A*A-Met2D.Mod2()*Lep2D.Mod2())
        MetZ1 = (A*vector1.Pz()+math.sqrt(Delta))/Lep2D.Mod2()
        MetZ2 = (A*vector1.Pz()-math.sqrt(Delta))/Lep2D.Mod2()
    else:
        MetZ1 =vector1.Pz()*vector2.Pt()/vector1.Pt()
        MetZ2 =vector1.Pz()*vector2.Pt()/vector1.Pt()
    return [MT,MetZ1,MetZ2]

def minValueForIdxList(values,idxlist):
    cleanedValueList = [val for i,val in enumerate(values) if i in idxlist]
    if len(cleanedValueList)>0: return min(cleanedValueList)
    else: return -999
#  print cleanedValueList, min(cleanedValueList)#d, key=d.get)


class susyTOP(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree

		self.out.branch("minDPhiBMET","F")
		self.out.branch("minDPhiJMET","F")
		self.out.branch("idxMinDPhiBMET","I")
		self.out.branch("mTClBPlusMET","F")
		self.out.branch("mTBJetMET","F")
		self.out.branch("mTLepMET","F")
		self.out.branch("mLepBJet","F")
		self.out.branch("LepBMass","F",50,"nCentralJet30")
		self.out.branch("MTbnu","F",50,"nCentralJet30")
		self.out.branch("Mtop","F",50,"nCentralJet30")
		self.out.branch("MTtop","F",50,"nCentralJet30"),
		self.out.branch("METovTop","F",50,"nCentralJet30")
		self.out.branch("METTopPhi","F",50,"nCentralJet30")
		self.out.branch("MtopDecor","F",50,"nCentralJet30")
		self.out.branch("TopPt","F",50,"nCentralJet30")
		self.out.branch("TopEt","F",50,"nCentralJet30")
		self.out.branch("nBMinVariantsTopVars","I")
		self.out.branch("TopVarsMTbnuMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("TopVarsLepBMassMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("TopVarsMTtopMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("TopVarsMtopMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("TopVarsMETovTopMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("TopVarsMtopDecorMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("TopVarsTopPtMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("TopVarsTopEtMin","F",10,"nBMinVariantsTopVars")
		self.out.branch("MTW","F")
		self.out.branch("MW1","F")
		self.out.branch("MW2","F")
		self.out.branch("minDphiLepB","F")
		self.out.branch("minDphiLepBidx","F")
                
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
		electrons = Collection(event, "Electron")
		muons = Collection(event, "Muon")
		Jets = Collection(event, "Jet")
		met = Object(event, "MET")
		genmet = Object(event, "GenMET")
		# for all leptons (veto or tight)
		Elecs = [x for x in electrons if x.eta < 2.4 and x.cutBased == 4 and x.convVeto == True and x.pt > 10 and x.miniPFRelIso_all < 0.1]      
		Mus = [x for x in muons if x.eta < 2.4 and  x.pt > 10 and x.mediumId >= 1 and x.miniPFRelIso_all < 0.2 and x.sip3d < 4 ]
		goodLep = Elecs + Mus 
		
		leps = [l for l in goodLep]
		nlep = len(leps)
        ### LEPTONS
		
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
			
					passTightID = (eidCB == 4 and lep.convVeto)
					passMediumID = (eidCB >= 3 and lep.convVeto)
					#passLooseID = (eidCB >= 2)
					passVetoID = (eidCB >= 1 and lep.convVeto)
			
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
						passConv = True # cuts already included in POG_Cuts_ID_SPRING15_25ns_v1_ConvVetoDxyDz_X

					passPostICHEPHLTHOverE = True # comment out again if (lep.hOverE < 0.04 and abs(lep.eta)>1.479) or abs(lep.eta)<=1.479 else False
			
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

		otherleps = [l for l in goodLep]
        #otherleps = []

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
				#passID = lep.mediumMuonId == 1
				passIso = lep.miniPFRelIso_all > muo_miniIsoCut
				# cuts like for the LepGood muons
				#passIP = abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1
		
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
		
					passMediumID = (eidCB >= 3 and lep.convVeto)
					passVetoID = (eidCB >= 1 and lep.convVeto)
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
		
		# choose common lepton collection: select selected or anti lepton
		if len(selectedTightLeps) > 0:
			tightLeps = selectedTightLeps
			tightLepsIdx = selectedTightLepsIdx
			vetoLeps = selectedVetoLeps
            
            #############
            #############
		# retrive and fill branches 
            #############
        # choose common lepton collection: select selected or anti lepton
		if len(selectedTightLeps) > 0:
			tightLeps = selectedTightLeps
			tightLepsIdx = selectedTightLepsIdx
		
			vetoLeps = selectedVetoLeps
		
		elif len(antiTightLeps) > 0:
			tightLeps = antiTightLeps
			tightLepsIdx = antiTightLepsIdx
		
			vetoLeps = antiVetoLeps
		
		else:
			tightLeps = []
			tightLepsIdx = []
		
			vetoLeps = []
		nTightLeps = len(tightLeps) 
		# Met 
		metp4 = ROOT.TLorentzVector(0,0,0,0)
		if hasattr(event, 'metMuEGClean_pt'):
			metp4.SetPtEtaPhiM(event.metMuEGClean_pt,event.metMuEGClean_eta,event.metMuEGClean_phi,event.metMuEGClean_mass)
		else:
			metp4.SetPtEtaPhiM(met.pt,0,met.phi,0)
		if hasattr(event, 'metMuEGClean_pt'):
			pmiss  =array.array('d',[event.metMuEGClean_pt * math.cos(event.metMuEGClean_phi), event.metMuEGClean_pt * math.sin(event.metMuEGClean_phi)] )
		else:
			pmiss  =array.array('d',[met.pt * math.cos(met.phi), met.pt * math.sin(met.phi)] )

#####################################################################
#####################################################################
#################Jets, BJets, METS and filters ######################
#####################################################################
#####################################################################
		########
		### Jets
		########
		jets = [j for j in Jets ]
		njet = len(jets)
		# it's not needed for nanoAOD there is a module to do the job for you 
		# Apply JEC up/down variations if needed (only MC!)
		'''if self.isMC == True:
			if corrJEC == "central":
				pass # don't do anything
				#for jet in jets: jet.pt = jet.rawPt * jet.corr
			elif corrJEC == "up":
				for jet in jets: jet.pt = jet.rawPt * jet.corr_JECUp
			elif corrJEC == "down":
				for jet in jets: jet.pt = jet.rawPt * jet.corr_JECDown
			else:
				pass
			if smearJER!= "None":
				for jet in jets: jet.pt = returnJERSmearedPt(jet.pt,abs(jet.eta),jet.mcPt,smearJER)'''
		
		centralJet30 = []; centralJet30idx = []
		centralJet40 = []
		cleanJets25 = []; cleanJets25idx = [] 
		cleanJets = []; cleanJetsidx = [] 
		# fill this flage but defults to 1 and then change it after the proper selection 
		for i,j in enumerate(jets):
			if j.pt>25 :
				cleanJets25.append(j)
				cleanJets25idx.append(j)
			if j.pt > 30 and abs(j.eta)<centralEta:
				centralJet30.append(j)
				centralJet30idx.append(i)
			if j.pt>40 and abs(j.eta)<centralEta:
				centralJet40.append(j)
		
		# jets 30 (cmg cleaning only)
		nCentralJet30 = len(centralJet30)
		
		btagWP = btag_MediumWP
		
		BJetMedium30 = []
		BJetMedium30idx = []
		
		BJetMedium40 = []
		
		nBJetDeep = 0
		
		cJet30Clean = []
		dRminCut = 0.4
		
		# Do cleaning a la CMG: clean max 1 jet for each lepton (the nearest)
		cJet30Clean = centralJet30
		cleanJets30 = centralJet30

		for i,j in enumerate(cJet30Clean):
			if j.btagCSVV2 > btagWP:
				BJetMedium30.append(j)
				BJetMedium30idx.append(j)
			if (j.btagDeepB) > btag_DeepMediumWP:
				nBJetDeep += 1
		
		for i,j in enumerate(centralJet40):
			if j.btagCSVV2 > btagWP:
				BJetMedium40.append(j)
		
		nBJetMedium30 = len(BJetMedium30)
		##################################################################
		# The following variables need to be double-checked for validity #
		##################################################################
		
		## B tagging WPs for CSVv2 (CSV-IVF)
		## L: 0.423, M: 0.814, T: 0.941
		## from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging#Preliminary_working_or_operating
		
		bTagWP = 0.814 # MediumWP for CSVv2
		#bTagWP = 0.732 # MediumWP for CMVA
		
		# min deltaPhi between a b-jet and MET; needs to be double-checked
		minDPhiBMET    = 100
		idxMinDPhiBMET = -999
		for i, jet in enumerate(centralJet30):
			if jet.btagCSVV2>bTagWP:
				dPhiBMET = abs(jet.p4().DeltaPhi(metp4))
				if dPhiBMET<minDPhiBMET:
					minDPhiBMET=dPhiBMET
					idxMinDPhiBMET = i
		
		self.out.fillBranch("idxMinDPhiBMET", idxMinDPhiBMET)
		self.out.fillBranch("minDPhiBMET", minDPhiBMET)
		
		# min deltaPhi between a jet (first three jets) and MET; needs to be double-checked
		minDPhiJMET    = 100
		for i, jet in enumerate(centralJet30[:3]):
			dPhiJMET = abs(jet.p4().DeltaPhi(metp4))
			if dPhiJMET<minDPhiJMET:
				minDPhiJMET=dPhiJMET
		
		self.out.fillBranch("minDPhiJMET",minDPhiJMET)
		
		# transverse mass of (closest (to MET) BJet, MET), (closest (to MET) BJet, lepton),
		# mass of (closest (to MET) BJet, lepton); need to be double-checked
		mTBJetMET      = -999
		mTLepMET       = -999
		mLepBJet       = -999
		if(idxMinDPhiBMET>=0):
			SumMetClosestBJet = jets[idxMinDPhiBMET].p4() + metp4
			self.out.fillBranch("mTClBPlusMET",SumMetClosestBJet.Mt())
			mTBJetMET = mt_2(centralJet30[idxMinDPhiBMET].p4(),metp4)
			if nTightLeps>=1:
				mLepBJet = (centralJet30[idxMinDPhiBMET].p4() + tightLeps[0].p4()).M()
				mTLepMET = mt_2(tightLeps[0].p4(),metp4)
		else:
			self.out.fillBranch("mTClBPlusMET", -999)
		
		self.out.fillBranch("mTBJetMET", mTBJetMET)
		self.out.fillBranch("mTLepMET", mTLepMET)
		self.out.fillBranch("mLepBJet", mLepBJet)
		
		# projection of MET along (MET + lepton + (closest (to MET) BJet)) sum; needs to be double-checked...
		MetZ1 = -9999
		MetZ2 = -9999
		MTW = -9999
		MW1 = -9999
		MW2 = -9999
		neutrino1 = ROOT.TLorentzVector(0,0,0,0)
		neutrino2 = ROOT.TLorentzVector(0,0,0,0)
		if(nTightLeps==1) :
			NeutrZList = GetZfromM(tightLeps[0].p4(),metp4,81)
			MTW = NeutrZList[0]
			MetZ1= NeutrZList[1]
			MetZ2= NeutrZList[2]
			neutrino1.SetXYZM(metp4.Px(),metp4.Py(), MetZ1, 0)
			neutrino2.SetXYZM(metp4.Px(),metp4.Py(), MetZ2, 0)
			MW1 = (neutrino1+tightLeps[0].p4()).M()
			MW2 = (neutrino2+tightLeps[0].p4()).M()
		self.out.fillBranch("MTW", MTW)
		self.out.fillBranch("MW1", MW1)
		self.out.fillBranch("MW2", MW2)
		# some extra plots
		
		MTbnu = []
		LepBMass = []
		MTtop = []
		Mtop = []
		METovTop = []
		METTopPhi = []
		MtopDecor = []
		
		TopEt = []
		TopPt = []
		
		if(nTightLeps==1) :
			for i,jet in  enumerate(centralJet30): #testing all jets as b-jet in top-reco
				if hasattr(event, 'metMuEGClean_pt'):
					ThisMTnub = math.sqrt(2*event.metMuEGClean_pt*jet.pt* (1-math.cos( metp4.DeltaPhi(jet.p4() ))))
				else:
					ThisMTnub = math.sqrt(2*met.pt*jet.pt* (1-math.cos( metp4.DeltaPhi(jet.p4() ))))
				MTbnu.append(ThisMTnub)
		
				# lep + jet vector
				lepJp4 = tightLeps[0].p4()+jet.p4()
		
				ThislepBMass = lepJp4.M()
				LepBMass.append(ThislepBMass )
		
				# top vector: MET + lep + jet
				topP4 = metp4+lepJp4
				TopEt.append(topP4.Et())
				TopPt.append(topP4.Pt())
		
				ThisMTtop =  math.sqrt( 81.*81. + ThislepBMass *ThislepBMass + ThisMTnub*ThisMTnub)
				MTtop.append(ThisMTtop)
		
				if hasattr(event, 'metMuEGClean_pt'):
					ThisMetovTop =  event.metMuEGClean_pt/ topP4.Pt()
				else:
					ThisMetovTop =  met.pt/ topP4.Pt()
				METovTop.append(ThisMetovTop)
		
				ThisMetTop = metp4.DeltaPhi(metp4+lepJp4)
				METTopPhi.append(ThisMetTop)
		
				TopMass1 = (neutrino1+lepJp4).M()
				TopMass2 = (neutrino2+lepJp4).M()
		
				#take smaller mtop of the two nu pz-solutions
				if TopMass1 > TopMass2:
					Mtop.append(TopMass2)
				else:
					Mtop.append(TopMass1)
		
				ThisMtopDecor1  = math.sqrt(lepJp4.M()*lepJp4.M()+ (neutrino1+jet.p4()).M()*(neutrino1+jet.p4()).M()+81*81)
				ThisMtopDecor2  = math.sqrt(lepJp4.M()*lepJp4.M()+ (neutrino2+jet.p4()).M()*(neutrino2+jet.p4()).M()+81*81)
		
				if ThisMtopDecor1 > ThisMtopDecor2:
					MtopDecor.append(ThisMtopDecor2)
				else:
					MtopDecor.append(ThisMtopDecor1)
		
		# fill them
		self.out.fillBranch("MTbnu",MTbnu)
		self.out.fillBranch("LepBMass",LepBMass)
		self.out.fillBranch("MTtop",MTtop)
		self.out.fillBranch("Mtop",Mtop)
		self.out.fillBranch("METovTop",METovTop)
		self.out.fillBranch("METTopPhi",METTopPhi)
		self.out.fillBranch("MtopDecor",MtopDecor)
		
		self.out.fillBranch("TopPt", TopPt)
		self.out.fillBranch("TopEt", TopEt)
		
		# nearest b jet to lead lepton
		minDphiLepB = 100
		minDphiLepBidx = -1
		
		if nTightLeps == 1:
			for i, jet in enumerate(centralJet30):
				if jet.btagCSVV2>bTagWP:
					dPhiLepB = abs(jet.p4().DeltaPhi(tightLeps[0].p4()))
					if dPhiLepB < minDphiLepB:
						minDphiLepB = dPhiLepB
						minDphiLepBidx = i
		
		self.out.fillBranch("minDphiLepB",minDphiLepB)
		self.out.fillBranch("minDphiLepBidx", minDphiLepBidx)
		
		#        TopVarsJetIdx = []
		TopVarsMTbnuMin = []
		TopVarsLepBMassMin = []
		TopVarsMTtopMin = []
		TopVarsMtopMin = []
		TopVarsMETovTopMin = []
		TopVarsMtopDecorMin = []
		
		TopVarsTopPtMin = []
		TopVarsTopEtMin = []
		
		iBTagDict = {i: jets[idx].btagCSVV2 for i, idx in enumerate(centralJet30idx)}
		sortIdsByBTag = sorted(iBTagDict.items(), key=operator.itemgetter(1), reverse=True)
		bTaggedJetsSorted = sortIdsByBTag[:nBJetMedium30]
		#        print bTaggedJetsSorted
		bTaggedJetsPPSorted = sortIdsByBTag[:nBJetMedium30+1]
		#        print bTaggedJetsPPSorted
		ThreeBestBTags = sortIdsByBTag[:3]
		#        print ThreeBestBTags
		#        print sortIdsByBTag
		
		if nTightLeps == 1 and nCentralJet30 > 0:
			#0: minimum of b-tagged jets
			TopVarsMTbnuMin      .append(minValueForIdxList(MTbnu     , [ids[0] for ids in bTaggedJetsSorted]))
			TopVarsLepBMassMin   .append(minValueForIdxList(LepBMass  , [ids[0] for ids in bTaggedJetsSorted]))
			TopVarsMTtopMin      .append(minValueForIdxList(MTtop     , [ids[0] for ids in bTaggedJetsSorted]))
			TopVarsMtopMin       .append(minValueForIdxList(Mtop      , [ids[0] for ids in bTaggedJetsSorted]))
			TopVarsMETovTopMin   .append(minValueForIdxList(METovTop  , [ids[0] for ids in bTaggedJetsSorted]))
			TopVarsMtopDecorMin  .append(minValueForIdxList(MtopDecor , [ids[0] for ids in bTaggedJetsSorted]))
			TopVarsTopPtMin      .append(minValueForIdxList(TopPt     , [ids[0] for ids in bTaggedJetsSorted]))
			TopVarsTopEtMin      .append(minValueForIdxList(TopEt     , [ids[0] for ids in bTaggedJetsSorted]))
		
			#1: minimum of b-tagged jets (+ adding next-best b-disc. jet)
			TopVarsMTbnuMin      .append(minValueForIdxList(MTbnu     , [ids[0] for ids in bTaggedJetsPPSorted]))
			TopVarsLepBMassMin   .append(minValueForIdxList(LepBMass  , [ids[0] for ids in bTaggedJetsPPSorted]))
			TopVarsMTtopMin      .append(minValueForIdxList(MTtop     , [ids[0] for ids in bTaggedJetsPPSorted]))
			TopVarsMtopMin       .append(minValueForIdxList(Mtop      , [ids[0] for ids in bTaggedJetsPPSorted]))
			TopVarsMETovTopMin   .append(minValueForIdxList(METovTop  , [ids[0] for ids in bTaggedJetsPPSorted]))
			TopVarsMtopDecorMin  .append(minValueForIdxList(MtopDecor , [ids[0] for ids in bTaggedJetsPPSorted]))
			TopVarsTopPtMin      .append(minValueForIdxList(TopPt     , [ids[0] for ids in bTaggedJetsPPSorted]))
			TopVarsTopEtMin      .append(minValueForIdxList(TopEt     , [ids[0] for ids in bTaggedJetsPPSorted]))
		
			#2: always consider the three jets with the best b-tag values (no pass of any working-point required)
			TopVarsMTbnuMin      .append(minValueForIdxList(MTbnu     , [ids[0] for ids in ThreeBestBTags]))
			TopVarsLepBMassMin   .append(minValueForIdxList(LepBMass  , [ids[0] for ids in ThreeBestBTags]))
			TopVarsMTtopMin      .append(minValueForIdxList(MTtop     , [ids[0] for ids in ThreeBestBTags]))
			TopVarsMtopMin       .append(minValueForIdxList(Mtop      , [ids[0] for ids in ThreeBestBTags]))
			TopVarsMETovTopMin   .append(minValueForIdxList(METovTop  , [ids[0] for ids in ThreeBestBTags]))
			TopVarsMtopDecorMin  .append(minValueForIdxList(MtopDecor , [ids[0] for ids in ThreeBestBTags]))
			TopVarsTopPtMin      .append(minValueForIdxList(TopPt     , [ids[0] for ids in ThreeBestBTags]))
			TopVarsTopEtMin      .append(minValueForIdxList(TopEt     , [ids[0] for ids in ThreeBestBTags]))
			if hasattr(tightLeps[0],'genPartIdx'):
				mcMatchIdLep = tightLeps[0].genPartIdx
				iCorrectJet=-999
				correctJetBTagged = False
				if abs(mcMatchIdLep)==6:
					for i,jet in  enumerate(centralJet30):
						if abs(jet.partonFlavour)==5 and jet.genJetIdx==mcMatchIdLep:
							iCorrectJet=i
							if jet.btagCSVV2>bTagWP: correctJetBTagged=True
		
			#3: value for the correct b-jet (i.e. the one matching the lepton)
				TopVarsMTbnuMin      .append(MTbnu     [iCorrectJet] if iCorrectJet>-999 else -999)
				TopVarsLepBMassMin   .append(LepBMass  [iCorrectJet] if iCorrectJet>-999 else -999)
				TopVarsMTtopMin      .append(MTtop     [iCorrectJet] if iCorrectJet>-999 else -999)
				TopVarsMtopMin       .append(Mtop      [iCorrectJet] if iCorrectJet>-999 else -999)
				TopVarsMETovTopMin   .append(METovTop  [iCorrectJet] if iCorrectJet>-999 else -999)
				TopVarsMtopDecorMin  .append(MtopDecor [iCorrectJet] if iCorrectJet>-999 else -999)
				TopVarsTopPtMin      .append(TopPt     [iCorrectJet] if iCorrectJet>-999 else -999)
				TopVarsTopEtMin      .append(TopEt     [iCorrectJet] if iCorrectJet>-999 else -999)
				
				foundCorrectBJetAndIsTagged = iCorrectJet>-999 and correctJetBTagged
				
				#4: value for the correct b-jet (if actually identified as b-jet (tagged))
				TopVarsMTbnuMin      .append(MTbnu     [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				TopVarsLepBMassMin   .append(LepBMass  [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				TopVarsMTtopMin      .append(MTtop     [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				TopVarsMtopMin       .append(Mtop      [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				TopVarsMETovTopMin   .append(METovTop  [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				TopVarsMtopDecorMin  .append(MtopDecor [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				TopVarsTopPtMin      .append(TopPt     [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				TopVarsTopEtMin      .append(TopEt     [iCorrectJet] if foundCorrectBJetAndIsTagged else -999)
				
				#5: consider the jet closest in dphi to MET
				TopVarsMTbnuMin      .append(MTbnu    [idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				TopVarsLepBMassMin   .append(LepBMass [idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				TopVarsMTtopMin      .append(MTtop    [idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				TopVarsMtopMin       .append(Mtop     [idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				TopVarsMETovTopMin   .append(METovTop [idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				TopVarsMtopDecorMin  .append(MtopDecor[idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				TopVarsTopPtMin      .append(TopPt    [idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				TopVarsTopEtMin      .append(TopEt    [idxMinDPhiBMET] if idxMinDPhiBMET!=-999 else -999)
				
				#6: nearest to lepton b jet
				TopVarsMTbnuMin      .append(MTbnu[minDphiLepBidx]    if minDphiLepBidx>-1 else -999)
				TopVarsLepBMassMin   .append(LepBMass[minDphiLepBidx] if minDphiLepBidx>-1 else -999)
				TopVarsMTtopMin      .append(MTtop[minDphiLepBidx]    if minDphiLepBidx>-1 else -999)
				TopVarsMtopMin       .append(Mtop[minDphiLepBidx]     if minDphiLepBidx>-1 else -999)
				TopVarsMETovTopMin   .append(METovTop[minDphiLepBidx] if minDphiLepBidx>-1 else -999)
				TopVarsMtopDecorMin  .append(MtopDecor[minDphiLepBidx]if minDphiLepBidx>-1 else -999)
				TopVarsTopPtMin      .append(TopPt[minDphiLepBidx]    if minDphiLepBidx>-1 else -999)
				TopVarsTopEtMin      .append(TopEt[minDphiLepBidx]    if minDphiLepBidx>-1 else -999)
				
		else:
			for i in range(7):
				TopVarsMTbnuMin      .append(-999)
				TopVarsLepBMassMin   .append(-999)
				TopVarsMTtopMin      .append(-999)
				TopVarsMtopMin       .append(-999)
				TopVarsMETovTopMin   .append(-999)
				TopVarsMtopDecorMin  .append(-999)
				TopVarsTopPtMin      .append(-999)
				TopVarsTopEtMin      .append(-999)
		
		
		self.out.fillBranch("nBMinVariantsTopVars",7)
		self.out.fillBranch("TopVarsMTbnuMin",TopVarsMTbnuMin)
		self.out.fillBranch("TopVarsLepBMassMin",TopVarsLepBMassMin)
		self.out.fillBranch("TopVarsMTtopMin",TopVarsMTtopMin)
		self.out.fillBranch("TopVarsMtopMin",TopVarsMtopMin)
		self.out.fillBranch("TopVarsMETovTopMin",TopVarsMETovTopMin)
		self.out.fillBranch("TopVarsMtopDecorMin",TopVarsMtopDecorMin)
		self.out.fillBranch("TopVarsTopPtMin",TopVarsTopPtMin)
		self.out.fillBranch("TopVarsTopEtMin",TopVarsTopEtMin)


		
		return True
		

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy1lepTOP = lambda : susyTOP()
 
