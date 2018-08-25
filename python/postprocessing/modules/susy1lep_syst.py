import ROOT
import math, os
import array
import operator
from leptonselector import leptonSel
ROOT.PyConfig.IgnoreCommandLineOptions = True
# for met object 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple

import susy1lep_base

import itertools
from ROOT import TLorentzVector, TVector2, std
mt2obj = ROOT.heppy.Davismt2.Davismt2()
from deltar import bestMatch

def getLV(p4):
    if p4 != None: return ROOT.LorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E())
    else: return p4

def getGenWandLepton(event):

    wP4 = None
    lepP4 = None

    genParts = [p for p in Collection(event,"GenPart","nGenPart")]
    genLeps = filter(lambda l:abs(l.pdgId) in [11,13,15], genParts)

    if len(genLeps) == 0:
        print "no gen lepton found!"
        return wP4, lepP4

    lFromW = filter(lambda w:abs(w[w.genPartIdxMother].pdgId)==24, genLeps)

    if len(lFromW) == 0:
        print "no gen W found!", genLeps
        return wP4, lepP4

    elif len(lFromW)>1:
        print 'this should not have happened'
        return wP4, lepP4

    elif len(lFromW) == 1:

        genLep = lFromW[0]
        genW = genParts[genLep.genPartIdxMother]

        wP4 = getLV(genW.p4())
        lepP4 = getLV(genLep.p4())

        #print genW.p4().M(),genLep.p4().M()
        #print wP4.M(),lepP4.M()

    return wP4, lepP4

def getGenTopWLepton(event):

    topP4 = None
    wP4 = None
    lepP4 = None

    genParts = [p for p in Collection(event,"GenPart","nGenPart")]
    genLeps = filter(lambda l:abs(l.pdgId) in [11,13,15], genParts)

    if len(genLeps) == 0:
        #print "no gen lepton found!" # happens in TTJets ;)
        return topP4, wP4, lepP4

    lFromW = filter(lambda w:abs(genParts[w.genPartIdxMother].pdgId)==24, genLeps)
    

    if len(lFromW) == 0:
        print "no gen W found!", genLeps
        return topP4, wP4, lepP4

    elif len(lFromW) > 2:
        print "More than 2 W's found!"
        return topP4, wP4, lepP4

    elif len(lFromW) == 1:

        genLep = lFromW[0]
        genW = genParts[genLep.genPartIdxMother]
        genTop = genParts[genW.genPartIdxMother]

        topP4 = getLV(genTop.p4())
        wP4 = getLV(genW.p4())
        lepP4 = getLV(genLep.p4())

        return topP4, wP4, lepP4

    elif len(lFromW) == 2:
        match = False
        goodLep = [l for l in event.selectedLeptons]
        LepOther = [l for l in event.otherLeptons]
        nLepGood = len(goodLep)
        if nLepGood > 0:
        	leadLep = event.selectedLeptons[0]
        	for genLep in lFromW:
				# we do not have the GenPart.charge in nanoAOD (bypass it for now) --> to be reported 
				# update (we can use the pdgid to identify the charge)
        		if (leadLep.pdgId > 0 and genLep.pdgId > 0) or (leadLep.pdgId < 0 and genLep.pdgId < 0):
        			match == True
        
        			genW = genParts[genLep.genPartIdxMother]
        			genTop = genParts[genW.genPartIdxMother]
        
        			topP4 = getLV(genTop.p4())
        			wP4 = getLV(genW.p4())
        			lepP4 = getLV(genLep.p4())
        
        			return topP4, wP4, lepP4
        
        if not match:
        	#print 'No match at all!'
        	return topP4, wP4, lepP4
        
    return topP4, wP4, lepP4

def getWPolWeights(event, sample):

    wUp = 1
    wDown = 1

    if "TTJets" in sample: #W polarization in TTbar
        topP4, wP4, lepP4 = getGenTopWLepton(event)

        if topP4 != None:
            #print topP4.M(), wP4.M(), lepP4.M()
            cosTheta = ROOT.ttbarPolarizationAngle(topP4, wP4, lepP4)
            #print cosTheta
            wUp = (1. + 0.05*(1.-cosTheta)**2) * 1./(1.+0.05*2./3.) * (1./1.0323239521945559)
            wDown = (1. - 0.05*(1.-cosTheta)**2) * 1./(1.-0.05*2./3.) * (1.034553190276963956)

    elif "WJets" in sample: #W polarization in WJets
        wP4, lepP4 = getGenWandLepton(event)

        if wP4 != None:
            #print wP4.M(), lepP4.M()
            cosTheta = ROOT.WjetPolarizationAngle(wP4, lepP4)
            #print cosTheta
            wUp = (1. + 0.1*(1.-cosTheta)**2) * 1./(1.+0.1*2./3.) * (1./1.04923678332724659)
            wDown = (1. - 0.1*(1.-cosTheta)**2) * 1./(1.-0.1*2./3.) * (1.05627060952003952)

    #print wUp, wDown

    return wUp, wDown

class susysinglelep_syst(Module):
	def __init__(self, isMC , isSig, era, isTTJets,isWJets):#, HTFilter, LTFilter):#, muonSelection, electronSelection):
		self.isMC = isMC
		self.isSig = isSig
		self.era = era
		self.isTTJets = isTTJets
		self.isWJets = isWJets
		if "/WPolarizationVariation_C.so" not in ROOT.gSystem.GetLibraries():
			print "Load C++ Worker"
			ROOT.gROOT.ProcessLine(".L %s/src/tthAnalysis/NanoAODTools/python/postprocessing/modules/WPolarizationVariation.C+" % os.environ['CMSSW_BASE'])
		if "/TTbarPolarization_C.so" not in ROOT.gSystem.GetLibraries():
			print "Load C++ Worker"
			ROOT.gROOT.ProcessLine(".L %s/src/tthAnalysis/NanoAODTools/python/postprocessing/modules/TTbarPolarization.C+" % os.environ['CMSSW_BASE'])
		pass
	def beginJob(self):
		pass
	def endJob(self):
		pass
	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		self.branchlist = [
			# Top related
			"GenTopPt", "GenAntiTopPt", "TopPtWeight","TopPtWeightII", "GenTTBarPt", "GenTTBarWeight",
			# ISR
			"ISRTTBarWeight", "GenGGPt", "ISRSigUp", "ISRSigDown",
			# DiLepton
			"DilepNJetCorr", #to shift central value 
			"DilepNJetWeightConstUp", "DilepNJetWeightSlopeUp", "DilepNJetWeightConstDn", "DilepNJetWeightSlopeDn",
			# W polarisation
			"WpolWup","WpolWdown",
			# TTJets ISR 
			'nISRtt','nISRttweight','nISRttweightsyst_up', 'nISRttweightsyst_down',
			]
		for branch in self.branchlist : 
			#print branch
			self.out.branch(str(branch),  "F");
		self.sample = inputFile
	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass
	def analyze(self, event):
 		"""process event, return True (go to next module) or False (fail, go to next event)"""
		leptonSel(event)
		goodLep = [l for l in event.selectedLeptons]
		LepOther = [l for l in event.otherLeptons]
				
		#LepOther = goodLep
		leps = goodLep 
		nlep = len(leps)

 		#print "from anayzer sample is : ", str(self.sample)
 		#### W polarisation
 		if self.isWJets : 
			wPolWup, wPolWdown = getWPolWeights(event, "WJets")
		elif self.isTTJets : 
			wPolWup, wPolWdown = getWPolWeights(event, "TTJets")
		else : wPolWup, wPolWdown = getWPolWeights(event, "anythingelse")
		
		self.out.fillBranch("WpolWup",wPolWup)
		self.out.fillBranch("WpolWdown",wPolWdown)
		### TOP RELATED VARS
		gParts = Collection(event,"GenPart")
		genParts = [g for g in gParts]
		GenTopPt = -999
		GenTopIdx = -999
		GenAntiTopPt = -999
		GenAntiTopIdx = -999
		TopPtWeight = 1.
		TopPtWeightII = 1.
		GenTTBarPt = -999
		GenTTBarWeight = 1.
		ISRTTBarWeight = 1.
		GenGGPt = -999
		ISRSigUp = 1.
		ISRSigDown = 1.
		nGenTops = 0
		GluinoIdx = []
		
		for i_part, genPart in enumerate(genParts):
			if genPart.pdgId ==  6:
				# protect double counting the GenTop as mother and as daughters
				if genParts[genPart.genPartIdxMother].pdgId == 6 : continue  
				GenTopPt = genPart.pt
				GenTopIdx = i_part
				nGenTops+=1
			if genPart.pdgId == -6:
				# protect double counting the GenTop as mother and as daughters
				if genParts[genPart.genPartIdxMother].pdgId == -6 : continue  
				GenAntiTopPt = genPart.pt
				GenAntiTopIdx = i_part
				nGenTops+=1
			if genPart.pdgId == 1000021:
				GluinoIdx.append(i_part)

		if len(GluinoIdx)==2:
			GenGluinoGluinop4 = genParts[GluinoIdx[0]].p4()+ genParts[GluinoIdx[1]].p4()
			GenGluinoGluinoPt = GenGluinoGluinop4.Pt()
			GenGGPt = GenGluinoGluinoPt
			if GenGluinoGluinoPt > 400: ISRSigUp = 1.15; ISRSigDown = 0.85
			if GenGluinoGluinoPt > 600: ISRSigUp = 1.30; ISRSigDown = 0.70
		#print GenTopPt," ",GenAntiTopPt," ",nGenTops
		if GenTopPt!=-999 and GenAntiTopPt!=-999 and nGenTops==2:
			SFTop     = math.exp(0.156    -0.00137*GenTopPt    )
			SFAntiTop = math.exp(0.156    -0.00137*GenAntiTopPt)
			SFTopII     = math.exp(0.0615    -0.0005*GenTopPt    )
			SFAntiTopII = math.exp(0.0615    -0.0005*GenAntiTopPt)
			TopPtWeight = math.sqrt(SFTop*SFAntiTop)
			TopPtWeightII = math.sqrt(SFTopII*SFAntiTopII)
			if TopPtWeight<0.5: TopPtWeight=0.5
			if TopPtWeightII<0.5: TopPtWeightII=0.5
			if GenAntiTopIdx!=-999 and GenTopIdx!=-999:
				GenTTBarp4 = genParts[GenTopIdx].p4()+ genParts[GenAntiTopIdx].p4()
				GenTTBarPt = GenTTBarp4.Pt()
				if GenTTBarPt>120: GenTTBarWeight= 0.95
				if GenTTBarPt>150: GenTTBarWeight= 0.90
				if GenTTBarPt>250: GenTTBarWeight= 0.80
				if GenTTBarPt>400: GenTTBarWeight= 0.70
				if GenTTBarPt>400: ISRTTBarWeight = 0.85
				if GenTTBarPt>600: ISRTTBarWeight = 0.7
			#print  "ISRTTBarWeight,GenTTBarPt ",ISRTTBarWeight,GenTTBarPt 
			#print  "TopPtWeight, TopPtWeightII ",TopPtWeight, TopPtWeightII
		####################################
		### For DiLepton systematics
		# values in sync with AN2015_207_v3
		#        Const weight
		# const: 0.85 +-0.06
		#        16%
		# 2016 update with 35.9/fb
		#https://indico.cern.ch/event/611061/contributions/2464202/attachments/1406419/2149049/diLepStudy.pdf
		constVariation= math.sqrt(0.030*0.030 +0.023*0.023)
		slopevariation = math.sqrt(0.017*0.017 +0.014*0.014)
		wmean = 6.93 - 0.5
		Jets = Collection(event, "Jet")
		jets = [j for j in Jets if j.pt > 20 and abs(j.eta) < 2.4]
		njet = len(jets)
		centralJet30 = []; centralJet30idx = []
		for i,j in enumerate(jets):
			if j.pt > 30 and abs(j.eta)<2.4:
				centralJet30.append(j)
				centralJet30idx.append(i)
		
		cJet30Clean = []
		dRminCut = 0.4
		# Do cleaning a la CMG: clean max 1 jet for each lepton (the nearest)
		cJet30Clean = centralJet30
		cleanJets30 = centralJet30
		#clean selected leptons at First 
		for lep in goodLep:
			if lep.pt < 20 : continue 
			jNear, dRmin = None, 99
			# find nearest jet
			for jet in centralJet30:
				dR = jet.p4().DeltaR(lep.p4())
				if dR < dRmin:
					jNear, dRmin = jet, dR
			# remove nearest jet
			if dRmin < dRminCut:
				cJet30Clean.remove(jNear)
		nJets30Clean = len(cJet30Clean)
		
		genTau = [l for l in Collection(event, "GenPart") if abs(l.pdgId) == 15]
		ngenTau = len(genTau)
		genLep = [l for l in Collection(event, "GenPart") if abs(l.pdgId) in [11,13]]
		ngenLep = len(genLep)
		
		sumnGenLepTau=0
		genTau_Mothidx = -1 
		genTau_gMothidx = -1 
		genLep_Mothidx = -1 
		genLep_gMothidx = -1
		 
		for gT in genTau:
			genTau_Mothidx = gT.genPartIdxMother
			genTau_gMothidx = genParts[genTau_Mothidx].genPartIdxMother
			if genTau_Mothidx >= 0 or genTau_gMothidx >= 0 : 
				if abs(genParts[genTau_gMothidx].pdgId)==6 and abs(genParts[genTau_Mothidx].pdgId)==24: sumnGenLepTau+=1
		for gL in genLep:
			genLep_Mothidx = gL.genPartIdxMother
			genLep_gMothidx = genParts[genLep_Mothidx].genPartIdxMother
			if genLep_Mothidx >= 0 or genLep_gMothidx >= 0 :
				if abs(genParts[genLep_gMothidx].pdgId)==6 and abs(genParts[genLep_Mothidx].pdgId)==24: sumnGenLepTau+=1
		#if (event.ngenLep+event.ngenTau)==2: #would like to restore this behavior...
		if sumnGenLepTau==2:
			self.out.fillBranch("DilepNJetCorr",1.030-0.017*(nJets30Clean-wmean))
			self.out.fillBranch("DilepNJetWeightConstUp",1-constVariation)
			self.out.fillBranch("DilepNJetWeightSlopeUp",1+ (nJets30Clean-wmean)*slopevariation)
			self.out.fillBranch("DilepNJetWeightConstDn",1+constVariation)
			self.out.fillBranch("DilepNJetWeightSlopeDn",1- (nJets30Clean-wmean)*slopevariation)
		else:
			self.out.fillBranch("DilepNJetCorr",1.)
			self.out.fillBranch("DilepNJetWeightConstUp",1.)
			self.out.fillBranch("DilepNJetWeightSlopeUp",1.)
			self.out.fillBranch("DilepNJetWeightConstDn",1.)
			self.out.fillBranch("DilepNJetWeightSlopeDn",1.)
		self.out.fillBranch("GenTopPt",GenTopPt)
		self.out.fillBranch("GenAntiTopPt",GenAntiTopPt)
		self.out.fillBranch("TopPtWeight",TopPtWeight)
		self.out.fillBranch("TopPtWeightII",TopPtWeightII)
		self.out.fillBranch("GenTTBarPt",GenTTBarPt)
		self.out.fillBranch("GenTTBarWeight",GenTTBarWeight)
		self.out.fillBranch("ISRTTBarWeight",ISRTTBarWeight)
		self.out.fillBranch("GenGGPt",GenGGPt)
		self.out.fillBranch("ISRSigUp",ISRSigUp)
		self.out.fillBranch("ISRSigDown",ISRSigDown)

		
		
		return True
susy1lepSIG_syst  = lambda : susysinglelep_syst(True , True,  "2016",True , False)
susy1lepMC_syst   = lambda : susysinglelep_syst(True , False, "2016",True , False)
