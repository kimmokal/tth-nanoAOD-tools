import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from math import floor
import math
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

# import gluino xsec table
xsecGlu = {} # dict for xsecs
xsecFile = "/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src/tthAnalysis/NanoAODTools/python/postprocessing/modules/glu_xsecs_13TeV.txt"
cntsSusy = {} # dict for signal counts
C_ISRweightsSusy = {}
cntFile = "/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src/tthAnalysis/NanoAODTools/python/postprocessing/modules/counts_T1tttt_2016.txt"
ISRweightFile = "/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src/tthAnalysis/NanoAODTools/python/postprocessing/modules/ISRnormWeights_T1tttt2016.txt"

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

def loadSUSYparams():
	
	global xsecGlu
	global cntsSusy
	global C_ISRweightsSusy
	
	print 80*"#"
	print "Loading SUSY parameters"
	
	with open(xsecFile,"r") as xfile:
		lines = xfile.readlines()
		print 'Found %i lines in %s' %(len(lines),xsecFile)
		for line in lines:
			if line[0] == '#': continue
			(mGo,xsec,err) = line.split()
			#print 'Importet', mGo, xsec, err, 'from', line
			xsecGlu[int(mGo)] = (float(xsec),float(err))
	
		print 'Filled %i items to dict' % (len(xsecGlu))
		#print sorted(xsecGlu.keys())
	
	with open(cntFile,"r") as cfile:
		lines = cfile.readlines()
		print 'Found %i lines in %s' %(len(lines),cntFile)
	
		for line in lines:
			if line[0] == '#': continue
			else:
				(mGo,mLSP,tot,totW,cnt,wgt) = line.split()
				#print 'Importet', mGo, mLSP, cnt, 'from', line
				#cntsSusy[(int(mGo),int(mLSP))] = (int(tot),int(cnt),float(wgt))
				cntsSusy[(int(mGo),int(mLSP))] = (float(totW),int(cnt),float(wgt))
	
		print 'Filled %i items to dict' % (len(cntsSusy))
		print "Finished signal parameter load"
	
	with open(ISRweightFile,"r") as cfile:
		lines = cfile.readlines()
		print 'Found %i lines in %s' %(len(lines),ISRweightFile)
	
		for line in lines:
			if line[0] == '#': continue
			else:
				(mGo,mLSP,C_ISRweight,C_ISRweight_up,C_ISRweight_down) = line.split()
				#print 'Importet', mGo, mLSP, cnt, 'from', line
				#cntsSusy[(int(mGo),int(mLSP))] = (int(tot),int(cnt),float(wgt))
				C_ISRweightsSusy[(int(mGo),int(mLSP))] = (C_ISRweight,C_ISRweight_up,C_ISRweight_down)
	
		print 'Filled ISR weights %i items to dict' % (len(C_ISRweightsSusy))
	
	return 1

class susy1lepSig(Module):
    def __init__(self,isScan):
		self.isScan = isScan 
		self.susyParticles = {
			100001 : 'Squark',
			(1000000 + 21) : 'Gluino',
			(1000000 + 39) : 'Gravitino',
			(1000000 + 5) : 'Sbottom',
			(2000000 + 5) : 'Sbottom2',
			(1000000 + 6) : 'Stop',
			(2000000 + 6) : 'Stop2',
			(1000000 + 15) : 'Stau',
			(2000000 + 15) : 'Stau2',
			(1000000 + 16) : 'SnuTau',
			(1000000 + 22) : 'Neutralino',
			(1000000 + 23) : 'Neutralino2',
			(1000000 + 25) : 'Neutralino3',
			(1000000 + 35) : 'Neutralino4',
			(1000000 + 24) : 'Chargino',
			(1000000 + 37) : 'Chargino2',
		} 
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		self.out.branch("nIsr","I");
		self.out.branch("mGo","F");
		self.out.branch("mLSP","F");
		self.out.branch("susyXsec","F");
		self.out.branch("susyNgen","F");
		self.out.branch("totalNgen","F");
		self.out.branch("susyWgen","F");
		self.out.branch("nISRweight","F");
		self.out.branch("nISRweightsyst_up","F");
		self.out.branch("nISRweightsyst_down","F");
		self.inFile = inputFile.GetName();
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
		masses = {}
		for p in Collection(event, "GenPart"):
			id = abs(p.pdgId)
			if (id / 1000000) % 10 in [1,2]:
				particle = None
				if id % 100 in [1,2,3,4]:
					particle = "Squark"
				elif id in self.susyParticles:
					particle = self.susyParticles[id]
				if particle != None:
					if particle not in masses: masses[particle] = []
					masses[particle].append(p.mass)
		for p,ms in masses.iteritems():
			avgmass = floor(sum(ms)/len(ms)+0.5)
			setattr(event, "GenSusyM"+p, avgmass)
		global xsecGlu
		global cntsSusy
		global C_ISRweightsSusy 
		
		if len(xsecGlu) == 0: loadSUSYparams()
		
		## MASS POINT
		mGo = 0
		mLSP = 0
		# Gluino Mass
		mGo = getattr(event, "GenSusyMGluino")
		# LSP Mass
		mLSP = getattr(event,"GenSusyMNeutralino")
		# set LSP mass of 1 to zero
		if mLSP == 1: mLSP = 0;

		self.out.fillBranch("mGo", mGo); self.out.fillBranch("mLSP", mLSP)
		
		# SUSY Xsec
		if mGo in xsecGlu:
			self.out.fillBranch("susyXsec", xsecGlu[mGo][0])
			#ret['susyXsecErr'] = xsecGlu[mGo][1]
		elif mGo > 0:
			print 'Xsec not found for mGo', mGo	
		nISR = 0
		Jets = Collection(event, "Jet")
		jets = [j for j in Jets if j.pt > 20 ]
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
		nISRforWeights = int(nIsr)
		if nIsr > 6:
			nISRforWeights = 6
			
			
		C_ISR = float(C_ISRweightsSusy[(mGo,mLSP)][0])
		C_ISR_up = float(C_ISRweightsSusy[(mGo,mLSP)][1])
		C_ISR_down = float(C_ISRweightsSusy[(mGo,mLSP)][2])
		
		nISRweight = C_ISR * ISRweights[nISRforWeights]
		nISRweightsyst_up   =  C_ISR_up   * (ISRweights[nISRforWeights] + ISRweightssyst[nISRforWeights])
		nISRweightsyst_down =  C_ISR_down * (ISRweights[nISRforWeights] - ISRweightssyst[nISRforWeights])
		
		self.out.fillBranch("nISRweight",nISRweight)
		self.out.fillBranch("nISRweightsyst_up",nISRweightsyst_up)
		self.out.fillBranch("nISRweightsyst_down",nISRweightsyst_down)		
		##############
		if (mGo,mLSP) in cntsSusy:
			#ret['totalNgen'] = cntsSusy[(mGo,mLSP)][0] # merged scan: 93743963
			if self.isScan == True : self.out.fillBranch("totalNgen", 93743963)
			else: self.out.fillBranch("totalNgen", cntsSusy[(mGo,mLSP)][0])
			self.out.fillBranch("susyNgen",cntsSusy[(mGo,mLSP)][1])
			self.out.fillBranch("susyWgen", cntsSusy[(mGo,mLSP)][2])
		else:
			self.out.fillBranch("totalNgen", 1)
			self.out.fillBranch("susyNgen",1)
			self.out.fillBranch("susyWgen",1)
		
		return True
        

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy_1l_Sig  = lambda : susy1lepSig(False)

