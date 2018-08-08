#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import os
f = sys.argv[-1]

TreeName="Events"
file0 = ROOT.TFile.Open(f)

t = file0.Get(TreeName)	
hcount = file0.Get("Count")
hcountPos = file0.Get("CountPosWeight")
hcountNeg = file0.Get("CountNegWeight")
t.SetBranchStatus("*", 1)
f = f.strip()
outname = f.split("/")[-1]
fout = ROOT.TFile(outname, 'RECREATE' )
fout.cd()
outTree = t.CloneTree(0)
for i in range(t.GetEntries()):
	t.GetEntry(i)
	print "analyzing entry number ",i, " out of ", t.GetEntries()
	if t.nLep == 1:
		if t.Lep_pt > 25 : 
			if t.Selected == 1 : 
				if t.nVeto == 0 : 
					if (not t.isData) or t.HLT_EleOR or t.HLT_MuOR or t.HLT_MetOR: 
						if (not t.isData) or (t.PD_SingleEle and t.HLT_EleOR) or (t.PD_SingleMu and t.HLT_MuOR and not t.HLT_EleOR) or (t.PD_MET and t.HLT_MetOR and not t.HLT_EleOR and not t.HLT_MuOR):
							if (not t.isData) or t.Flag_METFilters : 
								if  t.RA2_muJetFilter == 1 :
									if t.Flag_fastSimCorridorJetCleaning:
										if t.nJets30Clean >= 5 :
											if  t.Jet2_pt > 80 :
												if t.HT > 500 :
													if t.LT > 250 :
														if t.nBJet >= 1: 
															outTree.Fill()

outTree.AutoSave()
if file0.GetListOfKeys().Contains("Count"):
	hcount.Write()
	hcountPos.Write()
	hcountNeg.Write()
fout.Write()
fout.Close()
