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
		if t.Selected == 1 : 
			if t.passFilters : 
				if t.HT > 350 : 
					if t.LT > 150:
						outTree.Fill()
outTree.AutoSave()
fout.Write()
fout.Close()
