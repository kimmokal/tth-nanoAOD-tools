#!/usr/bin/env python
########### author : Ashraf Kasem Mohamed ########## 
########### institute : DESY #######################
########### Email : ashraf.mohamed@desy.de #########
########### Date : April 2018#######################
import sys,os, re, pprint
import re

import argparse
import commands
import subprocess
import shutil

## to be changed by user 
#X509name = 'x509up_u29118'
## to be kept 
crabTEMP = './templates/crab_cfg.py'
crabscriptTEMP = './templates/crab_script.py'


if  os.path.exists('submit_all_CRABs.sh'):
   os.remove('submit_all_CRABs.sh')

def moduletorun(sample):
	isMC = False 
	isSIG = False 
	isUSER = False 
	era = 'unknown'
	modu = 'unknown'
	if 'RunIISummer16' in str(sample) or "Run2016" in str(sample): era = 2016
	elif 'RunIIFall17' in str(sample) or "Run2017" in str(sample): era = 2017
	if '/NANOAODSIM' in sample:
		isMC = True 
	if "/SMS-T1tttt" in sample :
		isSIG = True
	if "//SMS-T1tttt" in sample : 
		isSIG = True
		isMC = True 
		era = 2016
	if isMC and not isSIG : 
		modu = "susy_1l_FiltersMC(),jecUncert(),susy1lepTOPMC(),susy_1l_gen()"#,xsec,genpartsusymod
		if era == 2016 : 
			modu+=",jetmetUncertainties16(),puWeight_2016(),btagSF_csvv2_2016(),btagSF_cmva_2016(),countHistogramAll_2016(),susy_1l_Trigg2016(),susy1lepMC()"
		if era == 2017 :
			# Temporarly use the jetmet uncertainty for 2016 
			modu+=",jetmetUncertainties16(),puWeight_2017(),btagSF_csvv2_2017(),btagSF_deep_2017(),countHistogramAll_2017(),susy_1l_Trigg2017(),susy1lepMC17()"
		if "TTJets" in str(sample) and era == 2016: modu+=",susy_1l_nISR16(),susy1lepTT_syst()"
		if "TTJets" in str(sample) and era == 2017: modu+=",susy_1l_nISR17(),susy1lepTT_syst()"
		if "WJets"  in str(sample): modu+=",susy1lepWJets_syst()"
	elif isMC and isSIG :
		modu = "susy_1l_FiltersMC(),jecUncert(),puWeight(),susy1lepTOPMC(),susy_1l_gen(),susy1lepSIG_syst()"#,xsec,genpartsusymod
			# Temporarly use the jetmet uncertainty for 2016 
		if era == 2016 : modu+=",jetmetUncertainties16(),puWeight_2016(),btagSF_csvv2_2016(),btagSF_cmva_2016(),susy_1l_Sig16(),countHistogramAll_2016(),susy1lepSIG()"
		if era == 2017 : modu+=",jetmetUncertainties16(),puWeight_2017(),btagSF_csvv2_2017(),btagSF_deep_2017(),susy_1l_Sig17(),countHistogramAll_2017(),susy1lepSIG17()"
	else :
		if era == 2016 : 
			modu = "susy1lepdata(),susy_1l_Trigg2016(),susy_1l_FiltersData(),susy1lepTOPData()"
		elif era == 2017:
			modu = "susy1lepdata17(),susy_1l_Trigg2017(),susy_1l_FiltersData(),susy1lepTOPData()"
	return modu , isMC

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Runs a NAF batch system for nanoAOD', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('dataList', help='List of datasets to process', metavar='datasets.txt')
	parser.add_argument('--out', help='out put directory',default=None, metavar='Dir')
	
	args = parser.parse_args()
	listtxt = open(args.dataList,"r")

	outdire = args.out

	if  os.path.exists(outdire):
		des = raw_input(" this dir is already exist : "+str(outdire)+" do you want to remove it [y/n]: ")
		if ( "y" in des or "Y" in des or "Yes" in des) : 
			shutil.rmtree(str(outdire))
			os.makedirs(str(outdire))
		elif ( "N" in des or  "n" in des or  "No" in des ): print str(outdire) , "will be ovewritten by the job output -- take care"  
		else :
			raise ValueError( "do not understand your potion")
	else : os.makedirs(str(outdire))
	flist = []
	modulelist = [] 
	for line in listtxt : 
		line = line.strip()
		if line.startswith('#') : continue 
		if len(line.strip()) == 0 : continue 
		print line 
		modulelist, MCSamples = moduletorun(line)
		print modulelist
		if MCSamples: 
			sp = 'FileBased'
			unts = 1
		else : 
			sp = 'FileBased'
			unts = 1
		getdir = line.split("/")
		ext = getdir[2].split("-")
		extension = ext[0]
		if "ver1" in getdir[2] : 
			ext1 = "ver1"
			dirname = outdire+"/"+getdir[1]+extension+ext1
			textname = getdir[1]+extension+ext1
		elif "ver2" in getdir[2] : 
			ext1 = "ver2"
			dirname = outdire+"/"+getdir[1]+extension+ext1
			textname = getdir[1]+extension+ext1
		elif "ext1" in getdir[2] : 
			ext1 = "ext1"
			dirname = outdire+"/"+getdir[1]+extension+ext1
			textname = getdir[1]+extension+ext1
		elif "ext2" in getdir[2] : 
			ext1 = "ext2"
			dirname = outdire+"/"+getdir[1]+extension+ext1
			textname = getdir[1]+extension+ext1
		elif "old_pmx" in getdir[2] : 
			ext1 = "old_pmx"
			dirname = outdire+"/"+getdir[1]+extension+ext1
			textname = getdir[1]+extension+ext1
		elif "new_pmx" in getdir[2] : 
			ext1 = "new_pmx"
			dirname = outdire+"/"+getdir[1]+extension+ext1
			textname = getdir[1]+extension+ext1			
		else :
			dirname = outdire+"/"+getdir[1]+extension
			textname = getdir[1]+extension
		
		if  os.path.exists(str(dirname)):
			des = raw_input(" this dir is already exist : "+str(dirname)+" do you want to remove it [y/n]: ")
			if (des == "y" or des =="Y" or des =="Yes") : 
				shutil.rmtree(str(dirname))
				os.makedirs(str(dirname))
			elif (des == "N" or des == "n" or des == "No"): print str(dirname) , "will be ovewritten by the job output -- take care" ; continue 
			else : raise ValueError( "do not understand your potion")
			
		else :os.makedirs(str(dirname))
		
		os.system("cp "+crabTEMP+" "+dirname+"/crab_cfg.py")
		os.system("cp "+crabscriptTEMP+" "+dirname+"/crab_script.py")
		
		s1 = open(dirname+"/crab_cfg.py").read()
		#print textname
		s1 = s1.replace('@REQUEST', str(textname)).replace('@CRAB_SCRIPT',str(dirname)+"/crab_script.py").replace('@UNITES',str(unts)).replace('@TAG',str(textname)).replace('@INPUTDATASET',line).replace('@SPLIT',sp)
		#if not MCSamples:
		#	s1 = s1.replace('#config.Data.lumiMask', 'config.Data.lumiMask')
		f1 = open(dirname+"/crab_cfg.py", 'w')
		f1.write(s1)
		f1.close()
		s2 = open(dirname+"/crab_script.py").read()
		s2 = s2.replace('@MODULES',modulelist)
		f2 = open(dirname+"/crab_script.py", 'w')
		f2.write(s2)
		f2.close()
		file = open('submit_all_CRABs.sh','a')
		file.write("\n")
		file.write("crab submit "+dirname+"/crab_cfg.py")
		file.close() 
	os.system('chmod a+x submit_all_CRABs.sh')
	print 'submit_all_CRABs.sh Created for you - you can run it now with ./'

