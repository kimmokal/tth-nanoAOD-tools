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
workarea_d = '%s/src' % os.environ['CMSSW_BASE']
exearea_d = '%s/src/tthAnalysis/NanoAODTools/batch'% os.environ['CMSSW_BASE']
X509_d = "%s/src/tthAnalysis/NanoAODTools/batch/x509up_u29118"% os.environ['CMSSW_BASE']

workarea_b = '%s/src'% os.environ['CMSSW_BASE']
exearea_b = '%s/src/tthAnalysis/NanoAODTools/batch'% os.environ['CMSSW_BASE']
X509_b = "%s/src/tthAnalysis/NanoAODTools/batch/x509up_u51021"% os.environ['CMSSW_BASE']

#X509name = 'x509up_u29118'
## to be kept 
condTEMP = './templates/submit.condor'
wrapTEMP = './templates/wrapnanoPost.sh'
SkimTEMP = './templates/Skim_tree.py'

if  os.path.exists('submit_all_to_batch_HTC.sh'):
   os.remove('submit_all_to_batch_HTC.sh')

def getlistoffiles(sample):
	filelist = []
	if sample == "//SMS-T1tttt-Scan":
		filelisttxt = open('SMS-T1tttt-Scan.txt','r')
		for line in filelisttxt : 
			line = line.strip()
			filelist.append(line)
	else : 
		CMD = 'dasgoclient -query=\"file dataset='+str(sample)+'\"'
		filelist = subprocess.check_output(CMD, shell=True)
		filelist = filelist.split('\n')
	fullpathlist = []
	#filelist = commands.getstatusoutput('dasgoclient -query=\"file dataset='+str(sample)+'\"'),"\n"
	for f in filelist :
		if f : 
			fi = 'root://cms-xrd-global.cern.ch/'+str(f)
			fullpathlist.append(fi)
		else: continue 
	return fullpathlist
	
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
		modu = "susy_1l_FiltersMC,jecUncert,susy1lepTOPMC,susy_1l_gen"#,xsec,genpartsusymod
		if era == 2016 : 
			modu+=",jetmetUncertainties16,puWeight_2016,btagSF_csvv2_2016,btagSF_cmva_2016,countHistogramAll_2016,susy_1l_Trigg2016,susy1lepMC"
		if era == 2017 :
			# Temporarly use the jetmet uncertainty for 2016 
			modu+=",jetmetUncertainties16,puWeight_2017,btagSF_csvv2_2017,btagSF_deep_2017,countHistogramAll_2017,susy_1l_Trigg2017,susy1lepMC17"
		if "TTJets" in str(sample) and era == 2016: modu+=",susy_1l_nISR16,susy1lepTT_syst"
		if "TTJets" in str(sample) and era == 2017: modu+=",susy_1l_nISR17,susy1lepTT_syst"
		if "WJets"  in str(sample): modu+=",susy1lepWJets_syst"
	elif isMC and isSIG :
		modu = "susy_1l_FiltersMC,jecUncert,puWeight,susy1lepTOPMC,susy_1l_gen,susy1lepSIG_syst"#,xsec,genpartsusymod
			# Temporarly use the jetmet uncertainty for 2016 
		if era == 2016 : modu+=",jetmetUncertainties16,puWeight_2016,btagSF_csvv2_2016,btagSF_cmva_2016,susy_1l_Sig16,countHistogramAll_2016,susy1lepSIG"
		if era == 2017 : modu+=",jetmetUncertainties16,puWeight_2017,btagSF_csvv2_2017,btagSF_deep_2017,susy_1l_Sig17,countHistogramAll_2017,susy1lepSIG17"
	else :
		if era == 2016 : 
			modu = "susy1lepdata,susy_1l_Trigg2016,susy_1l_FiltersData,susy1lepTOPData -J $CMSSW_BASE/src/tthAnalysis/NanoAODTools/data/JSONS/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt" 
		elif era == 2017:
			modu = "susy1lepdata17,susy_1l_Trigg2017,susy_1l_FiltersData,susy1lepTOPData -J $CMSSW_BASE/src/tthAnalysis/NanoAODTools/data/JSONS/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
	return modu

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Runs a NAF batch system for nanoAOD', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('dataList', help='List of datasets to process', metavar='datasets.txt')
	parser.add_argument('--out', help='out put directory',default=None, metavar='Dir')
	parser.add_argument('--site', help='Site DESY or Bari',default=None, metavar='site')
	args = parser.parse_args()
	listtxt = open(args.dataList,"r")
	if args.site == "DESY":
		workarea = workarea_d
		exearea = exearea_d
		X509 = X509_d
	elif args.site == "Bari":
		workarea = workarea_b
		exearea = exearea_b
		X509 = X509_b
		
	outdire = args.out
	
	os.environ['X509_USER_PROXY'] = X509
	os.system('voms-proxy-init --voms cms --valid 120:00')
	

		
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
		modulelist = moduletorun(line)
		print modulelist
		flist = getlistoffiles(line)
		#print flist
		getdir = line.split("/")
		if line == "//SMS-T1tttt-Scan":
			ext = ""
			extension = "Scan"
		else : 
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
		file = open(dirname+'/'+textname+'.txt',"w+")
		for f in flist : 
			file.write(str(f)+'\n')
		file.close()
		i = 0 
		logsdir = dirname+"/logs"
		os.makedirs(logsdir)
		os.system("cp "+SkimTEMP+" "+dirname)
		for f in flist:
			infile2 = f.split("/")
			f2 = infile2[-1].replace(".root","_Skim.root") 
			os.system("cp "+condTEMP+" "+dirname+"/Condor"+textname+str(i)+".submit")
			os.system("cp "+wrapTEMP+" "+dirname+"/Warp"+textname+str(i)+".sh")
			s1 = open(dirname+"/Condor"+textname+str(i)+".submit").read()
			#print textname
			s1 = s1.replace('@EXESH', dirname+"/Warp"+textname+str(i)+".sh").replace('@LOGS',logsdir).replace('@X509',X509).replace('@time','60*60*48')
			f1 = open(dirname+"/Condor"+textname+str(i)+".submit", 'w')
			f1.write(s1)
			f1.close()
			s2 = open(dirname+"/Warp"+textname+str(i)+".sh").read()
			s2 = s2.replace('@WORKDIR',workarea).replace('@EXEDIR',exearea).replace('@MODULES',modulelist).replace('@OUTDIR',dirname).replace('@INPUTFILE',f).replace('@X509',X509).replace("@STEP1",f2).replace("@TRIM",f2.replace("_Skim.root","_Trim.root"))
			f2 = open(dirname+"/Warp"+textname+str(i)+".sh", 'w')
			f2.write(s2)
			f2.close()
			file = open('submit_all_to_batch_HTC.sh','a')
			file.write("\n")
			if args.site == "DESY" : 
				file.write("condor_submit -name s02 "+dirname+"/Condor"+textname+str(i)+".submit")
			elif args.site == "Bari" :
				file.write("condor_submit -name ettore "+dirname+"/Condor"+textname+str(i)+".submit")
			file.close() 
			i+=1
	os.system('chmod a+x submit_all_to_batch_HTC.sh')
	print 'submit_all_to_batch_HTC.sh Created for you - you can run it now with ./'

