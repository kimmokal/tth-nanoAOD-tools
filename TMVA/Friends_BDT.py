#!/usr/bin/env python
########### author : Ashraf Kasem Mohamed ########## 
########### institute : DESY #######################
########### Email : ashraf.mohamed@desy.de #########
########### Date : April 2018#######################
import sys,os, re, pprint
import re
from os import listdir
from os.path import isfile, join
import argparse
import commands
import subprocess
import shutil
from ROOT import TFile

def find_all_matching(substring, path):
    result = []
    for root, dirs, files in os.walk(path):
        for thisfile in files:
            if substring in thisfile:
                result.append(os.path.join(root, thisfile ))
    return result

TMVA_Temp="./TMVA_1lep_App_Temp.C"

if  os.path.exists('submit_TMVA_HTC.sh'):
   os.remove('submit_TMVA_HTC.sh')

condTEMP = './templates_TMVA/submit.condor'
wrapTEMP = './templates_TMVA/wrapnanoPost.sh'

workarea = '%s/src' % os.environ['CMSSW_BASE']
exearea = '%s/src/tthAnalysis/NanoAODTools/TMVA'% os.environ['CMSSW_BASE']

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Runs a NAF batch system for nanoAOD', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--indir', help='List of datasets to process', metavar='indir')
	parser.add_argument('--outdir', help='output directory',default=None, metavar='outdir')
	parser.add_argument('--batchMode', '-b', help="Batch mode.", action='store_true')

	args = parser.parse_args()
	batch = args.batchMode
	indire = args.indir
	outdire = args.outdir
	
	if  os.path.exists(outdire):
		des = raw_input(" this dir is already exist : "+str(outdire)+" do you want to remove it [y/n]: ")
		if ( "y" in des or "Y" in des or "Yes" in des) : 
			shutil.rmtree(str(outdire))
			os.makedirs(str(outdire))
		elif ( "N" in des or  "n" in des or  "No" in des ): print str(outdire) , "will be ovewritten by the job output -- take care"  
		else :
			raise ValueError( "do not understand your potion")
	else : os.makedirs(str(outdire))
	commands = []
	rootfiles = find_all_matching('.root',indire)
	to_merge=''
	logsdir = outdire+"/logs"
	os.makedirs(logsdir)
	for rf in rootfiles:
		script_Name=outdire+"/"+rf.split("/")[-1].replace(".root",".C")
		os.system("cp "+TMVA_Temp+" "+script_Name)
		s1 = open(script_Name).read()
		#print textname
		s1 = s1.replace('@INFILE', rf).replace('@OUTFILE',outdire+"/"+rf.split("/")[-1]).replace("@SCRIPTNAME",rf.split("/")[-1].replace(".root",""))
		f1 = open(script_Name, 'w')
		f1.write(s1)
		f1.close()
		if not batch : 
			cmd = "root -l -q "+script_Name+'\\'+'('+'\\'+'\"'+"BDT"+'\\'+'\"'+'\\'+")"
		else : 
			cmd = "root -l -q -b "+script_Name+'\\'+'('+'\\'+'\"'+"BDT"+'\\'+'\"'+'\\'+")"
		#print "root -l "+script_Name+'\\'+'('+'\\'+'\"'+"BDT"+'\\'+'\"'+'\\'+")"
		if not batch : 
			os.system(cmd)
		else : 
			dirname = outdire+"/"+rf.split("/")[-1].replace(".root","")
			textname = rf.split("/")[-1].replace(".root","")
			os.makedirs(str(dirname))
			os.system("cp "+condTEMP+" "+dirname+"/Condor"+textname+".submit")
			os.system("cp "+wrapTEMP+" "+dirname+"/Warp"+textname+".sh")
			s1 = open(dirname+"/Condor"+textname+".submit").read()
			#print textname
			s1 = s1.replace('@EXESH', dirname+"/Warp"+textname+".sh").replace('@LOGS',logsdir).replace('@time','60*60*3')
			f1 = open(dirname+"/Condor"+textname+".submit", 'w')
			f1.write(s1)
			f1.close()
			s2 = open(dirname+"/Warp"+textname+".sh").read()
			s2 = s2.replace('@WORKDIR',workarea).replace('@EXEDIR',exearea).replace('@COMMAND',cmd)
			f2 = open(dirname+"/Warp"+textname+".sh", 'w')
			f2.write(s2)
			f2.close()
			file = open('submit_TMVA_HTC.sh','a')
			file.write("\n")
			file.write("condor_submit -name s02 "+dirname+"/Condor"+textname+".submit")
			file.close()		
	os.system('chmod a+x submit_TMVA_HTC.sh')
	print 'script is submit_TMVA_HTC is READY FOR BATCH'

