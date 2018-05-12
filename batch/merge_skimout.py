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
workarea = '/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src'
exearea = '/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src/tthAnalysis/NanoAODTools/batch'
condTEMP = './templates_merging/submit.condor'
wrapTEMP = './templates_merging/wrapnanoPost.sh'

if  os.path.exists('submit_MERGE_HTC.sh'):
   os.remove('submit_MERGE_HTC.sh')

def find_all_matching(substring, path):
    result = []
    for root, dirs, files in os.walk(path):
        for thisfile in files:
            if substring in thisfile:
                result.append(os.path.join(root, thisfile ))
    return result


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Runs a NAF batch system for nanoAOD', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--indir', help='List of datasets to process', metavar='indir')
	parser.add_argument('--outdir', help='output directory',default=None, metavar='outdir')
	parser.add_argument('--batchMode', '-b', help="Batch mode.", action='store_true')
	args = parser.parse_args()
	
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
	subdires = [f for f in listdir(indire) ]
	for subdir in subdires:
		if subdir.startswith('.'):  continue
		rootfiles = find_all_matching('root',indire+subdir)
		to_merge=''
		for rf in rootfiles:
			f = TFile(rf,'READ')
			if not f.Get("Events"): 
				print "corrupted file check if you want : ",f 
				continue
			if not f.IsZombie():
				to_merge += ' '+str(rf)
		print subdir+" is going to be merged"
		if not args.batchMode : 
			os.system('haddnano.py  '+outdire+'/'+subdir+'.root'+' '+str(to_merge))
			print subdir+" has been merged"
		else : 
			logsdir = outdire+'/'+subdir+"/logs"
			os.makedirs(logsdir)
			os.system("cp "+condTEMP+" "+outdire+'/'+subdir+"/Condor.submit")
			os.system("cp "+wrapTEMP+" "+outdire+'/'+subdir+"/Warp.sh")		
			s1 = open(outdire+'/'+subdir+"/Condor.submit").read()
			#print textname
			s1 = s1.replace('@EXESH', outdire+'/'+subdir+"/Warp.sh").replace('@LOGS',logsdir).replace('@time','60*60*6')
			f1 = open(outdire+'/'+subdir+"/Condor.submit", 'w')
			f1.write(s1)
			f1.close()
			s2 = open(outdire+'/'+subdir+"/Warp.sh").read()
			s2 = s2.replace('@WORKDIR',workarea).replace('@EXEDIR',exearea).replace('@COMMAND','haddnano.py  '+outdire+'/'+subdir+'.root'+' '+str(to_merge))
			f2 = open(outdire+'/'+subdir+"/Warp.sh", 'w')
			f2.write(s2)
			f2.close()
			file = open('submit_MERGE_HTC.sh','a')
			file.write("\n")
			file.write("condor_submit "+outdire+'/'+subdir+"/Condor.submit")
			file.close()		
	os.system('chmod a+x submit_MERGE_HTC.sh')
	print 'merging done OR READY FOR BATCH'

