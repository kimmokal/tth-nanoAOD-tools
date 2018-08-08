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
workarea_d = '/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src'
workarea_b = '/lustre/home/amohamed/susy-desy/nanoAOD/CMSSW_9_4_4/src'

condTEMP = './templates/submit.condor'
wrapTEMP = './templates/wrapnanoPost.sh'


if  os.path.exists('submit_filter_HTC.sh'):
   os.remove('submit_filter_HTC.sh')

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
	parser.add_argument('--site', help='Which site to submit DESY or Bari', metavar='site')
	parser.add_argument('--trig','-tr' ,help='use cut base for trigger study', action='store_true')
	args = parser.parse_args()
	
	if args.site == "DESY": 
		workarea = workarea_d 
	elif args.site == "Bari": 
		workarea = workarea_b
		
			
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
	
	rootfiles = find_all_matching('.root',indire)
	i = 0 
	for subdir in rootfiles:
		rf = subdir 
		subdir = subdir.replace(".root","").split("/")[-1]
		if subdir.startswith('.'):  continue
		os.makedirs(str(outdire)+"/"+str(subdir))
		if args.trig :
			os.system("cp ./templates/CutFlow_NANO_Trig.py "+outdire+"/"+subdir)
		else : 
			os.system("cp ./templates/CutFlow_NANO.py "+outdire+"/"+subdir)
		logsdir = outdire+'/'+subdir+"/logs"
		os.makedirs(logsdir)
		print " prepering ",subdir
		f = TFile(rf,'READ')
		if not f.Get("Events"): 
			print "corrupted file check if you want : ",f 
			continue
		if f.IsZombie():
			print f," IsZombie"
			continue
		if not args.batchMode : 
			retval = os.getcwd()
			os.chdir(outdire+"/"+subdir)
			if args.trig:
				os.system('./CutFlow_NANO_Trig.py  '+rf)
			else : os.system('./CutFlow_NANO.py  '+rf)
			os.chdir(retval)
			print subdir+" has been filtered"
		else : 
			os.system("cp "+condTEMP+" "+outdire+'/'+subdir+"/Condor"+str(i)+".submit")
			os.system("cp "+wrapTEMP+" "+outdire+'/'+subdir+"/Warp"+str(i)+".sh")		
			s1 = open(outdire+'/'+subdir+"/Condor"+str(i)+".submit").read()
			#print textname
			s1 = s1.replace('@EXESH', outdire+'/'+subdir+"/Warp"+str(i)+".sh").replace('@LOGS',logsdir).replace('@time','60*60*6')
			f1 = open(outdire+'/'+subdir+"/Condor"+str(i)+".submit", 'w')
			f1.write(s1)
			f1.close()
			s2 = open(outdire+'/'+subdir+"/Warp"+str(i)+".sh").read()
			s2 = s2.replace('@WORKDIR',workarea).replace('@EXEDIR',outdire+'/'+subdir).replace('@ROOTFILE', rf )
			if args.trig:
				s2 = s2.replace('@COMAND','./CutFlow_NANO_Trig.py ')
			else :s2 = s2.replace('@COMAND','./CutFlow_NANO.py ')
			f2 = open(outdire+'/'+subdir+"/Warp"+str(i)+".sh", 'w')
			f2.write(s2)
			f2.close()
			file = open('submit_filter_HTC.sh','a')
			file.write("\n")
			if args.site == "DESY":
				file.write("condor_submit -name s02 "+outdire+'/'+subdir+"/Condor"+str(i)+".submit")
			elif args.site == "Bari" : 
				file.write("condor_submit -name ettore "+outdire+'/'+subdir+"/Condor"+str(i)+".submit")
			file.close()		
		i+=1
	os.system('chmod a+x submit_filter_HTC.sh')
	print 'merging done OR READY FOR BATCH'

