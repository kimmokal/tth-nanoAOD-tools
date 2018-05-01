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


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Runs a NAF batch system for nanoAOD', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--indir', help='List of datasets to process', metavar='indir')
	parser.add_argument('--outdir', help='output directory',default=None, metavar='outdir')
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
		rootfiles = find_all_matching('root',indire+subdir)
		to_merge=''
		for rf in rootfiles:
			f = TFile(rf,'READ')
			if not f.IsZombie():
				to_merge += ' '+str(rf)
		print subdir+" is going to be merged"
		os.system('haddnano.py  '+outdire+'/'+subdir+'.root'+' '+str(to_merge))
		print subdir+" has been merged"
	print 'merging done'

