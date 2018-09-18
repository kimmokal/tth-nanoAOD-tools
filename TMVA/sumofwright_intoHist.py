#!/usr/bin/env python
import pickle , re , os
from ROOT import *
import sys,os, re, pprint
import re
from os import listdir
from os.path import isfile, join
import argparse
import commands
import subprocess
import shutil

def find_all_matching(substring, path):
    result = "None"
    for pck, dirs, files in os.walk(path):
        for thisfile in files:
            if substring in thisfile:
                result=os.path.join(pck, thisfile )
    return result

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Get Sum of weight into root file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--trees', help='List of datasets to process', metavar='indir')
	parser.add_argument('--outdir', help='output directory',default=None, metavar='outdir')
	parser.add_argument('--friends','--f', help='output directory',default=None, metavar='outdir')

	args = parser.parse_args()
	indir = args.trees
	outdir = args.outdir
	friends = args.friends
	newfriendsdir = outdir+'/newfriends'
	if  os.path.exists(outdir):
		des = raw_input(" this dir is already exist : "+str(outdir)+" do you want to remove it [y/n]: ")
		if ( "y" in des or "Y" in des or "Yes" in des) : 
			shutil.rmtree(str(outdir))
			os.makedirs(str(outdir))
		elif ( "N" in des or  "n" in des or  "No" in des ): print str(outdir) , "will be ovewritten by the job output -- take care"  
		else :
			raise ValueError( "do not understand your potion")
	else : os.makedirs(str(outdir))
	
	if  os.path.exists(newfriendsdir):
		des = raw_input(" this dir is already exist : "+str(newfriendsdir)+" do you want to remove it [y/n]: ")
		if ( "y" in des or "Y" in des or "Yes" in des) : 
			shutil.rmtree(str(newfriendsdir))
			os.makedirs(str(newfriendsdir))
		elif ( "N" in des or  "n" in des or  "No" in des ): print str(newfriendsdir) , "will be ovewritten by the job output -- take care"  
		else :
			raise ValueError( "do not understand your potion")
	else : os.makedirs(str(newfriendsdir))
	
	
	subdires = [f for f in listdir(indir) ]
	for subdir in subdires:
		if subdir.startswith('.'):  continue
		pckfile = find_all_matching("SkimReport.pck",indir+subdir)
		pckobj  = pickle.load(open(str(pckfile),'r'))
		counters = dict(pckobj)
		friendname = friends+'/evVarFriend_'+pckfile.split("/")[-3]+'.root'
		if  not os.path.exists(friendname): 
			print friendname, " is not there i will not produce weightfile for it chack please "
			continue 
		total_w = 0
		if ('Sum Weights' in counters) :
			total_w += counters['Sum Weights']
		else:
			total_w += counters['All Events']
		print total_w
		hfile = TFile( outdir+'/'+pckfile.split("/")[-3]+'_W.root', 'RECREATE' )
		hfile.cd()
		Sum_Weights = TH1F( 'Sum_Weights', 'Sum Weights', 1, 0, 2 )
		Sum_Weights.Fill(1,total_w)
		hfile.Write()
		hfile.Close()
		os.system("hadd "+newfriendsdir+'/evVarFriend_'+pckfile.split("/")[-3]+'.root '+outdir+'/'+pckfile.split("/")[-3]+'_W.root'+' '+friendname)
