# tth-nanoAOD-tools
Postprocessing scripts to add branches specific to ttH analysis to nanoAOD Ntuples.

## Setup

```bash
# set up the CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh # !! or .csh
export SCRAM_ARCH=slc6_amd64_gcc630 # !! or setenv SCRAM_ARCH slc6_amd64_gcc630
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src/
cmsenv

# clone necessary repositories
git cms-merge-topic cms-nanoAOD:master
git checkout -b nanoAOD cms-nanoAOD/master
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git CMSSW_BASE/src/PhysicsTools/NanoAODTools
git clone https://github.com/ashrafkasem/tth-nanoAOD-tools.git $CMSSW_BASE/src/tthAnalysis/NanoAODTools
# compile the thing
cd $CMSSW_BASE/src
scram b -j 16
```

## Generate an example nanoAOD file

```bash
cd $CMSSW_BASE/src/PhysicsTools/NanoAOD/test
cmsRun nano_cfg.py                           # probably change the paths to the input files
```

## How to run the post-processing steps:

The output file name is created by appending a suffix to the basename of the input file, specified by the option `-s`.
example for running the post-processing interactively
```bash
nano_postproc.py -s _I2 -I tthAnalysis.NanoAODTools.postprocessing.tthModules susy1lepMC ./nanoAODs_HT/ inputfile.root 

```
use the option `-b keep-and-drop.txt` for keeping or dropping branches after proccessing
use the option `-J $CMSSW_BASE/src/tthAnalysis/NanoAODTools/data/JSONS/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt` in case you want to use the lumiMask

in order to run over all out susy1L MC, Data or signals 
```bash
cd tthAnalysis/NanoAODTools/batch
python getfilelist.py datesets.txt --out /path/to/out/dir --site DESY/Bari
```
i've listed all of them in text files. I wrote this script to run over all the list you provide and it will recognaize automatically the modules to run and it will differentiate between data/MC/Signal also it will recognize automatically the era (2016/2017)
each time the `getfilelist.py` will produce a `.sh` script called `submit_all_to_batch_HTC.sh` all what you need is to run run it and it will submitt the jobs to HTC either at DEST or Bari

you can do the same with crab if you want to skim the samples from a cluster rather then DESY/BARI by : 
```bash
cd tthAnalysis/NanoAODTools/crab
python prepare_crabdir.py dataList --out /path/to/out/dir
./submit_all_CRABs.sh
```

after you get the output you can merge each sample with `merge_skimout.py`
it can run either interactively or on batch as : 
```bash
/merge_skimout.py --indir /path/to/in/dir --outdir /path/to/out/dir [in case of batch mode use] -b --site DESY/Bari
``` 
removing -b and --site options means you will run the merging interactively 

If you want to add more modules then you must add the relevant import statements to `$CMSSW_BASE/src/tthAnalysis/NanoAODTools/python/postprocessing/tthModules.py` and recompile the NanoAODTools packages in `PhysicsTools` and `tthAnalysis` in order for the changes to take effect.

## Links

1. Official tool for post-processing the nanoAOD Ntuples: https://github.com/cms-nanoAOD/nanoAOD-tools
1. nanoAOD fork of CMSSW (complementary): https://github.com/cms-nanoAOD/cmssw/tree/master/PhysicsTools/NanoAOD
1. Our own nanoAOD-tools fork: https://github.com/HEP-KBFI/nanoAOD-tools
1. Our own fork of the CMSSW FW, based off the official nanoAOD fork: https://github.com/HEP-KBFI/cmssw

## ToDo 
1. push the plotter from my local repo to the github along with the skimmer
1. i've modifed all the proper scripts for the analysis and limit setting to run with HTC and nanoAOD, i just need to see how to push them to the same repo as the dir is contains alot of root files, i need to exclude them. 
