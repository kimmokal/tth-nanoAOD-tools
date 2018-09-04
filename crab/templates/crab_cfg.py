from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = '@REQUEST'
config.General.transferLogs=True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['@CRAB_SCRIPT','../../../PhysicsTools/NanoAODTools/scripts/haddnano.py'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True
config.section_("Data")
#config.Data.inputDataset = '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PUMoriond17_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM'
config.Data.inputDataset = '@INPUTDATASET'
#config.Data.inputDBS = 'global'
#
config.Data.splitting = '@SPLIT'
#config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'LumiBased'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.unitsPerJob = @UNITES
#config.Data.totalUnits = 2000
#config.Data.inputDBS='phys03'
config.Data.outLFNDirBase = '/store/user/amohamed/NanoPost/'
config.Data.publication = False 
config.Data.outputDatasetTag = '@TAG'
config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"

