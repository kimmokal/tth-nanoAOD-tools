from tthAnalysis.NanoAODTools.postprocessing.modules.genParticleProducer import genAll
from tthAnalysis.NanoAODTools.postprocessing.modules.lepJetVarProducer import *
from tthAnalysis.NanoAODTools.postprocessing.modules.genHiggsDecayModeProducer import genHiggsDecayMode
from tthAnalysis.NanoAODTools.postprocessing.modules.tauIDLogProducer import tauIDLog_2016, tauIDLog_2017
from tthAnalysis.NanoAODTools.postprocessing.modules.absIsoProducer import absIso

from tthAnalysis.NanoAODTools.postprocessing.modules.jetSubstructureObservablesProducerHTTv2 import jetSubstructureObservablesHTTv2

############################################### we need these modules as a starting #######################
from tthAnalysis.NanoAODTools.postprocessing.modules.eventCountHistogramProducer import eventCountHistogram
from tthAnalysis.NanoAODTools.postprocessing.modules.countHistogramProducer import countHistogramAll_2016, countHistogramAll_2017
from tthAnalysis.NanoAODTools.postprocessing.modules.btagSFProducer_explicitBranchNames import btagSF_csvv2_2016, btagSF_cmva_2016, btagSF_csvv2_2017, btagSF_deep_2017
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight as puWeight_2016

from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import lepSF

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2017 as jetmetUncertainties17
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2016 as jetmetUncertainties16

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import jecUncert
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepdata
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepSIG
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepMC


from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepdata17
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepSIG17
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepMC17

from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_TOP import susy1lepTOPMC ,susy1lepTOPData
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_Signal import susy_1l_Sig16 ,susy_1l_Sig17
#correct the nJets 
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_nISR import susy_1l_nISR16
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_nISR import susy_1l_nISR17
# it has  an issue --> posted to HN --> now it works 16April 2018 --> Fixed 
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_gen import susy_1l_gen
# for Trig better to use Oleksii's for now 16April 2018 
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_trig import trigsusysusymod
#lepton SF already in base module
#from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_lepSF import susy_lepSF16, susy_lepSF17
###### Oleksii
from tthAnalysis.NanoAODTools.postprocessing.modules.susy_1l_triggers import susy_1l_Trigg2016,susy_1l_Trigg2017
from tthAnalysis.NanoAODTools.postprocessing.modules.susy_1l_filters import susy_1l_FiltersMC, susy_1l_FiltersData

from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puAutoWeight as puWeight_2017
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import muonScaleRes2016, muonScaleRes2017
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2016, jetmetUncertainties2017
