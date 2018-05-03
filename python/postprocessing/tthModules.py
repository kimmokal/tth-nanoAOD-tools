from tthAnalysis.NanoAODTools.postprocessing.modules.genParticleProducer import genAll
from tthAnalysis.NanoAODTools.postprocessing.modules.lepJetVarProducer import lepJetVarBTagAll
from tthAnalysis.NanoAODTools.postprocessing.modules.genHiggsDecayModeProducer import genHiggsDecayMode
from tthAnalysis.NanoAODTools.postprocessing.modules.tauIDLogProducer import tauIDLog
from tthAnalysis.NanoAODTools.postprocessing.modules.absIsoProducer import absIso


from tthAnalysis.NanoAODTools.postprocessing.modules.jetSubstructureObservablesProducerHTTv2 import jetSubstructureObservablesHTTv2

############################################### we need these modules as a starting #######################
from tthAnalysis.NanoAODTools.postprocessing.modules.eventCountHistogramProducer import eventCountHistogram
from tthAnalysis.NanoAODTools.postprocessing.modules.countHistogramProducer import countHistogramAll
from tthAnalysis.NanoAODTools.postprocessing.modules.btagSFProducer_explicitBranchNames import btagSF_csvv2_2016, btagSF_cmva_2016, btagSF_csvv2_2017, btagSF_deep_2017
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight as puWeight_2016
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puAutoWeight as puWeight_2017

from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import lepSF

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2017 as jetmetUncertainties17
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2016 as jetmetUncertainties16

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import jecUncert
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepdata
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepSIG
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepMC
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_TOP import susy1lepTOP
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_Signal import susy_1l_Sig
#correct the nJets 
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_nISR import susy_1l_nISR
# it has  an issue --> posted to HN --> now it works 16April 2018 --> Fixed 
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_gen import susy_1l_gen
# for Trig better to use Oleksii's for now 16April 2018 
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_trig import trigsusysusymod
###### Oleksii
from tthAnalysis.NanoAODTools.postprocessing.modules.susy_1l_triggers import susy_1l_Trigg
from tthAnalysis.NanoAODTools.postprocessing.modules.susy_1l_filters import susy_1l_FiltersMC, susy_1l_FiltersData
###### temp 
#from tthAnalysis.NanoAODTools.postprocessing.modules.new import susy1lepdata
