from tthAnalysis.NanoAODTools.postprocessing.modules.genParticleProducer import genAll
from tthAnalysis.NanoAODTools.postprocessing.modules.lepJetVarProducer import lepJetVarBTagAll
from tthAnalysis.NanoAODTools.postprocessing.modules.genHiggsDecayModeProducer import genHiggsDecayMode
from tthAnalysis.NanoAODTools.postprocessing.modules.tauIDLogProducer import tauIDLog
from tthAnalysis.NanoAODTools.postprocessing.modules.absIsoProducer import absIso


from tthAnalysis.NanoAODTools.postprocessing.modules.jetSubstructureObservablesProducerHTTv2 import jetSubstructureObservablesHTTv2

############################################### we need these modules as a starting #######################
from tthAnalysis.NanoAODTools.postprocessing.modules.eventCountHistogramProducer import eventCountHistogram
from tthAnalysis.NanoAODTools.postprocessing.modules.countHistogramProducer import countHistogramAll
from tthAnalysis.NanoAODTools.postprocessing.modules.btagSFProducer_explicitBranchNames import btagSF_csvv2, btagSF_deep
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2017 as jetmetUncertainties17
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2016 as jetmetUncertainties16
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import jecUncert
######## ِAShraf ٍ
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepdata
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepSIG
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_base import susy1lepMC
# it has  an issue 
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_gen import genpartsusymod
# for Trig
from tthAnalysis.NanoAODTools.postprocessing.modules.susy1lep_trig import trigsusysusymod

###### Oleksii
from tthAnalysis.NanoAODTools.postprocessing.modules.susy_1l_triggers import susy_1l_Trigg
from tthAnalysis.NanoAODTools.postprocessing.modules.susy_1l_filters import susy_1l_FiltersMC, susy_1l_FiltersData