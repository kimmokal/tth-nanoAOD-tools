/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This example shows the training of signal with three different backgrounds
/// Then in the application a tree is created with all signal and background
/// events where the true class ID and the three classifier outputs are added
/// finally with the application tree, the significance is maximized with the
/// help of the TMVA genetic algrorithm.
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Exectuable: TMVAGAexample
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker


#include <iostream> // Stream declarations
#include <vector>
#include <limits>

#include "TChain.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"

#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"//required to load dataset
#include "TMVA/Reader.h"

using namespace std;

using namespace TMVA;

// ----------------------------------------------------------------------------------------------
// Training
// ----------------------------------------------------------------------------------------------
//
void TMVA_1lep(){
		std::string factoryOptions( "!V:!Silent:Transformations=I;D;P;G,D:AnalysisType=Classification" );
		string inputFolder = "/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/Cut-NANO/all_trees_2016_Cut/";
		//TString DYJetsToLL_M_50_HT_100to200 = inputFolder+"DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		//TString DYJetsToLL_M_50_HT_100to200_1 	= inputFolder+"DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString DYJetsToLL_M_50_HT_1200to2500	= inputFolder+"DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		//TString DYJetsToLL_M_50_HT_200to400		= inputFolder+"DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		//TString DYJetsToLL_M_50_HT_200to400_1	= inputFolder+"DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString DYJetsToLL_M_50_HT_2500toInf	= inputFolder+"DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString DYJetsToLL_M_50_HT_400to600		= inputFolder+"DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString DYJetsToLL_M_50_HT_400to600_1	= inputFolder+"DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString DYJetsToLL_M_50_HT_600to800		= inputFolder+"DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString DYJetsToLL_M_50_HT_800to1200	= inputFolder+"DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString QCD_HT1000to1500				= inputFolder+"QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString QCD_HT1000to1500_1				= inputFolder+"QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString QCD_HT1500to2000				= inputFolder+"QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString QCD_HT1500to2000_1				= inputFolder+"QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString QCD_HT2000toInf					= inputFolder+"QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString QCD_HT2000toInf_1				= inputFolder+"QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString QCD_HT500to700					= inputFolder+"QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString QCD_HT500to700_1				= inputFolder+"QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString QCD_HT700to1000					= inputFolder+"QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString QCD_HT700to1000_1				= inputFolder+"QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString ST_s_channel_4f_leptonDecays	= inputFolder+"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1RunIISummer16NanoAOD.root";
		TString ST_t_channel_antitop_4f_inclusiveDecays	= inputFolder+"ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1RunIISummer16NanoAOD.root";
		TString ST_t_channel_top_4f_inclusiveDecays		= inputFolder+"ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1RunIISummer16NanoAOD.root";
		TString ST_tW_top_5f_NoFullyHadronicDecays		= inputFolder+"ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1RunIISummer16NanoAOD.root";
		TString TTJets_DiLept					= inputFolder+"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString TTJets_DiLept_ext1				= inputFolder+"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString TTJets_HT_1200to2500			= inputFolder+"TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString TTJets_HT_2500toInf				= inputFolder+"TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString TTJets_HT_600to800				= inputFolder+"TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString TTJets_HT_800to1200				= inputFolder+"TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString TTJets_SingleLeptFromT			= inputFolder+"TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString TTJets_SingleLeptFromTbar		= inputFolder+"TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString TTJets_SingleLeptFromTbar_1		= inputFolder+"TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString TTTo2L2Nu						= inputFolder+"TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8RunIISummer16NanoAOD.root";
		TString TTWJetsToLNu_1					= inputFolder+"TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8RunIISummer16NanoAODext1.root";
		TString TTWJetsToLNu_2					= inputFolder+"TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8RunIISummer16NanoAODext2.root";
		TString TTZToLLNuNu_1					= inputFolder+"TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8RunIISummer16NanoAOD.root";
		TString TTZToLLNuNu_2					= inputFolder+"TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8RunIISummer16NanoAODext1.root";
		TString TTZToLLNuNu_3					= inputFolder+"TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8RunIISummer16NanoAODext2.root";
		TString TTZToQQ							= inputFolder+"TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8RunIISummer16NanoAOD.root";
		TString TT								= inputFolder+"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_HT_1200To2500		= inputFolder+"WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_HT_2500ToInf			= inputFolder+"WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_HT_2500ToInf_1		= inputFolder+"WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString WJetsToLNu_HT_400To600			= inputFolder+"WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_HT_400To600_1		= inputFolder+"WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString WJetsToLNu_HT_600To800			= inputFolder+"WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_HT_600To800_1		= inputFolder+"WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString WJetsToLNu_HT_800To1200			= inputFolder+"WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_HT_800To1200_1		= inputFolder+"WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext1.root";
		TString WJetsToLNu_1					= inputFolder+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_2					= inputFolder+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8RunIISummer16NanoAODext2.root";
		TString WJetsToLNu_madgraphMLM			= inputFolder+"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		TString WJetsToLNu_madgraphMLM_1		= inputFolder+"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAODext2.root";
		TString WWTo2L2Nu						= inputFolder+"WWTo2L2Nu_13TeV-powhegRunIISummer16NanoAOD.root";
		TString WWToLNuQQ						= inputFolder+"WWToLNuQQ_13TeV-powhegRunIISummer16NanoAOD.root";
		TString WWToLNuQQ_1						= inputFolder+"WWToLNuQQ_13TeV-powhegRunIISummer16NanoAODext1.root";
		TString WZTo1L1Nu2Q_2					= inputFolder+"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8RunIISummer16NanoAOD.root";
		TString WZTo2L2Q						= inputFolder+"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8RunIISummer16NanoAOD.root";
		TString ZZTo2L2Nu						= inputFolder+"ZZTo2L2Nu_13TeV_powheg_pythia8RunIISummer16NanoAOD.root";
		TString ZZTo2L2Q						= inputFolder+"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8RunIISummer16NanoAOD.root";
		TString ZZTo2L2Q_1						= inputFolder+"ZZTo2L2Q_13TeV_powheg_pythia8RunIISummer16NanoAOD.root";
		
		TString sfname = inputFolder+"SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8RunIISummer16NanoAOD.root";
		
		TFile *sinput(0);
		//TFile *input_DYJetsToLL_M_50_HT_100to200(0);
		//TFile *input_DYJetsToLL_M_50_HT_100to200_1(0);
		TFile *input_DYJetsToLL_M_50_HT_1200to2500(0);
		//TFile *input_DYJetsToLL_M_50_HT_200to400(0);
		//TFile *input_DYJetsToLL_M_50_HT_200to400_1(0);
		TFile *input_DYJetsToLL_M_50_HT_2500toInf(0);	
		TFile *input_DYJetsToLL_M_50_HT_400to600(0);
		TFile *input_DYJetsToLL_M_50_HT_400to600_1(0);
		TFile *input_DYJetsToLL_M_50_HT_600to800(0);
		TFile *input_DYJetsToLL_M_50_HT_800to1200(0);
		TFile *input_QCD_HT1000to1500(0);
		TFile *input_QCD_HT1000to1500_1(0);
		TFile *input_QCD_HT1500to2000(0);
		TFile *input_QCD_HT1500to2000_1(0);
		TFile *input_QCD_HT2000toInf(0);
		TFile *input_QCD_HT2000toInf_1(0);
		TFile *input_QCD_HT500to700(0);	
		TFile *input_QCD_HT500to700_1(0);
		TFile *input_QCD_HT700to1000(0);	
		TFile *input_QCD_HT700to1000_1(0);	
		TFile *input_ST_s_channel_4f_leptonDecays(0);
		TFile *input_ST_t_channel_antitop_4f_inclusiveDecays(0);
		TFile *input_ST_t_channel_top_4f_inclusiveDecays(0);
		TFile *input_ST_tW_top_5f_NoFullyHadronicDecays(0);
		TFile *input_TTJets_DiLept(0);
		TFile *input_TTJets_DiLept_ext1(0);
		TFile *input_TTJets_HT_1200to2500(0);
		TFile *input_TTJets_HT_2500toInf(0);	
		TFile *input_TTJets_HT_600to800(0);	
		TFile *input_TTJets_HT_800to1200(0);	
		TFile *input_TTJets_SingleLeptFromT(0);	
		TFile *input_TTJets_SingleLeptFromTbar(0);	
		TFile *input_TTJets_SingleLeptFromTbar_1(0);	
		TFile *input_TTTo2L2Nu(0);	
		TFile *input_TTWJetsToLNu_1(0);	
		TFile *input_TTWJetsToLNu_2(0);
		TFile *input_TTZToLLNuNu_1(0);	
		TFile *input_TTZToLLNuNu_2(0);	
		TFile *input_TTZToLLNuNu_3(0);	
		TFile *input_TTZToQQ(0);	
		TFile *input_TT(0);	
		TFile *input_WJetsToLNu_HT_1200To2500(0);
		TFile *input_WJetsToLNu_HT_2500ToInf(0);	
		TFile *input_WJetsToLNu_HT_2500ToInf_1(0);	
		TFile *input_WJetsToLNu_HT_400To600(0);	
		TFile *input_WJetsToLNu_HT_400To600_1	(0);
		TFile *input_WJetsToLNu_HT_600To800(0);	
		TFile *input_WJetsToLNu_HT_600To800_1(0);
		TFile *input_WJetsToLNu_HT_800To1200(0);	
		TFile *input_WJetsToLNu_HT_800To1200_1(0);	
		TFile *input_WJetsToLNu_1(0);
		TFile *input_WJetsToLNu_2(0);
		TFile *input_WJetsToLNu_madgraphMLM(0);	
		TFile *input_WJetsToLNu_madgraphMLM_1(0);
		TFile *input_WWTo2L2Nu(0);	
		TFile *input_WWToLNuQQ(0);	
		TFile *input_WWToLNuQQ_1(0);	
		TFile *input_WZTo1L1Nu2Q_2(0);	
		TFile *input_WZTo2L2Q(0);
		TFile *input_ZZTo2L2Nu(0);	
		TFile *input_ZZTo2L2Q(0);
		TFile *input_ZZTo2L2Q_1(0);	
		
		sinput = TFile::Open( sfname );
		//input_DYJetsToLL_M_50_HT_100to200=TFile::Open(DYJetsToLL_M_50_HT_100to200);
		//input_DYJetsToLL_M_50_HT_100to200_1=TFile::Open(DYJetsToLL_M_50_HT_100to200_1);
		input_DYJetsToLL_M_50_HT_1200to2500=TFile::Open(DYJetsToLL_M_50_HT_1200to2500);
		//input_DYJetsToLL_M_50_HT_200to400=TFile::Open(DYJetsToLL_M_50_HT_200to400);
		//input_DYJetsToLL_M_50_HT_200to400_1=TFile::Open(DYJetsToLL_M_50_HT_200to400_1);
		input_DYJetsToLL_M_50_HT_2500toInf=TFile::Open(DYJetsToLL_M_50_HT_2500toInf);
		input_DYJetsToLL_M_50_HT_400to600=TFile::Open(DYJetsToLL_M_50_HT_400to600);
		input_DYJetsToLL_M_50_HT_400to600_1=TFile::Open(DYJetsToLL_M_50_HT_400to600_1);
		input_DYJetsToLL_M_50_HT_600to800=TFile::Open(DYJetsToLL_M_50_HT_600to800);
		input_DYJetsToLL_M_50_HT_800to1200=TFile::Open(DYJetsToLL_M_50_HT_800to1200);
		input_QCD_HT1000to1500=TFile::Open(QCD_HT1000to1500);
		input_QCD_HT1000to1500_1=TFile::Open(QCD_HT1000to1500_1);
		input_QCD_HT1500to2000=TFile::Open(QCD_HT1500to2000);
		input_QCD_HT1500to2000_1=TFile::Open(QCD_HT1500to2000_1);
		input_QCD_HT2000toInf=TFile::Open(QCD_HT2000toInf);
		input_QCD_HT2000toInf_1=TFile::Open(QCD_HT2000toInf_1);
		input_QCD_HT500to700=TFile::Open(QCD_HT500to700);
		input_QCD_HT500to700_1=TFile::Open(QCD_HT500to700_1);
		input_QCD_HT700to1000=TFile::Open(QCD_HT700to1000);
		input_QCD_HT700to1000_1=TFile::Open(QCD_HT700to1000_1);
		input_ST_s_channel_4f_leptonDecays=TFile::Open(ST_s_channel_4f_leptonDecays);
		input_ST_t_channel_antitop_4f_inclusiveDecays=TFile::Open(ST_t_channel_antitop_4f_inclusiveDecays);
		input_ST_t_channel_top_4f_inclusiveDecays=TFile::Open(ST_t_channel_top_4f_inclusiveDecays);
		input_ST_tW_top_5f_NoFullyHadronicDecays=TFile::Open(ST_tW_top_5f_NoFullyHadronicDecays);
		input_TTJets_DiLept=TFile::Open(TTJets_DiLept);
		input_TTJets_DiLept_ext1=TFile::Open(TTJets_DiLept_ext1);
		input_TTJets_HT_1200to2500=TFile::Open(TTJets_HT_1200to2500);
		input_TTJets_HT_2500toInf=TFile::Open(TTJets_HT_2500toInf);
		input_TTJets_HT_600to800=TFile::Open(TTJets_HT_600to800);
		input_TTJets_HT_800to1200=TFile::Open(TTJets_HT_800to1200);
		input_TTJets_SingleLeptFromT=TFile::Open(TTJets_SingleLeptFromT);
		input_TTJets_SingleLeptFromTbar=TFile::Open(TTJets_SingleLeptFromTbar);
		input_TTJets_SingleLeptFromTbar_1=TFile::Open(TTJets_SingleLeptFromTbar_1);
		input_TTTo2L2Nu=TFile::Open(TTTo2L2Nu);
		input_TTWJetsToLNu_1=TFile::Open(TTWJetsToLNu_1);
		input_TTWJetsToLNu_2=TFile::Open(TTWJetsToLNu_2);
		input_TTZToLLNuNu_1=TFile::Open(TTZToLLNuNu_1);
		input_TTZToLLNuNu_2=TFile::Open(TTZToLLNuNu_2);
		input_TTZToLLNuNu_3=TFile::Open(TTZToLLNuNu_3);
		input_TTZToQQ=TFile::Open(TTZToQQ);
		input_TT=TFile::Open(TT);
		input_WJetsToLNu_HT_1200To2500=TFile::Open(WJetsToLNu_HT_1200To2500);
		input_WJetsToLNu_HT_2500ToInf=TFile::Open(WJetsToLNu_HT_2500ToInf);
		input_WJetsToLNu_HT_2500ToInf_1=TFile::Open(WJetsToLNu_HT_2500ToInf_1);
		input_WJetsToLNu_HT_400To600=TFile::Open(WJetsToLNu_HT_400To600);
		input_WJetsToLNu_HT_400To600_1=TFile::Open(WJetsToLNu_HT_400To600_1);
		input_WJetsToLNu_HT_600To800=TFile::Open(WJetsToLNu_HT_600To800);
		input_WJetsToLNu_HT_600To800_1=TFile::Open(WJetsToLNu_HT_600To800_1);
		input_WJetsToLNu_HT_800To1200=TFile::Open(WJetsToLNu_HT_800To1200);
		input_WJetsToLNu_HT_800To1200_1=TFile::Open(WJetsToLNu_HT_800To1200_1);
		input_WJetsToLNu_1=TFile::Open(WJetsToLNu_1);
		input_WJetsToLNu_2=TFile::Open(WJetsToLNu_2);
		input_WJetsToLNu_madgraphMLM=TFile::Open(WJetsToLNu_madgraphMLM);
		input_WJetsToLNu_madgraphMLM_1=TFile::Open(WJetsToLNu_madgraphMLM_1);
		input_WWTo2L2Nu=TFile::Open(WWTo2L2Nu);
		input_WWToLNuQQ=TFile::Open(WWToLNuQQ);
		input_WWToLNuQQ_1=TFile::Open(WWToLNuQQ_1);
		input_WZTo1L1Nu2Q_2=TFile::Open(WZTo1L1Nu2Q_2);
		input_WZTo2L2Q=TFile::Open(WZTo2L2Q);
		input_ZZTo2L2Nu=TFile::Open(ZZTo2L2Nu);
		input_ZZTo2L2Q=TFile::Open(ZZTo2L2Q);
		input_ZZTo2L2Q_1=TFile::Open(ZZTo2L2Q_1);   
   // global event weights per tree (see below for setting event-wise weights)
		TTree *signal = (TTree*)sinput->Get("Events");
		//TTree *tree_DYJetsToLL_M_50_HT_100to200=(TTree*)input_DYJetsToLL_M_50_HT_100to200->Get("Events");
		//TTree *tree_DYJetsToLL_M_50_HT_100to200_1 =(TTree*)input_DYJetsToLL_M_50_HT_100to200_1->Get("Events");
		TTree *tree_DYJetsToLL_M_50_HT_1200to2500 =(TTree*)input_DYJetsToLL_M_50_HT_1200to2500->Get("Events");
		//TTree *tree_DYJetsToLL_M_50_HT_200to400 =(TTree*)input_DYJetsToLL_M_50_HT_200to400->Get("Events");
		//TTree *tree_DYJetsToLL_M_50_HT_200to400_1 =(TTree*)input_DYJetsToLL_M_50_HT_200to400_1->Get("Events");
		TTree *tree_DYJetsToLL_M_50_HT_2500toInf =(TTree*)input_DYJetsToLL_M_50_HT_2500toInf->Get("Events");
		TTree *tree_DYJetsToLL_M_50_HT_400to600 =(TTree*)input_DYJetsToLL_M_50_HT_400to600->Get("Events");
		TTree *tree_DYJetsToLL_M_50_HT_400to600_1 =(TTree*)input_DYJetsToLL_M_50_HT_400to600_1->Get("Events");
		TTree *tree_DYJetsToLL_M_50_HT_600to800 =(TTree*)input_DYJetsToLL_M_50_HT_600to800->Get("Events");
		TTree *tree_DYJetsToLL_M_50_HT_800to1200 =(TTree*)input_DYJetsToLL_M_50_HT_800to1200->Get("Events");
		TTree *tree_QCD_HT1000to1500 =(TTree*)input_QCD_HT1000to1500->Get("Events");
		TTree *tree_QCD_HT1000to1500_1 =(TTree*)input_QCD_HT1000to1500_1->Get("Events");
		TTree *tree_QCD_HT1500to2000 =(TTree*)input_QCD_HT1500to2000->Get("Events");
		TTree *tree_QCD_HT1500to2000_1 =(TTree*)input_QCD_HT1500to2000_1->Get("Events");
		TTree *tree_QCD_HT2000toInf =(TTree*)input_QCD_HT2000toInf->Get("Events");
		TTree *tree_QCD_HT2000toInf_1 =(TTree*)input_QCD_HT2000toInf_1->Get("Events");
		TTree *tree_QCD_HT500to700 =(TTree*)input_QCD_HT500to700->Get("Events");
		TTree *tree_QCD_HT500to700_1 =(TTree*)input_QCD_HT500to700_1->Get("Events");
		TTree *tree_QCD_HT700to1000 =(TTree*)input_QCD_HT700to1000->Get("Events");
		TTree *tree_QCD_HT700to1000_1 =(TTree*)input_QCD_HT700to1000_1->Get("Events");
		TTree *tree_ST_s_channel_4f_leptonDecays =(TTree*)input_ST_s_channel_4f_leptonDecays->Get("Events");
		TTree *tree_ST_t_channel_antitop_4f_inclusiveDecays =(TTree*)input_ST_t_channel_antitop_4f_inclusiveDecays->Get("Events");
		TTree *tree_ST_t_channel_top_4f_inclusiveDecays =(TTree*)input_ST_t_channel_top_4f_inclusiveDecays->Get("Events");
		TTree *tree_ST_tW_top_5f_NoFullyHadronicDecays =(TTree*)input_ST_tW_top_5f_NoFullyHadronicDecays->Get("Events");
		TTree *tree_TTJets_DiLept =(TTree*)input_TTJets_DiLept->Get("Events");
		TTree *tree_TTJets_DiLept_ext1 =(TTree*)input_TTJets_DiLept_ext1->Get("Events");
		TTree *tree_TTJets_HT_1200to2500 =(TTree*)input_TTJets_HT_1200to2500->Get("Events");
		TTree *tree_TTJets_HT_2500toInf =(TTree*)input_TTJets_HT_2500toInf->Get("Events");
		TTree *tree_TTJets_HT_600to800 =(TTree*)input_TTJets_HT_600to800->Get("Events");
		TTree *tree_TTJets_HT_800to1200 =(TTree*)input_TTJets_HT_800to1200->Get("Events");
		TTree *tree_TTJets_SingleLeptFromT =(TTree*)input_TTJets_SingleLeptFromT->Get("Events");
		TTree *tree_TTJets_SingleLeptFromTbar =(TTree*)input_TTJets_SingleLeptFromTbar->Get("Events");
		TTree *tree_TTJets_SingleLeptFromTbar_1 =(TTree*)input_TTJets_SingleLeptFromTbar->Get("Events");
		TTree *tree_TTTo2L2Nu =(TTree*)input_TTTo2L2Nu->Get("Events");
		TTree *tree_TTWJetsToLNu_1 =(TTree*)input_TTWJetsToLNu_1->Get("Events");
		TTree *tree_TTWJetsToLNu_2 =(TTree*)input_TTWJetsToLNu_2->Get("Events");
		TTree *tree_TTZToLLNuNu_1 =(TTree*)input_TTZToLLNuNu_1->Get("Events");
		TTree *tree_TTZToLLNuNu_2 =(TTree*)input_TTZToLLNuNu_2->Get("Events");
		TTree *tree_TTZToLLNuNu_3 =(TTree*)input_TTZToLLNuNu_3->Get("Events");
		TTree *tree_TTZToQQ =(TTree*)input_TTZToQQ->Get("Events");
		TTree *tree_TT =(TTree*)input_TT->Get("Events");
		TTree *tree_WJetsToLNu_HT_1200To2500=(TTree*)input_WJetsToLNu_HT_1200To2500->Get("Events");
		TTree *tree_WJetsToLNu_HT_2500ToInf=(TTree*)input_WJetsToLNu_HT_2500ToInf->Get("Events");
		TTree *tree_WJetsToLNu_HT_2500ToInf_1 =(TTree*)input_WJetsToLNu_HT_2500ToInf_1->Get("Events");
		TTree *tree_WJetsToLNu_HT_400To600 =(TTree*)input_WJetsToLNu_HT_400To600->Get("Events");
		TTree *tree_WJetsToLNu_HT_400To600_1=(TTree*)input_WJetsToLNu_HT_400To600_1->Get("Events");
		TTree *tree_WJetsToLNu_HT_600To800 =(TTree*)input_WJetsToLNu_HT_600To800->Get("Events");
		TTree *tree_WJetsToLNu_HT_600To800_1 =(TTree*)input_WJetsToLNu_HT_600To800_1->Get("Events");
		TTree *tree_WJetsToLNu_HT_800To1200 =(TTree*)input_WJetsToLNu_HT_800To1200->Get("Events");
		TTree *tree_WJetsToLNu_HT_800To1200_1 =(TTree*)input_WJetsToLNu_HT_800To1200_1->Get("Events");
		TTree *tree_WJetsToLNu_1 =(TTree*)input_WJetsToLNu_1->Get("Events");
		TTree *tree_WJetsToLNu_2 =(TTree*)input_WJetsToLNu_2->Get("Events");
		TTree *tree_WJetsToLNu_madgraphMLM =(TTree*)input_WJetsToLNu_madgraphMLM->Get("Events");
		TTree *tree_WJetsToLNu_madgraphMLM_1 =(TTree*)input_WJetsToLNu_madgraphMLM_1->Get("Events");
		TTree *tree_WWTo2L2Nu =(TTree*)input_WWTo2L2Nu->Get("Events");
		TTree *tree_WWToLNuQQ =(TTree*)input_WWToLNuQQ->Get("Events");
		TTree *tree_WWToLNuQQ_1 =(TTree*)input_WWToLNuQQ_1->Get("Events");
		TTree *tree_WZTo1L1Nu2Q_2 =(TTree*)input_WZTo1L1Nu2Q_2->Get("Events");
		TTree *tree_WZTo2L2Q =(TTree*)input_WZTo2L2Q->Get("Events");
		TTree *tree_ZZTo2L2Nu =(TTree*)input_ZZTo2L2Nu->Get("Events");
		TTree *tree_ZZTo2L2Q =(TTree*)input_ZZTo2L2Q->Get("Events");
		TTree *tree_ZZTo2L2Q_1 =(TTree*)input_ZZTo2L2Q_1->Get("Events");

   
   Double_t signalWeight      = 1.0;
   Double_t backgroundWeight = 1.0;
   Double_t lumi=35.9 ;
   // Create a new root output file.
   TString outfileName( "TMVASignalBackground.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   // ____________
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, factoryOptions );
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("datasetBkg");
   
	dataloader->AddVariable( "Lep_pt"		,"Lep_pt"		, "", 'F' );
	dataloader->AddVariable( "Lep_eta"		,"Lep_eta"		, "", 'F' );
	dataloader->AddVariable( "Lep_miniIso"	,"Lep_miniIso"	, "", 'F' );
	dataloader->AddVariable( "MT"			,"MT"			, "", 'F' );
	dataloader->AddVariable( "Met"			,"Met"			, "", 'F' );
	dataloader->AddVariable( "Jet1_pt"		,"Jet1_pt"		, "", 'F' );
	dataloader->AddVariable( "Jet2_pt"		,"Jet2_pt"		, "", 'F' );
	dataloader->AddVariable( "Jet1_eta"		,"Jet1_eta"		, "", 'F' );
	dataloader->AddVariable( "Jet2_eta"		,"Jet2_eta"		, "", 'F' );
	dataloader->AddVariable( "HT"			, "HT"			, "", 'F' );
	dataloader->AddVariable( "LT"			, "LT"			, "", 'F' );
	dataloader->AddVariable( "nJet"			, "nJet"		, "", 'F' );
	dataloader->AddVariable( "nBJet"		, "nBJet"		, "", 'F' );
	dataloader->AddVariable( "dPhi"			, "dPhi"		, "", 'F' );
	
	dataloader->AddSpectator("Jet_btagSF_csvv2","Jet_btagSF_csvv2","", 'D');
	dataloader->AddSpectator("puWeight","puWeight","", 'D');
	dataloader->AddSpectator("lepSF","lepSF","", 'D');
	//dataloader->AddSpectator("nISRweight","nISRweight","", 'D');
	//dataloader->AddSpectator("Xsec","Xsec","", 'D');
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables

   //dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );
	
	TH1F *weight_signal = (TH1F*)sinput->Get("Count");
	TH1F *weight_DYJetsToLL_M_50_HT_1200to2500=(TH1F*)input_DYJetsToLL_M_50_HT_1200to2500->Get("Count");
	TH1F *weight_DYJetsToLL_M_50_HT_2500toInf=(TH1F*)input_DYJetsToLL_M_50_HT_2500toInf->Get("Count");
	TH1F *weight_DYJetsToLL_M_50_HT_400to600=(TH1F*)input_DYJetsToLL_M_50_HT_400to600->Get("Count");
	TH1F *weight_DYJetsToLL_M_50_HT_400to600_1=(TH1F*)input_DYJetsToLL_M_50_HT_400to600_1->Get("Count");
	TH1F *weight_DYJetsToLL_M_50_HT_600to800=(TH1F*)input_DYJetsToLL_M_50_HT_600to800->Get("Count");
	TH1F *weight_DYJetsToLL_M_50_HT_800to1200=(TH1F*)input_DYJetsToLL_M_50_HT_800to1200->Get("Count");
	TH1F *weight_QCD_HT1000to1500=(TH1F*)input_QCD_HT1000to1500->Get("Count");
	TH1F *weight_QCD_HT1000to1500_1=(TH1F*)input_QCD_HT1000to1500_1->Get("Count");
	TH1F *weight_QCD_HT1500to2000=(TH1F*)input_QCD_HT1500to2000->Get("Count");
	TH1F *weight_QCD_HT1500to2000_1=(TH1F*)input_QCD_HT1500to2000_1->Get("Count");
	TH1F *weight_QCD_HT2000toInf=(TH1F*)input_QCD_HT2000toInf->Get("Count");
	TH1F *weight_QCD_HT2000toInf_1=(TH1F*)input_QCD_HT2000toInf_1->Get("Count");
	TH1F *weight_QCD_HT500to700=(TH1F*)input_QCD_HT500to700->Get("Count");
	TH1F *weight_QCD_HT500to700_1=(TH1F*)input_QCD_HT500to700_1->Get("Count");
	TH1F *weight_QCD_HT700to1000=(TH1F*)input_QCD_HT700to1000->Get("Count");
	TH1F *weight_QCD_HT700to1000_1=(TH1F*)input_QCD_HT700to1000_1->Get("Count");
	TH1F *weight_ST_s_channel_4f_leptonDecays=(TH1F*)input_ST_s_channel_4f_leptonDecays->Get("Count");
	TH1F *weight_ST_t_channel_antitop_4f_inclusiveDecays=(TH1F*)input_ST_t_channel_antitop_4f_inclusiveDecays->Get("Count");
	TH1F *weight_ST_t_channel_top_4f_inclusiveDecays=(TH1F*)input_ST_t_channel_top_4f_inclusiveDecays->Get("Count");
	TH1F *weight_ST_tW_top_5f_NoFullyHadronicDecays=(TH1F*)input_ST_tW_top_5f_NoFullyHadronicDecays->Get("Count");
	TH1F *weight_TTJets_DiLept=(TH1F*)input_TTJets_DiLept->Get("Count");
	TH1F *weight_TTJets_DiLept_ext1=(TH1F*)input_TTJets_DiLept_ext1->Get("Count");
	TH1F *weight_TTJets_HT_1200to2500=(TH1F*)input_TTJets_HT_1200to2500->Get("Count");
	TH1F *weight_TTJets_HT_2500toInf=(TH1F*)input_TTJets_HT_2500toInf->Get("Count");
	TH1F *weight_TTJets_HT_600to800=(TH1F*)input_TTJets_HT_600to800->Get("Count");
	TH1F *weight_TTJets_HT_800to1200=(TH1F*)input_TTJets_HT_800to1200->Get("Count");
	TH1F *weight_TTJets_SingleLeptFromT=(TH1F*)input_TTJets_SingleLeptFromT->Get("Count");
	TH1F *weight_TTJets_SingleLeptFromTbar=(TH1F*)input_TTJets_SingleLeptFromTbar->Get("Count");
	TH1F *weight_TTJets_SingleLeptFromTbar_1=(TH1F*)input_TTJets_SingleLeptFromTbar_1->Get("Count");
	TH1F *weight_TTTo2L2Nu=(TH1F*)input_TTTo2L2Nu->Get("Count");
	TH1F *weight_TTWJetsToLNu_1=(TH1F*)input_TTWJetsToLNu_1->Get("Count");
	TH1F *weight_TTWJetsToLNu_2=(TH1F*)input_TTWJetsToLNu_2->Get("Count");
	TH1F *weight_TTZToLLNuNu_1=(TH1F*)input_TTZToLLNuNu_1->Get("Count");
	TH1F *weight_TTZToLLNuNu_2=(TH1F*)input_TTZToLLNuNu_2->Get("Count");
	TH1F *weight_TTZToLLNuNu_3=(TH1F*)input_TTZToLLNuNu_3->Get("Count");
	TH1F *weight_TTZToQQ=(TH1F*)input_TTZToQQ->Get("Count");
	TH1F *weight_TT=(TH1F*)input_TT->Get("Count");
	TH1F *weight_WJetsToLNu_HT_1200To2500=(TH1F*)input_WJetsToLNu_HT_1200To2500->Get("Count");
	TH1F *weight_WJetsToLNu_HT_2500ToInf=(TH1F*)input_WJetsToLNu_HT_2500ToInf->Get("Count");
	TH1F *weight_WJetsToLNu_HT_2500ToInf_1=(TH1F*)input_WJetsToLNu_HT_2500ToInf_1->Get("Count");
	TH1F *weight_WJetsToLNu_HT_400To600=(TH1F*)input_WJetsToLNu_HT_400To600->Get("Count");
	TH1F *weight_WJetsToLNu_HT_400To600_1=(TH1F*)input_WJetsToLNu_HT_400To600_1->Get("Count");
	TH1F *weight_WJetsToLNu_HT_600To800=(TH1F*)input_WJetsToLNu_HT_600To800->Get("Count");
	TH1F *weight_WJetsToLNu_HT_600To800_1=(TH1F*)input_WJetsToLNu_HT_600To800_1->Get("Count");
	TH1F *weight_WJetsToLNu_HT_800To1200=(TH1F*)input_WJetsToLNu_HT_800To1200->Get("Count");
	TH1F *weight_WJetsToLNu_HT_800To1200_1=(TH1F*)input_WJetsToLNu_HT_800To1200_1->Get("Count");
	TH1F *weight_WJetsToLNu_1=(TH1F*)input_WJetsToLNu_1->Get("Count");
	TH1F *weight_WJetsToLNu_2=(TH1F*)input_WJetsToLNu_2->Get("Count");
	TH1F *weight_WJetsToLNu_madgraphMLM=(TH1F*)input_WJetsToLNu_madgraphMLM->Get("Count");
	TH1F *weight_WJetsToLNu_madgraphMLM_1=(TH1F*)input_WJetsToLNu_madgraphMLM_1->Get("Count");
	TH1F *weight_WWTo2L2Nu=(TH1F*)input_WWTo2L2Nu->Get("Count");
	TH1F *weight_WWToLNuQQ=(TH1F*)input_WWToLNuQQ->Get("Count");
	TH1F *weight_WWToLNuQQ_1=(TH1F*)input_WWToLNuQQ_1->Get("Count");
	TH1F *weight_WZTo1L1Nu2Q_2=(TH1F*)input_WZTo1L1Nu2Q_2->Get("Count");
	TH1F *weight_WZTo2L2Q=(TH1F*)input_WZTo2L2Q->Get("Count");
	TH1F *weight_ZZTo2L2Nu=(TH1F*)input_ZZTo2L2Nu->Get("Count");
	TH1F *weight_ZZTo2L2Q=(TH1F*)input_ZZTo2L2Q->Get("Count");
	TH1F *weight_ZZTo2L2Q_1=(TH1F*)input_ZZTo2L2Q_1->Get("Count");

	dataloader->AddSignalTree    ( signal,     weight_signal->GetBinContent(1));
	/*dataloader->AddBackgroundTree(tree_DYJetsToLL_M_50_HT_1200to2500	, 0.1514*1.23*1*lumi*1000/weight_DYJetsToLL_M_50_HT_1200to2500->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M_50_HT_2500toInf	, 0.003565*1.23*1*lumi*1000/weight_DYJetsToLL_M_50_HT_2500toInf->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M_50_HT_400to600,   5.678*1.23*1*lumi*1000/weight_DYJetsToLL_M_50_HT_400to600->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M_50_HT_400to600_1	, 5.678*1.23*1*lumi*1000/weight_DYJetsToLL_M_50_HT_400to600_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M_50_HT_600to800	,  1.367*1.23*1*lumi*1000/weight_DYJetsToLL_M_50_HT_600to800->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M_50_HT_800to1200, 0.6304*1.23*1*lumi*1000/weight_DYJetsToLL_M_50_HT_800to1200->GetBinContent(1));
	
	dataloader->AddBackgroundTree(tree_QCD_HT1000to1500, 1206.0*1*lumi*1000/weight_QCD_HT1000to1500->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT1000to1500_1, 1206.0*1*lumi*1000/weight_QCD_HT1000to1500_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT1500to2000, 120.4*1*lumi*1000/weight_QCD_HT1500to2000->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT1500to2000_1, 120.4*1*lumi*1000/weight_QCD_HT1500to2000_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT2000toInf, 25.25*1*lumi*1000/weight_QCD_HT2000toInf->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT2000toInf_1, 25.25*1*lumi*1000/weight_QCD_HT2000toInf_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT500to700, 31630.0*1*lumi*1000/weight_QCD_HT500to700->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT500to700_1, 31630.0*1*lumi*1000/weight_QCD_HT500to700_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT700to1000, 6802.0*1*lumi*1000/weight_QCD_HT700to1000->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_QCD_HT700to1000_1,6802.0*1*lumi*1000/weight_QCD_HT700to1000_1->GetBinContent(1));
	
	dataloader->AddBackgroundTree(tree_ST_s_channel_4f_leptonDecays, 3.365*1*lumi*1000/weight_ST_s_channel_4f_leptonDecays->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_ST_t_channel_antitop_4f_inclusiveDecays, 80.95*1*lumi*1000/weight_ST_t_channel_antitop_4f_inclusiveDecays->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_ST_t_channel_top_4f_inclusiveDecays, 80.95*1*lumi*1000/weight_ST_t_channel_top_4f_inclusiveDecays->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_ST_tW_top_5f_NoFullyHadronicDecays, 19.55*1*lumi*1000/weight_ST_tW_top_5f_NoFullyHadronicDecays->GetBinContent(1));
	*/
	dataloader->AddBackgroundTree(tree_TTJets_DiLept, 831.76*((3*0.108)*(3*0.108))*1*lumi*1000/weight_TTJets_DiLept->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTJets_DiLept_ext1, 831.76*((3*0.108)*(3*0.108))*1*lumi*1000/weight_TTJets_DiLept_ext1->GetBinContent(1));
	/*
	dataloader->AddBackgroundTree(tree_TTJets_HT_1200to2500, 0.12*831.76/502.2*1*lumi*1000/weight_TTJets_HT_1200to2500->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTJets_HT_2500toInf, 0.001430*831.76/502.2*1*lumi*1000/weight_TTJets_HT_2500toInf->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTJets_HT_600to800, 1.610*831.76/502.2*1*lumi*1000/weight_TTJets_HT_600to800->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTJets_HT_800to1200, 0.663*831.76/502.2*1*lumi*1000/weight_TTJets_HT_800to1200->GetBinContent(1));
	
	dataloader->AddBackgroundTree(tree_TTJets_SingleLeptFromT, 831.76*(3*0.108)*(1-3*0.108)*lumi*1000/weight_TTJets_SingleLeptFromT->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTJets_SingleLeptFromTbar, 831.76*(3*0.108)*(1-3*0.108)*lumi*1000/weight_TTJets_SingleLeptFromTbar->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTJets_SingleLeptFromTbar_1, 831.76*(3*0.108)*(1-3*0.108)*lumi*1000/weight_TTJets_SingleLeptFromTbar_1->GetBinContent(1));
	
	dataloader->AddBackgroundTree(tree_TTTo2L2Nu, 0.2529*lumi*1000/weight_TTTo2L2Nu->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTWJetsToLNu_1, 0.2043*1*lumi*1000/weight_TTWJetsToLNu_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTWJetsToLNu_2, 0.2043*1*lumi*1000/weight_TTWJetsToLNu_2->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTZToLLNuNu_1, 0.2529*lumi*1000/weight_TTZToLLNuNu_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTZToLLNuNu_2, 0.2529*lumi*1000/weight_TTZToLLNuNu_2->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTZToLLNuNu_3, 0.2529*lumi*1000/weight_TTZToLLNuNu_3->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_TTZToQQ, 0.5297*lumi*1000/weight_TTZToQQ->GetBinContent(1));
	
	//dataloader->AddBackgroundTree(tree_TT, *lumi*1000/weight_TT->GetBinContent(1));
	
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_1200To2500	, 1.329*1.21*1*lumi*1000/weight_WJetsToLNu_HT_1200To2500->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_2500ToInf	, 0.03216*1.21*lumi*1000/weight_WJetsToLNu_HT_2500ToInf->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_2500ToInf_1, 0.03216*1.21*lumi*1000/weight_WJetsToLNu_HT_2500ToInf_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_400To600, 48.91*1.21*1*lumi*1000/weight_WJetsToLNu_HT_400To600->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_400To600_1	,48.91*1.21*1*lumi*1000/weight_WJetsToLNu_HT_400To600_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_600To800, 12.05*1*lumi*1000/weight_WJetsToLNu_HT_600To800->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_600To800_1, 12.05*1*lumi*1000/weight_WJetsToLNu_HT_600To800_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_800To1200, 5.501*1*lumi*1000/weight_WJetsToLNu_HT_800To1200->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT_800To1200_1, 5.501*1*lumi*1000/weight_WJetsToLNu_HT_800To1200_1->GetBinContent(1));
	
	//dataloader->AddBackgroundTree(tree_WJetsToLNu_1, *lumi*1000/weight_WJetsToLNu_1->GetBinContent(1));
	//dataloader->AddBackgroundTree(tree_WJetsToLNu_2, *lumi*1000/weight_WJetsToLNu_2->GetBinContent(1));
	//dataloader->AddBackgroundTree(tree_WJetsToLNu_madgraphMLM, *lumi*1000/weight_WJetsToLNu_madgraphMLM->GetBinContent(1));
	//dataloader->AddBackgroundTree(tree_WJetsToLNu_madgraphMLM_1, *lumi*1000/weight_WJetsToLNu_madgraphMLM_1->GetBinContent(1));
	
	dataloader->AddBackgroundTree(tree_WWTo2L2Nu, 10.481*1*lumi*1000/weight_WWTo2L2Nu->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WWToLNuQQ, 43.53*1*lumi*1000/weight_WWToLNuQQ->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WWToLNuQQ_1, 43.53*1*lumi*1000/weight_WWToLNuQQ_1->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WZTo1L1Nu2Q_2, 10.71*1*lumi*1000/weight_WZTo1L1Nu2Q_2->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_WZTo2L2Q, 5.60*1*lumi*1000/weight_WZTo2L2Q->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_ZZTo2L2Nu, 0.564*1*lumi*1000/weight_ZZTo2L2Nu->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_ZZTo2L2Q, 3.28*1*lumi*1000/weight_ZZTo2L2Q->GetBinContent(1));
	dataloader->AddBackgroundTree(tree_ZZTo2L2Q_1, 3.28*1*lumi*1000/weight_ZZTo2L2Q_1->GetBinContent(1));*/
   // Set individual event weights (the variables must exist in the original TTree)
   // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
   // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
   dataloader->SetBackgroundWeightExpression( "Jet_btagSF_csvv2*puWeight*lepSF" );
   dataloader->SetSignalWeightExpression( "Jet_btagSF_csvv2*puWeight*lepSF*nISRweight*Xsec" );
   //     factory->SetBackgroundWeightExpression("weight");
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // tell the factory to use all remaining events in the trees after training for testing:
   // dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
   //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None:!V" );
   // Boosted Decision Trees
   factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
			"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=2" );
   // For an example of the category classifier usage, see: TMVAClassificationCategory
   //
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   //     factory->OptimizeAllMethods("SigEffAt001","Scan");
   //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------


   outputFile->Close();

   delete factory;
   delete dataloader;
   TMVA::TMVAGui( outfileName );

   return 0;
}

int main()
{
   TMVA_1lep();
}



