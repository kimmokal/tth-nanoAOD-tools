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
void TMVA_1lep_CMG(){
		std::string factoryOptions( "!V:!Silent:Transformations=I;D;P;G,D:AnalysisType=Classification" );
		string inputFolder = "/nfs/dust/cms/user/amohamed/susy-desy/CMGSamples/Friends_for_Limits/";
		TString DYJetsToLL_M50_HT1200to2500=inputFolder+"evVarFriend_DYJetsToLL_M50_HT1200to2500.root";//
		TString DYJetsToLL_M50_HT2500toInf=inputFolder+"evVarFriend_DYJetsToLL_M50_HT2500toInf.root";//
		TString DYJetsToLL_M50_HT400to600=inputFolder+"evVarFriend_DYJetsToLL_M50_HT400to600.root";//
		TString DYJetsToLL_M50_HT600to800=inputFolder+"evVarFriend_DYJetsToLL_M50_HT600to800.root";//
		TString DYJetsToLL_M50_HT800to1200=inputFolder+"evVarFriend_DYJetsToLL_M50_HT800to1200.root";//
		TString QCD_HT1000to1500=inputFolder+"evVarFriend_QCD_HT1000to1500.root";//
		TString QCD_HT1500to2000=inputFolder+"evVarFriend_QCD_HT1500to2000.root";//
		TString QCD_HT2000toInf=inputFolder+"evVarFriend_QCD_HT2000toInf.root";//
		TString QCD_HT500to700=inputFolder+"evVarFriend_QCD_HT500to700.root";//
		TString QCD_HT700to1000=inputFolder+"evVarFriend_QCD_HT700to1000.root";//
		TString SMS_T1tttt_TuneCUETP8M1_1=inputFolder+"evVarFriend_SMS_T1tttt_TuneCUETP8M1_1.root";//
		TString SMS_T1tttt_TuneCUETP8M1_2=inputFolder+"evVarFriend_SMS_T1tttt_TuneCUETP8M1_2.root";//
		TString SMS_T1tttt_TuneCUETP8M1_3=inputFolder+"evVarFriend_SMS_T1tttt_TuneCUETP8M1_3.root";//
		TString SMS_T1tttt_TuneCUETP8M1_4=inputFolder+"evVarFriend_SMS_T1tttt_TuneCUETP8M1_4.root";//
		TString SMS_T1tttt_TuneCUETP8M1_5=inputFolder+"evVarFriend_SMS_T1tttt_TuneCUETP8M1_5.root";//
		TString SMS_T1tttt_TuneCUETP8M1_6=inputFolder+"evVarFriend_SMS_T1tttt_TuneCUETP8M1_6.root";//
		TString TBar_tWch_ext1=inputFolder+"evVarFriend_TBar_tWch_ext1.root";//
		TString TBar_tch_powheg=inputFolder+"evVarFriend_TBar_tch_powheg.root";//
		TString TTJets_DiLepton=inputFolder+"evVarFriend_TTJets_DiLepton.root";//
		TString TTJets_LO_HT1200to2500_ext=inputFolder+"evVarFriend_TTJets_LO_HT1200to2500_ext.root";//
		TString TTJets_LO_HT2500toInf_ext=inputFolder+"evVarFriend_TTJets_LO_HT2500toInf_ext.root";//
		TString TTJets_LO_HT600to800=inputFolder+"evVarFriend_TTJets_LO_HT600to800.root";//
		TString TTJets_LO_HT800to1200_ext=inputFolder+"evVarFriend_TTJets_LO_HT800to1200_ext.root";//
		TString TTJets_SingleLeptonFromT=inputFolder+"evVarFriend_TTJets_SingleLeptonFromT.root";//
		TString TTJets_SingleLeptonFromTbar=inputFolder+"evVarFriend_TTJets_SingleLeptonFromTbar.root";//
		TString TTWToLNu=inputFolder+"evVarFriend_TTWToLNu.root";//
		TString TTWToQQ=inputFolder+"evVarFriend_TTWToQQ.root";//
		TString TTZToLLNuNu=inputFolder+"evVarFriend_TTZToLLNuNu.root";//
		TString TTZToQQ=inputFolder+"evVarFriend_TTZToQQ.root";//
		TString TToLeptons_sch=inputFolder+"evVarFriend_TToLeptons_sch.root";//
		TString T_tWch_ext1=inputFolder+"evVarFriend_T_tWch_ext1.root";//
		TString T_tch_powheg=inputFolder+"evVarFriend_T_tch_powheg.root";//
		TString WJetsToLNu_HT1200to2500=inputFolder+"evVarFriend_WJetsToLNu_HT1200to2500.root";//
		TString WJetsToLNu_HT2500toInf=inputFolder+"evVarFriend_WJetsToLNu_HT2500toInf.root";//
		TString WJetsToLNu_HT400to600=inputFolder+"evVarFriend_WJetsToLNu_HT400to600.root";//
		TString WJetsToLNu_HT600to800=inputFolder+"evVarFriend_WJetsToLNu_HT600to800.root";//
		TString WJetsToLNu_HT800to1200_ext=inputFolder+"evVarFriend_WJetsToLNu_HT800to1200_ext.root";//
		TString WWTo2L2Nu=inputFolder+"evVarFriend_WWTo2L2Nu.root";//
		TString WWToLNuQQ=inputFolder+"evVarFriend_WWToLNuQQ.root";//
		TString WZTo1L1Nu2Q=inputFolder+"evVarFriend_WZTo1L1Nu2Q.root";//
		TString WZTo1L3Nu=inputFolder+"evVarFriend_WZTo1L3Nu.root";//
		TString WZTo2L2Q=inputFolder+"evVarFriend_WZTo2L2Q.root";//
		TString ZZTo2L2Nu=inputFolder+"evVarFriend_ZZTo2L2Nu.root";//
		TString ZZTo2L2Q=inputFolder+"evVarFriend_ZZTo2L2Q.root";//
		
		TFile *input_DYJetsToLL_M50_HT1200to2500;// 
		TFile *input_DYJetsToLL_M50_HT2500toInf;//  
		TFile *input_DYJetsToLL_M50_HT400to600;//   
		TFile *input_DYJetsToLL_M50_HT600to800;//   
		TFile *input_DYJetsToLL_M50_HT800to1200;//  
		TFile *input_QCD_HT1000to1500;//            
		TFile *input_QCD_HT1500to2000;//            
		TFile *input_QCD_HT2000toInf;//             
		TFile *input_QCD_HT500to700;//              
		TFile *input_QCD_HT700to1000;//             
		TFile *input_SMS_T1tttt_TuneCUETP8M1_1;//   
		TFile *input_SMS_T1tttt_TuneCUETP8M1_2;//   
		TFile *input_SMS_T1tttt_TuneCUETP8M1_3;//   
		TFile *input_SMS_T1tttt_TuneCUETP8M1_4;//   
		TFile *input_SMS_T1tttt_TuneCUETP8M1_5;//   
		TFile *input_SMS_T1tttt_TuneCUETP8M1_6;//   
		TFile *input_TBar_tWch_ext1;//              
		TFile *input_TBar_tch_powheg;//             
		TFile *input_TTJets_DiLepton;//             
		TFile *input_TTJets_LO_HT1200to2500_ext;//  
		TFile *input_TTJets_LO_HT2500toInf_ext;//   
		TFile *input_TTJets_LO_HT600to800;//        
		TFile *input_TTJets_LO_HT800to1200_ext;//   
		TFile *input_TTJets_SingleLeptonFromT;//    
		TFile *input_TTJets_SingleLeptonFromTbar;// 
		TFile *input_TTWToLNu;//                    
		TFile *input_TTWToQQ;//                     
		TFile *input_TTZToLLNuNu;//                 
		TFile *input_TTZToQQ;//                     
		TFile *input_TToLeptons_sch;//              
		TFile *input_T_tWch_ext1;//                 
		TFile *input_T_tch_powheg;//                
		TFile *input_WJetsToLNu_HT1200to2500;//     
		TFile *input_WJetsToLNu_HT2500toInf;//      
		TFile *input_WJetsToLNu_HT400to600;//       
		TFile *input_WJetsToLNu_HT600to800;//       
		TFile *input_WJetsToLNu_HT800to1200_ext;//  
		TFile *input_WWTo2L2Nu;//                   
		TFile *input_WWToLNuQQ;//                   
		TFile *input_WZTo1L1Nu2Q;//                 
		TFile *input_WZTo1L3Nu;//                   
		TFile *input_WZTo2L2Q;//                    
		TFile *input_ZZTo2L2Nu;//                   
		TFile *input_ZZTo2L2Q;//                    
		
		input_DYJetsToLL_M50_HT1200to2500=TFile::Open(DYJetsToLL_M50_HT1200to2500);
		input_DYJetsToLL_M50_HT2500toInf=TFile::Open(DYJetsToLL_M50_HT2500toInf);
		input_DYJetsToLL_M50_HT400to600=TFile::Open(DYJetsToLL_M50_HT400to600);
		input_DYJetsToLL_M50_HT600to800=TFile::Open(DYJetsToLL_M50_HT600to800);
		input_DYJetsToLL_M50_HT800to1200=TFile::Open(DYJetsToLL_M50_HT800to1200);
		input_QCD_HT1000to1500=TFile::Open(QCD_HT1000to1500);
		input_QCD_HT1500to2000=TFile::Open(QCD_HT1500to2000);
		input_QCD_HT2000toInf=TFile::Open(QCD_HT2000toInf);
		input_QCD_HT500to700=TFile::Open(QCD_HT500to700);
		input_QCD_HT700to1000=TFile::Open(QCD_HT700to1000);
		input_SMS_T1tttt_TuneCUETP8M1_1=TFile::Open(SMS_T1tttt_TuneCUETP8M1_1);
		input_SMS_T1tttt_TuneCUETP8M1_2=TFile::Open(SMS_T1tttt_TuneCUETP8M1_2);
		input_SMS_T1tttt_TuneCUETP8M1_3=TFile::Open(SMS_T1tttt_TuneCUETP8M1_3);
		input_SMS_T1tttt_TuneCUETP8M1_4=TFile::Open(SMS_T1tttt_TuneCUETP8M1_4);
		input_SMS_T1tttt_TuneCUETP8M1_5=TFile::Open(SMS_T1tttt_TuneCUETP8M1_5);
		input_SMS_T1tttt_TuneCUETP8M1_6=TFile::Open(SMS_T1tttt_TuneCUETP8M1_6);
		input_TBar_tWch_ext1=TFile::Open(TBar_tWch_ext1);
		input_TBar_tch_powheg=TFile::Open(TBar_tch_powheg);
		input_TTJets_DiLepton=TFile::Open(TTJets_DiLepton);
		input_TTJets_LO_HT1200to2500_ext=TFile::Open(TTJets_LO_HT1200to2500_ext);
		input_TTJets_LO_HT2500toInf_ext=TFile::Open(TTJets_LO_HT2500toInf_ext);
		input_TTJets_LO_HT600to800=TFile::Open(TTJets_LO_HT600to800);
		input_TTJets_LO_HT800to1200_ext=TFile::Open(TTJets_LO_HT800to1200_ext);
		input_TTJets_SingleLeptonFromT=TFile::Open(TTJets_SingleLeptonFromT);
		input_TTJets_SingleLeptonFromTbar=TFile::Open(TTJets_SingleLeptonFromTbar);
		input_TTWToLNu=TFile::Open(TTWToLNu);
		input_TTWToQQ=TFile::Open(TTWToQQ);
		input_TTZToLLNuNu=TFile::Open(TTZToLLNuNu);
		input_TTZToQQ=TFile::Open(TTZToQQ);
		input_TToLeptons_sch=TFile::Open(TToLeptons_sch);
		input_T_tWch_ext1=TFile::Open(T_tWch_ext1);
		input_T_tch_powheg=TFile::Open(T_tch_powheg);
		input_WJetsToLNu_HT1200to2500=TFile::Open(WJetsToLNu_HT1200to2500);
		input_WJetsToLNu_HT2500toInf=TFile::Open(WJetsToLNu_HT2500toInf);
		input_WJetsToLNu_HT400to600=TFile::Open(WJetsToLNu_HT400to600);
		input_WJetsToLNu_HT600to800=TFile::Open(WJetsToLNu_HT600to800);
		input_WJetsToLNu_HT800to1200_ext=TFile::Open(WJetsToLNu_HT800to1200_ext);
		input_WWTo2L2Nu=TFile::Open(WWTo2L2Nu);
		input_WWToLNuQQ=TFile::Open(WWToLNuQQ);
		input_WZTo1L1Nu2Q=TFile::Open(WZTo1L1Nu2Q);
		input_WZTo1L3Nu=TFile::Open(WZTo1L3Nu);
		input_WZTo2L2Q=TFile::Open(WZTo2L2Q);
		input_ZZTo2L2Nu=TFile::Open(ZZTo2L2Nu);
		input_ZZTo2L2Q=TFile::Open(ZZTo2L2Q);

		TTree *tree_DYJetsToLL_M50_HT1200to2500=(TTree*)input_DYJetsToLL_M50_HT1200to2500->Get("sf/t");
		TTree *tree_DYJetsToLL_M50_HT2500toInf =(TTree*)input_DYJetsToLL_M50_HT2500toInf->Get("sf/t");
		TTree *tree_DYJetsToLL_M50_HT400to600  =(TTree*)input_DYJetsToLL_M50_HT400to600->Get("sf/t");
		TTree *tree_DYJetsToLL_M50_HT600to800  =(TTree*)input_DYJetsToLL_M50_HT600to800->Get("sf/t");
		TTree *tree_DYJetsToLL_M50_HT800to1200 =(TTree*)input_DYJetsToLL_M50_HT800to1200->Get("sf/t");
		TTree *tree_QCD_HT1000to1500           =(TTree*)input_QCD_HT1000to1500->Get("sf/t");
		TTree *tree_QCD_HT1500to2000           =(TTree*)input_QCD_HT1500to2000->Get("sf/t");
		TTree *tree_QCD_HT2000toInf            =(TTree*)input_QCD_HT2000toInf->Get("sf/t");
		TTree *tree_QCD_HT500to700             =(TTree*)input_QCD_HT500to700->Get("sf/t");
		TTree *tree_QCD_HT700to1000            =(TTree*)input_QCD_HT700to1000->Get("sf/t");
		TTree *tree_SMS_T1tttt_TuneCUETP8M1_1  =(TTree*)input_SMS_T1tttt_TuneCUETP8M1_1->Get("sf/t");
		TTree *tree_SMS_T1tttt_TuneCUETP8M1_2  =(TTree*)input_SMS_T1tttt_TuneCUETP8M1_2->Get("sf/t");
		TTree *tree_SMS_T1tttt_TuneCUETP8M1_3  =(TTree*)input_SMS_T1tttt_TuneCUETP8M1_3->Get("sf/t");
		TTree *tree_SMS_T1tttt_TuneCUETP8M1_4  =(TTree*)input_SMS_T1tttt_TuneCUETP8M1_4->Get("sf/t");
		TTree *tree_SMS_T1tttt_TuneCUETP8M1_5  =(TTree*)input_SMS_T1tttt_TuneCUETP8M1_5->Get("sf/t");
		TTree *tree_SMS_T1tttt_TuneCUETP8M1_6  =(TTree*)input_SMS_T1tttt_TuneCUETP8M1_6->Get("sf/t");
		TTree *tree_TBar_tWch_ext1             =(TTree*)input_TBar_tWch_ext1->Get("sf/t");
		TTree *tree_TBar_tch_powheg            =(TTree*)input_TBar_tch_powheg->Get("sf/t");
		TTree *tree_TTJets_DiLepton            =(TTree*)input_TTJets_DiLepton->Get("sf/t");
		TTree *tree_TTJets_LO_HT1200to2500_ext =(TTree*)input_TTJets_LO_HT1200to2500_ext->Get("sf/t");
		TTree *tree_TTJets_LO_HT2500toInf_ext  =(TTree*)input_TTJets_LO_HT2500toInf_ext->Get("sf/t");
		TTree *tree_TTJets_LO_HT600to800       =(TTree*)input_TTJets_LO_HT600to800->Get("sf/t");
		TTree *tree_TTJets_LO_HT800to1200_ext  =(TTree*)input_TTJets_LO_HT800to1200_ext->Get("sf/t");
		TTree *tree_TTJets_SingleLeptonFromT   =(TTree*)input_TTJets_SingleLeptonFromT->Get("sf/t");
		TTree *tree_TTJets_SingleLeptonFromTbar=(TTree*)input_TTJets_SingleLeptonFromTbar->Get("sf/t");
		TTree *tree_TTWToLNu                   =(TTree*)input_TTWToLNu->Get("sf/t");
		TTree *tree_TTWToQQ                    =(TTree*)input_TTWToQQ->Get("sf/t");
		TTree *tree_TTZToLLNuNu                =(TTree*)input_TTZToLLNuNu->Get("sf/t");
		TTree *tree_TTZToQQ                    =(TTree*)input_TTZToQQ->Get("sf/t");
		TTree *tree_TToLeptons_sch             =(TTree*)input_TToLeptons_sch->Get("sf/t");
		TTree *tree_T_tWch_ext1                =(TTree*)input_T_tWch_ext1->Get("sf/t");
		TTree *tree_T_tch_powheg               =(TTree*)input_T_tch_powheg->Get("sf/t");
		TTree *tree_WJetsToLNu_HT1200to2500    =(TTree*)input_WJetsToLNu_HT1200to2500->Get("sf/t");
		TTree *tree_WJetsToLNu_HT2500toInf     =(TTree*)input_WJetsToLNu_HT2500toInf->Get("sf/t");
		TTree *tree_WJetsToLNu_HT400to600      =(TTree*)input_WJetsToLNu_HT400to600->Get("sf/t");
		TTree *tree_WJetsToLNu_HT600to800      =(TTree*)input_WJetsToLNu_HT600to800->Get("sf/t");
		TTree *tree_WJetsToLNu_HT800to1200_ext =(TTree*)input_WJetsToLNu_HT800to1200_ext->Get("sf/t");
		TTree *tree_WWTo2L2Nu                  =(TTree*)input_WWTo2L2Nu->Get("sf/t");
		TTree *tree_WWToLNuQQ                  =(TTree*)input_WWToLNuQQ->Get("sf/t");
		TTree *tree_WZTo1L1Nu2Q                =(TTree*)input_WZTo1L1Nu2Q->Get("sf/t");
		TTree *tree_WZTo1L3Nu                  =(TTree*)input_WZTo1L3Nu->Get("sf/t");
		TTree *tree_WZTo2L2Q                   =(TTree*)input_WZTo2L2Q->Get("sf/t");
		TTree *tree_ZZTo2L2Nu                  =(TTree*)input_ZZTo2L2Nu->Get("sf/t");
		TTree *tree_ZZTo2L2Q                   =(TTree*)input_ZZTo2L2Q->Get("sf/t"); 
   
   Double_t signalWeight      = 1.0;
   Double_t backgroundWeight = 1.0;
   Double_t lumi=35.9 ;
   // Create a new root output file.
   TString outfileName( "TMVASignalBackground_CMG.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   // ____________
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, factoryOptions );
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("datasetBkg");
   
	dataloader->AddVariable( "Lep_pt"		,"Lep_pt"		, "", 'F' );
	dataloader->AddVariable( "Lep_eta"		,"Lep_eta"		, "", 'F' );
	dataloader->AddVariable( "Lep_miniIso"	,"Lep_miniIso"	, "", 'F' );
	dataloader->AddVariable( "MT"			,"MT"			, "", 'F' );
	dataloader->AddVariable( "MET"			,"MET"			, "", 'F' );
	dataloader->AddVariable( "Jet1_pt"		,"Jet1_pt"		, "", 'F' );
	dataloader->AddVariable( "Jet2_pt"		,"Jet2_pt"		, "", 'F' );
	//dataloader->AddVariable( "Jet1_eta"		,"Jet1_eta"		, "", 'F' );
	//dataloader->AddVariable( "Jet2_eta"		,"Jet2_eta"		, "", 'F' );
	dataloader->AddVariable( "HT"			, "HT"			, "", 'F' );
	dataloader->AddVariable( "LT"			, "LT"			, "", 'F' );
	//dataloader->AddVariable( "nJet"			, "nJet"		, "", 'F' );
	dataloader->AddVariable( "nBJet"		, "nBJet"		, "", 'F' );
	dataloader->AddVariable( "dPhi"			, "dPhi"		, "", 'F' );
	dataloader->AddVariable( "nEl"			,"nEl"			,""	,'F');
	dataloader->AddVariable( "nMu"			,"nMu"			,""	,'F');
	dataloader->AddVariable( "Selected"		,"Selected"		,""	,'F');
	dataloader->AddVariable( "nLep"		,"nLep"			,""	,'F');
	dataloader->AddVariable( "nVeto"		,"nVeto"		,""	,'F');
	//dataloader->AddVariable( "iso_Veto"		,"iso_Veto"		,""	,'O');
	//dataloader->AddVariable( "met_caloPt"	,"met_caloPt"	,""	,'F');
	

	dataloader->AddSpectator("RA2_muJetFilter","RA2_muJetFilter","", 'O');
	dataloader->AddSpectator("Flag_fastSimCorridorJetCleaning","Flag_fastSimCorridorJetCleaning","", 'O');
	dataloader->AddSpectator("btagSF","btagSF","", 'D');
	dataloader->AddSpectator("puRatio","puRatio","", 'D');
	dataloader->AddSpectator("lepSF","lepSF","", 'D');
	dataloader->AddSpectator("Xsec","Xsec","", 'D');
	dataloader->AddSpectator("genWeight","genWeight","", 'D');
	dataloader->AddSpectator("susyXsec","susyXsec","", 'D');
	dataloader->AddSpectator("susyNgen","susyNgen","", 'D');
	dataloader->AddSpectator("nISRweight","nISRweight","", 'D');
	dataloader->AddSpectator("nISRttweight","nISRttweight","", 'D');
	//dataloader->AddSpectator( "nJet", "nJet", 'D' );
	dataloader->AddSpectator( "nJets30Clean"	,"nJets30Clean"	,""	,'F');
	dataloader->AddSpectator( "mGo"	, "mGo"	, 'D' );
	dataloader->AddSpectator( "mLSP", "mLSP", 'D' );
	dataloader->AddSpectator( "nTrueInt", "nTrueInt", 'D' );
	//dataloader->AddSpectator( "lheHTIncoming","lheHTIncoming" , "F");
	//dataloader->AddSpectator( "genTau_grandmotherId","genTau_grandmotherId" , "I");
	//dataloader->AddSpectator( "genLep_grandmotherId","genLep_grandmotherId" , "I");
	
	//dataloader->AddSpectator("nISRweight","nISRweight","", 'D');
	//dataloader->AddSpectator("Xsec","Xsec","", 'D');
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables

   //dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );
												 
		/*TH1F *weight_DYJetsToLL_M50_HT1200to2500 =(TH1F*)input_DYJetsToLL_M50_HT1200to2500->Get("Sum_Weights");
		TH1F *weight_DYJetsToLL_M50_HT2500toInf  =(TH1F*)input_DYJetsToLL_M50_HT2500toInf->Get("Sum_Weights");
		TH1F *weight_DYJetsToLL_M50_HT400to600   =(TH1F*)input_DYJetsToLL_M50_HT400to600->Get("Sum_Weights");
		TH1F *weight_DYJetsToLL_M50_HT600to800   =(TH1F*)input_DYJetsToLL_M50_HT600to800->Get("Sum_Weights");
		TH1F *weight_DYJetsToLL_M50_HT800to1200  =(TH1F*)input_DYJetsToLL_M50_HT800to1200->Get("Sum_Weights");
		TH1F *weight_QCD_HT1000to1500            =(TH1F*)input_QCD_HT1000to1500->Get("Sum_Weights");
		TH1F *weight_QCD_HT1500to2000            =(TH1F*)input_QCD_HT1500to2000->Get("Sum_Weights");
		TH1F *weight_QCD_HT2000toInf             =(TH1F*)input_QCD_HT2000toInf->Get("Sum_Weights");
		TH1F *weight_QCD_HT500to700              =(TH1F*)input_QCD_HT500to700->Get("Sum_Weights");
		TH1F *weight_QCD_HT700to1000             =(TH1F*)input_QCD_HT700to1000->Get("Sum_Weights");
		TH1F *weight_SMS_T1tttt_TuneCUETP8M1_1   =(TH1F*)input_SMS_T1tttt_TuneCUETP8M1_1->Get("Sum_Weights");
		TH1F *weight_SMS_T1tttt_TuneCUETP8M1_2   =(TH1F*)input_SMS_T1tttt_TuneCUETP8M1_2->Get("Sum_Weights");
		TH1F *weight_SMS_T1tttt_TuneCUETP8M1_3   =(TH1F*)input_SMS_T1tttt_TuneCUETP8M1_3->Get("Sum_Weights");
		TH1F *weight_SMS_T1tttt_TuneCUETP8M1_4   =(TH1F*)input_SMS_T1tttt_TuneCUETP8M1_4->Get("Sum_Weights");
		TH1F *weight_SMS_T1tttt_TuneCUETP8M1_5   =(TH1F*)input_SMS_T1tttt_TuneCUETP8M1_5->Get("Sum_Weights");
		TH1F *weight_SMS_T1tttt_TuneCUETP8M1_6   =(TH1F*)input_SMS_T1tttt_TuneCUETP8M1_6->Get("Sum_Weights");
		TH1F *weight_TBar_tWch_ext1              =(TH1F*)input_TBar_tWch_ext1->Get("Sum_Weights");
		TH1F *weight_TBar_tch_powheg             =(TH1F*)input_TBar_tch_powheg->Get("Sum_Weights");
		TH1F *weight_TTJets_DiLepton             =(TH1F*)input_TTJets_DiLepton->Get("Sum_Weights");
		TH1F *weight_TTJets_LO_HT1200to2500_ext  =(TH1F*)input_TTJets_LO_HT1200to2500_ext->Get("Sum_Weights");
		TH1F *weight_TTJets_LO_HT2500toInf_ext   =(TH1F*)input_TTJets_LO_HT2500toInf_ext->Get("Sum_Weights");
		TH1F *weight_TTJets_LO_HT600to800        =(TH1F*)input_TTJets_LO_HT600to800->Get("Sum_Weights");
		TH1F *weight_TTJets_LO_HT800to1200_ext   =(TH1F*)input_TTJets_LO_HT800to1200_ext->Get("Sum_Weights");
		TH1F *weight_TTJets_SingleLeptonFromT    =(TH1F*)input_TTJets_SingleLeptonFromT->Get("Sum_Weights");
		TH1F *weight_TTJets_SingleLeptonFromTbar =(TH1F*)input_TTJets_SingleLeptonFromTbar->Get("Sum_Weights");
		TH1F *weight_TTWToLNu                    =(TH1F*)input_TTWToLNu->Get("Sum_Weights");
		TH1F *weight_TTWToQQ                     =(TH1F*)input_TTWToQQ->Get("Sum_Weights");
		TH1F *weight_TTZToLLNuNu                 =(TH1F*)input_TTZToLLNuNu->Get("Sum_Weights");
		TH1F *weight_TTZToQQ                     =(TH1F*)input_TTZToQQ->Get("Sum_Weights");
		TH1F *weight_TToLeptons_sch              =(TH1F*)input_TToLeptons_sch->Get("Sum_Weights");
		TH1F *weight_T_tWch_ext1                 =(TH1F*)input_T_tWch_ext1->Get("Sum_Weights");
		TH1F *weight_T_tch_powheg                =(TH1F*)input_T_tch_powheg->Get("Sum_Weights");
		TH1F *weight_WJetsToLNu_HT1200to2500     =(TH1F*)input_WJetsToLNu_HT1200to2500->Get("Sum_Weights");
		TH1F *weight_WJetsToLNu_HT2500toInf      =(TH1F*)input_WJetsToLNu_HT2500toInf->Get("Sum_Weights");
		TH1F *weight_WJetsToLNu_HT400to600       =(TH1F*)input_WJetsToLNu_HT400to600->Get("Sum_Weights");
		TH1F *weight_WJetsToLNu_HT600to800       =(TH1F*)input_WJetsToLNu_HT600to800->Get("Sum_Weights");
		TH1F *weight_WJetsToLNu_HT800to1200_ext  =(TH1F*)input_WJetsToLNu_HT800to1200_ext->Get("Sum_Weights");
		TH1F *weight_WWTo2L2Nu                   =(TH1F*)input_WWTo2L2Nu->Get("Sum_Weights");
		TH1F *weight_WWToLNuQQ                   =(TH1F*)input_WWToLNuQQ->Get("Sum_Weights");
		TH1F *weight_WZTo1L1Nu2Q                 =(TH1F*)input_WZTo1L1Nu2Q->Get("Sum_Weights");
		TH1F *weight_WZTo1L3Nu                   =(TH1F*)input_WZTo1L3Nu->Get("Sum_Weights");
		TH1F *weight_WZTo2L2Q                    =(TH1F*)input_WZTo2L2Q->Get("Sum_Weights");
		TH1F *weight_ZZTo2L2Nu                   =(TH1F*)input_ZZTo2L2Nu->Get("Sum_Weights");
		TH1F *weight_ZZTo2L2Q                    =(TH1F*)input_ZZTo2L2Q->Get("Sum_Weights"); */

		
	dataloader->AddSignalTree    ( tree_SMS_T1tttt_TuneCUETP8M1_1,     1);
	dataloader->AddSignalTree    ( tree_SMS_T1tttt_TuneCUETP8M1_2,     1);
	dataloader->AddSignalTree    ( tree_SMS_T1tttt_TuneCUETP8M1_3,     1);
	dataloader->AddSignalTree    ( tree_SMS_T1tttt_TuneCUETP8M1_4,     1);
	dataloader->AddSignalTree    ( tree_SMS_T1tttt_TuneCUETP8M1_5,     1);
	dataloader->AddSignalTree    ( tree_SMS_T1tttt_TuneCUETP8M1_6,     1);
	
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M50_HT1200to2500	,	1);
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M50_HT2500toInf 	,	1);
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M50_HT400to600  	,	1);
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M50_HT600to800  	,	1);
	dataloader->AddBackgroundTree(tree_DYJetsToLL_M50_HT800to1200 	,	1);
	dataloader->AddBackgroundTree(tree_QCD_HT1000to1500           	,	1);
	dataloader->AddBackgroundTree(tree_QCD_HT1500to2000           	,	1);
	dataloader->AddBackgroundTree(tree_QCD_HT2000toInf            	,	1);
	dataloader->AddBackgroundTree(tree_QCD_HT500to700             	,	1);
	dataloader->AddBackgroundTree(tree_QCD_HT700to1000            	,	1);
	dataloader->AddBackgroundTree(tree_TBar_tWch_ext1             	,	1);
	dataloader->AddBackgroundTree(tree_TBar_tch_powheg            	,	1);
	dataloader->AddBackgroundTree(tree_TTJets_DiLepton            	,	1);
	dataloader->AddBackgroundTree(tree_TTJets_LO_HT1200to2500_ext 	,	1);
	dataloader->AddBackgroundTree(tree_TTJets_LO_HT2500toInf_ext  	,	1);
	dataloader->AddBackgroundTree(tree_TTJets_LO_HT600to800       	,	1);
	dataloader->AddBackgroundTree(tree_TTJets_LO_HT800to1200_ext  	,	1);
	dataloader->AddBackgroundTree(tree_TTJets_SingleLeptonFromT   	,	1);
	dataloader->AddBackgroundTree(tree_TTJets_SingleLeptonFromTbar	,	1);
	dataloader->AddBackgroundTree(tree_TTWToLNu                   	,	1);
	dataloader->AddBackgroundTree(tree_TTWToQQ                    	,	1);
	dataloader->AddBackgroundTree(tree_TTZToLLNuNu                	,	1);
	dataloader->AddBackgroundTree(tree_TTZToQQ                    	,	1);
	dataloader->AddBackgroundTree(tree_TToLeptons_sch             	,	1);
	dataloader->AddBackgroundTree(tree_T_tWch_ext1                	,	1);
	dataloader->AddBackgroundTree(tree_T_tch_powheg               	,	1);
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT1200to2500    	,	1);
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT2500toInf     	,	1);
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT400to600      	,	1);
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT600to800      	,	1);
	dataloader->AddBackgroundTree(tree_WJetsToLNu_HT800to1200_ext 	,	1);
	dataloader->AddBackgroundTree(tree_WWTo2L2Nu                  	,	1);
	dataloader->AddBackgroundTree(tree_WWToLNuQQ                  	,	1);
	dataloader->AddBackgroundTree(tree_WZTo1L1Nu2Q                	,	1);
	dataloader->AddBackgroundTree(tree_WZTo1L3Nu                  	,	1);
	dataloader->AddBackgroundTree(tree_WZTo2L2Q                   	,	1);
	dataloader->AddBackgroundTree(tree_ZZTo2L2Nu                  	,	1);
	dataloader->AddBackgroundTree(tree_ZZTo2L2Q 					,	1);
	  
   // Set individual event weights (the variables must exist in the original TTree)
   // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
   // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
   dataloader->SetBackgroundWeightExpression( "" );
   dataloader->SetSignalWeightExpression( "" );
   //     factory->SetBackgroundWeightExpression("weight");
   TCut mycuts = ""; //"mGo == 1600 && mLSP == 100"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // tell the factory to use all remaining events in the trees after training for testing:
   // dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
   //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
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
   TMVA_1lep_CMG();
}



