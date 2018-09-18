/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Exectuable: TMVAClassificationApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void @SCRIPTNAME( TString myMethodList = ""  )
{

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();
   TString inFile = "@INFILE";
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 1;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 1; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 1;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object
   
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   Float_t Lep_pt,Lep_eta,Lep_miniIso,MT,MET,Jet1_pt,Jet2_pt,Run,Lumi,isData,HLT_EleOR,HLT_MuOR,HLT_LepOR,HLT_MetOR,iso_Veto,PD_JetHT,PD_SingleEle,PD_SingleMu,PD_MET;
   Long64_t Event;
   Float_t HT,LT,nBJet,dPhi,nEl,nMu,Selected,nLep,nVeto,RA2_muJetFilter,Flag_fastSimCorridorJetCleaning,btagSF,met_caloPt,isDPhiSignal,METfilters,Lp;
   Float_t puRatio,lepSF,Xsec,genWeight,susyXsec,susyNgen,nISRweight,nISRttweight,nJets30Clean,mGo,mLSP,nTrueInt;


   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   	reader->AddVariable( "Lep_pt"		,&Lep_pt		 );
	reader->AddVariable( "Lep_eta"		,&Lep_eta		 );
	reader->AddVariable( "Lep_miniIso"	,&Lep_miniIso	 );
	reader->AddVariable( "MT"			,&MT			 );
	reader->AddVariable( "MET"			,&MET			 );
	reader->AddVariable( "Jet1_pt"		,&Jet1_pt		 );
	reader->AddVariable( "Jet2_pt"		,&Jet2_pt		 );
	reader->AddVariable( "HT"			,&HT			 );
	reader->AddVariable( "LT"			,&LT			 );
	reader->AddVariable( "nBJet"		,&nBJet 		 );
	reader->AddVariable( "dPhi"			,&dPhi		  	 );
	reader->AddVariable( "nEl"			,&nEl			 );
	reader->AddVariable( "nMu"			,&nMu			 );
	reader->AddVariable( "Selected"		,&Selected		 );
	reader->AddVariable( "nLep "		,&nLep			 );
	reader->AddVariable( "nVeto"		,&nVeto  		 );

	reader->AddSpectator("RA2_muJetFilter"					,&RA2_muJetFilter					);
	reader->AddSpectator("Flag_fastSimCorridorJetCleaning"	,&Flag_fastSimCorridorJetCleaning	);
	reader->AddSpectator("btagSF"							,&btagSF							);
	reader->AddSpectator("puRatio"							,&puRatio							);
	reader->AddSpectator("lepSF"							,&lepSF								);
	reader->AddSpectator("Xsec"								,&Xsec								);
	reader->AddSpectator("genWeight"						,&genWeight							);
	reader->AddSpectator("susyXsec"							,&susyXsec							);
	reader->AddSpectator("susyNgen"							,&susyNgen							);
	reader->AddSpectator("nISRweight"						,&nISRweight						);
	reader->AddSpectator("nISRttweight"						,&nISRttweight						);
	reader->AddSpectator( "nJets30Clean"					,&nJets30Clean						);
	reader->AddSpectator( "mGo"								,&mGo								);
	reader->AddSpectator( "mLSP"							,&mLSP								);
	reader->AddSpectator( "nTrueInt"						,&nTrueInt							);
   // Book the MVA methods
   TString dir    = "dataset/weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = "/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src/tthAnalysis/NanoAODTools/TMVA/datasetBkg/weights/TMVAClassification_BDT.weights.xml";
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 100;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
   if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input(0);
   //TString fname = "./tmva_class_example.root";
   if (!gSystem->AccessPathName( inFile )) {
      input = TFile::Open( inFile ); // check if file in local directory exists
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("sf/t");
   theTree->SetBranchAddress("Run",&Run);
   theTree->SetBranchAddress("Lumi",&Lumi);
   theTree->SetBranchAddress("Event",&Event);
   theTree->SetBranchAddress("Lep_pt",&Lep_pt);
   theTree->SetBranchAddress("Lep_eta",&Lep_eta);
   theTree->SetBranchAddress("Lep_miniIso",&Lep_miniIso);
   theTree->SetBranchAddress("MT",&MT);
   theTree->SetBranchAddress("MET",&MET);
   theTree->SetBranchAddress("Jet1_pt",&Jet1_pt);
   theTree->SetBranchAddress("Jet2_pt",&Jet2_pt);
   theTree->SetBranchAddress("HT",&HT);
   theTree->SetBranchAddress("LT",&LT);
   theTree->SetBranchAddress("nBJet",&nBJet);
   theTree->SetBranchAddress("dPhi",&dPhi);
   theTree->SetBranchAddress("nEl",&nEl);
   theTree->SetBranchAddress("nMu",&nMu);
   theTree->SetBranchAddress("Selected",&Selected);
   theTree->SetBranchAddress("nLep",&nLep );
   theTree->SetBranchAddress("nVeto",&nVeto);
   theTree->SetBranchAddress("RA2_muJetFilter",&RA2_muJetFilter);
   theTree->SetBranchAddress("Flag_fastSimCorridorJetCleaning",&Flag_fastSimCorridorJetCleaning);
   theTree->SetBranchAddress("btagSF",&btagSF);
   theTree->SetBranchAddress("puRatio",&puRatio);
   theTree->SetBranchAddress("lepSF",&lepSF);
   theTree->SetBranchAddress("Xsec",&Xsec);
   theTree->SetBranchAddress("genWeight",&genWeight);
   theTree->SetBranchAddress("susyXsec",&susyXsec);
   theTree->SetBranchAddress("susyNgen",&susyNgen);
   theTree->SetBranchAddress("nISRweight",&nISRweight);
   theTree->SetBranchAddress("nISRttweight",&nISRttweight);
   theTree->SetBranchAddress("nJets30Clean",&nJets30Clean);
   theTree->SetBranchAddress("mGo",&mGo);
   theTree->SetBranchAddress("mLSP",&mLSP);
   theTree->SetBranchAddress("nTrueInt",&nTrueInt);
   theTree->SetBranchAddress("met_caloPt",&met_caloPt);
   
   theTree->SetBranchAddress("isData",&isData);
   theTree->SetBranchAddress("HLT_EleOR",&HLT_EleOR);
   theTree->SetBranchAddress("HLT_MuOR",&HLT_MuOR);
   theTree->SetBranchAddress("HLT_LepOR",&HLT_LepOR);
   theTree->SetBranchAddress("HLT_MetOR",&HLT_MetOR);
   theTree->SetBranchAddress("iso_Veto",&iso_Veto);
   theTree->SetBranchAddress("PD_JetHT",&PD_JetHT);
   theTree->SetBranchAddress("PD_SingleEle",&PD_SingleEle);
   theTree->SetBranchAddress("PD_SingleMu",&PD_SingleMu);
   theTree->SetBranchAddress("PD_MET",&PD_MET);
   theTree->SetBranchAddress("isDPhiSignal",&isDPhiSignal);
   theTree->SetBranchAddress("METfilters",&METfilters);
   theTree->SetBranchAddress("Lp",&Lp);


   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   TFile *target  = new TFile( "@OUTFILE","RECREATE" );
   Float_t b_BDT;
   TDirectory *sfdir;
   TTree *t=new TTree("t","t");
   t->Branch("Run",  &Run); 
   t->Branch("Lumi",  &Lumi); 
   t->Branch("Event",  &Event); 
   t->Branch("b_BDT",  &b_BDT); 
   t->Branch("Lep_pt",&Lep_pt);
   t->Branch("Lep_eta",&Lep_eta);
   t->Branch("Lep_miniIso",&Lep_miniIso);
   t->Branch("MT",&MT);
   t->Branch("MET",&MET);
   t->Branch("Jet1_pt",&Jet1_pt);
   t->Branch("Jet2_pt",&Jet2_pt);
   t->Branch("HT",&HT);
   t->Branch("LT",&LT);
   t->Branch("nBJet",&nBJet);
   t->Branch("dPhi",&dPhi);
   t->Branch("nEl",&nEl);
   t->Branch("nMu",&nMu);
   t->Branch("Selected",&Selected);
   t->Branch("nLep",&nLep );
   t->Branch("nVeto",&nVeto);
   t->Branch("RA2_muJetFilter",&RA2_muJetFilter);
   t->Branch("Flag_fastSimCorridorJetCleaning",&Flag_fastSimCorridorJetCleaning);
   t->Branch("btagSF",&btagSF);
   t->Branch("puRatio",&puRatio);
   t->Branch("lepSF",&lepSF);
   t->Branch("Xsec",&Xsec);
   t->Branch("genWeight",&genWeight);
   t->Branch("susyXsec",&susyXsec);
   t->Branch("susyNgen",&susyNgen);
   t->Branch("nISRweight",&nISRweight);
   t->Branch("nISRttweight",&nISRttweight);
   t->Branch("nJets30Clean",&nJets30Clean);
   t->Branch("mGo",&mGo);
   t->Branch("mLSP",&mLSP);
   t->Branch("nTrueInt",&nTrueInt);
   t->Branch("met_caloPt",&met_caloPt);
   t->Branch("isData",&isData);
   t->Branch("HLT_EleOR",&HLT_EleOR);
   t->Branch("HLT_MuOR",&HLT_MuOR);
   t->Branch("HLT_LepOR",&HLT_LepOR);
   t->Branch("HLT_MetOR",&HLT_MetOR);
   t->Branch("iso_Veto",&iso_Veto);
   t->Branch("PD_JetHT",&PD_JetHT);
   t->Branch("PD_SingleEle",&PD_SingleEle);
   t->Branch("PD_SingleMu",&PD_SingleMu);
   t->Branch("PD_MET",&PD_MET);
   t->Branch("isDPhiSignal",&isDPhiSignal);
   t->Branch("METfilters",&METfilters);
   t->Branch("Lp",&Lp);
   
   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

      //var1 = userVar1 + userVar2;
      //var2 = userVar1 - userVar2;

      // Return the MVA outputs and fill into histograms

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"));
      if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"));
      if (Use["BDT"          ]) {
		    histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
		    b_BDT= reader->EvaluateMVA( "BDT method");
		    t->Fill();		}
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );
      
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: "
                      << cutsMin[ivar]
                      << " < \""
                      << mcuts->GetInputVar(ivar)
                      << "\" <= "
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   // Write also error and significance histos
   sfdir = target->mkdir("sf");
   sfdir->cd();
   t->Write();
   target->Close();

   std::cout << "--- Created root file: " << target<<  "containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}

int main( int argc, char** argv )
{
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   @SCRIPTNAME(methodList);
   return 0;
}
