/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: ClassApplication                                                   *
 *                                                                                *
 * Test suit for comparison of Reader and standalone class outputs                *
 **********************************************************************************/
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <fstream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TBranch.h"
#include <vector>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
//#include "/nfs/dust/cms/user/bobovnii/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/MVA/weights/myTMVA_BDTupNTrees.class.C"
//#include "/nfs/dust/cms/user/bobovnii/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/functionsSUSY.h"
//#include "/nfs/dust/cms/user/bobovnii/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/interface/lester_mt2_bisect.h"
#include "/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/pureNANOAOD/CMSSW_9_4_4/src/tthAnalysis/NanoAODTools/TMVA/datasetBkg/weights/TMVAClassification_BDT.class.C"
void ClassApplication( TString myMethodList = "BDT" ) 
{

   // Declaration of leaf types
   Float_t         Lep_pt;
   Float_t         Lep_eta;
   Float_t         Lep_miniIso;
   Float_t         MT;
   Float_t         Met;
   Float_t         Jet1_pt;
   Float_t         Jet2_pt;
   Float_t         Jet1_eta;
   Float_t         Jet2_eta;
   Float_t         HT;
   Float_t         LT;
   Float_t         nJet;
   Float_t         nBJet;
   Float_t         dPhi;
   // List of branches
   TBranch        *b_classID;   //!
   TBranch        *b_className;   //!
   TBranch        *b_Lep_pt;   //!
   TBranch        *b_Lep_eta;   //!
   TBranch        *b_Lep_miniIso;   //!
   TBranch        *b_MT;   //!
   TBranch        *b_Met;   //!
   TBranch        *b_Jet1_pt;   //!
   TBranch        *b_Jet2_pt;   //!
   TBranch        *b_Jet1_eta;   //!
   TBranch        *b_Jet2_eta;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_LT;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_nBJet;   //!
   TBranch        *b_dPhi;   //!

///////////////////////////////////////// END of stupid definition!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111111

   cout << endl;
   cout << "==> start ClassApplication" << endl;
   const int Nmvas = 16;

   const char* bulkname[Nmvas] = { "BDT",/*"SVM"*/};

   bool iuse[Nmvas] = { Nmvas*kFALSE };

   // interpret input list
   if (myMethodList != "") {
      TList* mlist = TMVA::gTools().ParseFormatLine( myMethodList, " :," );
      for (int imva=0; imva<Nmvas; imva++) if (mlist->FindObject( bulkname[imva] )) iuse[imva] = kTRUE;
      delete mlist;
   }

   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to 
   // those given in the weight file(s) that you use
   std::vector<std::string> inputVars;

        inputVars.push_back( "Lep_pt" );
        inputVars.push_back( "Lep_eta" );
        inputVars.push_back( "Lep_miniIso" );
        inputVars.push_back( "MT" );
        inputVars.push_back( "Met" );
        inputVars.push_back( "Jet1_pt" );
        inputVars.push_back( "Jet2_pt" );
        inputVars.push_back( "Jet1_eta" );
        inputVars.push_back( "Jet2_eta" );
        inputVars.push_back( "HT" );
        inputVars.push_back( "LT" );
        inputVars.push_back( "nJet" );
        inputVars.push_back( "nBJet" );
        inputVars.push_back( "dPhi" );


      
   // preload standalone class(es)
   string dir    = "weights/";
   string prefix = "TMVAClassification";
	TH1* histBDTupNTrees;
	IClassifierReader* BDTupNTreesResponse = new ReadBDT( inputVars );
   	int nbin = 100;
        histBDTupNTrees = new TH1F( "MVA_BDTupNTrees",    "MVA_BDTupNTrees",    nbin,  -1, 1 );


   cout << "=== Macro        : Class creation was successful" << endl;

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input(0);

   if (!gSystem->AccessPathName("/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/Cut-NANO/all_trees_2016_Cut/all_2016_Data.root")) {
      // first we try to find tmva_example.root in the local directory
      cout << "=== Macro        : Accessing /nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/Cut-NANO/all_trees_2016_Cut/all_2016_Data.root" << endl;
      input = TFile::Open("/nfs/dust/cms/user/amohamed/susy-desy/nanoAOD/Cut-NANO/all_trees_2016_Cut/all_2016_Data.root");
   } 
  
/*
   if (!gSystem->AccessPathName("/nfs/dust/cms/user/bobovnii/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/mutau/stau-stau_100_LSP10_B_OS.root")) {
      // first we try to find tmva_example.root in the local directory
      cout << "=== Macro        : Accessing ./tmva_example.root" << endl;
      input = TFile::Open("/nfs/dust/cms/user/bobovnii/CMSSW_8_0_12/src/DesyTauAnalyses/NTupleMaker/test/mutau/stau-stau_100_LSP10_B_OS.root");
   }
*/
  
   if (!input) {
      cout << "ERROR: could not open data file" << endl;
      exit(1);
   }

   Float_t         v_Lep_pt;
   Float_t         v_Lep_eta;
   Float_t         v_Lep_miniIso;
   Float_t         v_MT;
   Float_t         v_Met;
   Float_t         v_Jet1_pt;
   Float_t         v_Jet2_pt;
   Float_t         v_Jet1_eta;
   Float_t         v_Jet2_eta;
   Float_t         v_HT;
   Float_t         v_LT;
   UInt_t         v_nJet;
   Int_t         v_nBJet;
   Float_t         v_dPhi;
   //
   // prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   TTree* theTree = (TTree*)input->Get("Events");
   cout << "=== Macro        : Loop over signal sample" << endl;

   // the references to the variables
		//////////////////////////////////////////////////////////////////////////////////////BAD waY!!!!!!!!!!!!!!!!!!!!!!!1

			 //theTree->SetBranchAddress("met_ex", &met_ex, &b_met_ex);
             theTree->SetBranchAddress( "Lep_pt"		, &v_Lep_pt	,&b_Lep_pt);
             theTree->SetBranchAddress( "Lep_eta"		, &v_Lep_eta,&b_Lep_eta);
             theTree->SetBranchAddress( "Lep_miniIso"	, &v_Lep_miniIso,&b_Lep_miniIso);
             theTree->SetBranchAddress( "MT"			, &v_MT,&b_MT);
             theTree->SetBranchAddress( "Met"			, &v_Met,&b_Met);
             theTree->SetBranchAddress( "Jet1_pt"		, &v_Jet1_pt,&b_Jet1_pt);
             theTree->SetBranchAddress( "Jet2_pt"		, &v_Jet2_pt,&b_Jet2_pt);
             theTree->SetBranchAddress( "Jet1_eta"		, &v_Jet1_eta,&b_Jet1_eta);
             theTree->SetBranchAddress( "Jet2_eta"		, &v_Jet2_eta,&b_Jet2_eta);
             theTree->SetBranchAddress( "HT"			, &v_HT,&b_HT);
             theTree->SetBranchAddress( "LT"			, &v_LT,&b_LT);
             theTree->SetBranchAddress( "nJet"			, &v_nJet,&b_nJet);
             theTree->SetBranchAddress( "nBJet"			, &v_nBJet,&b_nBJet);
             theTree->SetBranchAddress( "dPhi"			, &v_dPhi,&b_dPhi);

		//////////////////////////////////////////////////////////////////////////////////////End of BAD waY!!!!!!!!!!!!!!!!!!!!!!!1

   cout << "=== Macro        : Processing total of " << theTree->GetEntries() << " events ... " << endl;

   std::vector<double>* inputVec = new std::vector<double>( 14 );
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) cout << "=== Macro        : ... processing event: " << ievt << endl;

      theTree->GetEntry(ievt);
      
			(*inputVec)[0] = v_Lep_pt;
			(*inputVec)[1] = v_Lep_eta;
			(*inputVec)[2] = v_Lep_miniIso;
			(*inputVec)[3] = v_MT;
			(*inputVec)[4] = v_Met;
			(*inputVec)[5] = v_Jet1_pt;
			(*inputVec)[6] = v_Jet2_pt;
			(*inputVec)[7] = v_Jet1_eta;
			(*inputVec)[8] = v_Jet2_eta;
			(*inputVec)[9] = v_HT;
			(*inputVec)[10] = v_LT;
			(*inputVec)[11] = v_nJet;
			(*inputVec)[12] = v_nBJet;
			(*inputVec)[13] = v_dPhi;


			double retval = BDTupNTreesResponse->GetMvaValue( *inputVec );
        	histBDTupNTrees->Fill( retval, 1.0 );

   }
   
   cout << "=== Macro        : Event loop done! " << endl;

   TFile *target  = new TFile( "ClassAppData.root","RECREATE" );

         histBDTupNTrees->Write();
   cout << "=== Macro        : Created target file: " << target->GetName() << endl;
   target->Close();

   delete target;
   delete inputVec;
    
   cout << "==> ClassApplication is done!" << endl << endl;
} 
