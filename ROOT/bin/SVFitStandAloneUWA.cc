#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include <math.h> 
#include "TMath.h" 
#include <limits>
#include "TSystem.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

//If recoilType 0 then don't do recoil
//              FIXME amc@nlo is not ready yet!!! 1 then aMC@NLO DY and W+Jets MC samples
//                1 is not longer an option
//              2 MG5 DY and W+Jets MC samples or Higgs MC samples
//
//If doES       0 does not apply any ES shifts
//              1 applies ES shifts to TT channel, no effect on other channels
//
//If isWJets    0 no shift in number of jets used for recoil corrections
//              1 shifts njets + 1 for recoil corrections
//
//If metType    1 use mvamet
//        -1 use pf met

void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType, int doES, int isWJets, int metType) ;
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);
void runSVFit(std::vector<svFitStandalone::MeasuredTauLepton> & measuredTauLeptons, TFile * inputFile_visPtResolution, double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass);

int main (int argc, char* argv[]) 
{
   optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
   parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__svFit__");
   parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
   parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
   parser.addOption("newOutputFile",optutl::CommandLineParser::kDouble,"New Output File",0.0);
   parser.addOption("recoilType",optutl::CommandLineParser::kDouble,"recoilType",0.0);
   parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
   parser.addOption("isWJets",optutl::CommandLineParser::kDouble,"isWJets",0.0);
   parser.addOption("metType",optutl::CommandLineParser::kDouble,"metType",-1.0); // 1 = mvamet, -1 = pf met

   parser.parseArguments (argc, argv);

   std::cout << "EXTRA COMMANDS:"
    << "\n --- recoilType: " << parser.doubleValue("recoilType")
    << "\n --- doES: " << parser.doubleValue("doES")
    << "\n --- isWJets: " << parser.doubleValue("isWJets")
    << "\n --- metType: " << parser.doubleValue("metType") << std::endl;

   // Make sure a proper Met Type is chosen
   assert (parser.doubleValue("metType") == 1.0 || parser.doubleValue("metType") == -1.0);
   
   char TreeToUse[80]="first" ;

   TFile *fProduce;//= new TFile(parser.stringValue("newFile").c_str(),"UPDATE");

   if(parser.doubleValue("newOutputFile")>0){
   TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"READ");
     std::cout<<"Creating new outputfile"<<std::endl;
     std::string newFileName = parser.stringValue("newFile");

     fProduce = new TFile(newFileName.c_str(),"RECREATE");
     copyFiles(parser, f, fProduce);//new TFile(parser.stringValue("inputFile").c_str()+"SVFit","UPDATE");
     fProduce = new TFile(newFileName.c_str(),"UPDATE");
     std::cout<<"listing the directories================="<<std::endl;
     fProduce->ls();
     readdir(fProduce,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
            parser.doubleValue("isWJets"),parser.doubleValue("metType"));

     fProduce->Close();
     f->Close();
   }
   else{
     TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"UPDATE");
     readdir(f,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
            parser.doubleValue("isWJets"),parser.doubleValue("metType"));
     f->Close();
   }


} 


void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType, int doES, int isWJets, int metType) 
{
  std::string recoilFileName = "HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root";
  if(recoilType == 1) { //amc@nlo
    std::cout << "Alexei no long specified MG vs. AMC@NLO, so use recoilType = 2" << std::endl;
    return; }
  if(recoilType == 2 && metType == 1) { // mva met (Alexei no long specified MG vs. AMC@NLO)
    std::cout << "Alexei does not provide full 2016 data recoil corrections for Mva Met\n\n" << std::endl;
    std::cout << "Using ICHEP Mva Met corrections\n\n" << std::endl;
    recoilFileName = "HTT-utilities/RecoilCorrections/data/MvaMET_2016BCD.root";}
  if(recoilType == 2 && metType == -1) { // pf met (Alexei no long specified MG vs. AMC@NLO)
    recoilFileName = "HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root";}

  svFitStandalone::kDecayType decayType1 = svFitStandalone::kUndefinedDecayType; //svFitStandalone::kTauToElecDecay
  svFitStandalone::kDecayType decayType2 = svFitStandalone::kUndefinedDecayType; //svFitStandalone::kTauToElecDecay
  // Both masses should depend on decay mode and particle?
  float mass1=0.;
  float mass2;
  std::string channel = "x";
  std::cout << "TreeToUse: " << TreeToUse << std::endl;

  if((std::string(TreeToUse).find("TauNom")!= std::string::npos)||(std::string(TreeToUse).find("TauUp")!= std::string::npos)||(std::string(TreeToUse).find("TauDown")!= std::string::npos)){
      std::cout<<"TreeToUse "<< std::string(TreeToUse)<<" Is TauNom, TauUp, TauDown... Skipping!!"<<std::endl;
  }
  else if((std::string(TreeToUse).find("muTauEvent")!= std::string::npos) ||
            (std::string(TreeToUse).find("first")!= std::string::npos) ||
            (std::string(TreeToUse).find("mutau_tree")!= std::string::npos)){
           decayType1 = svFitStandalone::kTauToMuDecay;
           decayType2 = svFitStandalone::kTauToHadDecay;
           mass1 = 0.105658;
           mass2 = 0;
           channel = "mt";
  }
  else if(std::string(TreeToUse).find("eleTauEvent")!= std::string::npos){
           std::cout<<"eleTauTree"<<std::endl;
           decayType1 = svFitStandalone::kTauToElecDecay;
           decayType2 = svFitStandalone::kTauToHadDecay;
           mass1 = 0.00051100;
           mass2 = 0;
           channel = "et";
  }
  else if(std::string(TreeToUse).find("em")!= std::string::npos){
        std::cout<< "EMu sample" <<std::endl;
           decayType1 = svFitStandalone::kTauToElecDecay;
           decayType2 = svFitStandalone::kTauToMuDecay;
           mass1 = 0.00051100;
           mass2 = 0.105658;
           channel = "em";
  }
  else if((std::string(TreeToUse).find("tt")!= std::string::npos) ||
                (parser.stringValue("inputFile").find("_tt.root") != std::string::npos)){
        std::cout<< "Double Hadronic sample" <<std::endl;
           decayType1 = svFitStandalone::kTauToHadDecay;
           decayType2 = svFitStandalone::kTauToHadDecay;
           mass1 = 0.13957;
           mass2 = 0.13957;
           channel = "tt";
  }
  else{
      std::cout<<"TreeToUse "<< std::string(TreeToUse)<<" does not match muTauEvent or eleTauEvent... Skipping!!"<<std::endl;
  }

  TDirectory *dirsav = gDirectory;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  char stringA[80]="first";
  dir->cd();      
  while ((key = (TKey*)next())) {
      printf("Found key=%s \n",key->GetName());

      TObject *obj = key->ReadObj();
      if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
          dir->cd(key->GetName());
          TDirectory *subdir = gDirectory;
          sprintf(TreeToUse,"%s",key->GetName());
          readdir(subdir,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
                  parser.doubleValue("isWJets"),parser.doubleValue("metType"));

          dirsav->cd();
      }
      else if(obj->IsA()->InheritsFrom(TTree::Class())) {

          TTree *t = (TTree*)obj;
          float svFitMass = -10;
          float svFitPt = -10;
          float svFitEta = -10;
          float svFitPhi = -10;
          float svFitMET = -10;
          float svFitTransverseMass = -10;

          float metcorr_ex = -10; // corrected met px (float)
          float metcorr_ey = -10;  // corrected met py (float)
          float metcor = -10; // corrected metcor
          float metcorphi = -10; // corrected metcorphi
          float mt1 = -10; // corrected mvamet
          float mt1UP = -10; // corrected mvamet
          float mt1DOWN = -10; // corrected mvamet

          //TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&svFitMass,(parser.stringValue("branch")+"/F").c_str());
          TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
          TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
          TBranch *newBranch3 = t->Branch("eta_sv", &svFitEta, "eta_sv/F");
          TBranch *newBranch4 = t->Branch("phi_sv", &svFitPhi, "phi_sv/F");
          TBranch *newBranch5 = t->Branch("met_sv", &svFitMET, "met_sv/F");
          TBranch *newBranch6 = t->Branch("mt_sv", &svFitTransverseMass, "mt_sv/F");

          TBranch *newBranch7 = t->Branch("metcorr_ex", &metcorr_ex, "metcorr_ex/F");
          TBranch *newBranch8 = t->Branch("metcorr_ey", &metcorr_ey, "metcorr_ey/F");
          TBranch *newBranch9 = t->Branch("metcor", &metcor, "metcor/F");
          TBranch *newBranch10 = t->Branch("metcorphi", &metcorphi, "metcorphi/F");
          TBranch *newBranch23 = t->Branch("mtRecoil_1", &mt1, "mtRecoil_1/F");
          TBranch *newBranch24 = t->Branch("mtRecoil_1UP", &mt1UP, "mtRecoil_1UP/F");
          TBranch *newBranch25 = t->Branch("mtRecoil_1DOWN", &mt1DOWN, "mtRecoil_1DOWN/F");

          // If doing ES shifts, we need extra ouput branches
          //
          float svFitMass_UP = -10;
          float svFitPt_UP = -10;
          float svFitEta_UP = -10;
          float svFitPhi_UP = -10;
          float svFitMET_UP = -10;
          float svFitTransverseMass_UP = -10;
          float svFitMass_DOWN = -10;
          float svFitPt_DOWN = -10;
          float svFitEta_DOWN = -10;
          float svFitPhi_DOWN = -10;
          float svFitMET_DOWN = -10;
          float svFitTransverseMass_DOWN = -10;
          TBranch *newBranch11 = t->Branch("m_sv_UP", &svFitMass_UP, "m_sv_UP/F");
          TBranch *newBranch12 = t->Branch("pt_sv_UP", &svFitPt_UP, "pt_sv_UP/F");
          TBranch *newBranch13 = t->Branch("eta_sv_UP", &svFitEta_UP, "eta_sv_UP/F");
          TBranch *newBranch14 = t->Branch("phi_sv_UP", &svFitPhi_UP, "phi_sv_UP/F");
          TBranch *newBranch15 = t->Branch("met_sv_UP", &svFitMET_UP, "met_sv_UP/F");
          TBranch *newBranch16 = t->Branch("mt_sv_UP", &svFitTransverseMass_UP, "mt_sv_UP/F");

          TBranch *newBranch17 = t->Branch("m_sv_DOWN", &svFitMass_DOWN, "m_sv_DOWN/F");
          TBranch *newBranch18 = t->Branch("pt_sv_DOWN", &svFitPt_DOWN, "pt_sv_DOWN/F");
          TBranch *newBranch19 = t->Branch("eta_sv_DOWN", &svFitEta_DOWN, "eta_sv_DOWN/F");
          TBranch *newBranch20 = t->Branch("phi_sv_DOWN", &svFitPhi_DOWN, "phi_sv_DOWN/F");
          TBranch *newBranch21 = t->Branch("met_sv_DOWN", &svFitMET_DOWN, "met_sv_DOWN/F");
          TBranch *newBranch22 = t->Branch("mt_sv_DOWN", &svFitTransverseMass_DOWN, "mt_sv_DOWN/F");

          //unsigned long long evt;
          //     int run, lumi;
          unsigned int evt, run, lumi;
          float pt1;
          float eta1;
          float phi1;
          float gen_match_1;
          float pt2;
          float eta2;
          float phi2;
          float m2;
          float gen_match_2;
          float decayMode;
          float covMatrix00;
          float covMatrix10;
          float covMatrix01;
          float covMatrix11;
          //float mvamet_ex, // uncorrected mva met px (float)
          //  mvamet_ey, // uncorrected mva met py (float)
          float  genPx=-999.    , // generator Z/W/Higgs px (float)
                 genPy =-999.   , // generator Z/W/Higgs py (float)
                 visPx =-999.   , // generator visible Z/W/Higgs px (float)
                 visPy =-999.   ; // generator visible Z/W/Higgs py (float)

          // define MET
          double measuredMETx = 0.;
          double measuredMETy = 0.;
          int njets;
          float pfmet;
          float pfmetphi;
          TLorentzVector TMet(0,0,0,0);
          // define MET covariance
          TMatrixD covMET(2, 2);
          //ele/mu variables
          TBranch *pt1branch;

          t->SetBranchAddress("evt",&evt);
          t->SetBranchAddress("run",&run);
          t->SetBranchAddress("lumi",&lumi);
          t->SetBranchAddress("pt_1",&pt1,&pt1branch);
          t->SetBranchAddress("eta_1",&eta1);
          t->SetBranchAddress("phi_1",&phi1);
          t->SetBranchAddress("pt_2",&pt2);
          t->SetBranchAddress("eta_2",&eta2);
          t->SetBranchAddress("phi_2",&phi2);
          t->SetBranchAddress("m_2",&m2);
          t->SetBranchAddress("tauDecayMode",&decayMode);
          t->SetBranchAddress("gen_match_2",&gen_match_1);
          t->SetBranchAddress("gen_match_1",&gen_match_2);
          t->SetBranchAddress("njets", &njets);
          t->SetBranchAddress("pfmet",&pfmet);
          t->SetBranchAddress("pfmetphi",&pfmetphi);
          // Recoil variables below
          t->SetBranchAddress( "genpX", &genPx);
          t->SetBranchAddress( "genpY", &genPy);
          t->SetBranchAddress( "vispX", &visPx);
          t->SetBranchAddress( "vispY", &visPy);
          // FOR PF MET ANALYSIS
          t->SetBranchAddress("cov00",&covMatrix00);
          t->SetBranchAddress("cov01",&covMatrix01);
          t->SetBranchAddress("cov10",&covMatrix10);
          t->SetBranchAddress("cov11",&covMatrix11);

          // use this RooT file when running on aMC@NLO DY and W+Jets MC samples
          RecoilCorrector* recoilCorrector = new RecoilCorrector(recoilFileName);
          if (metType == 1) std::cout<<"MetType: MvaMet"<<std::endl;
          if (metType == -1) std::cout<<"MetType: PF Met"<<std::endl;
          std::cout<<"recoiltype "<<recoilType<<" recoilFileName "<<recoilFileName<<std::endl;

          printf("Found tree -> weighting\n");

          // Only open this once!
          edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
          TH1::AddDirectory(false);  
          TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());  

          for(Int_t i=0;i<t->GetEntries();++i){
              t->GetEntry(i);

              //Recoil Correction time
              // Correct WJets recoil for faked lepton / extra jet
              int recoilNJets;
              if(isWJets) {
                  recoilNJets = njets + 1;
                  std::cout << " - njets: " << njets << " recoilNJets: " << recoilNJets << std::endl;
              }
              else recoilNJets = njets;

              // Using PF Met or Mva Met?
              if (metType == 1) { // 1 = Mva Met
                  std::cout<<"THIS IS BROKEN FIX ME"<<std::endl;
              } // mva met

              if (metType == -1) { // -1 = PF Met
                  TMet.SetPtEtaPhiM(pfmet,0,pfmetphi,0);
                  measuredMETx = pfmet*TMath::Cos(pfmetphi);
                  measuredMETy = pfmet*TMath::Sin(pfmetphi);

                  covMET[0][0] =  covMatrix00;
                  covMET[1][0] =  covMatrix10;
                  covMET[0][1] =  covMatrix01;
                  covMET[1][1] =  covMatrix11;
              } // pf met

              // Do recoil corrections if requested
              if(recoilType != 0){
                  // Alexie shows that corrections via quantile mapping provide best results
                  // for mva met and pf met, so
                  // use that as the defaul.  People can switch if they want below
                  //
                  // RecoilCorrector::Correct == Quantile Mapping
                  // RecoilCorrector::CorrectByMeanResolution == By Mean Res

                  //recoilCorrector->Correct(
                  recoilCorrector->CorrectByMeanResolution( // This method is faster (Alexei)
                          measuredMETx, // uncorrected mva met px (float)
                          measuredMETy, // uncorrected mva met py (float)
                          genPx, // generator Z/W/Higgs px (float)
                          genPy, // generator Z/W/Higgs py (float)
                          visPx, // generator visible Z/W/Higgs px (float)
                          visPy, // generator visible Z/W/Higgs py (float)
                          recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                          metcorr_ex, // corrected met px (float)
                          metcorr_ey  // corrected met py (float)
                          );
                  std::cout << " - MEASURED:  met_ex: " << measuredMETx << "  met_ey: " << measuredMETy << std::endl;
                  std::cout << " - CORRECTED: met_ex: " << metcorr_ex << "  met_ey: " << metcorr_ey << std::endl;
                  std::cout << " - INPUT VAR:  genPx: " << genPx << "  genPy: " << genPy << "   visPx: "<<visPx<<"   visPy: "<<visPy<<std::endl;
                  std::cout << " - INPUT VAR: Recoil NJets: "<<recoilNJets<<std::endl;
              }
              else{
                  metcorr_ex = measuredMETx;
                  metcorr_ey = measuredMETy;
              }
              if(channel=="mt"||channel=="et"){

                  float shiftTau=1.0;
                  if (decayMode==0) shiftTau=0.982;
                  if (decayMode==1) shiftTau=1.010;
                  if (decayMode==10) shiftTau=1.004;
                  std::cout<<"gen_match: "<<gen_match_2<<"  DM: "<<decayMode<<"  shiftTau: "<<shiftTau<<std::endl;
                  pt2 = pt2 * shiftTau;
                  double dx2, dy2;
                  dx2 = pt2 * TMath::Cos( phi2 ) * (( 1. / shiftTau ) - 1.);
                  dy2 = pt2 * TMath::Sin( phi2 ) * (( 1. / shiftTau ) - 1.);
                  metcorr_ex = metcorr_ex + dx2;
                  metcorr_ey = metcorr_ey + dy2;

                  metcor = TMath::Sqrt( metcorr_ex*metcorr_ex + metcorr_ey*metcorr_ey);
                  metcorphi = TMath::ATan2( metcorr_ey, metcorr_ex );

                  float px1 = pt1*TMath::Cos(phi1) + metcorr_ex;	
                  float py1 = pt1*TMath::Sin(phi1) + metcorr_ey;	
                  float et = pt1+metcor;	
                  mt1= TMath::Sqrt(et*et - (px1*px1 + py1*py1));

                  std::cout << " INPUT FOR SVFIT: - metcor "<<metcor<<" metcorphi "<<metcorphi<<std::endl;
                  std::cout<< "met00: "<<covMatrix00<<" met11: "<<covMatrix11<<" met10: "<<covMatrix10<<" met01: "<<covMatrix01<<std::endl;



                  mass2 = m2*shiftTau;

                  //mods needed for tau tau here
                  if(channel=="et" || channel=="mt"){
                      mass2 = m2;
                      if(decayMode==0)
                          mass2 = 0.13957;
                  }


                  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
                  measuredTauLeptons.push_back(
                          svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));

                  measuredTauLeptons.push_back(
                          svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode)); 
                  std::cout<< "evt: "<<evt<<"  run: "<<run<<"  lumi: "<<lumi<< "  pt1 " << pt1 << "  mass1 " << mass1 << "  pt2: "<< pt2<< "  mass2: "<< mass2 <<std::endl;        
                  //modified
                  runSVFit(measuredTauLeptons, inputFile_visPtResolution, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
                  std::cout<<"finished runningSVFit"<<std::endl;

                  if(doES) {

                      std::cout<<"runningSVFit for SHIFTS UP"<<std::endl;
                      if (gen_match_2<=5){
                          float ES_UP_scale=1.000; // this value is for jet -> tau fakes
                          if (gen_match_2<5) ES_UP_scale=1.01; // for gen matched ele/muon
                          if (gen_match_2==5) ES_UP_scale=1.006; // for real taus
                          double pt2_UP;
                          double mass2_UP=mass2;
                          pt2_UP = pt2 * ES_UP_scale;
                          if (decayMode!=0) mass2_UP = mass2*ES_UP_scale;
                          double metcorr_ex_UP, metcorr_ey_UP;
                          double dx2_UP, dy2_UP;
                          dx2_UP = pt2_UP * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
                          dy2_UP = pt2_UP * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
                          metcorr_ex_UP = metcorr_ex + dx2_UP;
                          metcorr_ey_UP = metcorr_ey + dy2_UP;

                          float metcorUP = TMath::Sqrt( metcorr_ex_UP*metcorr_ex_UP + metcorr_ey_UP*metcorr_ey_UP);
                          float px1 = pt1*TMath::Cos(phi1) + metcorr_ex_UP;	
                          float py1 = pt1*TMath::Sin(phi1) + metcorr_ey_UP;	
                          float et = pt1+metcorUP;	
                          mt1UP= TMath::Sqrt(et*et - (px1*px1 + py1*py1));

                          std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsUP;
                          measuredTauLeptonsUP.push_back(
                                  svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
                          measuredTauLeptonsUP.push_back(
                                  svFitStandalone::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2_UP, decayMode));

                          runSVFit(measuredTauLeptonsUP, inputFile_visPtResolution, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_UP, svFitPt_UP, svFitEta_UP, svFitPhi_UP, svFitMET_UP, svFitTransverseMass_UP);
                      }
                      else {
                          svFitMass_UP=svFitMass;
                          svFitPt_UP=svFitPt;
                          svFitEta_UP=svFitEta;
                          svFitPhi_UP=svFitPhi;
                          svFitMET_UP=svFitMET;
                          svFitTransverseMass_UP=svFitTransverseMass;
                      }

                      // Second ES Down, x 0.97
                      if (gen_match_2<=5){
                          float ES_DOWN_scale=1.000; // jet
                          if (gen_match_2==5) ES_DOWN_scale=0.994; // tau
                          if (gen_match_2<5) ES_DOWN_scale=0.990;  // elec/mu
                          double pt2_DOWN;
                          double mass2_DOWN=mass2;
                          if (decayMode!=0) mass2_DOWN = mass2*ES_DOWN_scale;
                          pt2_DOWN = pt2 * ES_DOWN_scale;
                          double metcorr_ex_DOWN, metcorr_ey_DOWN;
                          double dx2_DOWN, dy2_DOWN;
                          dx2_DOWN = pt2_DOWN * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
                          dy2_DOWN = pt2_DOWN * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
                          metcorr_ex_DOWN = metcorr_ex + dx2_DOWN;
                          metcorr_ey_DOWN = metcorr_ey + dy2_DOWN;

                          float metcorDOWN = TMath::Sqrt( metcorr_ex_DOWN*metcorr_ex_DOWN + metcorr_ey_DOWN*metcorr_ey_DOWN);
                          float px1 = pt1*TMath::Cos(phi1) + metcorr_ex_DOWN;
                          float py1 = pt1*TMath::Sin(phi1) + metcorr_ey_DOWN;
                          float et = pt1+metcorDOWN;
                          mt1DOWN= TMath::Sqrt(et*et - (px1*px1 + py1*py1));

                          std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsDOWN;          
                          measuredTauLeptonsDOWN.push_back(
                                  svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
                          measuredTauLeptonsDOWN.push_back(
                                  svFitStandalone::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2_DOWN, decayMode));

                          runSVFit(measuredTauLeptonsDOWN, inputFile_visPtResolution, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DOWN, svFitPt_DOWN, svFitEta_DOWN, svFitPhi_DOWN, svFitMET_DOWN, svFitTransverseMass_DOWN);
                      }
                      else {
                          svFitMass_DOWN=svFitMass;
                          svFitPt_DOWN=svFitPt;
                          svFitEta_DOWN=svFitEta;
                          svFitPhi_DOWN=svFitPhi;
                          svFitMET_DOWN=svFitMET;
                          svFitTransverseMass_DOWN=svFitTransverseMass;
                      }
                  } // end doES
              } // eTau / muTau

              else {
                  svFitMass = -100;
                  svFitPt = -100;
                  svFitEta = -100;
                  svFitPhi = -100;
                  svFitMET = -100;
                  svFitTransverseMass = -100;
              }
              std::cout << "\n\n" << std::endl;
              //std::cout << "\n\nex: " << metcorr_ex << "   ey: " << metcorr_ey <<  " phi: " << metcorphi<<"\n"<<std::endl; 
              newBranch1->Fill();
              newBranch2->Fill();
              newBranch3->Fill();
              newBranch4->Fill();
              newBranch5->Fill();
              newBranch6->Fill();
              newBranch7->Fill();
              newBranch8->Fill();
              newBranch9->Fill();
              newBranch10->Fill();

              newBranch11->Fill();
              newBranch12->Fill();
              newBranch13->Fill();
              newBranch14->Fill();
              newBranch15->Fill();
              newBranch16->Fill();
              newBranch17->Fill();
              newBranch18->Fill();
              newBranch19->Fill();
              newBranch20->Fill();
              newBranch21->Fill();
              newBranch22->Fill();

              newBranch23->Fill(); //Recoil
              newBranch24->Fill(); //Recoil
              newBranch25->Fill(); //Recoil


          }
          delete inputFile_visPtResolution;
          dir->cd();
          t->Write("",TObject::kOverwrite);
          strcpy(TreeToUse,stringA) ;

      }
  }
}

void runSVFit(std::vector<svFitStandalone::MeasuredTauLepton> & measuredTauLeptons, TFile * inputFile_visPtResolution, double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass){

    //edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
    //TH1::AddDirectory(false);  
    //TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());  
    //float svFitMass = 0;
    SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
    algo.addLogM(false);  
    algo.shiftVisPt(true, inputFile_visPtResolution);
    algo.integrateMarkovChain();
    svFitMass = algo.getMass(); // return value is in units of GeV
    svFitPt = algo.pt();
    svFitEta = algo.eta();
    svFitPhi = algo.phi();
    svFitMET = algo.fittedMET().Rho();
    svFitTransverseMass = algo.transverseMass();
    if ( algo.isValidSolution() ) {
        std::cout << "found mass = " << svFitMass << std::endl;
    } else {
        std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
    }
    //delete inputFile_visPtResolution;

}

//Thank you Renee Brun :)
void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
    //copy all objects and subdirs of directory source as a subdir of the current directory
    TDirectory *savdir = gDirectory;
    TDirectory *adir = savdir; 
    if(source->GetName()!=parser.stringValue("inputFile")){
        adir = savdir->mkdir(source->GetName());
        std::cout<<"Source name is not outputfile name"<<std::endl;
        adir->cd();    
    }
    else{
        //adir = savdir->mkdir("input");
        adir->cd();    
    }

    //loop on all entries of this directory
    TKey *key;
    TIter nextkey(source->GetListOfKeys());
    while ((key = (TKey*)nextkey())) {
        const char *classname = key->GetClassName();
        TClass *cl = gROOT->GetClass(classname);
        if (!cl) continue;
        if (cl->InheritsFrom(TDirectory::Class())) {
            source->cd(key->GetName());
            TDirectory *subdir = gDirectory;
            adir->cd();
            CopyDir(subdir,parser);
            adir->cd();
        } else if (cl->InheritsFrom(TTree::Class())) {
            TTree *T = (TTree*)source->Get(key->GetName());
            adir->cd();
            TTree *newT = T->CloneTree(-1,"fast");
            newT->Write();
        } else {
            source->cd();
            TObject *obj = key->ReadObj();
            adir->cd();
            obj->Write();
            delete obj;
        }
    }
    adir->SaveSelf(kTRUE);
    savdir->cd();
}
void CopyFile(const char *fname, optutl::CommandLineParser parser) {
    //Copy all objects and subdirs of file fname as a subdir of the current directory
    TDirectory *target = gDirectory;
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
        printf("Cannot copy file: %s\n",fname);
        target->cd();
        return;
    }
    target->cd();
    CopyDir(f,parser);
    delete f;
    target->cd();
}
void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) 
{
    //prepare files to be copied
    if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
        gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
    }

    fNew->cd();
    CopyFile(parser.stringValue("inputFile").c_str(),parser);
    fNew->ls();
    fNew->Close();

}

