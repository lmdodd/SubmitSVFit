#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1.h"
#include <math.h> 
#include "TMath.h" 
#include <limits>
#include "TSystem.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"


//If recoilType 0 then don't do recoil
//              1 then aMC@NLO DY and W+Jets MC samples
//              2 MG5 DY and W+Jets MC samples or Higgs MC samples

void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
//void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType, int doES, int isWJets, int metType) ;
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);
void runSVFit(std::vector<svFitStandalone::MeasuredTauLepton> & measuredTauLeptons, TFile * inputFile_visPtResolution, double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass);

int main (int argc, char* argv[]) 
{
    optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
    parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__svFit__");
    parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile");
    parser.addOption("newOutputFile",optutl::CommandLineParser::kDouble,"New Output File",0.0);
    parser.addOption("recoilType",optutl::CommandLineParser::kDouble,"recoilType",0.0);
    parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
    parser.addOption("isWJets",optutl::CommandLineParser::kDouble,"isWJets",0.0);
    parser.addOption("metType",optutl::CommandLineParser::kDouble,"metType",-1.0); // -1 = pf met default


    parser.parseArguments (argc, argv);

    assert (parser.doubleValue("metType") == 1.0 || parser.doubleValue("metType") == -1.0);
    std::cout << "EXTRA COMMANDS:"
        << "\n --- recoilType: " << parser.doubleValue("recoilType")
        << "\n --- doES: " << parser.doubleValue("doES")
        << "\n --- isWJets: " << parser.doubleValue("isWJets")
        << "\n --- metType: " << parser.doubleValue("metType") << std::endl;


    char TreeToUse[80]="first" ;

    TFile *fProduce;//= new TFile(parser.stringValue("newFile").c_str(),"UPDATE");

    if(parser.doubleValue("newOutputFile")>0){
        TFile *f = new TFile(parser.stringValue("outputFile").c_str(),"READ");
        std::cout<<"Creating new outputfile"<<std::endl;
        std::string newFileName = parser.stringValue("newFile");

        fProduce = new TFile(newFileName.c_str(),"RECREATE");
        copyFiles(parser, f, fProduce);//new TFile(parser.stringValue("outputFile").c_str()+"SVFit","UPDATE");
        fProduce = new TFile(newFileName.c_str(),"UPDATE");
        std::cout<<"listing the directories================="<<std::endl;
        fProduce->ls();
        readdir(fProduce,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
                parser.doubleValue("isWJets"),parser.doubleValue("metType"));
        fProduce->Close();
        f->Close();
    }
    else{
        TFile *f = new TFile(parser.stringValue("outputFile").c_str(),"UPDATE");
        readdir(f,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
                parser.doubleValue("isWJets"),parser.doubleValue("metType"));
        f->Close();
    }


} 


void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType,int doES, int isWJets, int metType ) 
{
    std::string recoilFileName = "HTT-utilities/RecoilCorrections/data/MvaMET_2016BCD.root";
    if(recoilType == 1) { //amc@nlo
        std::cout << "Recoil Corrections for amc@nlo are not ready yet, use MadGraph samples!" << std::endl;
        return; }
    if(recoilType == 2 && metType == 1) //MG5 mva met
        recoilFileName = "HTT-utilities/RecoilCorrections/data/MvaMET_2016BCD.root";
    if(recoilType == 2 && metType == -1) //MG5 pf met
        recoilFileName = "HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root";



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

            float mvametcorr_ex = -10; // corrected met px (float)
            float mvametcorr_ey = -10;  // corrected met py (float)
            float mt1 = -10; // corrected mvamet
            float mvamet = -10; // corrected mvamet
            float mvametphi = -10; // corrected mvametphi

            //TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&svFitMass,(parser.stringValue("branch")+"/F").c_str());
            TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
            TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
            TBranch *newBranch3 = t->Branch("eta_sv", &svFitEta, "eta_sv/F");
            TBranch *newBranch4 = t->Branch("phi_sv", &svFitPhi, "phi_sv/F");
            TBranch *newBranch5 = t->Branch("met_sv", &svFitMET, "met_sv/F");
            TBranch *newBranch6 = t->Branch("mt_sv", &svFitTransverseMass, "mt_sv/F");

            TBranch *newBranch7 = t->Branch("metcorr_ex", &mvametcorr_ex, "metcorr_ex/F");
            TBranch *newBranch8 = t->Branch("metcorr_ey", &mvametcorr_ey, "metcorr_ey/F");
            TBranch *newBranch9 = t->Branch("metcorr", &mvamet, "metcorr/F");
            TBranch *newBranch10 = t->Branch("metphicorr", &mvametphi, "metphicorr/F");
            TBranch *newBranch11 = t->Branch("mtRecoil_1", &mt1, "mtRecoil_1/F");


            unsigned int evt, run, lumi;
            int evt2, run2, lumi2;
            evt2=0; run2=0; lumi2=0;
            float px1;
            float py1;
            float pt1;
            float eta1;
            float phi1;
            float pt2;
            float eta2;
            float phi2;
            float m2;
            float decayMode;
            float genMatch2;
            float covMatrix00;
            float covMatrix10;
            float covMatrix01;
            float covMatrix11;
            //float mvamet_ex, // uncorrected mva met px (float)
            //  mvamet_ey, // uncorrected mva met py (float)
            float  genPx    , // generator Z/W/Higgs px (float)
                   genPy    , // generator Z/W/Higgs py (float)
                   visPx    , // generator visible Z/W/Higgs px (float)
                   visPy    ; // generator visible Z/W/Higgs py (float)
            int    njets    ;  // number of jets (hadronic jet multiplicity) (int)
            // define MET
            double measuredMETx;
            double measuredMETy;
            float met;
            float metphi;
            TLorentzVector TMet(0,0,0,0);
            // define MET covariance
            TMatrixD covMET(2, 2);
            //ele/mu variables
            svFitStandalone::kDecayType decayType1 = svFitStandalone::kUndefinedDecayType; //svFitStandalone::kTauToElecDecay
            svFitStandalone::kDecayType decayType2 = svFitStandalone::kUndefinedDecayType; //svFitStandalone::kTauToElecDecay
            // Both masses should depend on decay mode and particle?
            float mass1;
            float mass2;
            std::string channel = "x";
            std::string shift = "x";

            if(std::string(TreeToUse).find("TauNom")!=std::string::npos){
                shift="nom";
            }
            else if(std::string(TreeToUse).find("TauUp")!=std::string::npos){
                shift="up";
            }
            else if(std::string(TreeToUse).find("TauDown")!=std::string::npos){
                shift="down";
            }

            if(std::string(TreeToUse).find("muTauEvent")!= std::string::npos){
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
            //else if(parser.stringValue("outputFile").find("_em.root") != std::string::npos){
            else if(std::string(TreeToUse).find("em")!= std::string::npos){
                //std::cout<< parser.stringValue("outputFile").c_str() << std::endl;
                std::cout<< "EMu sample" <<std::endl;
                decayType1 = svFitStandalone::kTauToElecDecay;
                decayType2 = svFitStandalone::kTauToMuDecay;
                mass1 = 0.00051100;
                mass2 = 0;
                channel = "em";
            }
            //else if(parser.stringValue("outputFile").find("_tt.root") != std::string::npos){
            else if(std::string(TreeToUse).find("tt")!= std::string::npos){
                //std::cout<< parser.stringValue("outputFile").c_str() << std::endl;
                std::cout<< "Double Hadronic sample" <<std::endl;
                decayType1 = svFitStandalone::kTauToHadDecay;
                decayType2 = svFitStandalone::kTauToHadDecay;
                mass1 = 0.13957;
                mass2 = 0.13957;
                channel = "tt";
            }
            else{
                std::cout<<"TreeToUse "<< std::string(TreeToUse)<<" does not match muTauEvent or eleTauEvent... Skipping!!"<<std::endl;
                continue;
            }

            //MET, MET Covariance, lepton objects,
            TBranch *pt1branch;
            //TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&svFitMass,(parser.stringValue("branch")+"/F").c_str());
            //newBranch->SetFile("svfitout.root:/input/muTauEventTree");
            t->SetBranchAddress("EVENT",&evt);
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
            t->SetBranchAddress("gen_match_2",&genMatch2);
            if (metType==1.0){
                t->SetBranchAddress("mvacov00",&covMatrix00);
                t->SetBranchAddress("mvacov01",&covMatrix01);
                t->SetBranchAddress("mvacov10",&covMatrix10);
                t->SetBranchAddress("mvacov11",&covMatrix11);
                t->SetBranchAddress("mvamet",&met);
                t->SetBranchAddress("mvametphi",&metphi);
            }
            if (metType==-1.0){
                t->SetBranchAddress("cov00",&covMatrix00);
                t->SetBranchAddress("cov01",&covMatrix01);
                t->SetBranchAddress("cov10",&covMatrix10);
                t->SetBranchAddress("cov11",&covMatrix11);
                t->SetBranchAddress("pfmet",&met);
                t->SetBranchAddress("pfmetphi",&metphi);
            }

            t->SetBranchAddress("njets", &njets);
            t->SetBranchAddress( "genpX", &genPx);
            t->SetBranchAddress( "genpY", &genPy);
            t->SetBranchAddress( "vispX", &visPx);
            t->SetBranchAddress( "vispY", &visPy);

            // use this RooT file when running on aMC@NLO DY and W+Jets MC samples
            RecoilCorrector* recoilMvaMetCorrector = new RecoilCorrector(recoilFileName);
            std::cout<<"recoiltype "<<recoilType<<" recoilFileName "<<recoilFileName<<std::endl;

            printf("Found tree -> weighting\n");

            // Only open this once!
            edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
            TH1::AddDirectory(false);  
            TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());  

            for(Int_t i=0;i<t->GetEntries();++i)
            {
                t->GetEntry(i);

                if(channel=="tt" || channel=="em"){
                    evt = (unsigned int)evt2;
                    run = (unsigned int)run2;
                    lumi = (unsigned int)lumi2;
                }
                //Recoil Correction time
                // Correct WJets recoil for faked lepton / extra jet

                int recoilNJets;
                if(isWJets) {
                    recoilNJets = njets + 1;
                    std::cout << " - njets: " << njets << " recoilNJets: " << recoilNJets << std::endl;
                }
                else recoilNJets = njets;


                TMet.SetPtEtaPhiM(met,0,metphi,0);
                measuredMETx = met*TMath::Cos(metphi);
                measuredMETy = met*TMath::Sin(metphi);
                //Recoil Correction time
                if(recoilType != 0){

                    //RecoilCorrector recoilMvaMetCorrector(recoilFileName);
                    recoilMvaMetCorrector->CorrectByMeanResolution(measuredMETx, // uncorrected mva met px (float)
                            measuredMETy, // uncorrected mva met py (float)
                            genPx, // generator Z/W/Higgs px (float)
                            genPy, // generator Z/W/Higgs py (float)
                            visPx, // generator visible Z/W/Higgs px (float)
                            visPy, // generator visible Z/W/Higgs py (float)
                            recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                            mvametcorr_ex, // corrected met px (float)
                            mvametcorr_ey  // corrected met py (float)
                            );}
                else{
                    mvametcorr_ex = measuredMETx;
                    mvametcorr_ey = measuredMETy;
                }

                //std::cout<<"metcorr_ex: "<<mvametcorr_ex<<std::endl;
                //std::cout<<"metcorr_ey: "<<mvametcorr_ey<<std::endl;
                mvamet = TMath::Sqrt( mvametcorr_ex*mvametcorr_ex + mvametcorr_ey*mvametcorr_ey);
                mvametphi = TMath::ATan2( mvametcorr_ey, mvametcorr_ex );

                px1 = pt1*TMath::Cos(phi1) + mvametcorr_ex;	
                py1 = pt1*TMath::Sin(phi1) + mvametcorr_ey;	
                float et = pt1+mvamet;	
                mt1= TMath::Sqrt(et*et - (px1*px1 + py1*py1));


                //mods needed for tau tau here
                if(channel=="et" || channel=="mt"){
                    mass2 = m2;
                    if(decayMode==0)
                        mass2 = 0.13957;
                }

                covMET[0][0] =  covMatrix00;
                covMET[1][0] =  covMatrix10;
                covMET[0][1] =  covMatrix01;
                covMET[1][1] =  covMatrix11;
                //std::cout<<"getting decay mode "<< decayMode<<std::endl;
                if((channel!="tt"&&channel!="em")&&(decayMode==0||decayMode==1||decayMode==10)){
                    // define lepton four vectors
                    std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
                    // check if electron or muon
                    measuredTauLeptons.push_back(
                            svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1)
                            ); // tau -> electron decay (Pt, eta, phi, mass)

                    measuredTauLeptons.push_back(
                            svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode)
                            ); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())
                    //std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
                    //std::cout<< "tauDM: "<<decayMode<<std::endl;
                    //std::cout<< "met: "<<met<<" recoilmet: "<<mvamet<<std::endl;
                    //std::cout<< "met00: "<<covMatrix00<<" met11: "<<covMatrix11<<" met10: "<<covMatrix10<<" met01: "<<covMatrix01<<std::endl;



                    //modified
                    runSVFit(measuredTauLeptons, inputFile_visPtResolution, mvametcorr_ex, mvametcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
                    //std::cout<<"finished runningSVFit"<<std::endl;
                } // eTau / muTau
                else {
                    svFitMass = -100;
                    svFitPt = -100;
                    svFitEta = -100;
                    svFitPhi = -100;
                    svFitMET = -100;
                    svFitTransverseMass = -100;
                }
                //std::cout << "\n\nex: " << mvametcorr_ex << "   ey: " << mvametcorr_ey <<  " phi: " << mvametphi<<"\n"<<std::endl; 
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
            if(source->GetName()!=parser.stringValue("outputFile")){
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
            if(gSystem->AccessPathName(parser.stringValue("outputFile").c_str())) {
                gSystem->CopyFile("hsimple.root", parser.stringValue("outputFile").c_str());
            }

            fNew->cd();
            CopyFile(parser.stringValue("outputFile").c_str(),parser);
            fNew->ls();
            fNew->Close();

        }

