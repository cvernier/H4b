#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <TMath.h>

#include "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/bJetRegression/TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

#include "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/bJetRegression/HelperFunctions.h"


bool rescaleEnergy=true;

double sigmaJECUnc=0; // (-1, 0, 1)
double sigmaJERUnc=0; // (-1, 0, 1)
double sigmaTrigSF=1; // (-1, 0, 1)
double pi=3.14159265358979;

double bTagCSV_tightCut=0.898;
double bTagCSV_mediumCut=0.679;
double bTagCSV_looseCut=0.244;
double bTagCSV_noCut=0.;

double jetpT_cut=40.;
double jeteta_cut=2.5;
double H_mass=125.0;
double mH_diff_cut=50.;
double mH_mean_cut=15.;
double HpT_cut=0.;
double jetCSV_cut=bTagCSV_mediumCut;

typedef struct
{
  float et;
  float sumet;
  float sig;
  float phi;
} METInfo;

typedef struct
{
  float CSV;
  float E;
  float pT;
  float eta;
  float phi;
} JetInfo;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} GenParticleInfo;

JetInfo jet_CSV1={-1,-1,-1,-1,-1};
JetInfo jet_CSV2={-1,-1,-1,-1,-1};
JetInfo jet_CSV3={-1,-1,-1,-1,-1};
JetInfo jet_CSV4={-1,-1,-1,-1,-1};
JetInfo jet_CSV5={-1,-1,-1,-1,-1};

typedef std::map<double, int> JetList;

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector jet_p4;
  jet_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return jet_p4;
} 
/*
double deltaPhi(double phi1, double phi2)
{
  double dphi=phi1-phi2;
  if (dphi<-pi) dphi=2*pi+dphi;
  if (dphi>pi) dphi=2*pi-dphi;
  return dphi;
}
*/
double azimuthalAngle(TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, TLorentzVector v4)
{
  TLorentzVector v12=v1+v2;
  TVector3 v12_boost=v12.BoostVector();
  TLorentzVector v1_unboosted=v1; v1_unboosted.Boost(-v12_boost);
  TLorentzVector v2_unboosted=v2; v2_unboosted.Boost(-v12_boost);
  TVector3 v1v2Plane=v1_unboosted.Vect().Cross(v12_boost);
  
  TLorentzVector v34=v3+v4;
  TVector3 v34_boost=v34.BoostVector();
  TLorentzVector v3_unboosted=v3; v3_unboosted.Boost(-v34_boost);
  TLorentzVector v4_unboosted=v4; v4_unboosted.Boost(-v34_boost);
  TVector3 v3v4Plane=v3_unboosted.Vect().Cross(v34_boost);
  
  double angle=v1v2Plane.Angle(v3v4Plane);
  return angle;
}


int helper_TrigBin(double csv)
{
  int i=-1;
  if (csv<0.3) i=1;
  else if (csv<0.65) i=2;
  else if (csv<0.9) i=3;
  else if (csv<0.95) i=4;
  else if (csv>=0.95) i=5;
  return i;
}

int helper_TrigBin_pT2(double pT)
{
  int i=-1;
  if (pT<40) i=1;
  else if (pT<60) i=2;
  else if (pT<80) i=3;
  else if (pT<100) i=4;
  else if (pT<120) i=5;
  else if (pT<150) i=6;
  else i=7;
  return i;
}

int helper_TrigBin_pT4(double pT)
{
  int i=-1;
  if (pT<30) i=1;
  else if (pT<40) i=2;
  else if (pT<60) i=3;
  else if (pT<80) i=4;
  else if (pT<100) i=5;
  else if (pT<150) i=6;
  else i=7;
  return i;
}


TH2F *h_udsg_JetpTeta_den, *h_udsg_JetpTeta_num;
TH2F *h_c_JetpTeta_den, *h_c_JetpTeta_num;
TH2F *h_b_JetpTeta_den, *h_b_JetpTeta_num;

double CSVEfficiency(std::string j, double jetpT, double jeteta, double jetCSV)
{
  double weight, num, den;
  
  if (j=="j")
  {
    if (jetCSV>jetCSV_cut) weight=1.; else weight=0;
  }
  else if (j=="Q")
  {
    num=h_udsg_JetpTeta_num->GetBinContent(h_udsg_JetpTeta_num->FindBin(jetpT, jeteta));
    den=h_udsg_JetpTeta_den->GetBinContent(h_udsg_JetpTeta_den->FindBin(jetpT, jeteta));
    weight=num/den;
  }
  else if (j=="C")
  {
    num=h_c_JetpTeta_num->GetBinContent(h_c_JetpTeta_num->FindBin(jetpT, jeteta));
    den=h_c_JetpTeta_den->GetBinContent(h_c_JetpTeta_den->FindBin(jetpT, jeteta));
    weight=num/den;
  }
  else if (j=="B")
  {
    num=h_b_JetpTeta_num->GetBinContent(h_b_JetpTeta_num->FindBin(jetpT, jeteta));
    den=h_b_JetpTeta_den->GetBinContent(h_b_JetpTeta_den->FindBin(jetpT, jeteta));
    weight=num/den;
  }
  
  return weight;
}

bool returnRandom()
{
  double rand1=(double)rand()/double(RAND_MAX);
  bool returnValue;
  if (rand1<0.5) returnValue=false; else returnValue=true;
  return returnValue;
}

void fillStringsRandomly(std::string templ, std::string &j1, std::string &j2, std::string &j3, std::string &j4)
{
  if (templ=="QQjj") {j1="Q"; j2="Q"; j3="j"; j4="j";}
  if (templ=="QCjj") {j1="Q"; j2="C"; j3="j"; j4="j";}
  if (templ=="QBjj") {j1="Q"; j2="B"; j3="j"; j4="j";}
  if (templ=="CCjj") {j1="C"; j2="C"; j3="j"; j4="j";}
  if (templ=="CBjj") {j1="C"; j2="B"; j3="j"; j4="j";}
  if (templ=="BBjj") {j1="B"; j2="B"; j3="j"; j4="j";}
  if (templ=="QjQj") {j1="Q"; j2="j"; j3="Q"; j4="j";}
  if (templ=="QjCj") {j1="Q"; j2="j"; j3="C"; j4="j";}
  if (templ=="QjBj") {j1="Q"; j2="j"; j3="B"; j4="j";}
  if (templ=="CjCj") {j1="C"; j2="j"; j3="C"; j4="j";}
  if (templ=="CjBj") {j1="C"; j2="j"; j3="B"; j4="j";}
  if (templ=="BjBj") {j1="B"; j2="j"; j3="B"; j4="j";}
  
  // Now randomize between Higgs and between jets
  if (returnRandom())
  {
    swap(j1, j3);
    swap(j2, j4);
  }
  if (returnRandom())
  {
    swap(j1, j2);
  }
  if (returnRandom())
  {
    swap(j3, j4);
  }
}

std::string returnFlavorString(int flavor)
{
  std::string flav="j";
  if (fabs(flavor)==1 || fabs(flavor)==2 || fabs(flavor)==3 || fabs(flavor)==21) flav="Q";
  if (fabs(flavor)==4) flav="C";
  if (fabs(flavor)==5) flav="B";
  return flav;
}

int withinRegion(double mH1, double mH2, double r1=15., double r2=30., double mH1_c=H_mass, double mH2_c=H_mass)
{
  double r=pow(pow(mH1-mH1_c, 2)+pow(mH2-mH2_c, 2), 0.5);
  double angle=atan2(mH2-mH2_c, mH1-mH1_c);
  // std::cout<<"(mH1, mH2) = ("<<mH1<<", "<<mH2<<") lies in region ";
  int ret=-1;
  if (r<r1) ret=0;
  else if (r>r1 && r<r2)
  {
    if (angle>=0 && angle<pi/2.) ret=1;
    else if (angle>=pi/2. && angle<pi) ret=4;
    else if (angle<0 && angle>=-pi/2.) ret=2;
    else if (angle<pi/2.&& angle>=-pi) ret=3;
    else std::cout<<"This is within annulus but not within any CR!"<<std::endl;
  }
  else ret=5;
  // std::cout<<ret<<std::endl;
  return ret;
}

int HbbHbb_bTagging(std::string templ, std::string sample, std::string tagString, std::string csvReshape="", std::string PUWeight="", int bJetRegression=0)
{
  std::cout<<"H mass set at "<<H_mass<<std::endl;
  
  std::cout<<"templ = "<<templ<<std::endl;
  if (templ=="QCD")
  {
    std::string csvEfficienciesFilename="/Users/souvik/HbbHbb/Analysis/CSVFlavorEfficiencies_M.root";
    TFile *csvFlavorEfficiencies=new TFile(csvEfficienciesFilename.c_str());
    std::cout<<"Opened "<<csvEfficienciesFilename<<" for CSV efficiencies."<<std::endl;
    h_udsg_JetpTeta_den=(TH2F*)csvFlavorEfficiencies->Get("h_udsg_JetpTeta_den");
    h_udsg_JetpTeta_num=(TH2F*)csvFlavorEfficiencies->Get("h_udsg_JetpTeta_num");
    h_c_JetpTeta_den=(TH2F*)csvFlavorEfficiencies->Get("h_c_JetpTeta_den");
    h_c_JetpTeta_num=(TH2F*)csvFlavorEfficiencies->Get("h_c_JetpTeta_num");
    h_b_JetpTeta_den=(TH2F*)csvFlavorEfficiencies->Get("h_b_JetpTeta_den");
    h_b_JetpTeta_num=(TH2F*)csvFlavorEfficiencies->Get("h_b_JetpTeta_num");
    h_udsg_JetpTeta_den->Rebin2D(2,2);
    h_udsg_JetpTeta_num->Rebin2D(2,2);
    h_c_JetpTeta_den->Rebin2D(2,2);
    h_c_JetpTeta_num->Rebin2D(2,2);
    h_b_JetpTeta_den->Rebin2D(2,2);
    h_b_JetpTeta_num->Rebin2D(2,2);
  }
  std::string inputfilename="../"+sample+"_selected.root"; 
  //std::string inputfilename="/Users/souvik/HbbHbb/Analysis/OfficialStep2_KinematicallySelected_"+sample+".root";
  TChain *tree=new TChain("tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
  
  TFile *file_PUWeight;
  TH1F *h_PUWeight;
  if (PUWeight!="")
  {
    file_PUWeight=new TFile("/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/PUWeight.root");
    std::cout<<"Opened PU weight file = /Users/souvik/HbbHbb/Analysis/PUWeight.root"<<std::endl;
    h_PUWeight=(TH1F*)gDirectory->Get("h_PUWeight");
  }
 


   TFile *file_triggerSF=new TFile("/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/ratio_2D.root");

 
  std::cout<<"tagString = "<<tagString<<std::endl;
  
  // Book variables
  int vType;
  bool triggerFlags[500];
  int nPV;
  int nhJets, naJets;
  int nJets, nCJets;
  float hJetE[100], hJetpT[100], hJeteta[100], hJetphi[100], hJetCSV[100], hJetflavor[100], hJetpTRaw[100], hJet_ptLeadTrack[100], hJet_vtx3dL[100], hJet_vtx3deL[100], hJet_vtxMass[100], hJet_vtxPt[100], hJet_cef[100], hJet_nconstituents[100], hJet_JECUnc[100], hJet_genpT[100], hJet_SoftLeptptRel[100], hJet_SoftLeptPt[100], hJet_SoftLeptdR[100];
  float aJetE[100], aJetpT[100], aJeteta[100], aJetphi[100], aJetCSV[100], aJetflavor[100], aJetpTRaw[100], aJet_ptLeadTrack[100], aJet_vtx3dL[100], aJet_vtx3deL[100], aJet_vtxMass[100], aJet_vtxPt[100], aJet_cef[100], aJet_nconstituents[100], aJet_JECUnc[100], aJet_genpT[100], aJet_SoftLeptptRel[100], aJet_SoftLeptPt[100], aJet_SoftLeptdR[100];
  float jetE[100], jetpT[100], jeteta[100], jetphi[100], jetCSV[100], jetflavor[100], jetpTRaw[100], jet_ptLeadTrack[100], jet_vtx3dL[100], jet_vtx3deL[100], jet_vtxMass[100], jet_vtxPt[100], jet_cef[100], jet_nconstituents[100], jet_JECUnc[100], jet_genpT[100], jet_SoftLeptptRel[100], jet_SoftLeptPt[100], jet_SoftLeptdR[100];
  int hJet_SoftLeptIdlooseMu[100], hJet_SoftLeptId95[100];
  int aJet_SoftLeptIdlooseMu[100], aJet_SoftLeptId95[100];
  int jet_SoftLeptIdlooseMu[100], jet_SoftLeptId95[100];
  METInfo metObj;
  float met, metSig, ht;
  double weightPU;
  GenParticleInfo genX, genH1, genH2;
  int H1jet1_i, H1jet2_i;
  int H2jet1_i, H2jet2_i;
  int nJetCSV=0;
  float JetE[100], JetPt[100], JetEta[100], JetPhi[100], JetCSV[100];
  int  PtInd[4]; 
  float this_jetpT, this_jetpTRaw, this_et, this_mt, this_jet_eta, this_jet_phi, this_jet_e, 
        this_jet_ptLeadTrack, this_jet_vtx3dL, this_jet_vtx3deL, this_jet_vtxMass, this_jet_vtxPt, 
        this_jet_cef, this_jet_nconstituents, this_jet_JECUnc, 
        this_jet_SoftLeptptRel, this_jet_SoftLeptPt, this_jet_SoftLeptdR,
        this_jet_genpT;
  
  TMVA::Reader *reader=new TMVA::Reader("!Color:!Silent");
  if (bJetRegression==1)
  {
    reader->AddVariable("breg_pt := jetpT", &this_jetpT);
    reader->AddVariable("breg_rawptJER := smear_pt_res(jetpTRaw, jet_genpT, jet_eta)", &this_jetpTRaw);
    reader->AddVariable("breg_et := evalEt(jetpT, jet_eta, jet_phi, jet_e)", &this_et);
    reader->AddVariable("breg_mt := evalMt(jetpT, jet_eta, jet_phi, jet_e)", &this_mt);
    reader->AddVariable("breg_leadtrackpt := max(0, jet_ptLeadTrack)",&this_jet_ptLeadTrack);
    reader->AddVariable("breg_vtx3dL := max(0, jet_vtx3dL)", &this_jet_vtx3dL);
    reader->AddVariable("breg_vtx3deL := max(0, jet_vtx3deL)", &this_jet_vtx3deL);
    reader->AddVariable("breg_vtxMass := max(0, jet_vtxMass)", &this_jet_vtxMass);
    reader->AddVariable("breg_vtxPt := max(0, jet_vtxPt)", &this_jet_vtxPt);
    reader->AddVariable("breg_cef := jet_cef", &this_jet_cef);
    reader->AddVariable("breg_nconstituents := jet_nconstituents", &this_jet_nconstituents);
    reader->AddVariable("breg_JECUnc := jet_JECUnc", &this_jet_JECUnc);
    reader->AddVariable("breg_softlepptrel := max(0, jet_SoftLeptptRel*(jet_SoftLeptIdlooseMu==1 || jet_SoftLeptId95==1))", &this_jet_SoftLeptptRel);
    reader->AddVariable("breg_softleppt := max(0, jet_SoftLeptPt*(jet_SoftLeptIdlooseMu==1 || jet_SoftLeptId95==1))", &this_jet_SoftLeptPt);
    reader->AddVariable("breg_softlepdR := max(0, jet_SoftLeptdR*(jet_SoftLeptIdlooseMu==1 || jet_SoftLeptId95==1))", &this_jet_SoftLeptdR);
    // reader->BookMVA("BDTG method", "/Users/souvik/HbbHbb/Analysis/bJetRegression/weights/TMVARegression_BDTG.weights.xml");
    reader->BookMVA("BDTG method", "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/bJetRegression/weights/TMVARegression_BDTG.weights.xml");
  }
  
  // Retrieve variables
  tree->SetBranchAddress("Vtype", &(vType));
  tree->SetBranchAddress("triggerFlags", &(triggerFlags));
  tree->SetBranchAddress("nPVs", &(nPV));
  tree->SetBranchAddress("nhJets", &(nhJets));
  tree->SetBranchAddress("hJet_e", &(hJetE)); 
  tree->SetBranchAddress("hJet_pt", &(hJetpT));
  tree->SetBranchAddress("hJet_eta", &(hJeteta));
  tree->SetBranchAddress("hJet_phi", &(hJetphi));
  tree->SetBranchAddress(("hJet_csv_"+csvReshape).c_str(), &(hJetCSV));
  tree->SetBranchAddress("hJet_flavour", &(hJetflavor));
  tree->SetBranchAddress("hJet_ptRaw", &(hJetpTRaw));
  tree->SetBranchAddress("hJet_ptLeadTrack", &(hJet_ptLeadTrack));
  tree->SetBranchAddress("hJet_vtx3dL", &(hJet_vtx3dL));
  tree->SetBranchAddress("hJet_vtx3deL", &(hJet_vtx3deL));
  tree->SetBranchAddress("hJet_vtxMass", &(hJet_vtxMass));
  tree->SetBranchAddress("hJet_vtxPt", &(hJet_vtxPt));
  tree->SetBranchAddress("hJet_cef", &(hJet_cef));
  tree->SetBranchAddress("hJet_nconstituents", &(hJet_nconstituents));
  tree->SetBranchAddress("hJet_JECUnc", &(hJet_JECUnc));
  tree->SetBranchAddress("hJet_genPt", &(hJet_genpT));
  tree->SetBranchAddress("naJets", &(naJets));
  tree->SetBranchAddress("aJet_e", &(aJetE));
  tree->SetBranchAddress("aJet_pt", &(aJetpT));
  tree->SetBranchAddress("aJet_eta", &(aJeteta));
  tree->SetBranchAddress("aJet_phi", &(aJetphi));
  tree->SetBranchAddress(("aJet_csv_"+csvReshape).c_str(), &(aJetCSV));
  tree->SetBranchAddress("aJet_flavour", &(aJetflavor));
  tree->SetBranchAddress("aJet_ptRaw", &(aJetpTRaw));
  tree->SetBranchAddress("aJet_ptLeadTrack", &(aJet_ptLeadTrack));
  tree->SetBranchAddress("aJet_vtx3dL", &(aJet_vtx3dL));
  tree->SetBranchAddress("aJet_vtx3deL", &(aJet_vtx3deL));
  tree->SetBranchAddress("aJet_vtxMass", &(aJet_vtxMass));
  tree->SetBranchAddress("aJet_vtxPt", &(aJet_vtxPt));
  tree->SetBranchAddress("aJet_cef", &(aJet_cef));
  tree->SetBranchAddress("aJet_nconstituents", &(aJet_nconstituents));
  tree->SetBranchAddress("aJet_JECUnc", &(aJet_JECUnc));
  tree->SetBranchAddress("aJet_genPt", &(aJet_genpT));
  tree->SetBranchAddress("MET", &(metObj));
  tree->SetBranchAddress("genX", &(genX));
  tree->SetBranchAddress("genH1", &(genH1));
  tree->SetBranchAddress("genH2", &(genH2));
  tree->SetBranchAddress("H1jet1_i", &H1jet1_i);
  tree->SetBranchAddress("H1jet2_i", &H1jet2_i);
  tree->SetBranchAddress("H2jet1_i", &H2jet1_i);
  tree->SetBranchAddress("H2jet2_i", &H2jet2_i);
  tree->SetBranchAddress("nJetCSV", &(nJetCSV));
  tree->SetBranchAddress("JetE", &(JetE));
  tree->SetBranchAddress("JetPt", &(JetPt));
  tree->SetBranchAddress("JetEta", &(JetEta));
  tree->SetBranchAddress("JetPhi", &(JetPhi));
  tree->SetBranchAddress("JetCSV", &(JetCSV));
  tree->SetBranchAddress("PtInd",&(PtInd));

  
  TH1F *h_nJets_3Tags=new TH1F("h_nJets_3Tags", "# Cleaned PAT Jets; n", 20, 0., 20.);
  
  TH1F *h_CSV1=new TH1F("h_CSV1", "CSV1", 50, 0., 1.);
  TH1F *h_CSV2=new TH1F("h_CSV2", "CSV2", 50, 0., 1.);
  TH1F *h_CSV3=new TH1F("h_CSV3", "CSV3", 50, 0., 1.);
  TH1F *h_CSV4=new TH1F("h_CSV4", "CSV4", 50, 0., 1.);
  
  TH1F *h_MET=new TH1F("h_MET", "MET; MET (GeV)", 100, 0, 200.);
  TH1F *h_MET_sig=new TH1F("h_MET_sig", "MET Significance; Sig", 20, 0., 20.);
  
  TH1F *h_H1_mass_reg=new TH1F("h_H1_mass_reg", "H1 mass (after regression); mass (GeV)", 50, 50., 200.);
  TH1F *h_H2_mass_reg=new TH1F("h_H2_mass_reg", "H2 mass (after regression); mass (GeV)", 50, 50., 200.);
  TH1F *h_H1_mass_bTagged=new TH1F("h_H1_mass_bTagged", "H1 mass (after b-tagging); mass (GeV)", 50, 50., 200.);
  TH1F *h_H2_mass_bTagged=new TH1F("h_H2_mass_bTagged", "H2 mass (after b-tagging); mass (GeV)", 50, 50., 200.);
  
  TH1F *h_H1_pT_bTagged=new TH1F("h_H1_pT_bTagged", "H1 p_{T} (after b-tagging); p_{T} (GeV/c)", 50, 0., 800.);
  TH1F *h_H2_pT_bTagged=new TH1F("h_H2_pT_bTagged", "H2 p_{T} (after b-tagging); p_{T} (GeV/c)", 50, 0., 800.);
  TH1F *h_H1mH2_pT_bTagged=new TH1F("h_H1mH2_pT_bTagged", "| H1 p_{T} - H2 p_{T}| (after b-tagging); p_{T} (GeV/c)", 50, 0., 800.);
  
  TH1F *h_H1_CSV1=new TH1F("h_H1_CSV1", "H1 CSV1", 50, 0., 1.);
  TH1F *h_H1_CSV2=new TH1F("h_H1_CSV2", "H1 CSV2", 50, 0., 1.);
  TH1F *h_H2_CSV1=new TH1F("h_H2_CSV1", "H2 CSV1", 50, 0., 1.);
  TH1F *h_H2_CSV2=new TH1F("h_H2_CSV2", "H2 CSV2", 50, 0., 1.);
  
  TH1F *h_X_mass=new TH1F("h_X_mass", "X mass; mass (GeV)", 100, 0., 1200.);
  TH1F *h_X_pT=new TH1F("h_X_pT", "X p_{T}; p_{T} (GeV)", 100, 0., 500.);
  TH1F *h_X_mass_bTagged=new TH1F("h_X_mass_bTagged", "X mass (after b-tagging); mass (GeV)", 100, 0., 1200.);
  TH1F *h_X_pT_bTagged=new TH1F("h_X_pT_bTagged", "X p_{T} (after b-tagging); p_{T} (GeV)", 100, 0., 500.);
  
  TH1F *h_HH_dPhi=new TH1F("h_HH_dPhi", "#Delta#Phi(H,H); #Delta#Phi(H,H)", 20, 0., pi);
  TH1F *h_HH_deta=new TH1F("h_HH_deta", "#Delta#eta(H,H); #Delta#eta(H,H)", 20, 0., pi);
  TH1F *h_Azimuth=new TH1F("h_Azimuth", "Azimuthal angle between H(bb) and H(bb) decay planes", 10, 0., 1.);
  
  TH2F *h_mH_mH=new TH2F("h_mH_mH", "h_mH_mH", 50, 0., 400., 50, 0., 400.);
  
  TH1F *h_H1Jet1pT=new TH1F("h_H1Jet1pT", "h_H1Jet1pT", 50, 0., 200.);
  TH1F *h_H1Jet2pT=new TH1F("h_H1Jet2pT", "h_H1Jet2pT", 50, 0., 200.);
  TH1F *h_H1sumJetpT=new TH1F("h_H1sumJetpT", "h_H1sumJetpT", 50, 0., 500.);
  TH1F *h_H2sumJetpT=new TH1F("h_H2sumJetpT", "h_H2sumJetpT", 50, 0., 500.);
  TH1F *h_H1sumJetpTMinusH2sumJetpT=new TH1F("h_H1sumJetpTMinusH2sumJetpT", "h_H1sumJetpTMinusH2sumJetpT", 50, 0., 200.);
  
  TH1F *h_nTags=new TH1F("h_nTags", "h_nTags", 10, 0., 10.);
  
  TH3F *h_mH_mH_mX=new TH3F("h_mH_mH_mX", "mH vs mH", 50, 0., 200., 50, 0., 200., 20, 200., 800.);
  
  TH2F *h_regpT_genpT_1=new TH2F("h_regpT_genpT_1", "Reg pT vs gen pT", 50, 0., 500., 50, 0., 500.);
  TH2F *h_regpT_genpT_2=new TH2F("h_regpT_genpT_2", "Reg pT vs gen pT", 50, 0., 500., 50, 0., 500.);
  TH2F *h_regpT_genpT_3=new TH2F("h_regpT_genpT_3", "Reg pT vs gen pT", 50, 0., 500., 50, 0., 500.);
  TH2F *h_regpT_genpT_4=new TH2F("h_regpT_genpT_4", "Reg pT vs gen pT", 50, 0., 500., 50, 0., 500.);
  
  TH2F *h_RecoGenpT_pT_1_b=new TH2F("h_RecoGenpT_pT_1_b", "(Reco pT - Gen pT) / GenpT ~ Jet 1", 50, 0., 500., 50, -1., 1.);
  TH2F *h_RecoGenpT_pT_2_b=new TH2F("h_RecoGenpT_pT_2_b", "(Reco pT - Gen pT) / GenpT ~ Jet 2", 50, 0., 500., 50, -1., 1.);
  TH2F *h_RecoGenpT_pT_3_b=new TH2F("h_RecoGenpT_pT_3_b", "(Reco pT - Gen pT) / GenpT ~ Jet 3", 50, 0., 500., 50, -1., 1.);
  TH2F *h_RecoGenpT_pT_4_b=new TH2F("h_RecoGenpT_pT_4_b", "(Reco pT - Gen pT) / GenpT ~ Jet 4", 50, 0., 500., 50, -1., 1.);
  
  TH2F *h_RegGenpT_pT_1_b=new TH2F("h_RegGenpT_pT_1_b", "(Reg pT - Gen pT) / GenpT ~ Jet 1; Gen pT (GeV); (Reg pT - Gen pT) / GenpT", 50, 0., 500., 50, -1., 1.);
  TH2F *h_RegGenpT_pT_2_b=new TH2F("h_RegGenpT_pT_2_b", "(Reg pT - Gen pT) / GenpT ~ Jet 2; Gen pT (GeV); (Reg pT - Gen pT) / GenpT", 50, 0., 500., 50, -1., 1.);
  TH2F *h_RegGenpT_pT_3_b=new TH2F("h_RegGenpT_pT_3_b", "(Reg pT - Gen pT) / GenpT ~ Jet 3; Gen pT (GeV); (Reg pT - Gen pT) / GenpT", 50, 0., 500., 50, -1., 1.);
  TH2F *h_RegGenpT_pT_4_b=new TH2F("h_RegGenpT_pT_4_b", "(Reg pT - Gen pT) / GenpT ~ Jet 4; Gen pT (GeV); (Reg pT - Gen pT) / GenpT", 50, 0., 500., 50, -1., 1.);
  
  TH2F *h_RegmGenpT_pT_1_b=new TH2F("h_RegmGenpT_pT_1_b", "(Reg pT - Gen pT)~ Jet 1", 50, 0., 500., 20, -200., 200.);
  TH2F *h_RegmGenpT_pT_2_b=new TH2F("h_RegmGenpT_pT_2_b", "(Reg pT - Gen pT)~ Jet 2", 50, 0., 500., 20, -200., 200.);
  TH2F *h_RegmGenpT_pT_3_b=new TH2F("h_RegmGenpT_pT_3_b", "(Reg pT - Gen pT)~ Jet 3", 50, 0., 500., 20, -200., 200.);
  TH2F *h_RegmGenpT_pT_4_b=new TH2F("h_RegmGenpT_pT_4_b", "(Reg pT - Gen pT)~ Jet 4", 50, 0., 500., 20, -200., 200.);
  
  TH2F *h_diffpT_pT_1=new TH2F("h_diffpT_pT_1", "Reg pT - Reco pT vs pT for H1 Jet 1; Reco pT; Reg pT - Reco pT", 50, 0., 500., 20, -50., 50.);
  TH2F *h_diffpT_pT_2=new TH2F("h_diffpT_pT_2", "Reg pT - Reco pT vs pT for H1 Jet 2; Reco pT; Reg pT - Reco pT", 50, 0., 500., 20, -50., 50.);
  TH2F *h_diffpT_pT_3=new TH2F("h_diffpT_pT_3", "Reg pT - Reco pT vs pT for H2 Jet 1; Reco pT; Reg pT - Reco pT", 50, 0., 500., 20, -50., 50.);
  TH2F *h_diffpT_pT_4=new TH2F("h_diffpT_pT_4", "Reg pT - Reco pT vs pT for H2 Jet 2; Reco pT; Reg pT - Reco pT", 50, 0., 500., 20, -50., 50.);
  
  TH2F *h_diffpT_pT_1_b=new TH2F("h_diffpT_pT_1_b", "Reg pT - Reco pT vs pT for H1 Jet 1 b-tagged; Reco pT; Reg pT - Reco pT", 50, 0., 500., 50, -1., 1.);
  TH2F *h_diffpT_pT_2_b=new TH2F("h_diffpT_pT_2_b", "Reg pT - Reco pT vs pT for H1 Jet 2 b-tagged; Reco pT; Reg pT - Reco pT", 50, 0., 500., 50, -1., 1.);
  TH2F *h_diffpT_pT_3_b=new TH2F("h_diffpT_pT_3_b", "Reg pT - Reco pT vs pT for H2 Jet 1 b-tagged; Reco pT; Reg pT - Reco pT", 50, 0., 500., 50, -1., 1.);
  TH2F *h_diffpT_pT_4_b=new TH2F("h_diffpT_pT_4_b", "Reg pT - Reco pT vs pT for H2 Jet 2 b-tagged; Reco pT; Reg pT - Reco pT", 50, 0., 500., 50, -1., 1.);
  
  TH2F *h_dpT_HpT=new TH2F("h_dpT_HpT", "#Deltap_{T}(H_{1}, H_{2}) vs <H p_{T}>; <H p_{T}> (GeV); #Deltap_{T}(H_{1}, H_{2}) (GeV)", 50, 0., 500., 50, -2., 2.);
  
  
// Get the trigger efficiency data/MC scale factor
  TH1F *h_trigSF_CSV=(TH1F*)(((TCanvas*)file_triggerSF->Get("c5"))->GetPrimitive("Ratio"));
  std::cout<<"h_trigSF_CSV->GetBinContent(helper_TrigBin(0.679), helper_TrigBin(0.679)) = "<<h_trigSF_CSV->GetBinContent(helper_TrigBin(0.679), helper_TrigBin(0.679))<<std::endl;

  std::string histfilename="Histograms_"+sample+".root";
  gSystem->Exec(("cp ../"+histfilename+" "+histfilename).c_str());
  TFile *tFile1=new TFile(("../"+histfilename).c_str(), "READ");
  TH1F h_Cuts=*((TH1F*)((TH1F*)tFile1->Get("h_Cuts"))->Clone("h_Cuts"));
  tFile1->Close();
 
  std::string outfilename=sample+"_selected"+"_bTagged_.root"; 
  //std::string outfilename="OfficialStep2_KinematicallySelected_bTagged_"+sample+".root";
  TFile *outfile=new TFile(outfilename.c_str(), "recreate");
  TTree *outtree=tree->CloneTree(0);
  float regJet1E, regJet2E, regJet3E, regJet4E;
  float regJet1pT, regJet2pT, regJet3pT, regJet4pT;
  float eventWeight;	
  outtree->Branch("regJet1E", &regJet1E);
  outtree->Branch("regJet2E", &regJet2E);
  outtree->Branch("regJet3E", &regJet3E);
  outtree->Branch("regJet4E", &regJet4E); 
  outtree->Branch("regJet1pT", &regJet1pT);
  outtree->Branch("regJet2pT", &regJet2pT);
  outtree->Branch("regJet3pT", &regJet3pT);
  outtree->Branch("regJet4pT", &regJet4pT); 
  outtree->Branch("eventWeight", &eventWeight); 
  // Loop over events
  double n4Tags=0;
  int nEvents=tree->GetEntries();
  std::cout<<"Number of events in input file = "<<nEvents<<std::endl;
  double nCut6=0;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);
    
    weightPU=1.;
    if (PUWeight!="") weightPU=h_PUWeight->GetBinContent(h_PUWeight->FindBin(nPV));
    
    // Collate hJets and aJets into Jets
    nJets=nhJets+naJets;
    for (int j=0; j<nhJets; ++j)
    {
      jetE[j]=hJetE[j];
      jetpT[j]=hJetpT[j];
      jeteta[j]=hJeteta[j];
      jetphi[j]=hJetphi[j];
      jetCSV[j]=hJetCSV[j];
      jetflavor[j]=hJetflavor[j];
      jetpTRaw[j]=hJetpTRaw[j];
      jet_ptLeadTrack[j]=hJet_ptLeadTrack[j];
      jet_vtx3dL[j]=hJet_vtx3dL[j];
      jet_vtx3deL[j]=hJet_vtx3deL[j];
      jet_vtxMass[j]=hJet_vtxMass[j];
      jet_vtxPt[j]=hJet_vtxPt[j];
      jet_cef[j]=hJet_cef[j];
      jet_nconstituents[j]=hJet_nconstituents[j];
      jet_JECUnc[j]=hJet_JECUnc[j];
      jet_SoftLeptptRel[j]=hJet_SoftLeptptRel[j];
      jet_SoftLeptPt[j]=hJet_SoftLeptPt[j];
      jet_SoftLeptdR[j]=hJet_SoftLeptdR[j];
      jet_SoftLeptIdlooseMu[j]=hJet_SoftLeptIdlooseMu[j];
      jet_SoftLeptId95[j]=hJet_SoftLeptId95[j];
      jet_genpT[j]=hJet_genpT[j];
    }
    for (int j=0; j<naJets; ++j)
    {
      jetE[j+nhJets]=aJetE[j];
      jetpT[j+nhJets]=aJetpT[j];
      jeteta[j+nhJets]=aJeteta[j];
      jetphi[j+nhJets]=aJetphi[j];
      jetCSV[j+nhJets]=aJetCSV[j];
      jetflavor[j+nhJets]=aJetflavor[j];
      jetpTRaw[j+nhJets]=aJetpTRaw[j];
      jet_ptLeadTrack[j+nhJets]=aJet_ptLeadTrack[j];
      jet_vtx3dL[j+nhJets]=aJet_vtx3dL[j];
      jet_vtx3deL[j+nhJets]=aJet_vtx3deL[j];
      jet_vtxMass[j+nhJets]=aJet_vtxMass[j];
      jet_vtxPt[j+nhJets]=aJet_vtxPt[j];
      jet_cef[j+nhJets]=aJet_cef[j];
      jet_nconstituents[j+nhJets]=aJet_nconstituents[j];
      jet_JECUnc[j+nhJets]=aJet_JECUnc[j];
      jet_SoftLeptptRel[j+nhJets]=aJet_SoftLeptptRel[j];
      jet_SoftLeptPt[j+nhJets]=aJet_SoftLeptPt[j];
      jet_SoftLeptdR[j+nhJets]=aJet_SoftLeptdR[j];
      jet_SoftLeptIdlooseMu[j+nhJets]=aJet_SoftLeptIdlooseMu[j];
      jet_SoftLeptId95[j+nhJets]=aJet_SoftLeptId95[j];
      jet_genpT[j+nhJets]=aJet_genpT[j];
    }
    
    nCJets=0;
    for (int j=0; j<nJets; ++j)
    {
      if (fabs(jeteta[j])<jeteta_cut && jetpT[j]>jetpT_cut) ++nCJets;
    }
          
    TLorentzVector jet1_p4;      
    TLorentzVector jet2_p4;      
    TLorentzVector jet3_p4;      
    TLorentzVector jet4_p4;      
          
    double old_jet1pT=0, old_jet2pT=0, old_jet3pT=0, old_jet4pT=0;
    if (bJetRegression==0)
    {
      jet1_p4=fillTLorentzVector(JetPt[H1jet1_i], JetEta[H1jet1_i], JetPhi[H1jet1_i], JetE[H1jet1_i]);
      jet2_p4=fillTLorentzVector(JetPt[H1jet2_i], JetEta[H1jet2_i], JetPhi[H1jet2_i], JetE[H1jet2_i]);            
      jet3_p4=fillTLorentzVector(JetPt[H2jet1_i], JetEta[H2jet1_i], JetPhi[H2jet1_i], JetE[H2jet1_i]);             
      jet4_p4=fillTLorentzVector(JetPt[H2jet2_i], JetEta[H2jet2_i], JetPhi[H2jet2_i], JetE[H2jet2_i]);
    } 
    else
    {
      // -- Jet 1
      this_jetpT=jetpT[H1jet1_i]; old_jet1pT=this_jetpT;
      this_jet_ptLeadTrack=TMath::Max(float(0.), jet_ptLeadTrack[H1jet1_i]);
      this_jet_vtx3dL=TMath::Max(float(0.), jet_vtx3dL[H1jet1_i]);
      this_jet_vtx3deL=TMath::Max(float(0.), jet_vtx3deL[H1jet1_i]);
      this_jet_vtxMass=TMath::Max(float(0.), jet_vtxMass[H1jet1_i]);
      this_jet_vtxPt=TMath::Max(float(0.), jet_vtxPt[H1jet1_i]);
      this_jet_cef=jet_cef[H1jet1_i];
      this_jet_nconstituents=jet_nconstituents[H1jet1_i];
      this_jet_JECUnc=jet_JECUnc[H1jet1_i];
      this_jet_SoftLeptptRel=TMath::Max(float(0.), jet_SoftLeptptRel[H1jet1_i]*(jet_SoftLeptIdlooseMu[H1jet1_i]==1 || jet_SoftLeptId95[H1jet1_i]==1));
      this_jet_SoftLeptPt=TMath::Max(float(0.), jet_SoftLeptPt[H1jet1_i]*(jet_SoftLeptIdlooseMu[H1jet1_i]==1 || jet_SoftLeptId95[H1jet1_i]==1));
      this_jet_SoftLeptdR=TMath::Max(float(0.), jet_SoftLeptdR[H1jet1_i]*(jet_SoftLeptIdlooseMu[H1jet1_i]==1 || jet_SoftLeptId95[H1jet1_i]==1));
      if (templ=="Data")
      {
        this_jetpTRaw=jetpTRaw[H1jet1_i];
      } 
      else
      {
        this_jetpTRaw=smear_pt_res(jetpTRaw[H1jet1_i], jet_genpT[H1jet1_i], jeteta[H1jet1_i]);
      }
      this_et=evalEt(jetpT[H1jet1_i], jeteta[H1jet1_i], jetphi[H1jet1_i], jetE[H1jet1_i]);
      this_mt=evalMt(jetpT[H1jet1_i], jeteta[H1jet1_i], jetphi[H1jet1_i], jetE[H1jet1_i]);
      jetpT[H1jet1_i]=(reader->EvaluateRegression("BDTG method"))[0];
      // std::cout<<"i = "<<i<<", jetpT[H1jet1_i] = "<<jetpT[H1jet1_i]<<std::endl;
      jet1_p4=fillTLorentzVector(jetpT[H1jet1_i], jeteta[H1jet1_i], jetphi[H1jet1_i], (rescaleEnergy ? jetE[H1jet1_i]*jetpT[H1jet1_i]/this_jetpT : jetE[H1jet1_i]));
      
      
      // --- Jet 2
      this_jetpT=jetpT[H1jet2_i]; old_jet2pT=this_jetpT;
      this_jet_ptLeadTrack=TMath::Max(float(0.), jet_ptLeadTrack[H1jet2_i]);
      this_jet_vtx3dL=TMath::Max(float(0.), jet_vtx3dL[H1jet2_i]);
      this_jet_vtx3deL=TMath::Max(float(0.), jet_vtx3deL[H1jet2_i]);
      this_jet_vtxMass=TMath::Max(float(0.), jet_vtxMass[H1jet2_i]);
      this_jet_vtxPt=TMath::Max(float(0.), jet_vtxPt[H1jet2_i]);
      this_jet_cef=jet_cef[H1jet2_i];
      this_jet_nconstituents=jet_nconstituents[H1jet2_i];
      this_jet_JECUnc=jet_JECUnc[H1jet2_i];
      this_jet_SoftLeptptRel=TMath::Max(float(0.), jet_SoftLeptptRel[H1jet2_i]*(jet_SoftLeptIdlooseMu[H1jet2_i]==1 || jet_SoftLeptId95[H1jet2_i]==1));
      this_jet_SoftLeptPt=TMath::Max(float(0.), jet_SoftLeptPt[H1jet2_i]*(jet_SoftLeptIdlooseMu[H1jet2_i]==1 || jet_SoftLeptId95[H1jet2_i]==1));
      this_jet_SoftLeptdR=TMath::Max(float(0.), jet_SoftLeptdR[H1jet2_i]*(jet_SoftLeptIdlooseMu[H1jet2_i]==1 || jet_SoftLeptId95[H1jet2_i]==1));
      if (templ=="Data")
      {
        this_jetpTRaw=jetpTRaw[H1jet2_i];
      } 
      else
      {
        this_jetpTRaw=smear_pt_res(jetpTRaw[H1jet2_i], jet_genpT[H1jet2_i], jeteta[H1jet2_i]);
      }
      this_et=evalEt(jetpT[H1jet2_i], jeteta[H1jet2_i], jetphi[H1jet2_i], jetE[H1jet2_i]);
      this_mt=evalMt(jetpT[H1jet2_i], jeteta[H1jet2_i], jetphi[H1jet2_i], jetE[H1jet2_i]);
      jetpT[H1jet2_i]=(reader->EvaluateRegression("BDTG method"))[0];
      jet2_p4=fillTLorentzVector(jetpT[H1jet2_i], jeteta[H1jet2_i], jetphi[H1jet2_i], (rescaleEnergy ? jetE[H1jet2_i]*jetpT[H1jet2_i]/this_jetpT : jetE[H1jet2_i]));
            
      // --- Jet 3
      this_jetpT=jetpT[H2jet1_i]; old_jet3pT=this_jetpT;
      this_jet_ptLeadTrack=TMath::Max(float(0.), jet_ptLeadTrack[H2jet1_i]);
      this_jet_vtx3dL=TMath::Max(float(0.), jet_vtx3dL[H2jet1_i]);
      this_jet_vtx3deL=TMath::Max(float(0.), jet_vtx3deL[H2jet1_i]);
      this_jet_vtxMass=TMath::Max(float(0.), jet_vtxMass[H2jet1_i]);
      this_jet_vtxPt=TMath::Max(float(0.), jet_vtxPt[H2jet1_i]);
      this_jet_cef=jet_cef[H2jet1_i];
      this_jet_nconstituents=jet_nconstituents[H2jet1_i];
      this_jet_JECUnc=jet_JECUnc[H2jet1_i];
      this_jet_SoftLeptptRel=TMath::Max(float(0.), jet_SoftLeptptRel[H2jet1_i]*(jet_SoftLeptIdlooseMu[H2jet1_i]==1 || jet_SoftLeptId95[H2jet1_i]==1));
      this_jet_SoftLeptPt=TMath::Max(float(0.), jet_SoftLeptPt[H2jet1_i]*(jet_SoftLeptIdlooseMu[H2jet1_i]==1 || jet_SoftLeptId95[H2jet1_i]==1));
      this_jet_SoftLeptdR=TMath::Max(float(0.), jet_SoftLeptdR[H2jet1_i]*(jet_SoftLeptIdlooseMu[H2jet1_i]==1 || jet_SoftLeptId95[H2jet1_i]==1));
      if (templ=="Data")
      {
        this_jetpTRaw=jetpTRaw[H2jet1_i];
      } 
      else
      {
        this_jetpTRaw=smear_pt_res(jetpTRaw[H2jet1_i], jet_genpT[H2jet1_i], jeteta[H2jet1_i]);
      }
      this_et=evalEt(jetpT[H2jet1_i], jeteta[H2jet1_i], jetphi[H2jet1_i], jetE[H2jet1_i]);
      this_mt=evalMt(jetpT[H2jet1_i], jeteta[H2jet1_i], jetphi[H2jet1_i], jetE[H2jet1_i]);
      jetpT[H2jet1_i]=(reader->EvaluateRegression("BDTG method"))[0];
      jet3_p4=fillTLorentzVector(jetpT[H2jet1_i], jeteta[H2jet1_i], jetphi[H2jet1_i], (rescaleEnergy ? jetE[H2jet1_i]*jetpT[H2jet1_i]/this_jetpT : jetE[H2jet1_i]));
      
      // --- Jet 4
      this_jetpT=jetpT[H2jet2_i]; old_jet4pT=this_jetpT;
      this_jet_ptLeadTrack=TMath::Max(float(0.), jet_ptLeadTrack[H2jet2_i]);
      this_jet_vtx3dL=TMath::Max(float(0.), jet_vtx3dL[H2jet2_i]);
      this_jet_vtx3deL=TMath::Max(float(0.), jet_vtx3deL[H2jet2_i]);
      this_jet_vtxMass=TMath::Max(float(0.), jet_vtxMass[H2jet2_i]);
      this_jet_vtxPt=TMath::Max(float(0.), jet_vtxPt[H2jet2_i]);
      this_jet_cef=jet_cef[H2jet2_i];
      this_jet_nconstituents=jet_nconstituents[H2jet2_i];
      this_jet_JECUnc=jet_JECUnc[H2jet2_i];
      this_jet_SoftLeptptRel=TMath::Max(float(0.), jet_SoftLeptptRel[H2jet2_i]*(jet_SoftLeptIdlooseMu[H2jet2_i]==1 || jet_SoftLeptId95[H2jet2_i]==1));
      this_jet_SoftLeptPt=TMath::Max(float(0.), jet_SoftLeptPt[H2jet2_i]*(jet_SoftLeptIdlooseMu[H2jet2_i]==1 || jet_SoftLeptId95[H2jet2_i]==1));
      this_jet_SoftLeptdR=TMath::Max(float(0.), jet_SoftLeptdR[H2jet2_i]*(jet_SoftLeptIdlooseMu[H2jet2_i]==1 || jet_SoftLeptId95[H2jet2_i]==1));
      if (templ=="Data")
      {
        this_jetpTRaw=jetpTRaw[H2jet2_i];
      } 
      else
      {
        this_jetpTRaw=smear_pt_res(jetpTRaw[H2jet2_i], jet_genpT[H2jet2_i], jeteta[H2jet2_i]);
      }
      this_et=evalEt(jetpT[H2jet2_i], jeteta[H2jet2_i], jetphi[H2jet2_i], jetE[H2jet2_i]);
      this_mt=evalMt(jetpT[H2jet2_i], jeteta[H2jet2_i], jetphi[H2jet2_i], jetE[H2jet2_i]);
      jetpT[H2jet2_i]=(reader->EvaluateRegression("BDTG method"))[0];
      jet4_p4=fillTLorentzVector(jetpT[H2jet2_i], jeteta[H2jet2_i], jetphi[H2jet2_i], (rescaleEnergy ? jetE[H2jet2_i]*jetpT[H2jet2_i]/this_jetpT : jetE[H2jet2_i]));
    }     
        
    // Now evaluate the weight from the trigger efficiency data/MC scale factor 
    double H1CSVhi=TMath::Max(JetCSV[H1jet1_i], JetCSV[H1jet2_i]);
    double H1CSVlo=TMath::Min(JetCSV[H1jet1_i], JetCSV[H1jet2_i]);
    double H2CSVhi=TMath::Max(JetCSV[H2jet1_i], JetCSV[H2jet2_i]);
    double H2CSVlo=TMath::Min(JetCSV[H2jet1_i], JetCSV[H2jet2_i]);
    double HCSVhi=(H1CSVhi+H2CSVhi)/2.;
    double HCSVlo=(H1CSVlo+H2CSVlo)/2.;
    double triggerWeight=h_trigSF_CSV->GetBinContent(helper_TrigBin(HCSVhi), helper_TrigBin(HCSVlo));
    triggerWeight=triggerWeight+sigmaTrigSF*(triggerWeight-1.)/2.;

    if (templ!="Data")
    {
      weightPU=weightPU*triggerWeight;
    }	



  
    h_regpT_genpT_1->Fill(jet_genpT[H1jet1_i], jetpT[H1jet1_i], weightPU);
    h_regpT_genpT_2->Fill(jet_genpT[H1jet2_i], jetpT[H1jet2_i], weightPU);
    h_regpT_genpT_3->Fill(jet_genpT[H2jet1_i], jetpT[H2jet1_i], weightPU);
    h_regpT_genpT_4->Fill(jet_genpT[H2jet2_i], jetpT[H2jet2_i], weightPU);
    
    h_diffpT_pT_1->Fill(old_jet1pT, jetpT[H1jet1_i]-old_jet1pT, weightPU);
    h_diffpT_pT_2->Fill(old_jet2pT, jetpT[H1jet2_i]-old_jet2pT, weightPU);
    h_diffpT_pT_3->Fill(old_jet3pT, jetpT[H2jet1_i]-old_jet3pT, weightPU);
    h_diffpT_pT_4->Fill(old_jet4pT, jetpT[H2jet2_i]-old_jet4pT, weightPU);
    
    TLorentzVector H1_p4=jet1_p4+jet2_p4;
    TLorentzVector H2_p4=jet3_p4+jet4_p4;
    
    h_H1_mass_reg->Fill(H1_p4.M(), weightPU);
    h_H2_mass_reg->Fill(H2_p4.M(), weightPU);
    
    TLorentzVector X_p4=H1_p4+H2_p4;
    double X_mass=X_p4.M();
    double X_pT=X_p4.Pt();
    h_X_mass->Fill(X_mass, weightPU);
    
	h_X_pT->Fill(X_pT, weightPU);
    
    // Jet pT study
    h_H1_CSV1->Fill(jetCSV[H1jet1_i], weightPU);
    // if (jetCSV[H1jet1_i]>jetCSV_cut) 
    {
      h_H1_CSV2->Fill(jetCSV[H1jet2_i], weightPU);
      // if (jetCSV[H1jet2_i]>jetCSV_cut) 
      {
        h_H2_CSV1->Fill(jetCSV[H2jet1_i], weightPU);
        // if (jetCSV[H2jet1_i]>jetCSV_cut) 
        {
          h_H2_CSV2->Fill(jetCSV[H2jet2_i], weightPU);
          // if (jetCSV[H2jet2_i]>jetCSV_cut) 
          {
            h_H1sumJetpT->Fill(jet1_p4.Pt()+jet2_p4.Pt(), weightPU);
            h_H2sumJetpT->Fill(jet3_p4.Pt()+jet4_p4.Pt(), weightPU);
            h_H1sumJetpTMinusH2sumJetpT->Fill((jet1_p4.Pt()+jet2_p4.Pt())-(jet3_p4.Pt()+jet4_p4.Pt()), weightPU);
            h_H1Jet1pT->Fill(jet1_p4.Pt(), weightPU);
            h_H1Jet2pT->Fill(jet2_p4.Pt(), weightPU);
          }
        }
      }
    }
            
    // Cut on nTags
    int nTags=0; int jetTag[4]={0,0,0,0};
    int nTagsM=0;
    int nTagsL=0;
    
    if (JetCSV[H1jet1_i]>bTagCSV_looseCut) {++nTagsL;jetTag[0]=1;} 
    if (JetCSV[H1jet2_i]>bTagCSV_looseCut) {++nTagsL;jetTag[1]=1;} 
    if (JetCSV[H2jet1_i]>bTagCSV_looseCut) {++nTagsL;jetTag[2]=1;} 
    if (JetCSV[H2jet2_i]>bTagCSV_looseCut) {++nTagsL;jetTag[3]=1;} 
    
    if (JetCSV[H1jet1_i]>bTagCSV_mediumCut) {++nTagsM; jetTag[0]=2;}
    if (JetCSV[H1jet2_i]>bTagCSV_mediumCut) {++nTagsM; jetTag[1]=2;}
    if (JetCSV[H2jet1_i]>bTagCSV_mediumCut) {++nTagsM; jetTag[2]=2;}
    if (JetCSV[H2jet2_i]>bTagCSV_mediumCut) {++nTagsM; jetTag[3]=2;}
    
    if (JetCSV[H1jet1_i]>bTagCSV_tightCut) {++nTags; jetTag[0]=3;}
    if (JetCSV[H1jet2_i]>bTagCSV_tightCut) {++nTags; jetTag[1]=3;}
    if (JetCSV[H2jet1_i]>bTagCSV_tightCut) {++nTags; jetTag[2]=3;}
    if (JetCSV[H2jet2_i]>bTagCSV_tightCut) {++nTags; jetTag[3]=3;}
    
    h_nTags->Fill(nTagsM, weightPU);
          
    // CSV weight events if QCD
    if (templ=="QCD")
    {
      double csvweight1=CSVEfficiency(returnFlavorString(jetflavor[H1jet1_i]), jetpT[H1jet1_i], jeteta[H1jet1_i], jetCSV[H1jet1_i]);
      double csvweight2=CSVEfficiency(returnFlavorString(jetflavor[H1jet2_i]), jetpT[H1jet2_i], jeteta[H1jet2_i], jetCSV[H1jet2_i]);
      double csvweight3=CSVEfficiency(returnFlavorString(jetflavor[H2jet1_i]), jetpT[H2jet1_i], jeteta[H2jet1_i], jetCSV[H2jet1_i]);
      double csvweight4=CSVEfficiency(returnFlavorString(jetflavor[H2jet2_i]), jetpT[H2jet2_i], jeteta[H2jet2_i], jetCSV[H2jet2_i]);
      weightPU=weightPU*csvweight1*csvweight2*csvweight3*csvweight4;
    }
            
    if ( 
         (tagString=="MMM" && nTagsM>=3)                                              ||
         (tagString=="MMMM" && (nTags>=2&&nTagsL>=4))                                             ||
         (tagString=="TTM" && (nTagsM>=3&&nTags>=2))                                               ||
         (tagString=="TMTM" && ((jetTag[0]+jetTag[1])>4 && (jetTag[2]+jetTag[3])>4))  ||
         (tagString=="MJMJ" && ((jetTag[0]+jetTag[1])>1 && (jetTag[2]+jetTag[3])>1))
       )
    {
      ++nCut6;
      
      // The most tagged is tagged as H1
      // if (jetTag[0]==0 || jetTag[1]==0) swap(H1_p4, H2_p4);
      
      h_RecoGenpT_pT_1_b->Fill(jet_genpT[H1jet1_i], (old_jet1pT-jet_genpT[H1jet1_i])/jet_genpT[H1jet1_i], weightPU);
      h_RecoGenpT_pT_2_b->Fill(jet_genpT[H1jet2_i], (old_jet2pT-jet_genpT[H1jet2_i])/jet_genpT[H1jet2_i], weightPU);
      h_RecoGenpT_pT_3_b->Fill(jet_genpT[H2jet1_i], (old_jet3pT-jet_genpT[H2jet1_i])/jet_genpT[H2jet1_i], weightPU);
      h_RecoGenpT_pT_4_b->Fill(jet_genpT[H2jet2_i], (old_jet4pT-jet_genpT[H2jet2_i])/jet_genpT[H2jet2_i], weightPU);
      
      h_RegGenpT_pT_1_b->Fill(jet_genpT[H1jet1_i], (jetpT[H1jet1_i]-jet_genpT[H1jet1_i])/jet_genpT[H1jet1_i], weightPU);
      h_RegGenpT_pT_2_b->Fill(jet_genpT[H1jet2_i], (jetpT[H1jet2_i]-jet_genpT[H1jet2_i])/jet_genpT[H1jet2_i], weightPU);
      h_RegGenpT_pT_3_b->Fill(jet_genpT[H2jet1_i], (jetpT[H2jet1_i]-jet_genpT[H2jet1_i])/jet_genpT[H2jet1_i], weightPU);
      h_RegGenpT_pT_4_b->Fill(jet_genpT[H2jet2_i], (jetpT[H2jet2_i]-jet_genpT[H2jet2_i])/jet_genpT[H2jet2_i], weightPU);
      
      h_RegmGenpT_pT_1_b->Fill(jet_genpT[H1jet1_i], (jetpT[H1jet1_i]-jet_genpT[H1jet1_i]), weightPU);
      h_RegmGenpT_pT_2_b->Fill(jet_genpT[H1jet2_i], (jetpT[H1jet2_i]-jet_genpT[H1jet2_i]), weightPU);
      h_RegmGenpT_pT_3_b->Fill(jet_genpT[H2jet1_i], (jetpT[H2jet1_i]-jet_genpT[H2jet1_i]), weightPU);
      h_RegmGenpT_pT_4_b->Fill(jet_genpT[H2jet2_i], (jetpT[H2jet2_i]-jet_genpT[H2jet2_i]), weightPU);
      
      h_diffpT_pT_1_b->Fill(old_jet1pT, (jetpT[H1jet1_i]-old_jet1pT)/old_jet1pT, weightPU);
      h_diffpT_pT_2_b->Fill(old_jet2pT, (jetpT[H1jet2_i]-old_jet2pT)/old_jet2pT, weightPU);
      h_diffpT_pT_3_b->Fill(old_jet3pT, (jetpT[H2jet1_i]-old_jet3pT)/old_jet3pT, weightPU);
      h_diffpT_pT_4_b->Fill(old_jet4pT, (jetpT[H2jet2_i]-old_jet4pT)/old_jet4pT, weightPU);
              
      double H1_mass=H1_p4.M();
      double H2_mass=H2_p4.M();
      
      h_X_mass_bTagged->Fill(X_mass, weightPU);
      h_X_pT_bTagged->Fill(X_pT, weightPU);
      
      h_H1_mass_bTagged->Fill(H1_mass, weightPU);       
      h_H2_mass_bTagged->Fill(H2_mass, weightPU);
      h_H1_pT_bTagged->Fill(H1_p4.Pt(), weightPU);
      h_H2_pT_bTagged->Fill(H2_p4.Pt(), weightPU);
      h_H1mH2_pT_bTagged->Fill(fabs(H1_p4.Pt()-H2_p4.Pt()), weightPU);
      
      // Needed for b-jet regression validation
      if (nCJets==4)
      {
        h_dpT_HpT->Fill((H1_p4.Pt()+H2_p4.Pt())/2., 2.*(H1_p4.Pt()-H2_p4.Pt())/(H1_p4.Pt()+H2_p4.Pt()), weightPU);
      }
      
      double dPhi_H1H2=deltaPhi(H1_p4.Phi(), H2_p4.Phi());
      h_HH_dPhi->Fill(dPhi_H1H2, weightPU);
      double dEta_H1H2=fabs(H1_p4.Eta()-H2_p4.Eta());
      h_HH_deta->Fill(dEta_H1H2, weightPU);
      // double dR_H1H2=H1_p4.DeltaR(H2_p4);
      double azimuth_H1H2=fabs(cos(azimuthalAngle(jet1_p4, jet2_p4, jet3_p4, jet4_p4)));
      h_Azimuth->Fill(azimuth_H1H2, weightPU);
              
      h_mH_mH->Fill(H1_mass, H2_mass, weightPU);
      h_mH_mH_mX->Fill(H1_mass, H2_mass, X_mass, weightPU);
      
      // h_dRJets_mH->Fill(jet1_p4.DeltaR(jet2_p4), H1_mass, weightPU);
      
      // Fill extra variables in tree
      regJet1E=jet1_p4.E();
      regJet2E=jet2_p4.E();
      regJet3E=jet3_p4.E();
      regJet4E=jet4_p4.E();
      regJet1pT=jet1_p4.Pt();
      regJet2pT=jet2_p4.Pt();
      regJet3pT=jet3_p4.Pt();
      regJet4pT=jet4_p4.Pt();
      eventWeight=weightPU;
	outtree->Fill();
      
      
    } // nTags
      
    // std::cout<<"i = "<<i<<std::endl;
  } // Event loop
  
  std::cout<<"out of event loop"<<std::endl;
  std::cout<<"nCut6 = "<<nCut6<<std::endl;
  h_Cuts.Fill(13, nCut6);
  
  outtree->Write();
  outfile->Close();
  std::cout<<"Wrote output file "<<outfilename<<std::endl;
  
  TFile *tFile2=new TFile(histfilename.c_str(), "UPDATE");
  tFile2->Delete("h_Cuts;1");
  h_CSV1->Write();
  h_CSV2->Write();
  h_CSV3->Write();
  h_CSV4->Write();
  h_H1_mass_reg->Write();
  h_H2_mass_reg->Write();
  h_H1_mass_bTagged->Write();
  h_H2_mass_bTagged->Write();
  h_H1_pT_bTagged->Write();
  h_H2_pT_bTagged->Write();
  h_H1mH2_pT_bTagged->Write();
  h_H1_CSV1->Write();
  h_H1_CSV2->Write();
  h_H2_CSV1->Write();
  h_H2_CSV2->Write();
  h_X_mass->Write();
  h_X_pT->Write();
  h_X_mass_bTagged->Write();
  h_X_pT_bTagged->Write();
  h_HH_dPhi->Write();
  h_HH_deta->Write();
  h_dpT_HpT->Write();
  h_Azimuth->Write();
  h_mH_mH->Write();
  h_H1sumJetpT->Write();
  h_H2sumJetpT->Write();
  h_H1sumJetpTMinusH2sumJetpT->Write();
  h_H1Jet1pT->Write();
  h_H1Jet2pT->Write();
  h_nTags->Write();
  h_Cuts.Write();
  h_mH_mH_mX->Write();
  h_regpT_genpT_1->Write();
  h_regpT_genpT_2->Write();
  h_regpT_genpT_3->Write();
  h_regpT_genpT_4->Write();
  h_RecoGenpT_pT_1_b->Write();
  h_RecoGenpT_pT_2_b->Write();
  h_RecoGenpT_pT_3_b->Write();
  h_RecoGenpT_pT_4_b->Write();
  h_RegGenpT_pT_1_b->Write();
  h_RegGenpT_pT_2_b->Write();
  h_RegGenpT_pT_3_b->Write();
  h_RegGenpT_pT_4_b->Write();
  h_RegmGenpT_pT_1_b->Write();
  h_RegmGenpT_pT_2_b->Write();
  h_RegmGenpT_pT_3_b->Write();
  h_RegmGenpT_pT_4_b->Write();
  h_diffpT_pT_1->Write();
  h_diffpT_pT_2->Write();
  h_diffpT_pT_3->Write();
  h_diffpT_pT_4->Write();
  h_diffpT_pT_1_b->Write();
  h_diffpT_pT_2_b->Write();
  h_diffpT_pT_3_b->Write();
  h_diffpT_pT_4_b->Write();
  tFile2->Write();
  tFile2->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  
  
  std::cout<<"=== Cut Efficiencies === "<<std::endl;
  std::cout<<"Number of events after b-tagging = "<<(float)nCut6/nEvents<<std::endl;
  std::cout<<"========================"<<std::endl;
  
  return 0;
}
            
      
  
  
  
  
