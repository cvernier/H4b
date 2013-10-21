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

#include "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas//HbbHbb/Analysis/bJetRegression/HelperFunctions.h"

#include "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas//HbbHbb/Analysis/kinFit4b.h"

double sigmaJECUnc=0; // (-1, 0, 1)
double sigmaJERUnc=0; // (-1, 0, 1)

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
 

int HbbHbb_Partitioning(std::string sample, double h_mass, int kinConstraint=0, std::string PUWeight="")
{
  H_mass=h_mass;
  std::cout<<"H mass set at "<<H_mass<<std::endl;
	
	 gSystem->Load("libPhysicsToolsKinFitter.so");
 	gSystem->Load("/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/kinFit4b_h.so"); 
  std::string inputfilename="../"+sample+"_selected_bTagged_.root";//"../OfficialStep2_KinematicallySelected_bTagged_"+sample+".root";
  TChain *tree=new TChain("tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
  
  TFile *file_PUWeight;
  TH1F *h_PUWeight;
  if (PUWeight!="")
  {
    file_PUWeight=new TFile("/uscms_data/d2/souvik/HbbHbb/CMSSW_5_3_3_patch2/src/Analysis/PUWeight.root");
    std::cout<<"Opened PU weight file = /uscms_data/d2/souvik/HbbHbb/CMSSW_5_3_3_patch2/src/Analysis/PUWeight.root"<<std::endl;
    h_PUWeight=(TH1F*)gDirectory->Get("h_PUWeight");
  }
  
  // Book variables
  int vType;
  bool triggerFlags[500],QuadJetFilterFlag;
  int nPV;
  int nhJets, naJets;
  int nJets, nCJets;
  float hJetE[100], hJetpT[100], hJeteta[100], hJetphi[100], hJetCSV[100], hJetflavor[100], hJetpTRaw[100], hJet_ptLeadTrack[100], hJet_vtx3dL[100], hJet_vtx3deL[100], hJet_vtxMass[100], hJet_vtxPt[100], hJet_cef[100], hJet_nconstituents[100], hJet_JECUnc[100], hJet_genpT[100];
  float aJetE[100], aJetpT[100], aJeteta[100], aJetphi[100], aJetCSV[100], aJetflavor[100], aJetpTRaw[100], aJet_ptLeadTrack[100], aJet_vtx3dL[100], aJet_vtx3deL[100], aJet_vtxMass[100], aJet_vtxPt[100], aJet_cef[100], aJet_nconstituents[100], aJet_JECUnc[100], aJet_genpT[100];
  float jetE[100], jetpT[100], jeteta[100], jetphi[100], jetCSV[100], jetflavor[100], jetpTRaw[100], jet_ptLeadTrack[100], jet_vtx3dL[100], jet_vtx3deL[100], jet_vtxMass[100], jet_vtxPt[100], jet_cef[100], jet_nconstituents[100], jet_JECUnc[100], jet_genpT[100];
  METInfo metObj;
  float met, metSig, ht; float dR_min_5thJet;
  GenParticleInfo genX, genH1, genH2;
  int H1jet1_i, H1jet2_i;
  int H2jet1_i, H2jet2_i;
  int pass125125;
   float Ht_less4;	
  int pass9090; int matched;
  float regJet1E, regJet2E, regJet3E, regJet4E;
  float regJet1pT, regJet2pT, regJet3pT, regJet4pT;
  float eventWeight; float corr[4];
  // Retrieve variables
  tree->SetBranchAddress("dR_min_5thJet", &(dR_min_5thJet));
  tree->SetBranchAddress("QuadJetFilterFlag", &(QuadJetFilterFlag));  
  tree->SetBranchAddress("Vtype", &(vType));
  tree->SetBranchAddress("triggerFlags", &(triggerFlags));
  tree->SetBranchAddress("nPVs", &(nPV));
  tree->SetBranchAddress("corr", &(corr));
  tree->SetBranchAddress("nhJets", &(nhJets));
  tree->SetBranchAddress("Ht_less4", &(Ht_less4));
  tree->SetBranchAddress("hJet_e", &(hJetE)); 
  tree->SetBranchAddress("hJet_pt", &(hJetpT));
  tree->SetBranchAddress("hJet_eta", &(hJeteta));
  tree->SetBranchAddress("hJet_phi", &(hJetphi));
  tree->SetBranchAddress("hJet_csv", &(hJetCSV));
  tree->SetBranchAddress("hJet_flavour", &(hJetflavor));
  tree->SetBranchAddress("naJets", &(naJets));
  tree->SetBranchAddress("aJet_e", &(aJetE));
  tree->SetBranchAddress("aJet_pt", &(aJetpT));
  tree->SetBranchAddress("aJet_eta", &(aJeteta));
  tree->SetBranchAddress("aJet_phi", &(aJetphi));
  tree->SetBranchAddress("aJet_csv", &(aJetCSV));
  tree->SetBranchAddress("MET", &(metObj));
  tree->SetBranchAddress("genX", &(genX));
  tree->SetBranchAddress("genH1", &(genH1));
  tree->SetBranchAddress("genH2", &(genH2));
  tree->SetBranchAddress("H1jet1_i", &H1jet1_i);
  tree->SetBranchAddress("H1jet2_i", &H1jet2_i);
  tree->SetBranchAddress("H2jet1_i", &H2jet1_i);
  tree->SetBranchAddress("H2jet2_i", &H2jet2_i);
  tree->SetBranchAddress("regJet1E", &regJet1E);
  tree->SetBranchAddress("regJet2E", &regJet2E);
  tree->SetBranchAddress("regJet3E", &regJet3E);
  tree->SetBranchAddress("regJet4E", &regJet4E); 
  tree->SetBranchAddress("regJet1pT", &regJet1pT);
  tree->SetBranchAddress("regJet2pT", &regJet2pT);
  tree->SetBranchAddress("regJet3pT", &regJet3pT);
  tree->SetBranchAddress("regJet4pT", &regJet4pT); 
  tree->SetBranchAddress("pass125125", &pass125125);
  tree->SetBranchAddress("pass9090", &pass9090);
  tree->SetBranchAddress("matched", &matched);
   float JetE[100], JetPt[100], JetEta[100], JetPhi[100], JetCSV[100]; 
int nJetCSV;
tree->SetBranchAddress("nJetCSV", &(nJetCSV));
  tree->SetBranchAddress("JetE", &(JetE));
  tree->SetBranchAddress("JetPt", &(JetPt));
  tree->SetBranchAddress("JetEta", &(JetEta));
  tree->SetBranchAddress("JetPhi", &(JetPhi));
  tree->SetBranchAddress("JetCSV", &(JetCSV));
 tree->SetBranchAddress("eventWeight", &(eventWeight)); 
  TH1F *h_mX_CR1=new TH1F("h_mX_CR1", "h_mX_CR1", 200, 0., 2000.); h_mX_CR1->Sumw2();
  TH1F *h_mX_CR2=new TH1F("h_mX_CR2", "h_mX_CR2", 200, 0., 2000.); h_mX_CR2->Sumw2();
  TH1F *h_mX_CR3=new TH1F("h_mX_CR3", "h_mX_CR3", 200, 0., 2000.); h_mX_CR3->Sumw2();
  TH1F *h_mX_CR4=new TH1F("h_mX_CR4", "h_mX_CR4", 200, 0., 2000.); h_mX_CR4->Sumw2();
  TH1F *h_mX_CR5=new TH1F("h_mX_CR5", "h_mX_CR5", 200, 0., 2000.); h_mX_CR5->Sumw2();
  TH1F *h_mX_SR=new TH1F("h_mX_SR", "h_mX_SR",200, 0., 2000.);     h_mX_SR->Sumw2();
  TH1F * h_jet_pt_right[4];
  h_jet_pt_right[0]=new TH1F("h_jet_pt_right1","jet pt_{1} - gen jet pt_{1}; p_{T1} [GeV]", 150, -50.,50.);
  h_jet_pt_right[1]=new TH1F("h_jet_pt_right2","jet pt_{2} - gen jet pt_{2}; p_{T2} [GeV]", 150, -50.,50.);
  h_jet_pt_right[2]=new TH1F("h_jet_pt_right3","jet pt_{3} - gen jet pt_{3}; p_{T3} [GeV]", 150, -50.,50.);
  h_jet_pt_right[3]=new TH1F("h_jet_pt_right4","jet pt_{4} - gen jet pt_{4}; p_{T4} [GeV]", 150, -50.,50.);
  TH1F * h_jet_pt_wrong[4];
  h_jet_pt_wrong[0]=new TH1F("h_jet_pt_wrong1","jet pt_{1} - gen jet pt_{1}; p_{T1} [GeV]", 150, -50.,50.);
  h_jet_pt_wrong[1]=new TH1F("h_jet_pt_wrong2","jet pt_{2} - gen jet pt_{2}; p_{T2} [GeV]", 150, -50.,50.);
  h_jet_pt_wrong[2]=new TH1F("h_jet_pt_wrong3","jet pt_{3} - gen jet pt_{3}; p_{T3} [GeV]", 150, -50.,50.);
  h_jet_pt_wrong[3]=new TH1F("h_jet_pt_wrong4","jet pt_{4} - gen jet pt_{4}; p_{T4} [GeV]", 150, -50.,50.);
  TH1F *h_CSV2_CR1=new TH1F("h_CSV2_CR1", "h_CSV2_CR1", 25, 0., 1.); 
  TH1F *h_CSV2_CR2=new TH1F("h_CSV2_CR2", "h_CSV2_CR2", 25, 0., 1.); 
  TH1F *h_CSV2_CR3=new TH1F("h_CSV2_CR3", "h_CSV2_CR3", 25, 0., 1.); 
  TH1F *h_CSV2_CR4=new TH1F("h_CSV2_CR4", "h_CSV2_CR4", 25, 0., 1.); 
  TH1F *h_CSV2_CR5=new TH1F("h_CSV2_CR5", "h_CSV2_CR5", 25, 0., 1.); 
  TH1F *h_CSV2_SR=new TH1F("h_CSV2_SR", "h_CSV2_SR", 25, 0., 1.); 
  TH1F *h_mX_SR_right=new TH1F("h_mX_SR_right", "X mass (after KinFit) matched",200, 0., 2000.);  
  TH1F *h_mX_SR_wrong=new TH1F("h_mX_SR_wrong", "X mass (after KinFit) not matched",200, 0., 2000.);  
  TH1F *h_pt1= new TH1F("h_pt1", "pt jet 1 ; p_{T} [GeV]", 200, 0., 300.);
  TH1F *h_pt2= new TH1F("h_pt2", "pt jet 2 ; p_{T} [GeV]", 200, 0., 300.);
  TH1F *h_pt3= new TH1F("h_pt3", "pt jet 3 ; p_{T} [GeV]", 200, 0., 300.);
  TH1F *h_pt4= new TH1F("h_pt4", "pt jet 4 ; p_{T} [GeV]", 200, 0., 300.); 
  TH1F *h_pt1_corr= new TH1F("h_pt1_corr", "pt jet 1 (after KinFit); p_{T} [GeV]", 200, 0., 300.);
  TH1F *h_pt2_corr= new TH1F("h_pt2_corr", "pt jet 2 (after KinFit); p_{T} [GeV]", 200, 0., 300.);
  TH1F *h_pt3_corr= new TH1F("h_pt3_corr", "pt jet 3 (after KinFit); p_{T} [GeV]", 200, 0., 300.);
  TH1F *h_pt4_corr= new TH1F("h_pt4_corr", "pt jet 4 (after KinFit); p_{T} [GeV]", 200, 0., 300.);
  TH1F *h_mh1 = new TH1F("h_mh1", " mh1 ; m_{H1} [GeV]", 200, 0., 400.);
  TH1F *h_mh2 = new TH1F("h_mh2", " mh2 ; m_{H2} [GeV]", 200, 0., 400.);
  TH1F *h_mh1_corr = new TH1F("h_mh1_corr", " mh1 (after KinFit); m_{H1} [GeV]", 200, 0., 400.);
  TH1F *h_mh2_corr = new TH1F("h_mh2_corr", " mh2 (after KinFit); m_{H2} [GeV]", 200, 0., 400.); 
  TH1F *h_Xpt_right = new TH1F("h_Xpt_right", " Xpt matched; [GeV]", 1000, 0.,1000.);
  TH1F *h_Xpt_wrong = new TH1F("h_Xpt_wrong", " Xpt not matched; [GeV]", 1000, 0.,1000.);
  TH1F *h_chi_right = new TH1F("h_chi_right", " #chi_{2} matched; #chi_{2}", 1000, 0.,100.);
  TH1F *h_chi_wrong = new TH1F("h_chi_wrong", " #chi_{2} not matched; #chi_{2}", 1000, 0.,100.);
  TH1F *h_najets_right = new TH1F("h_najets_right", " njets matched; njets", 7, 3.,10.);
  TH1F *h_najets_wrong = new TH1F("h_najets_wrong", " njets not matched; njets", 7, 3.,10.); 
  TH1F *h_dR_right = new TH1F("h_dR_right", " min dR 5th jet matched; dR(5th)", 30, 0.,3.);
  TH1F *h_dR_wrong = new TH1F("h_dR_wrong", " min dR 5th jet not matched; dR(5th)", 30, 0.,3.);
 
  TH1F *h_mX_VB1=new TH1F("h_mX_VB1", "h_mX_VB1", 200, 0., 2000.); h_mX_VB1->Sumw2();
  TH1F *h_mX_VB2=new TH1F("h_mX_VB2", "h_mX_VB2", 200, 0., 2000.); h_mX_VB2->Sumw2();
  TH1F *h_mX_VB3=new TH1F("h_mX_VB3", "h_mX_VB3", 200, 0., 2000.); h_mX_VB3->Sumw2();
  TH1F *h_mX_VB4=new TH1F("h_mX_VB4", "h_mX_VB4", 200, 0., 2000.); h_mX_VB4->Sumw2();
  TH1F *h_mX_VB5=new TH1F("h_mX_VB5", "h_mX_VB5", 200, 0., 2000.); h_mX_VB5->Sumw2();
  TH1F *h_mX_VR=new TH1F("h_mX_VR", "h_mX_VR",200, 0., 2000.);     h_mX_VR->Sumw2(); 
  TH2F *h_deltaMX_MX = new TH2F("h_deltaMX_MX", "#Delta(MX) vs MX ", 100, 200., 600., 60, -60., 60. );  
  TH2F *h_deltaMX_MXfit = new TH2F("h_deltaMX_MXfit", "#Delta(MX) vs MX (after KinFit)", 100, 200., 600., 60, -60., 60. );
  TH2F *h_deltaMX_MX_najet = new TH2F("h_deltaMX_MX_najet", "#Delta(MX) vs MX if najet>4", 100, 200., 600., 60, -60., 60. );
  TH2F *h_deltaMX_MXfit_najet = new TH2F("h_deltaMX_MXfit_najet", "#Delta(MX) vs MX (after KinFit) if najet>4", 100, 200., 600., 60, -60., 60. );
  TH1F *h_X_mass_SR=new TH1F("h_X_mass_SR", "X mass; mass (GeV)", 200, 0., 2000.);
  TH1F *h_X_mass_SR_right=new TH1F("h_X_mass_SR_right", "X mass (in SR) matched; mass (GeV)", 200, 0., 2000.);
  TH1F *h_X_mass_SR_wrong=new TH1F("h_X_mass_SR_wrong", "X mass (in SR) not matched; mass (GeV)", 200, 0., 2000.);
  TH1F *h_X_mass_VR=new TH1F("h_X_mass_VR", "X mass; mass (GeV)", 200, 0., 2000.);
	
  TH1F *h_mX_SR_chi2=new TH1F("h_mX_SR_chi2", "h_mX_SR_chi2", 100, 0., 1200.); //h_mX_SR_chi2->Sumw2();
  
  std::string histfilename="Histograms_"+sample+".root";
  gSystem->Exec(("cp ../"+histfilename+" "+histfilename).c_str());
  TFile *tFile1=new TFile(("../"+histfilename).c_str(), "READ");
  TH1F h_Cuts=*((TH1F*)((TH1F*)tFile1->Get("h_Cuts"))->Clone("h_Cuts"));
  tFile1->Close();
  
  // Loop over events
  int nEvents=tree->GetEntries();
  std::cout<<"Number of events in input file = "<<nEvents<<std::endl;
  double nCut7=0;
double nSb=0;
double nSb_norm=0;
double nSc=0;
double nSc_norm=0;
double nSrest=0;
double nSrest_norm=0;

  double nCut8=0;	
double nCut9=0;
double nCut10=0;
double nCut11=0;  
for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);
    
    
    // Collate hJets and aJets into Jets
    nJets=nhJets+naJets;
    for (int j=0; j<nhJets; ++j)
    {
      jetE[j]=hJetE[j];
      jetpT[j]=smear_pt_resErr(hJetpT[j], hJet_genpT[j], hJeteta[j], sigmaJERUnc)+sigmaJECUnc*hJet_JECUnc[j];
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
      jet_genpT[j]=hJet_genpT[j];
    }
    for (int j=0; j<naJets; ++j)
    {
      jetE[j+nhJets]=aJetE[j];
      jetpT[j+nhJets]=smear_pt_resErr(aJetpT[j], aJet_genpT[j], aJeteta[j], sigmaJERUnc)+sigmaJECUnc*aJet_JECUnc[j];
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
      jet_genpT[j+nhJets]=aJet_genpT[j];
    }
    
    /*
    std::cout<<jetE[H1jet1_i]-regJet1E<<std::endl;
    std::cout<<jetE[H1jet2_i]-regJet2E<<std::endl;
    std::cout<<jetE[H2jet1_i]-regJet3E<<std::endl;
    std::cout<<jetE[H2jet2_i]-regJet4E<<std::endl;
    
    std::cout<<jetpT[H1jet1_i]-regJet1pT<<std::endl;
    std::cout<<jetpT[H1jet2_i]-regJet2pT<<std::endl;
    std::cout<<jetpT[H2jet1_i]-regJet3pT<<std::endl;
    std::cout<<jetpT[H2jet2_i]-regJet4pT<<std::endl;
    */
    double csv[4];
    for(int k =0; k<4; k++){ csv[k]=-1;}

    csv[0]=JetCSV[H1jet1_i];
    csv[1]=JetCSV[H1jet2_i];
    csv[2]=JetCSV[H2jet1_i];
    csv[3]=JetCSV[H2jet2_i];		
    float sCSV=-1.;
    float fCSV=-1.;	
   
    for(int k=0; k<4; k++){
	    if(csv[k]>fCSV) fCSV=csv[k];
	    else if(csv[k]>sCSV) sCSV=csv[k];
		}
    JetE[H1jet1_i]=regJet1E;
    JetE[H1jet2_i]=regJet2E;
    JetE[H2jet1_i]=regJet3E;
    JetE[H2jet2_i]=regJet4E;
    JetPt[H1jet1_i]=regJet1pT;
    JetPt[H1jet2_i]=regJet2pT;
    JetPt[H2jet1_i]=regJet3pT;
    JetPt[H2jet2_i]=regJet4pT;
          
    TLorentzVector jet1_p4=fillTLorentzVector(JetPt[H1jet1_i], JetEta[H1jet1_i], JetPhi[H1jet1_i], JetE[H1jet1_i]);      
    TLorentzVector jet2_p4=fillTLorentzVector(JetPt[H1jet2_i], JetEta[H1jet2_i], JetPhi[H1jet2_i], JetE[H1jet2_i]);      
    TLorentzVector jet3_p4=fillTLorentzVector(JetPt[H2jet1_i], JetEta[H2jet1_i], JetPhi[H2jet1_i], JetE[H2jet1_i]);      
    TLorentzVector jet4_p4=fillTLorentzVector(JetPt[H2jet2_i], JetEta[H2jet2_i], JetPhi[H2jet2_i], JetE[H2jet2_i]);      
   
     if (int((jet1_p4+jet2_p4).Pt()*100.) % 2 == 1) {swap(H1jet1_i, H2jet1_i); swap(H1jet2_i, H2jet2_i);} 

    jet1_p4=fillTLorentzVector(JetPt[H1jet1_i], JetEta[H1jet1_i], JetPhi[H1jet1_i], JetE[H1jet1_i]);
    jet2_p4=fillTLorentzVector(JetPt[H1jet2_i], JetEta[H1jet2_i], JetPhi[H1jet2_i], JetE[H1jet2_i]);
    jet3_p4=fillTLorentzVector(JetPt[H2jet1_i], JetEta[H2jet1_i], JetPhi[H2jet1_i], JetE[H2jet1_i]);
    jet4_p4=fillTLorentzVector(JetPt[H2jet2_i], JetEta[H2jet2_i], JetPhi[H2jet2_i], JetE[H2jet2_i]);
    
    TLorentzVector H1_p4=jet1_p4+jet2_p4;              
    TLorentzVector H2_p4=jet3_p4+jet4_p4;
    TLorentzVector X_p4=H1_p4+H2_p4;
          
    double H1_mass=H1_p4.M();
    double H2_mass=H2_p4.M();
    double X_mass=X_p4.M();
    double X_mass_pre = X_mass;
    int region=withinRegion(H1_mass, H2_mass, 20., 35., H_mass, H_mass);
    int region2=withinRegion(H1_mass, H2_mass, 15., 35., 90, 90);	
    if (region==0 && pass125125) 
    {
      h_X_mass_SR->Fill(X_mass, eventWeight);	
			if(matched==1 )
      {
			  h_X_mass_SR_right->Fill(X_mass, eventWeight); 
        if(X_mass<340 && X_mass >310) 
        {  
			    h_mh1->Fill(H1_mass);
			    h_mh2->Fill(H2_mass);
			    h_pt1->Fill(jet1_p4.Pt());
          h_pt2->Fill(jet2_p4.Pt());
          h_pt3->Fill(jet3_p4.Pt());
          h_pt4->Fill(jet4_p4.Pt());
			  }
			  else
        {
			    h_mh1_corr->Fill(H1_mass);
          h_mh2_corr->Fill(H2_mass);
			  }	
			}
			else h_X_mass_SR_wrong->Fill(X_mass, eventWeight);
		}	
    double chi2=0; 		
    double X_pt=0;
		// Put in kinematic constraining here
    TLorentzVector jet1_new_p4, jet2_new_p4;
		if (kinConstraint==1)
		{
			TLorentzVector X_chi2_p4=H4b::Xchi2(jet1_p4, jet2_p4, jet3_p4, jet4_p4, H_mass, jet1_new_p4, jet2_new_p4);
			X_mass=X_chi2_p4.M();
			X_pt = X_chi2_p4.Pt();
			H1_p4=jet1_p4+jet2_p4;
      H2_p4=jet3_p4+jet4_p4;
			H1_mass=H1_p4.M();
   		H2_mass=H2_p4.M();
		}
    if(region2==0)//VR
	  {
	    h_X_mass_VR->Fill(X_mass);
	  }
    
    if(pass125125)
    {

    if (region==0) // SR                                                  
    { ++nCut7;                                                                    
    if(QuadJetFilterFlag)  ++nCut8;  
    if(triggerFlags[57])  ++nCut9;
    if(triggerFlags[54])  ++nCut10;
    if(triggerFlags[57]||triggerFlags[54])  ++nCut11;	
    h_mX_SR->Fill(X_mass,eventWeight);  	                                                            
    if(matched==1)                                             
    {                                                          
	    if(nJetCSV==4)                                           
      {		                                                     
		    h_deltaMX_MXfit->Fill(X_mass, (X_mass_pre- X_mass));   
		    h_deltaMX_MX->Fill(X_mass_pre,(X_mass_pre- X_mass));   
		}
	  else{
		h_deltaMX_MXfit_najet->Fill(X_mass, (X_mass_pre- X_mass));
                h_deltaMX_MX_najet->Fill(X_mass_pre,(X_mass_pre- X_mass));
		}	
		h_mX_SR_right->Fill(X_mass, eventWeight);
                // h_najets_right->Fill(nJetCSV);
		//h_chi_right->Fill(chi2);
   if(X_mass<340 && X_mass >310) {
		for(int i=0;i<4;i++) h_jet_pt_right[i]->Fill(corr[i]);
		//h_mh1_corr->Fill(H1_mass);
                //h_mh2_corr->Fill(H2_mass);
                h_pt1_corr->Fill(jet1_p4.Pt());
                h_pt2_corr->Fill(jet2_p4.Pt());
                h_pt3_corr->Fill(jet3_p4.Pt());
                h_pt4_corr->Fill(jet4_p4.Pt());
		//	}
//		//if(X_mass<340 && X_mass >310) {
		//h_chi_right->Fill(chi2);
		if(dR_min_5thJet<999) h_dR_right->Fill(dR_min_5thJet);
		h_najets_right->Fill(nJetCSV);	
		}
		else { 
		for(int i=0;i<4;i++) h_jet_pt_wrong[i]->Fill(corr[i]);
	//	h_chi_wrong->Fill(chi2);
                if(dR_min_5thJet<999) h_dR_wrong->Fill(dR_min_5thJet);
		h_najets_wrong->Fill(nJetCSV);	
			}
		}
      else {
	//	if(X_mass<340 && X_mass >310) i
	//	h_najets_wrong->Fill(nJetCSV);  
		h_mX_SR_wrong->Fill(X_mass, eventWeight);	 
		//h_chi_wrong->Fill(chi2);

		}
     if(X_mass<340 && X_mass >310) {
                if(matched==1) h_Xpt_right->Fill(Ht_less4+X_pt);
		else h_Xpt_wrong->Fill(Ht_less4+X_pt);
	}
      h_mX_SR_chi2->Fill(X_mass, eventWeight);
      h_CSV2_SR->Fill(sCSV);	
    }                                                                     
    else if (region==1) // CR1                                            
    {                                                                     
      h_mX_CR1->Fill(X_mass, eventWeight);                          
      h_CSV2_CR1->Fill(sCSV); 	
    }                                                                     
    else if (region==2) // CR2                                            
    {
      nSb_norm++;	                                                                     
      if(matched==1)	nSb++;
	
      h_mX_CR2->Fill(X_mass, eventWeight);                          
      h_CSV2_CR2->Fill(sCSV); 	
    }                                                                     
    else if (region==3) // CR3                                            
    {                                                                     
      nSc_norm++;
      if(matched==1)    nSc++;	
      h_mX_CR3->Fill(X_mass, eventWeight);                          
      h_CSV2_CR3->Fill(sCSV);	
    }                                                                     
    else if (region==4) // CR4                                            
    {            
      nSb_norm++;
      if(matched==1)    nSb++;	                                                         
      h_mX_CR4->Fill(X_mass, eventWeight);                           
      h_CSV2_CR4->Fill(sCSV);	 
    }                                                                     
    else if (region==5) // CR5                                            
    {                         
      nSc_norm++;
      if(matched==1)    nSc++;	                                            
      h_mX_CR5->Fill(X_mass, eventWeight);                           
      h_CSV2_CR5->Fill(sCSV);	
    }                                                                     
    else if (region==-1)                                                  
    { 
      nSrest_norm++;	                                                                    
      if(matched==1)    nSrest++;
      std::cout<<"Didn't fall in any region!"<<std::endl;                 
    }
   }

if(pass9090){
    if (region2==0) // VR                                                  
    { 
      h_mX_VR->Fill(X_mass);
    }
    else if (region2==1) // VB1                                            
    {
      h_mX_VB1->Fill(X_mass, eventWeight);
    }
    else if (region2==2) // VB2                                            
    {
      h_mX_VB2->Fill(X_mass, eventWeight);
    }
    else if (region2==3) // VB3                                            
    {
      h_mX_VB3->Fill(X_mass, eventWeight);
    }
    else if (region2==4) // VB4                                            
    {
      h_mX_VB4->Fill(X_mass, eventWeight);
    }
    else if (region2==5) // VB5                                            
    {
      h_mX_VB5->Fill(X_mass, eventWeight);
    }
    else if (region2==-1)
    {
      std::cout<<"Didn't fall in any region!"<<std::endl;
    }
    }



  
  } // Event loop
  
  std::cout<<"out of event loop"<<std::endl;
  std::cout<<"nCut7 Signal = "<<nCut7<<std::endl;
  std::cout<<"nCut8  Quad45= "<<nCut8<<std::endl;	
  std::cout<<"nCut9 Quad50= "<<nCut9<<std::endl;
  std::cout<<"nCut10 CSV= "<<nCut10<<std::endl;
  std::cout<<"nCut11 OR 45 CSV = "<<nCut11<<std::endl;
  h_Cuts.Fill(15, nCut7);
  
  TFile *tFile2=new TFile(histfilename.c_str(), "UPDATE");
  tFile2->Delete("h_Cuts;1");
  h_mX_CR1->Write();
  h_mX_CR2->Write();
  h_mX_CR3->Write();
  h_mX_CR4->Write();
  h_mX_CR5->Write();
  h_mX_SR->Write();
  h_mX_SR_right->Write();
  h_mX_SR_wrong->Write();
  h_chi_right->Write();
  h_chi_wrong->Write();
  h_Xpt_right->Write();
  h_Xpt_wrong->Write();

  h_pt1->Write();
  h_pt2->Write();
  h_pt3->Write();
  h_pt4->Write();
  h_pt1_corr->Write();
  h_pt2_corr->Write();
  h_pt3_corr->Write();
  h_pt4_corr->Write();
  h_mh1->Write();
  h_mh2->Write();
  h_najets_right->Write();
  h_najets_wrong->Write();
  h_dR_right->Write();
  h_dR_wrong->Write();	
  h_mh1_corr->Write();
  h_mh2_corr->Write();
  h_mX_VB1->Write();
  h_mX_VB2->Write();
  h_mX_VB3->Write();
  h_mX_VB4->Write();
  h_mX_VB5->Write();
  h_mX_VR->Write();
  h_CSV2_CR1->Write();
  h_CSV2_CR2->Write();
  h_CSV2_CR3->Write();
  h_CSV2_CR4->Write();
  h_CSV2_CR5->Write();
  h_CSV2_SR->Write();
  h_X_mass_VR->Write();
  h_X_mass_SR->Write();	
  h_X_mass_SR_right->Write();
  h_X_mass_SR_wrong->Write();
  h_mX_SR_chi2->Write();
  h_jet_pt_right[0]->Write();
  h_jet_pt_right[1]->Write();
  h_jet_pt_right[2]->Write();
  h_jet_pt_right[3]->Write();
  h_jet_pt_wrong[0]->Write();
  h_jet_pt_wrong[1]->Write();
  h_jet_pt_wrong[2]->Write();
  h_jet_pt_wrong[3]->Write(); 
  h_deltaMX_MXfit->Write();
  h_deltaMX_MX_najet->Write();
  h_deltaMX_MXfit_najet->Write();
  h_deltaMX_MX->Write();
  h_Cuts.Write();
  tFile2->Write();
  tFile2->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
 std::cout<< " sb fraction "<<nSb/nSb_norm<<"  matched "<< nSb << " norm "<< nSb_norm<<std::endl;
std::cout<< " sc fraction "<<nSc/nSc_norm<<"  matched "<< nSc << " norm "<< nSc_norm<<std::endl;
std::cout<< " srest fraction "<<nSrest/nSrest_norm<<"  matched "<< nSrest << " norm "<< nSrest_norm<<std::endl;
  std::cout<<"=== Cut Efficiencies === "<<std::endl;
  std::cout<<"Number of events in SR = "<<nCut7<<std::endl;
  std::cout<<"========================"<<std::endl;
  
  return 0;
}
            
      
  
  
  
  
