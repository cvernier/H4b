#include <TH1F.h>
#include <string>
#include <sstream>
#include <cmath>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/HbbHbb/Analysis/bJetRegression/HelperFunctions.h"


double sigmaJECUnc=-1; // (-1, 0, 1)
double sigmaJERUnc=0; // (-1, 0, 1)

double pi=3.14159265358979;

double jetpT_cut=30.;
double jeteta_cut=2.5;
double H_mass=160.0;
double mH_diff_cut=50.;
double mH_mean_cut=15.;

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


 int withinRegion(double mH1, double mH2, double r1=15., double r2=30., double mH1_c=H_mass, double mH2_c=H_mass)
{
  double r=pow(pow(mH1-mH1_c, 2)+pow(mH2-mH2_c, 2), 0.5);
  double angle=atan2(mH2-mH2_c, mH1-mH1_c);
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
  return ret;
}


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

int Hbb_one(std::string dir, std::string sample, std::string PUWeight="", std::string csvReshape="")
{
  
  std::string inputfilename=dir+sample+".root";
  TChain *tree=new TChain("tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
  
  TFile *file_PUWeight;
  TH1F *h_PUWeight;
  if (PUWeight!="")
  {
    file_PUWeight=new TFile(PUWeight.c_str());
    std::cout<<"Opened PU weight file = "<<PUWeight<<std::endl;
    h_PUWeight=(TH1F*)gDirectory->Get("h_PUWeight");
  }
  
  // Book variables
  int vType;
  int pass125125=0;
  int pass9090=0;	
  bool triggerFlags[500],QuadJetFilterFlag;
  int nPV;
  int nhJets, naJets;
  float dR_min_5thJet; 	
  float corr[4];	
  int nJets, nCJets, n8Jets, nBJets, nBMJets ; 
  float hJetE[100], hJetpT[100], hJeteta[100], hJetphi[100], hJetCSV[100], hJetflavor[100], hJetpTRaw[100], hJet_ptLeadTrack[100], hJet_vtx3dL[100], hJet_vtx3deL[100], hJet_vtxMass[100], hJet_vtxPt[100], hJet_cef[100], hJet_nconstituents[100], hJet_JECUnc[100], hJet_genpT[100];
  float aJetE[100], aJetpT[100], aJeteta[100], aJetphi[100], aJetCSV[100], aJetflavor[100], aJetpTRaw[100], aJet_ptLeadTrack[100], aJet_vtx3dL[100], aJet_vtx3deL[100], aJet_vtxMass[100], aJet_vtxPt[100], aJet_cef[100], aJet_nconstituents[100], aJet_JECUnc[100], aJet_genpT[100];
  float jetE[100], jetpT[100], jeteta[100], jetphi[100], jetCSV[100], jetflavor[100], jetpTRaw[100], jet_ptLeadTrack[100], jet_vtx3dL[100], jet_vtx3deL[100], jet_vtxMass[100], jet_vtxPt[100], jet_cef[100], jet_nconstituents[100], jet_JECUnc[100], jet_genpT[100];
  METInfo metObj;
 float ajetE[100], ajetpT[100], ajeteta[100], ajetphi[100], ajetCSV[100], in[100];
 int ind[4];
 for( int i =0; i<100 ; i++) {
		
		ajetE[i]=-999;
		ajetpT[i]=-999;
		ajeteta[i]=-999;
		ajetphi[i]=-999;
		ajetCSV[i]=-999;	
  }

for(int i=0;i<4;i++) {ind[i]=-999; corr[i]=-999;}
  float met, metSig, ht;
  double weightPU=1.;
 
  GenParticleInfo genX, genH1, genH2,genH1B, genH1Bbar, genH2B, genH2Bbar;

  //std::string a_csv = "aJet_csv_"+csvReshape;	
  //std::string h_csv = "hJet_csv_"+csvReshape;
 
  // Retrieve variables
  tree->SetBranchAddress("QuadJetFilterFlag", &(QuadJetFilterFlag));
  tree->SetBranchAddress("Vtype", &(vType));
  tree->SetBranchAddress("triggerFlags", &(triggerFlags));
  tree->SetBranchAddress("nPVs", &(nPV));
  tree->SetBranchAddress("nhJets", &(nhJets));
  tree->SetBranchAddress("hJet_e", &(hJetE)); 
  tree->SetBranchAddress("hJet_pt", &(hJetpT));
  tree->SetBranchAddress("hJet_eta", &(hJeteta));
  tree->SetBranchAddress("hJet_phi", &(hJetphi));
  tree->SetBranchAddress("hJet_csv_nominal", &(hJetCSV));
  tree->SetBranchAddress("hJet_flavour", &(hJetflavor));
  tree->SetBranchAddress("naJets", &(naJets));
  tree->SetBranchAddress("aJet_e", &(aJetE));
  tree->SetBranchAddress("aJet_pt", &(aJetpT));
  tree->SetBranchAddress("aJet_eta", &(aJeteta));
  tree->SetBranchAddress("aJet_phi", &(aJetphi));
  tree->SetBranchAddress("aJet_csv_nominal", &(aJetCSV));
  tree->SetBranchAddress("MET", &(metObj));
  tree->SetBranchAddress("aJet_JECUnc", &(aJet_JECUnc));
  tree->SetBranchAddress("hJet_JECUnc", &(hJet_JECUnc));
  tree->SetBranchAddress("genX", &(genX));
  tree->SetBranchAddress("genH1", &(genH1));
  tree->SetBranchAddress("genH2", &(genH2));
  tree->SetBranchAddress("genH1B", &(genH1B));
  tree->SetBranchAddress("genH1Bbar", &(genH1Bbar));
  tree->SetBranchAddress("genH2B", &(genH2B));
  tree->SetBranchAddress("genH2Bbar", &(genH2Bbar));    	
  tree->SetBranchAddress("PUweight", &(weightPU));
  tree->SetBranchAddress("aJet_genPt", &(aJet_genpT));
  
  TH1F *h_nJets=new TH1F("h_nJets", "# Central CSVL Jets with p_{T}> 30 GeV; n", 10, 0., 10.);
  TH1F *h_n8Jets=new TH1F("h_n8Jets", "# Central Jets with p_{T}> 80 GeV;  n", 10, 0., 10.); 
  TH1F *h_nCand=new TH1F("h_nCand", "# HH candidates", 20, 0., 20.);
  TH1F *h_nCand_true=new TH1F("h_nCand_true", "# HH candidates if true", 20, 0., 20.);	
  TH1F *h_nPV=new TH1F("h_nPV", "# of Primary Vertices; nPV", 51, 0., 50.); h_nPV->Sumw2();
  TH1F *h_nPV_weighted=new TH1F("h_nPV_weighted", "# of Primary Vertices after Reweighting; nPV", 51, 0., 50.); h_nPV_weighted->Sumw2();
  
  TH1F *h_JetpT1=new TH1F("h_JetpT1", "JetpT1", 50, 0., 800.);
  TH1F *h_JetpT2=new TH1F("h_JetpT2", "JetpT2", 50, 0., 500.);
  TH1F *h_JetpT3=new TH1F("h_JetpT3", "JetpT3", 50, 0., 350.);
  TH1F *h_JetpT4=new TH1F("h_JetpT4", "JetpT4", 50, 0., 250.);
  
  
 /* TH1F *h_JetpT1=new TH1F("h_JetpT1", "JetCSV1", 50, -1., 1.);
  TH1F *h_JetpT2=new TH1F("h_JetpT2", "JetCSV2", 50, -1., 1.);
  TH1F *h_JetpT3=new TH1F("h_JetpT3", "JetCSV3", 50, -1.,1.);
  TH1F *h_JetpT4=new TH1F("h_JetpT4", "JetCSV4", 50, -1., 1.);
*/
  TH1F *h_MET=new TH1F("h_MET", "MET; MET (GeV)", 100, 0, 200.);
  TH1F *h_MET_sig=new TH1F("h_MET_sig", "MET Significance; Sig", 20, 0., 20.);
  TH1F *h_HT=new TH1F("h_HT", "HT Distribution", 100, 0., 3000.);
  
  TH1F *h_H1_mass=new TH1F("h_H1_mass", "H mass; mass (GeV)", 100, 50., 300.);
  TH1F *h_H1_pT=new TH1F("h_H1_pT", "H1 p_{T}; p_{T} (GeV/c)", 50, 0., 800.);
  TH1F *h_H2_mass=new TH1F("h_H2_mass", "H2 mass; mass (GeV)", 50, 50., 200.);
  TH1F *h_H2_pT=new TH1F("h_H2_pT", "H2 p_{T}; p_{T} (GeV/c)", 50, 0., 800.);
  TH2F *h_H1_H2_mass = new TH2F("h_H1_H2_mass", " all comb if min(#sigma (#delta m)^{2}) ", 100, 30., 200., 100, 30., 200.);
  TH1F *h_HH_mass_diff_cand=new TH1F("h_HH_mass_diff_cand", "|#Deltam| between Higgs masses - all candidates", 50, 0., 200.);
  TH1F *h_HH_massNorm_diff=new TH1F("h_HH_massNorm_diff", "|#Deltam| between Higgs masses", 50, 0., 2.);
  TH1F *h_HH_CSV = new TH1F("h_HH_CSV"," Sum CSV | between the two higgs ", 70, -4., 4.); 
  TH1F *h_genX_mass=new TH1F("h_genX_mass", "Generator X mass; m_{X} GeV", 100, 0., 1000.);
  TH1F *h_genH1_mass=new TH1F("h_genH1_mass", "Generator H1 mass; m_{X} GeV", 100, 0., 1000.);
  TH1F *h_genH2_mass=new TH1F("h_genH2_mass", "Generator H2 mass; m_{X} GeV", 100, 0., 1000.);
  TH1F *h_Xmass = new TH1F("h_Xmass"," h_Xmass" , 100, 0., 1000.); 
  TH1F *h_Xmass_right = new TH1F("h_Xmass_right"," h_Xmass right" , 100, 0., 1000.); 
  TH1F *h_Xmass_wrong = new TH1F("h_Xmass_wrong"," h_Xmass wrong" , 100, 0., 1000.);
  TH1F *h_Xmass_first = new TH1F("h_Xmass_first"," h_Xmass first" , 100, 0., 1000.);
  TH1F *h_Xmass_other = new TH1F("h_Xmass_other"," h_Xmass other" , 100, 0., 1000.);		
  TH2F *h_Xmass2 = new TH2F("h_Xmass2"," h_Xmass vs njets " , 10, 4, 14., 100, 0., 1000.); 	
  TH1F *h_Xmass_comb = new TH1F("h_Xmass_comb"," h_Xmass SR" , 200, 0., 2000.);
  TH1F *h_Cuts=new TH1F("h_Cuts", "Cut flow", 16, 0, 16);
  TH1F * h_jet_pt[4];
  h_jet_pt[0]=new TH1F("h_jet_pt1","jet pt_{1} - gen jet pt_{1}; p_{T1} [GeV]", 150, -50.,50.);
  h_jet_pt[1]=new TH1F("h_jet_pt2","jet pt_{2} - gen jet pt_{2}; p_{T2} [GeV]", 150, -50.,50.);
  h_jet_pt[2]=new TH1F("h_jet_pt3","jet pt_{3} - gen jet pt_{3}; p_{T3} [GeV]", 150, -50.,50.);
  h_jet_pt[3]=new TH1F("h_jet_pt4","jet pt_{4} - gen jet pt_{4}; p_{T4} [GeV]", 150, -50.,50.);			
  TAxis *a_Cuts=h_Cuts->GetXaxis();
  a_Cuts->SetBinLabel(2, "Step2");
  a_Cuts->SetBinLabel(4, "Trigger");
  a_Cuts->SetBinLabel(6, "Vtype");
  a_Cuts->SetBinLabel(8, "nCJets>3");
  a_Cuts->SetBinLabel(10, "HH Cand");
  a_Cuts->SetBinLabel(12, "massDiff<50");
  a_Cuts->SetBinLabel(14, "b-tagging");
  a_Cuts->SetBinLabel(16, "SR");
  
  
  std::string outfilename=sample+"_selected"+".root";

 // std::string outfilename="OfficialStep2_KinematicallySelected_"+sample+".root";
  TFile *outfile=new TFile(outfilename.c_str(), "recreate");

  TTree *outtree=tree->CloneTree(0);
  int H1jet1_i, H1jet2_i;
  int nJetCSV;	int matched=-999;
  int H2jet1_i, H2jet2_i;
  outtree->Branch("H1jet1_i", &H1jet1_i, "H1jet1_i/I");
  outtree->Branch("H1jet2_i", &H1jet2_i, "H1jet2_i/I");
  outtree->Branch("H2jet1_i", &H2jet1_i, "H2jet1_i/I");
  outtree->Branch("H2jet2_i", &H2jet2_i, "H2jet2_i/I");
  outtree->Branch("JetCSV", &ajetCSV, "ajetCSV[100]/F");
  outtree->Branch("JetPt", &ajetpT, "ajetpT[100]/F");
  outtree->Branch("JetEta", &ajeteta, "ajeteta[100]/F");
  outtree->Branch("JetPhi", &ajetphi, "ajetphi[100]/F");	
  outtree->Branch("JetE", &ajetE, "ajetE[100]/F");	
  outtree->Branch("nJetCSV", &nJetCSV, "nJetCSV/I");   	
  outtree->Branch("PtInd", &ind, "ind[4]/I");	
  outtree->Branch("pass125125", &pass125125, "pass125125/I");
  outtree->Branch("pass9090", &pass9090, "pass9090/I");
  outtree->Branch("matched", &matched, "matched/I");
  outtree->Branch("dR_min_5thJet", &dR_min_5thJet, "dR_min_5thJet/F");
  outtree->Branch("corr", &corr, "corr[4]/F");   
  // Loop over events
  int nEvents=tree->GetEntries();
  double nCut0=0, nCut1=0, nCut2=0, nCut3=0, nCut4=0, nCut5=0, nCut6=0, nCut7=0, nCutGen=0, nCut50=0, nCutT=0, nCutOR=0, ncutFail=0, ncut4G=0, ncut5G=0, nCountS=0;
  double contComb=0;	
  int nCandidate =0;
  for (int i=0; i<nEvents; ++i)
    {
    
    ++nCut0;
    tree->GetEvent(i);
    
    h_nPV->Fill(nPV);
    if (PUWeight!="") weightPU=h_PUWeight->GetBinContent(h_PUWeight->FindBin(nPV));
    h_nPV_weighted->Fill(nPV, weightPU);
    TLorentzVector h1b, h1bar, h2b, h2bar;
    h1b.SetPtEtaPhiM(genH1B.pt, genH1B.eta, genH1B.phi, genH1B.mass);
    h2b.SetPtEtaPhiM(genH2B.pt, genH2B.eta, genH2B.phi, genH2B.mass);
    h1bar.SetPtEtaPhiM(genH1Bbar.pt, genH1Bbar.eta, genH1Bbar.phi, genH1Bbar.mass);
    h2bar.SetPtEtaPhiM(genH2Bbar.pt, genH2Bbar.eta, genH2Bbar.phi, genH2Bbar.mass);
	
   // Fill generator level info
    TLorentzVector h1, h2;
    h1.SetPtEtaPhiM(genH1.pt, genH1.eta, genH1.phi, genH1.mass);
    h2.SetPtEtaPhiM(genH2.pt, genH2.eta, genH2.phi, genH2.mass);
  //  h_genX_mass->Fill((h1+h2).M());
    // h_genX_mass->Fill(genX.mass);
   // h_genH1_mass->Fill(genH1.mass);
   // h_genH2_mass->Fill(genH2.mass);
    
    // Collate hJets and aJets into Jets
    nJets=nhJets+naJets;
    for (int j=0; j<nhJets; ++j)
    {
      jetpT[j]=smear_pt_resErr(hJetpT[j], hJet_genpT[j], hJeteta[j], sigmaJERUnc)+sigmaJECUnc*hJet_JECUnc[j];
      jetE[j]=hJetE[j]*jetpT[j]/hJetpT[j]; 
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
      jetpT[j+nhJets]=smear_pt_resErr(aJetpT[j], aJet_genpT[j], aJeteta[j], sigmaJERUnc)+sigmaJECUnc*aJet_JECUnc[j];
      jetE[j+nhJets]=aJetE[j]*jetpT[j+nhJets]/aJetpT[j];
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
    TLorentzVector ajet1_p4, ajet2_p4, ajet3_p4, ajet4_p4;
    JetList jetList_CSV, jetList_pT; 
    nCJets=0;
    n8Jets=0;	
    nBJets=0;	
    nBMJets=0;	
    ht=0;
    for (int j=0; j<nJets; ++j)
    {
      if (fabs(jeteta[j])<2.5) 
      {
        jetList_pT[jetpT[j]]=j;
        if (jetpT[j]>jetpT_cut&&jetCSV[j]>0.244)
	      {
          jetList_CSV[jetCSV[j]]=j;
	        ++nCJets;
          ht+=jetpT[j];
        }
			
      }
    }
    
    met=metObj.et;
    metSig=metObj.sig;
    h_MET_sig->Fill(metSig);
    h_HT->Fill(ht);
    
    // Analysis begins here
    if (triggerFlags[54])//||triggerFlags[57])
    {
      if(QuadJetFilterFlag) ++nCutT;
      if(triggerFlags[57]) ++nCut50;
      if(triggerFlags[54]) ++nCut1;
      if(QuadJetFilterFlag||triggerFlags[54]) ++nCutOR;		
      if (vType==4 || vType==8 || vType==9 || vType==10)
      {
        ++nCut2;
        if(nCJets>3){
	int c =0;  
          for (JetList::reverse_iterator iJet=jetList_CSV.rbegin(); iJet!=jetList_CSV.rend(); ++iJet)
          {
	  ajetpT[c]=jetpT[iJet->second];//iJet->first;
          ajeteta[c] = jeteta[iJet->second];
          ajetphi[c] = jetphi[iJet->second];
          ajetE[c] = jetE[iJet->second];	
	  ajetCSV[c]=iJet->first;//jetCSV[iJet->second];	
	  in[c]= iJet->second;
	  if(jetpT[iJet->second]>80.) ++n8Jets;
	  if(iJet->first > 0.244) ++nBJets;
	  if(iJet->first > 0.679) ++nBMJets;
           
	  ++c;

           } //end for jet iterator
	  nJetCSV =c;	
        }
        if (nCJets>3 ) 
        {
          ++nCut3;
	  h_nJets->Fill(nBJets, weightPU); 
          h_JetpT1->Fill(ajetCSV[0]);
          h_JetpT2->Fill(ajetCSV[1]);
          h_JetpT3->Fill(ajetCSV[2]);
          h_JetpT4->Fill(ajetCSV[3]); 
	  h_n8Jets->Fill(n8Jets, weightPU);
	  bool found4jet = false;
          bool foundHH=false; 
	  bool match=false;
	  int nCC=0;int cont=0;
	  TLorentzVector jet1_p4, jet2_p4, jet3_p4, jet4_p4, jet5;
          TLorentzVector H1_p4, H2_p4;
	  double min_diff=50;
	  for (int j=0; j<nCJets; ++j)
	  {
			  jet1_p4=fillTLorentzVector(ajetpT[j], ajeteta[j], ajetphi[j], ajetE[j]);  
			  for (int k=0; k<nCJets; ++k)
			  {
				  if (k!=j)
				  {
					  jet2_p4=fillTLorentzVector(ajetpT[k], ajeteta[k], ajetphi[k], ajetE[k]);
					  cont++;        
					  {

						  for (int l=0; l<nCJets; ++l)	
						  {     
							  if(l!=j && l!=k )
							  {
								  jet3_p4=fillTLorentzVector(ajetpT[l], ajeteta[l], ajetphi[l], ajetE[l]);
								  for (int m=0; m<nCJets; ++m)
								  {
									  if (m!=l && m!=j && m!=k )
									  {
										  jet4_p4=fillTLorentzVector(ajetpT[m], ajeteta[m], ajetphi[m], ajetE[m]);
										  TLorentzVector diJet2_p4, diJet1_p4;	

										  diJet2_p4=jet3_p4+jet4_p4;
										  diJet1_p4=jet1_p4+jet2_p4;
										double diff=fabs(diJet1_p4.M()-diJet2_p4.M());     
										  if((diff<min_diff)&&(((diJet2_p4.M()<160. && diJet2_p4.M()> 90.) && (diJet1_p4.M()<160. && diJet1_p4.M()> 90.))||((diJet2_p4.M()<125. && diJet2_p4.M()> 55.) && (diJet1_p4.M()<125. && diJet1_p4.M()>  55.))))
										  {


	
											  if(((jet1_p4.Pt()>80. && jet2_p4.Pt()>80.) || (jet1_p4.Pt()>80. && jet3_p4.Pt()>80.) || (jet1_p4.Pt()>80. && jet4_p4.Pt()>80.) ||(jet2_p4.Pt()>80.&&jet3_p4.Pt()>80.) || (jet2_p4.Pt()>80. && jet4_p4.Pt()>80.) || (jet3_p4.Pt()>80. && jet4_p4.Pt()>80.)) && ((ajetCSV[j]>0.898 && ajetCSV[k]>0.898)||(ajetCSV[k]>0.898 && ajetCSV[l]>0.898)||(ajetCSV[l]>0.898 && ajetCSV[j]>0.898)||(ajetCSV[m]>0.898 && ajetCSV[l]>0.898)||(ajetCSV[m]>0.898 && ajetCSV[k]>0.898)||(ajetCSV[j]>0.898 && ajetCSV[m]>0.898))) {	

											  if((diJet2_p4.M()<160 && diJet2_p4.M()> 90.) && (diJet1_p4.M()<160 && diJet1_p4.M()> 90.)) pass125125 =1;
                                                                                          if((diJet2_p4.M()<125 && diJet2_p4.M()> 55.) && (diJet1_p4.M()<125 && diJet1_p4.M()>  55.)) pass9090=1;

											  foundHH=true;	
         										   H1jet1_i=j;
                                                                                          H1jet2_i=k;
                                                                                          H2jet1_i=l;
                                                                                          H2jet2_i=m;
                                                                                          ind[2]=in[l];
                                                                                          ind[3]=in[m];


											  //std::cout<<"Event "<<i<<"  0 "<<ajetCSV[H1jet1_i]<<"  1  "<<ajetCSV[H1jet2_i]<< "   2   "<<ajetCSV[H2jet1_i]<<"  3  "<<ajetCSV[H2jet2_i]<<std::endl;	 
										          min_diff=diff;
												  nCC++;
												  JetList jet4pt;
												  ajet1_p4=fillTLorentzVector(ajetpT[H1jet1_i], ajeteta[H1jet1_i], ajetphi[H1jet1_i], ajetE[H1jet1_i]);
												  ajet2_p4=fillTLorentzVector(ajetpT[H1jet2_i], ajeteta[H1jet2_i], ajetphi[H1jet2_i], ajetE[H1jet2_i]);
												  ajet3_p4=fillTLorentzVector(ajetpT[H2jet1_i], ajeteta[H2jet1_i], ajetphi[H2jet1_i], ajetE[H2jet1_i]);
												  ajet4_p4=fillTLorentzVector(ajetpT[H2jet2_i], ajeteta[H2jet2_i], ajetphi[H2jet2_i], ajetE[H2jet2_i]);
											          jet4pt[ajetpT[H1jet1_i]]=0;
												  jet4pt[ajetpT[H1jet2_i]]=1;
												  jet4pt[ajetpT[H2jet1_i]]=2;
												  jet4pt[ajetpT[H2jet2_i]]=3;	 			
												  contComb++;
													  //*************** matching to gen jets   ***************************/
													  double dRhj[4][4];
													  for (int ij =0; ij<4; ij++){
														  for(int ji =0 ;ji<4; ji++){
															  dRhj[ij][ji]=999.;
														  }
													  }	
													  match = false;	
													  dRhj[0][0] = h1b.DeltaR(ajet1_p4);
													  dRhj[1][0] = h1bar.DeltaR(ajet1_p4);
													  dRhj[0][1] = h1b.DeltaR(ajet2_p4);
													  dRhj[1][1] = h1bar.DeltaR(ajet2_p4);
													  dRhj[0][2] = h1b.DeltaR(ajet3_p4);
													  dRhj[1][2] = h1bar.DeltaR(ajet3_p4);
													  dRhj[0][3] = h1b.DeltaR(ajet4_p4);
													  dRhj[1][3] = h1bar.DeltaR(ajet4_p4);
													  dRhj[2][2] = h2b.DeltaR(ajet3_p4);
													  dRhj[3][2] = h2bar.DeltaR(ajet3_p4);
													  dRhj[2][3] = h2b.DeltaR(ajet4_p4);
													  dRhj[3][3] = h2bar.DeltaR(ajet4_p4);	
													  dRhj[2][0] = h2b.DeltaR(ajet1_p4);
													  dRhj[3][0] = h2bar.DeltaR(ajet1_p4);
													  dRhj[2][1] = h2b.DeltaR(ajet2_p4);
													  dRhj[3][1] = h2bar.DeltaR(ajet2_p4);

													  if((((dRhj[0][1]<0.3 && dRhj[1][0]<0.3)||(dRhj[0][0]<0.3 && dRhj[1][1]<0.3))&&((dRhj[2][2]<0.3 && dRhj[3][3]<0.3)||(dRhj[3][2]<0.3 && dRhj[2][3]<0.3)))||(((dRhj[0][2]<0.3 && dRhj[1][3]<0.3)||(dRhj[0][3]<0.3 && dRhj[1][2]<0.3))&&((dRhj[2][0]<0.3 && dRhj[3][1]<0.3)||(dRhj[3][0]<0.3 && dRhj[2][1]<0.3)))){

														  match=true;int n=0;
														  h_nCand_true->Fill(nCC);
														  for (JetList::reverse_iterator iJet=jet4pt.rbegin(); iJet!=jet4pt.rend(); ++iJet)
														  {
															  //std::cout<<" n "<<n<<" iJet->second "<<iJet->second<< " iJet->first "<<iJet->first<<std::endl;
															  if(dRhj[0][iJet->second]<0.3) {h_jet_pt[n]->Fill(iJet->first-h1b.Pt()); corr[0] = iJet->first-h1b.Pt();}
															  if(dRhj[1][iJet->second]<0.3) {h_jet_pt[n]->Fill(iJet->first-h1bar.Pt()); corr[1] = iJet->first-h1bar.Pt();} 
															  if(dRhj[2][iJet->second]<0.3) {h_jet_pt[n]->Fill(iJet->first-h2b.Pt()); corr[2]=iJet->first-h2b.Pt();}
															  if(dRhj[3][iJet->second]<0.3) {h_jet_pt[n]->Fill(iJet->first-h2bar.Pt()); corr[3]=iJet->first-h2bar.Pt();}	
															  n++;			
														  }
													  }	
													  /**************** 5th jet ************/
													  if(nCJets==5){
														  //TLorentzVector jet5;	
														  float dRR[4];
														  for(int ii =0 ; ii<4; ii++) {dRR[ii]=999.;}
														  float dR_min = 999.;
														  for(int q=0; q<nCJets; q++){
															  if(q!=H1jet1_i && q!=H1jet2_i && q!=H2jet1_i && q!=H2jet2_i) 	{
																  jet5.SetPtEtaPhiE(ajetpT[q], ajeteta[q], ajetphi[q], ajetE[q]);
																  dRR[0]= jet5.DeltaR(ajet1_p4);
																  if(dRR[0]<dR_min) dR_min= dRR[0];
																  dRR[1]= jet5.DeltaR(ajet2_p4); if(dRR[1]<dR_min) dR_min= dRR[1];
																  dRR[2]= jet5.DeltaR(ajet3_p4); if(dRR[2]<dR_min) dR_min= dRR[2];
																  dRR[3]= jet5.DeltaR(ajet4_p4); if(dRR[3]<dR_min) dR_min= dRR[3]; 

															  }

														  }
														  dR_min_5thJet= dR_min;		
													  }
													  else dR_min_5thJet =-999;
											  } // Loop over 4th jet
										  } // Conditions on 3rd jet
									  } // Loop over 3rd jet
								  } // Conditions on 2nd jet
							  }
						  } // Loop over 1st jet
					  }
			  }
		  }

	  }  if (foundHH)
	  {											 TLorentzVector diJet2_p4, diJet1_p4;
		  diJet2_p4=ajet3_p4+ajet4_p4;
		  diJet1_p4=ajet1_p4+ajet2_p4;
		  if(match) nCutGen++;		
		  H1_p4=diJet1_p4;
		  H2_p4=diJet2_p4;
		  float h11 = H1_p4.M();
		  float h22 = H2_p4.M();
		  h_H1_mass->Fill(h11);
		  h_H2_mass->Fill(h22);
		  h_Xmass->Fill((H1_p4+H2_p4).M());
		  int region=withinRegion(h11, h22, 20., 35., 125, 125);
		  if(region==0) h_Xmass_comb->Fill((H1_p4+H2_p4).M());
		  h_H1_H2_mass->Fill(h11, h22);
		  if(match)  {
			  h_Xmass_right->Fill((H1_p4+H2_p4).M());
			  matched=1;
		  }
		  else {
			  h_Xmass_wrong->Fill((H1_p4+H2_p4).M()); 
			  matched=0;
		  }
		  if((!match) &&( !foundHH)) ncutFail++;




		  outtree->Fill();			

		  h_nCand->Fill(nCC);	
		  ++nCut4;
	  } // found HH

	} // nCJets>3
      } // vType==5
    } // Trigger
    } // Event loop

  h_Cuts->Fill(1, nCut0);
  h_Cuts->Fill(3, nCut1);
  h_Cuts->Fill(5, nCut2);
  h_Cuts->Fill(7, nCut3);
  h_Cuts->Fill(9, nCut4);
  h_Cuts->Fill(11, nCut5);
  h_Cuts->Fill(13, nCut6);
  h_Cuts->Fill(15, nCut7);

  outtree->Write();
  outfile->Close();


  std::string histfilename="Histograms_"+sample+".root";
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_nCand->Write();
  h_nCand_true->Write();
  h_Xmass_right->Write();
  h_Xmass_wrong->Write();	
  h_Xmass_first->Write();
  h_Xmass_other->Write();
  h_nJets->Write();
  h_n8Jets->Write();
  h_nPV->Write();
  h_nPV_weighted->Write();
  h_JetpT1->Write();
  h_JetpT2->Write();
  h_JetpT3->Write();
  h_JetpT4->Write();
  h_MET->Write();
  h_MET_sig->Write();
  h_HT->Write();
  h_H1_mass->Write();
  h_H1_pT->Write();
  h_H2_mass->Write();
  h_H2_pT->Write();
  h_HH_mass_diff_cand->Write();
  h_HH_massNorm_diff->Write();
  h_genX_mass->Write();
  h_genH1_mass->Write();
  h_genH2_mass->Write();
  h_H1_H2_mass->Write();
  h_HH_CSV->Write();
  h_Xmass_comb->Write();
  h_Xmass->Write();
  h_Xmass2->Write(); 
  h_Cuts->Write();
  h_jet_pt[0]->Write();	
  h_jet_pt[1]->Write();
  h_jet_pt[2]->Write();
  h_jet_pt[3]->Write();
  tFile->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  std::cout<<"Number of events at the end of step 2 = "<<nCut0<<std::endl;
  std::cout<<"Number of events after trigger = "<<(float)nCut1/nCut0<<std::endl;
  std::cout<<"Number of events after trigger 50 = "<<(float)nCut50/nCut0<<std::endl;
  std::cout<<"Number of events after trigger  csv = "<<(float)nCutT/nCut0<<std::endl;
  std::cout<<"Number of events after trigger  or = "<<(float)nCutOR/nCut0<<std::endl;	
  std::cout<<"Number of events after Vtype==5 = "<<(float)nCut2/nCut1<<std::endl;
  std::cout<<"Number of events after nCJets>3 = "<<(float)nCut3/nCut2<<std::endl;
  std::cout<<"Number of events after finding HH kinematic candidate = "<<(float)nCut4/nCut3<<" ncut3 "<<nCut3<<std::endl;
  std::cout<<"Number of matched events "<<(float)nCutGen/nCut4<<"   "<<contComb<<"  nCut3 : "<<nCut3<<"  fail "<<ncutFail<<" ratio "<<(float)ncutFail/nCutGen<< std::endl; 	
  std::cout<<"Number of candidates "<<(float)nCandidate<<"Numeber of event passing HHfound criterion "<<nCut4<<std::endl;	

  return 0;
}




