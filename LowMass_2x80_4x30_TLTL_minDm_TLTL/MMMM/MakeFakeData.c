#include <TH1F.h>
#include <TH2F.h>
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
#include <TFractionFitter.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TArrow.h>
#include <TArc.h>
#include <TRandom3.h>

#include "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas//CommonAnalysisTools/Significance.c"


using namespace RooFit ;

double H_mass=125.0;
double mH_diff_cut=40.;
double mH_mean_cut=20.;

double rebin=2;

double SR_lo=380.;
double SR_hi=1600.;

SR_lo=200.;

double r1=15., r2=35.;

int fit=2; // 0 = Fit 2&4, 1 = Fit 1, 24, 3, 2 = Kolmogorov with 2+4.

double Normalize(TH1* h)
{
  double nEntries=h->GetSumOfWeights();
  h->Scale(1./nEntries);
  return nEntries;
}


void fillUnweightedTH1(TH1F *h_in, TH1F *h_out, double weight=1.)
{
  TRandom3 *r3=new TRandom3();
  for (unsigned int i=1; i<=h_in->GetNbinsX(); ++i)
  {
    double binCenter=h_in->GetBinCenter(i);
    double entry=h_in->GetBinContent(i);
    // std::cout<<"entry = "<<entry<<std::endl;
    for (unsigned int j=0; j<entry; ++j)
    {
      double throw=r3->Rndm();
      // std::cout<<"throw = "<<throw<<std::endl;
      if (throw<weight) 
      {
        x=binCenter;
        // std::cout<<"x = "<<x<<std::endl;
        h_out->Fill(x);
      }
    }
  }
}

void fillByThrowTH1(TH1F *h_in, TH1F *h_out, double weight=1.)
{
  TRandom3 *r3=new TRandom3();
  
  // Decide total number of throws
  double max_h_in=h_in->GetMaximum();
  double new_max_h_in=max_h_in*weight;
  int throws=int(new_max_h_in*50.);
  std::cout<<"throws = "<<throws<<std::endl;
  std::cout<<"h_in->GetNbinsX() = "<<h_in->GetNbinsX()<<std::endl;
  for (unsigned int i=1; i<=h_in->GetNbinsX(); ++i)
  {
    double binCenter=h_in->GetBinCenter(i);
    double entry=h_in->GetBinContent(i);
    
    double newEntry=entry*weight;
    double prob=newEntry/double(throws);
    for (unsigned int j=0; j<throws; ++j)
    {
      double throw=r3->Rndm();
      if (throw<prob) h_out->Fill(binCenter);
    }
  }
}   

Double_t crystalBall(Double_t *x, Double_t *par)
{
  Double_t std=(x[0]-par[2])/par[3];
  Double_t A=pow(par[1]/par[0], par[1])*exp(-0.5*pow(par[0], 2));
  Double_t B=par[1]/par[0]-par[0];
  Double_t result=0.;
  
  if (std<par[0]) // Gaussian region
  {
    result=exp(-0.5*pow(std, 2));
  }
  else // Power Law region
  {
    result=A/pow(B+std, par[1]);
  }
  
  result=result*par[4];
  
  return result;
}

void MakeFakeData()
{
  std::string filename="../MMMM_nominal/a_KinFit/Histograms_DiJetPt_glugluToX300ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root";
  int rebin=4;
  
  TFile *signal=new TFile(filename.c_str());
  TH1F *h_mX_SR_signal=(TH1F*)signal->Get("h_mX_SR");
  
  TFile *data_8TeVData2012=new TFile("Histograms_BJetPlusX_Run2012BCD_Skim.root");
  TH1F *h_mX_CR2=(TH1F*)data_8TeVData2012->Get("h_mX_CR2");
  TH1F *h_mX_CR4=(TH1F*)data_8TeVData2012->Get("h_mX_CR4");
  TH1F *h_mX_CR24=(TH1F*)h_mX_CR2->Clone("h_mX_CR24");
  h_mX_CR24->Add(h_mX_CR4);
  
  h_mX_SR_signal->Rebin(rebin);
  h_mX_CR24->Rebin(rebin);
  
  // Normalize data control region shape with numbers from signal region to pretend it is the signal region
  TH1F *h_mX_SR_fakeData=(TH1F*)h_mX_CR24->Clone("h_mX_SR_fakeData"); 
  h_mX_SR_fakeData->Reset(); 
  h_mX_SR_fakeData->SetLineColor(kBlack);
  std::cout<<"Background CR->SR Scale, 277.542/h_mX_CR24->GetSumOfWeights() = "<<277.542/h_mX_CR24->GetSumOfWeights()<<std::endl;
  // fillUnweightedTH1(h_mX_CR24, h_mX_SR_fakeData, (277.542/h_mX_CR24->GetSumOfWeights()));
  fillByThrowTH1(h_mX_CR24, h_mX_SR_fakeData, 1);//(277.542/h_mX_CR24->GetSumOfWeights()));
  
  // Calculate nSignal events given production cross section, branching fractions and efficiency
  TH1F *h_CountWithPU=(TH1F*)signal->Get("CountWithPU");
  double nSignal_init=h_CountWithPU->GetBinContent(1);
  double nSignal_now=h_mX_SR_signal->GetSumOfWeights();
  double eff_signal=nSignal_now/nSignal_init;
  std::cout<<"Signal efficiency = "<<eff_signal<<std::endl;
  // double totalLumi=13241.968; // /pb
  double totalLumi=18600.0; // /pb
  double prodXsec_1=1.0; // pb
  
  // fillUnweightedTH1(h_mX_SR_signal, h_mX_SR_fakeData, prodXsec_1*totalLumi/nSignal_init);
  fillByThrowTH1(h_mX_SR_signal, h_mX_SR_fakeData, prodXsec_1*totalLumi/nSignal_init);
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);
  RooDataHist data_obs("data_obs", "Data", RooArgList(x), h_mX_SR_fakeData);
  RooPlot *plot=x.frame();
  data_obs.plotOn(plot);
  TCanvas *c_data=new TCanvas("c_data", "c_data", 500, 500);
  plot->Draw();
  RooWorkspace *w=new RooWorkspace("HbbHbb");
  w->import(data_obs);
  std::string mass=filename.substr(27,3);
  std::cout<<"mass = "<<mass<<std::endl;
  w->SaveAs(("w_data_"+mass+".root").c_str());
  std::cout<<"1 pb signal = "<<h_mX_SR_signal->Integral(h_mX_SR_signal->FindBin(SR_lo), h_mX_SR_signal->FindBin(SR_hi))<<" (RAW) -- in 18.6 /fb --> "<<h_mX_SR_signal->Integral(h_mX_SR_signal->FindBin(SR_lo), h_mX_SR_signal->FindBin(SR_hi))*prodXsec_1*totalLumi/nSignal_init<<std::endl;
  std::cout<<"1 pb yield = "<<h_mX_SR_fakeData->Integral(h_mX_SR_fakeData->FindBin(SR_lo), h_mX_SR_fakeData->FindBin(SR_hi))<<std::endl;
  c_data->SaveAs(("c_data_"+mass+".png").c_str());
} 
