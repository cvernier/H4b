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
#include "/gpfs/ddn/cms/user/cvernier/H4b/CMSSW_5_3_3_patch2/src/UserCode/SouvikDas/CommonAnalysisTools/Significance.c"


using namespace RooFit ;

double H_mass=125.0;
double mH_diff_cut=40.;
double mH_mean_cut=20.;

double rebin=1;

double SR_lo=200.;
double SR_hi=1600.;

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

void Significance_Kinematic(std::string filename="Histograms_DiJetPt_glugluToX400ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root", double &signal=0., double &bg=0., double &bg_error=0.)
{
  TFile *DiJetPt_glugluToX500=new TFile(filename.c_str());
  std::cout<<"Opened Histograms_DiJetPt_glugluToX800ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root"<<std::endl;
  
  if (filename=="Histograms_DiJetPt_glugluToX400ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root") rebin=1;
  else if (filename=="Histograms_DiJetPt_glugluToX300ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root") rebin=1;
  else rebin=4;
  
  TH1F *h_mjjjj_3Tag_CR1_signal=(TH1F*)DiJetPt_glugluToX500->Get("h_mX_CR1"); h_mjjjj_3Tag_CR1_signal->SetTitle("; m_{jjjj} (GeV)");
  TH1F *h_mjjjj_3Tag_CR2_signal=(TH1F*)DiJetPt_glugluToX500->Get("h_mX_CR2");
  TH1F *h_mjjjj_3Tag_CR3_signal=(TH1F*)DiJetPt_glugluToX500->Get("h_mX_CR3");
  TH1F *h_mjjjj_3Tag_CR4_signal=(TH1F*)DiJetPt_glugluToX500->Get("h_mX_CR4");
  TH1F *h_mjjjj_3Tag_SR_signal=(TH1F*)DiJetPt_glugluToX500->Get("h_mX_SR");
  TH2F *h_mH_mH_3Tag_signal=(TH2F*)DiJetPt_glugluToX500->Get("h_mH_mH_3Tag");
  
  TFile *data_8TeVData2012=new TFile("Histograms_BJetPlusX_Run2012BCD_Skim.root");
  if (filename=="Histograms_DiJetPt_glugluToX300ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root") data_8TeVData2012=new TFile("Histograms_BJetPlusX_Run2012BCD_Skim.root");
  // TFile *data_8TeVData2012=new TFile("Histograms_8TeVData2012BCD_Skim.root");
  std::cout<<"Opened Histograms_8TeVData2012BCD_Data_2.root"<<std::endl;
  
  TH1F *h_mjjjj_3Tag_CR1=(TH1F*)data_8TeVData2012->Get("h_mX_CR1"); h_mjjjj_3Tag_CR1->SetTitle("; m_{jjjj} (GeV)");
  TH1F *h_mjjjj_3Tag_CR2=(TH1F*)data_8TeVData2012->Get("h_mX_CR2");
  TH1F *h_mjjjj_3Tag_CR3=(TH1F*)data_8TeVData2012->Get("h_mX_CR3");
  TH1F *h_mjjjj_3Tag_CR4=(TH1F*)data_8TeVData2012->Get("h_mX_CR4");
  TH1F *h_mjjjj_3Tag_SR=(TH1F*)data_8TeVData2012->Get("h_mX_SR");
  TH2F *h_mH_mH_3Tag_data=(TH2F*)data_8TeVData2012->Get("h_mH_mH_3Tag");
  
  h_mjjjj_3Tag_CR1_signal->Rebin(rebin);
  h_mjjjj_3Tag_CR2_signal->Rebin(rebin);
  h_mjjjj_3Tag_CR3_signal->Rebin(rebin);
  h_mjjjj_3Tag_CR4_signal->Rebin(rebin);
  if (filename=="Histograms_DiJetPt_glugluToX300ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root") h_mjjjj_3Tag_SR_signal->Rebin(1);
  else h_mjjjj_3Tag_SR_signal->Rebin(rebin);
  
  h_mjjjj_3Tag_CR1->Rebin(rebin);
  h_mjjjj_3Tag_CR2->Rebin(rebin);
  h_mjjjj_3Tag_CR3->Rebin(rebin);
  h_mjjjj_3Tag_CR4->Rebin(rebin);
  h_mjjjj_3Tag_SR->Rebin(rebin);
  
  h_mjjjj_3Tag_CR1_signal->SetLineColor(kOrange);   h_mjjjj_3Tag_CR1_signal->SetLineWidth(2);
  h_mjjjj_3Tag_CR2_signal->SetLineColor(kRed);      h_mjjjj_3Tag_CR2_signal->SetLineWidth(2);
  h_mjjjj_3Tag_CR3_signal->SetLineColor(kMagenta);  h_mjjjj_3Tag_CR3_signal->SetLineWidth(2);
  h_mjjjj_3Tag_CR4_signal->SetLineColor(kBlue);     h_mjjjj_3Tag_CR4_signal->SetLineWidth(2);
  h_mjjjj_3Tag_SR_signal->SetLineWidth(2);
  
  h_mjjjj_3Tag_CR1->SetLineColor(kOrange);   h_mjjjj_3Tag_CR1->SetLineWidth(2);
  h_mjjjj_3Tag_CR2->SetLineColor(kRed);      h_mjjjj_3Tag_CR2->SetLineWidth(2);
  h_mjjjj_3Tag_CR3->SetLineColor(kMagenta);  h_mjjjj_3Tag_CR3->SetLineWidth(2);
  h_mjjjj_3Tag_CR4->SetLineColor(kBlue);     h_mjjjj_3Tag_CR4->SetLineWidth(2);
  h_mjjjj_3Tag_SR->SetLineColor(kCyan);      h_mjjjj_3Tag_SR->SetLineWidth(2);
  
  TH1F *h_mjjjj_3Tag_CR24_signal=(TH1F*)h_mjjjj_3Tag_CR2_signal->Clone("h_mjjjj_3Tag_CR24_signal");
  h_mjjjj_3Tag_CR24_signal->Add(h_mjjjj_3Tag_CR4_signal);
  
  TH1F *h_mjjjj_3Tag_CR24=(TH1F*)h_mjjjj_3Tag_CR2->Clone("h_mjjjj_3Tag_CR24");
  h_mjjjj_3Tag_CR24->Add(h_mjjjj_3Tag_CR4);
  
  // Normalize data control region shape with numbers from signal region to pretend it is the signal region
  TH1F *h_mjjjj_3Tag_SR_pretend=(TH1F*)h_mjjjj_3Tag_CR24->Clone("h_mjjjj_3Tag_SR_pretend");
  // h_mjjjj_3Tag_SR_pretend->Scale(158.459/h_mjjjj_3Tag_CR24->GetSumOfWeights());
  h_mjjjj_3Tag_SR_pretend->Scale(277.542/h_mjjjj_3Tag_CR24->GetSumOfWeights());
  h_mjjjj_3Tag_SR_pretend->SetLineColor(kBlack);
  
  // Fill the RooDataSet here for the fake data
  TH1F *h_mjjjj_3Tag_SR_fakeData=(TH1F*)h_mjjjj_3Tag_CR24->Clone("h_mjjjj_3Tag_SR_pretend"); h_mjjjj_3Tag_SR_fakeData->Reset(); h_mjjjj_3Tag_SR_fakeData->SetLineColor(kBlack);
  std::cout<<"277.542/h_mjjjj_3Tag_CR24->GetSumOfWeights() = "<<277.542/h_mjjjj_3Tag_CR24->GetSumOfWeights()<<std::endl;
  if (filename=="Histograms_DiJetPt_glugluToX400ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root") fillUnweightedTH1(h_mjjjj_3Tag_CR24, h_mjjjj_3Tag_SR_fakeData, (16119/h_mjjjj_3Tag_CR24->GetSumOfWeights()));
  else fillUnweightedTH1(h_mjjjj_3Tag_CR24, h_mjjjj_3Tag_SR_fakeData, (277.542/h_mjjjj_3Tag_CR24->GetSumOfWeights()));
  
  // Calculate nSignal events given production cross section, branching fractions and efficiency
  double nSignal_init=100000.;
  double nSignal_now=h_mjjjj_3Tag_SR_signal->GetSumOfWeights();
  double eff_signal=nSignal_now/nSignal_init;
  std::cout<<"Signal efficiency = "<<eff_signal<<std::endl;
  // double totalLumi=13241.968; // /pb
  double totalLumi=18600.0; // /pb
  double prodXsec_1=1.0; // pb
  double prodXsec_2=.5; // pb
  double prodXsec_3=2; // pb
  std::cout<<"With prod X sec of "<<prodXsec_1<<", we produce "<<totalLumi*prodXsec_1<<" events, of which remain "<<totalLumi*prodXsec_1*eff_signal<<std::endl;
  std::cout<<"With prod X sec of "<<prodXsec_2<<", we produce "<<totalLumi*prodXsec_2<<" events, of which remain "<<totalLumi*prodXsec_2*eff_signal<<std::endl;
  std::cout<<"With prod X sec of "<<prodXsec_3<<", we produce "<<totalLumi*prodXsec_3<<" events, of which remain "<<totalLumi*prodXsec_3*eff_signal<<std::endl;
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  // Draw modeled background with signal embossed 
  TH1F *h_mjjjj_3Tag_SR_plusSig1=(TH1F*)h_mjjjj_3Tag_SR_pretend->Clone("h_mjjjj_3Tag_SR_plusSig1");
  TH1F *h_mjjjj_3Tag_SR_plusSig2=(TH1F*)h_mjjjj_3Tag_SR_pretend->Clone("h_mjjjj_3Tag_SR_plusSig2");
  TH1F *h_mjjjj_3Tag_SR_plusSig3=(TH1F*)h_mjjjj_3Tag_SR_pretend->Clone("h_mjjjj_3Tag_SR_plusSig3");
  h_mjjjj_3Tag_SR_plusSig1->Add(h_mjjjj_3Tag_SR_signal, prodXsec_1*totalLumi/nSignal_init);
  h_mjjjj_3Tag_SR_plusSig2->Add(h_mjjjj_3Tag_SR_signal, prodXsec_2*totalLumi/nSignal_init);
  h_mjjjj_3Tag_SR_plusSig3->Add(h_mjjjj_3Tag_SR_signal, prodXsec_3*totalLumi/nSignal_init);
  h_mjjjj_3Tag_SR_plusSig1->SetLineColor(kBlue);
  h_mjjjj_3Tag_SR_plusSig2->SetLineColor(kRed);
  h_mjjjj_3Tag_SR_plusSig3->SetLineColor(kBlue);
  std::cout<<"weight = "<<prodXsec_1*totalLumi/nSignal_init<<std::endl;
  std::cout<<"h_mjjjj_3Tag_SR_signal->GetBinWidth() = "<<h_mjjjj_3Tag_SR_signal->GetBinWidth(1)<<", h_mjjjj_3Tag_SR_fakeData->GetBinWidth() = "<<h_mjjjj_3Tag_SR_fakeData->GetBinWidth(1)<<std::endl;
  fillUnweightedTH1(h_mjjjj_3Tag_SR_signal, h_mjjjj_3Tag_SR_fakeData, prodXsec_1*totalLumi/nSignal_init);
  TCanvas *c_mjjjj_3Tag=new TCanvas("c_mjjjj_3Tag", "c_mjjjj_3Tag", 700, 700);
  h_mjjjj_3Tag_SR_plusSig1->SetTitle("m_{X} for signal overlaid with data; m_{X} GeV");
  h_mjjjj_3Tag_SR_plusSig1->GetXaxis()->SetRangeUser(250., 1100.);
  // h_mjjjj_3Tag_SR_plusSig1->GetXaxis()->SetRangeUser(200., 1200.);
  h_mjjjj_3Tag_SR_plusSig1->Draw("E");
  h_mjjjj_3Tag_SR_pretend->Draw("E SAME"); // h_mjjjj_3Tag_SR_pretend->Draw("C HIST SAME")
  h_mjjjj_3Tag_SR_plusSig1->Draw("E SAME"); // h_mjjjj_3Tag_SR_plusSig1->Draw("C HIST SAME");
  h_mjjjj_3Tag_SR_plusSig2->Draw("E SAME"); // h_mjjjj_3Tag_SR_plusSig2->Draw("C HIST SAME");
  // h_mjjjj_3Tag_SR_plusSig3->Draw("E SAME"); // h_mjjjj_3Tag_SR_plusSig3->Draw("C HIST SAME");
  // h_mjjjj_3Tag_SR->Draw("E SAME");
  TLegend *leg=new TLegend(0.6, 0.6, 0.9, 0.9);
  leg->AddEntry(h_mjjjj_3Tag_SR_plusSig1, "Bg + Sg (1 pb)");
  leg->AddEntry(h_mjjjj_3Tag_SR_plusSig2, "Bg + Sg (0.5 pb)");
  // leg->AddEntry(h_mjjjj_3Tag_SR_plusSig3, "Bg + Sg (0.1 pb)");
  leg->Draw();
  // Fit to background function
  TF1 *f_mjjjj_3Tag_SR_pretend=new TF1("f_mjjjj_3Tag_SR_pretend", crystalBall, 250., 1100., 5);
  f_mjjjj_3Tag_SR_pretend->SetLineWidth(0);
  f_mjjjj_3Tag_SR_pretend->SetParameter(0, 3.49638e-01);
  f_mjjjj_3Tag_SR_pretend->SetParameter(1, 1.08296e+02);
  f_mjjjj_3Tag_SR_pretend->SetParameter(2, 4.66769e+02);
  f_mjjjj_3Tag_SR_pretend->SetParameter(3, 6.17764e+01);
  f_mjjjj_3Tag_SR_pretend->SetParLimits(4, 20, 100);
  h_mjjjj_3Tag_SR_pretend->Fit(f_mjjjj_3Tag_SR_pretend, "REFMN");
  f_mjjjj_3Tag_SR_pretend->Draw("SAME");
  c_mjjjj_3Tag->SaveAs((filename+".png").c_str());
  
  // Draw fake data
  TCanvas *c_fakeData=new TCanvas("c_fakeData", "Fake Data", 700, 700);
  h_mjjjj_3Tag_SR_fakeData->Draw();
  
  // Fit with RooFit
  /*
  TCanvas *c_Illustration=new TCanvas("c_Illustration", "c_Illustration", 700, 700);
  gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/RevCrystalBall_cxx.so");
  RooRealVar mx("mx", "m_{X} (Gev)", 250., 1100.);
  RooRealVar bg_p0("bg_p0", "bg_p0", 0.1, 1.0);
  RooRealVar bg_p1("bg_p1", "bg_p1", 1., 5.);
  RooRealVar bg_p2("bg_p2", "bg_p2", 350., 850.);
  RooRealVar bg_p3("bg_p3", "bg_p3", 10., 50.);
  RevCrystalBall bg_func("bg_func", "Background Prediction PDF", mx, bg_p0, bg_p1, bg_p2, bg_p3);
  RooRealVar sg_p0("sg_p0", "sg_p0", 1.0, 5.0);
  RooRealVar sg_p1("sg_p1", "sg_p1", 1.0., 5.0.);
  RooRealVar sg_p2("sg_p2", "sg_p2", 250., 850.);
  RooRealVar sg_p3("sg_p3", "sg_p3", 25., 35.);
  RevCrystalBall sg_func("sg_func", "Signal", mx, sg_p0, sg_p1, sg_p2, sg_p3);
  RooRealVar bgfrac("bgfrac", "Background Fraction", 0.5, 0.1, 1.0);
  RooAddPdf model("model", "Back+sig", RooArgList(bg_func, sg_func), bgfrac);
  RooDataHist bg_illustration("bg_illustration", "Bg Illustration", RooArgList(mx), h_mjjjj_3Tag_SR_pretend);
  RooDataHist illustration("illustration", "Data", RooArgList(mx), h_mjjjj_3Tag_SR_plusSig1);
  // bg_func.fitTo(bg_illustration);
  model.fitTo(illustration);
  RooPlot *plotIllus=mx.frame();
  illustration.plotOn(plotIllus);
  // bg_illustration.plotOn(plotIllus);
  // model.plotOn(plotIllus);
  // bg_func.plotOn(plotIllus);
  plotIllus->Draw();
  c_Illustration->SaveAs("c_Illustration.png");
  */
  
  // Count background and 1 pb-signal events between 400 - 800 GeV
  double signal_400_800=h_mjjjj_3Tag_SR_plusSig1->Integral(h_mjjjj_3Tag_SR_plusSig1->FindBin(250.), h_mjjjj_3Tag_SR_plusSig1->FindBin(1100.));
  double bg_error_400_800;
  double bg_400_800=h_mjjjj_3Tag_SR_pretend->IntegralAndError(h_mjjjj_3Tag_SR_pretend->FindBin(250.), h_mjjjj_3Tag_SR_pretend->FindBin(1100.), bg_error_400_800);
 std::cout<<" here "<<std::endl;
  
  
  std::cout<<"signal_400_800 = "<<signal_400_800<<std::endl;
  std::cout<<"bg_400_800 = "<<bg_400_800<<" +- "<<bg_error_400_800<<std::endl;
  /*
  signal=signal_400_800-bg_400_800;
   std::cout<<" here "<<std::endl;
  bg=bg_400_800;
   std::cout<<" here "<<std::endl;
  bg_error=bg_400_800*0.20; // bg_error_400_800;
  std::cout<<" here "<<std::endl;
*/
  /*
  TFile *outFile=new TFile("../SignalBackgroundData.root", "UPDATE");
  TH1F *h_fakeData=(TH1F*)h_mjjjj_3Tag_SR_plusSig1->Clone("h_fakeData");
  h_fakeData->SetTitle("Fake 8 TeV data with 13.42 /fb (signal injected with x-sec = 1 pb)");
  h_fakeData->Write();
  outFile->Write();
  outFile->Close();
  */
  std::cout<<" here "<<std::endl;

  if (filename=="Histograms_DiJetPt_glugluToX300ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root") {SR_lo=200; SR_hi=1600;}


   RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);


  // RooDataHist data_obs("data_obs", "Data", RooArgList(x), h_mjjjj_3Tag_SR_plusSig1);
  

  RooDataHist data_obs("data_obs", "Data", RooArgList(x), h_mjjjj_3Tag_SR_fakeData);
  RooPlot *plot=x.frame();
  data_obs.plotOn(plot);
  TCanvas *c_data=new TCanvas("c_data", "c_data", 500, 500);
  plot->Draw();
  RooWorkspace *w=new RooWorkspace("HbbHbb");
  w->import(data_obs);
  std::string mass=filename.substr(28,3);
  std::cout<<"mass = "<<mass<<std::endl;
  w->SaveAs(("w_data_"+mass+".root").c_str());
  std::cout<<"1 pb signal = "<<h_mjjjj_3Tag_SR_signal->Integral(h_mjjjj_3Tag_SR_signal->FindBin(SR_lo), h_mjjjj_3Tag_SR_signal->FindBin(SR_hi))<<" (RAW) -- in 18.6 /fb --> "<<h_mjjjj_3Tag_SR_signal->Integral(h_mjjjj_3Tag_SR_signal->FindBin(SR_lo), h_mjjjj_3Tag_SR_signal->FindBin(SR_hi))*prodXsec_1*totalLumi/nSignal_init<<std::endl;
  std::cout<<"1 pb yield = "<<h_mjjjj_3Tag_SR_fakeData->Integral(h_mjjjj_3Tag_SR_fakeData->FindBin(SR_lo), h_mjjjj_3Tag_SR_fakeData->FindBin(SR_hi))<<std::endl;
  c_data->SaveAs(("c_data_"+mass+".png").c_str());
  
}
    
  
  

