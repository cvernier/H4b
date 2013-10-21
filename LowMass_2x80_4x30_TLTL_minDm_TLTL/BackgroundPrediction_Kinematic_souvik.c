#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <TStyle.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraphErrors.h>

double H_mass=125.0;
double mH_diff_cut=40.;
double mH_mean_cut=20.;

double rebin=1;
bool bReg=false;

std::string tags="MMMM"; // MMMM

double VR_lo=250.;
double VR_hi=600.;
double SR_lo=250.;
double SR_hi=600.;

double r1=15., r2=35.;

double quad(double a, double b, double c=0, double d=0, double e=0, double f=0)
{
  return pow(a*a+b*b+c*c+d*d+e*e+f*f, 0.5);
}

void fillByThrowTH1(TH1 *h_in, TH1 *h_out, double weight=1.)
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



TCanvas* comparePlots(TH1F *data, TH1F *qcd, std::string title)
{
  TCanvas *c=new TCanvas(("c"+title).c_str(), "c", 600, 700);
  TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.3);
  TPad *p_1=new TPad("p_1", "p_1", 0, 0.3, 1, 1);
  p_1->Draw();
  p_2->Draw();
  p_1->cd();
  // p_1->SetLogy();
  qcd->SetTitle((title+"; m_{X} (GeV)").c_str());
  double s=data->Integral(data->FindBin(200.), data->FindBin(600.))/qcd->Integral(qcd->FindBin(200.), qcd->FindBin(600.));
  qcd->Scale(s);
  qcd->Draw("HIST");
  data->Draw("Ep9 SAME");
  TLegend *leg=new TLegend(0.5, 0.9, 0.9, 0.7);
  leg->AddEntry(qcd, "Sideband Region");
  leg->AddEntry(data, "Central Region");
  leg->Draw();
  p_2->cd();
  p_2->SetGridy();
  TH1F *h_ratio=(TH1F*)data->Clone("h_ratio");
  h_ratio->SetTitle(("Data/MC Ratio "+title+" ; data/MC").c_str());
  h_ratio->Divide(qcd);                          
  h_ratio->SetMinimum(-1.); h_ratio->SetMaximum(3.);                  
  h_ratio->Draw();                                                    
  qcd->Scale(1./s); 
  p_1->cd(); 
  return c;                           
}

TCanvas* comparePlots2(RooPlot *plot_bC, RooPlot *plot_bS, TH1F *data, TH1F *qcd, std::string title)
{
  TCanvas *c=new TCanvas(("c_RooFit_"+title).c_str(), "c", 1000, 1000);
  TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.3);
  TPad *p_1=new TPad("p_1", "p_1", 0, 0.3, 1, 1);
  p_1->Draw();
  p_2->Draw();
  p_1->cd();
  plot_bS->Draw();
  plot_bC->Draw("same");
  p_2->cd();
  p_2->SetGridy();
  TH1F *h_ratio=(TH1F*)data->Clone("h_ratio");
  h_ratio->SetTitle(("VR/VR-SB Ratio "+title+" ; VR/VR-SB Ratio").c_str());
  h_ratio->Divide(qcd);
  h_ratio->GetXaxis()->SetRangeUser(VR_lo, VR_hi);                         
  h_ratio->SetMinimum(-1.); h_ratio->SetMaximum(3.);                  
  h_ratio->Draw();
  p_1->cd(); 
  return c;                           
}

// 0 = cut (x sigma), 1 = power, 2 = center, 3 = sigma
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

double returnCurveSyst(TF1 *h_CR, TF1 *h_SR, double mass)
{
  double err_lo, err_hi;
  if (mass==300) {err_lo=200.; err_hi=600.;}
  else if (mass==400) {err_lo=200.; err_hi=600.;}
  else if (mass==500) {err_lo=300.; err_hi=700.;}
  else if (mass==600) {err_lo=300.; err_hi=800.;}
  else if (mass==700) {err_lo=500.; err_hi=900.;}
  else if (mass==800) {err_lo=600.; err_hi=1000.;}
  double int_mMMMMb_3Tag_CR24=h_CR->Integral(err_lo, err_hi);
  double int_mMMMMb_3Tag_SR=h_SR->Integral(err_lo, err_hi);
  double fracSyst=2.*(int_mMMMMb_3Tag_SR-int_mMMMMb_3Tag_CR24)/(int_mMMMMb_3Tag_CR24+int_mMMMMb_3Tag_SR);
  return fracSyst;
} 

void BackgroundPrediction_Kinematic_souvik()
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  if (bReg) tags=tags+"_bReg";
  
/*  const unsigned int nPoints=4;
  double mass[nPoints]={90., 107.5, 142.5, 160.};
  double n_SB[nPoints], n_SR[nPoints];
  double ratio[nPoints];
  double errorsY[nPoints], errorsX[nPoints];
 */ 
  // === MMMM/b ===
  // TFile *f_MMMM_b=new TFile((tags+"/b/Histograms_8TeVData2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_b=new TFile((tags+"/b/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_MMMM_b=new TFile((tags+"/b/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_b=new TFile((tags+"/b/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TH1F *h_mX_CR2_b=(TH1F*)f_MMMM_b->Get("h_mX_CR2");
  TH1F *h_mX_CR4_b=(TH1F*)f_MMMM_b->Get("h_mX_CR4");
  TH1F *h_mX_SR_b=(TH1F*)f_MMMM_b->Get("h_mX_SR");
/*  n_SB[0]=(h_mX_CR2_b->GetSumOfWeights()+h_mX_CR4_b->GetSumOfWeights());
  n_SR[0]=h_mX_SR_b->GetSumOfWeights();
  
  // === MMMM/d ===
  // TFile *f_MMMM_d=new TFile((tags+"/d/Histograms_8TeVData2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_d=new TFile((tags+"/d/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_MMMM_d=new TFile((tags+"/c/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TH1F *h_mX_CR2_d=(TH1F*)f_MMMM_d->Get("h_mX_CR2");
  TH1F *h_mX_CR4_d=(TH1F*)f_MMMM_d->Get("h_mX_CR4");
  TH1F *h_mX_SR_d=(TH1F*)f_MMMM_d->Get("h_mX_SR");
  n_SB[1]=(h_mX_CR2_d->GetSumOfWeights()+h_mX_CR4_d->GetSumOfWeights());
  n_SR[1]=h_mX_SR_d->GetSumOfWeights();
  
  // === MMMM/e ===
  // TFile *f_MMMM_e=new TFile((tags+"/e/Histograms_8TeVData2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_e=new TFile((tags+"/e/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_MMMM_e=new TFile((tags+"/d/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TH1F *h_mX_CR2_e=(TH1F*)f_MMMM_e->Get("h_mX_CR2");
  TH1F *h_mX_CR4_e=(TH1F*)f_MMMM_e->Get("h_mX_CR4");
  TH1F *h_mX_SR_e=(TH1F*)f_MMMM_e->Get("h_mX_SR");
  n_SB[2]=(h_mX_CR2_e->GetSumOfWeights()+h_mX_CR4_e->GetSumOfWeights());
  n_SR[2]=h_mX_SR_e->GetSumOfWeights();
  
  // === MMMM/c ===
  // TFile *f_MMMM_c=new TFile((tags+"/c/Histograms_8TeVData2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_c=new TFile((tags+"/c/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_MMMM_c=new TFile((tags+"/e/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TH1F *h_mX_CR2_c=(TH1F*)f_MMMM_c->Get("h_mX_CR2");
  TH1F *h_mX_CR4_c=(TH1F*)f_MMMM_c->Get("h_mX_CR4");
  TH1F *h_mX_SR_c=(TH1F*)f_MMMM_c->Get("h_mX_SR");
  n_SB[3]=(h_mX_CR2_c->GetSumOfWeights()+h_mX_CR4_c->GetSumOfWeights());
  n_SR[3]=h_mX_SR_c->GetSumOfWeights();
  */
  // TFile *f_MMMM_a=new TFile((tags+"/a/Histograms_8TeVData2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_a=new TFile((tags+"/a/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  double mass=125.;
  double n_SB, n_SR;
  double ratio;
  double errorsY, errorsX;

  TFile *f_MMMM_a=new TFile((tags+"/a/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TH1F *h_mX_CR2_a=(TH1F*)f_MMMM_a->Get("h_mX_CR2");
  TH1F *h_mX_CR4_a=(TH1F*)f_MMMM_a->Get("h_mX_CR4");
  TH1F *h_mX_SR_a=(TH1F*)f_MMMM_a->Get("h_mX_SR");
  n_SB=(h_mX_CR2_a->GetSumOfWeights()+h_mX_CR4_a->GetSumOfWeights());
  n_SR=h_mX_SR_a->GetSumOfWeights();
  ratio=n_SR/n_SB;
  errorsY=ratio*pow(1./n_SR+1./n_SB, 0.5);
  errorsX=0.;
  // TFile *f_MMMM_a=new TFile((tags+"/a/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  
  /*for (unsigned int i=0; i<nPoints; ++i)
  {
    // if (n_SB[i]>0)
    if (n_SB[i]>0 && i!=1 && i!=2)
    {
      ratio[i]=n_SR[i]/n_SB[i];
      errorsY[i]=ratio[i]*pow(1./n_SR[i]+1./n_SB[i], 0.5);
      errorsX[i]=0.;
    } else 
    {
      ratio[i]=1;
      errorsX[i]=1;
      errorsY[i]=0;
    }
    std::cout<<"mass["<<i<<"] = "<<mass[i]<<std::endl;
    std::cout<<"n_SR["<<i<<"] = "<<n_SR[i]<<std::endl;
    std::cout<<"n_SB["<<i<<"] = "<<n_SB[i]<<std::endl;
  }
  */
  /*TGraphErrors *g_ratio=new TGraphErrors(nPoints, mass, ratio, errorsX, errorsY);
  g_ratio->SetTitle("SR/SB ratio");
  TCanvas *c_ratio=new TCanvas("c_ratio", "c_ratio", 700, 700);
  g_ratio->SetMinimum(-1.); g_ratio->SetMaximum(2.);
  g_ratio->Draw("A*");
  TF1 *f_ratio=new TF1("f_ratio", "pol1");
  g_ratio->Fit(f_ratio);
  c_ratio->SaveAs(("c_ratio_"+tags+".png").c_str());
  */
  double ratioAt125 =ratio;
  double errorAt125=errorsY;
  std::cout<<"ratioAt125 = "<<ratioAt125<<" +- "<<errorAt125<<std::endl;
  
  std::cout<<" = MMMM b ======================================== "<<std::endl;
  TH1F *h_mMMMMb_3Tag_CR2=(TH1F*)f_MMMM_b->Get("h_mX_CR2");
  TH1F *h_mMMMMb_3Tag_CR4=(TH1F*)f_MMMM_b->Get("h_mX_CR4");
  TH1F *h_mMMMMb_3Tag_SR=(TH1F*)f_MMMM_b->Get("h_mX_SR");
  // TH1F *h_mMMMMb_3Tag_SR=(TH1F*)f_MMMM_b->Get("h_X_mass_SR");
  
  TH1F *h_mMMMMb_3Tag_CR24=(TH1F*)h_mMMMMb_3Tag_CR2->Clone("h_mX_CR24");
  h_mMMMMb_3Tag_CR24->Add(h_mMMMMb_3Tag_CR4);
  h_mMMMMb_3Tag_CR24->GetXaxis()->SetRangeUser(VR_lo, VR_hi);
  h_mMMMMb_3Tag_SR->GetXaxis()->SetRangeUser(VR_lo, VR_hi);
  h_mMMMMb_3Tag_CR24->Rebin(rebin);
  h_mMMMMb_3Tag_SR->Rebin(rebin);
  h_mMMMMb_3Tag_CR24->SetLineColor(kRed);
  h_mMMMMb_3Tag_SR->SetLineColor(kBlue);
  double bS=h_mMMMMb_3Tag_SR->GetSumOfWeights();
  double bC=h_mMMMMb_3Tag_CR24->GetSumOfWeights();
  std::cout<<"Number of events in MMMM b signal region = "<<bS<<std::endl;
  std::cout<<"bC = "<<bC<<", bS = "<<bS<<std::endl;
h_mMMMMb_3Tag_SR->Scale(bC/bS);
  // Do the fits using RooFit
  gSystem->Load("LowMassBackgroundPDF_cxx.so");
  gSystem->Load("ExpGaussExp_cxx.so");
 RooRealVar *x, *sg_p0, *sg_p1, *sg_p2, *sg_p3, *sg_p4, *sg_p5, *sg_p6;
  x=new RooRealVar("x", "m_{X} (GeV)", VR_lo, VR_hi);
  
//RooRealVar x("x", "m_{X} (GeV)", VR_lo, VR_hi);
  // bC
   RooRealVar bC_p0("bC_p0", "bC_p0", 0., 500.);
  RooRealVar bC_p1("bC_p1", "bC_p1", 0., 150.1);
  RooRealVar bC_p2("bC_p2", "bC_p2", 0., 15.);
  RooRealVar bC_p3("bC_p3", "bC_p3", 0., 10.);

  /*RooRealVar bC_p0("bC_p0", "bC_p0", 250., 500.);
  RooRealVar bC_p1("bC_p1", "bC_p1", 40., 150.1);
  RooRealVar bC_p2("bC_p2", "bC_p2", 0.1, 5.1);
  RooRealVar bC_p3("bC_p3", "bC_p3", 0., 1.);
  RooRealVar bC_p4("bC_p4", "bC_p4", 250., 350.);*/
    sg_p0=new RooRealVar("sg_p0", "sg_p0", 0., 1.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 0., 10.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 10.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 10.);
  sg_p4=new RooRealVar("sg_p4", "sg_p4", 0., 10.);
sg_p5=new RooRealVar("sg_p5", "sg_p5", 0., 5.);
sg_p6=new RooRealVar("sg_p6", "sg_p6", 0., 5.);

  ExpGaussExp bC_fit("bC_fit", "bC_fit", *x, bC_p0, bC_p1, bC_p2, bC_p3);//, bC_p4);
  RooDataHist bC_data("bC_data", "bC Data", RooArgList(*x), h_mMMMMb_3Tag_CR24);
  //   RooCBShape bC_fit("bC_fit", "Signal Prediction", *x, *m0, *sigma,*alpha, *n);
  RooBernstein polyn("polyn", "Combinatoric", *x, RooArgList(*sg_p1,*sg_p2, *sg_p3, *sg_p4, *sg_p5 , *sg_p6));
  // RooGaussian polyn("polyn", "Combinatoric",*x,*sg_p1,*sg_p2);
        RooAddPdf signal("signal", "signal", RooArgList(bC_fit, polyn), *sg_p0);

  bC_fit.fitTo(bC_data, RooFit::Range(VR_lo, VR_hi));
  RooPlot *bC_plot=x.frame();
  bC_data.plotOn(bC_plot, RooFit::MarkerColor(kRed));
  bC_fit.plotOn(bC_plot, RooFit::LineColor(kRed));
  
  bC_plot->Draw();
  //return;
    // bS
  RooRealVar *x, *sg_0, *sg_1, *sg_2, *sg_3, *sg_4, *sg_5, *sg_6; 
   RooRealVar bS_p0("bS_p0", "bS_p0", 0., 500.);
  RooRealVar bS_p1("bS_p1", "bS_p1", 0., 150.1);
  RooRealVar bS_p2("bS_p2", "bS_p2", 0., 15.);
  RooRealVar bS_p3("bS_p3", "bS_p3", 0., 10.);
  /*RooRealVar bS_p0("bS_p0", "bS_p0", 250., 500.);
  RooRealVar bS_p1("bS_p1", "bS_p1", 40., 150.1);
  RooRealVar bS_p2("bS_p2", "bS_p2", 0.1, 5.1);
  RooRealVar bS_p3("bS_p3", "bS_p3", 0., 1.);
  RooRealVar bS_p4("bS_p4", "bS_p4", 250., 350.);
  sg_0=new RooRealVar("sg_0", "sg_0", 0., 1.);
  sg_1=new RooRealVar("sg_1", "sg_1", 50., 150.);
  sg_2=new RooRealVar("sg_2", "sg_2", 0., 60.);
  sg_3=new RooRealVar("sg_3", "sg_3", 0., 10.);
*/

  sg_0=new RooRealVar("sg_0", "sg_0", 0., 1.);
  sg_1=new RooRealVar("sg_1", "sg_1", 0., 10.);
  sg_2=new RooRealVar("sg_2", "sg_2", 0., 10.);
  sg_3=new RooRealVar("sg_3", "sg_3", 0., 10.);
  sg_4=new RooRealVar("sg_4", "sg_4", 0., 10.);
  sg_5=new RooRealVar("sg_5", "sg_5", 0., 5.); 
 sg_6=new RooRealVar("sg_6", "sg_6", 0., 5.);

  //LowMassBackgroundPDF bS_fit("bS_fit", "bS_fit", x, bS_p0, bS_p1, bS_p2, bS_p3, bS_p4);
  ExpGaussExp bS_fit("bS_fit", "bS_fit", *x, bS_p0, bS_p1, bS_p2, bS_p3);//, bS_p4);
  RooDataHist bS_data("bS_data", "bS Data", RooArgList(*x), h_mMMMMb_3Tag_SR);
  RooBernstein polyn_("polyn_", "Combinatoric", *x, RooArgList(*sg_1,*sg_2, *sg_3, *sg_4, *sg_5, *sg_6));
  RooAddPdf signal_("signal", "signal", RooArgList(bS_fit, polyn_), *sg_0);
  //RooDataHist bS_data("bS_data", "bS Data", RooArgList(x), h_mMMMMb_3Tag_SR);
  bS_fit.fitTo(bS_data, RooFit::Range(VR_lo, VR_hi));
  RooPlot *bS_plot=x.frame();
  bS_data.plotOn(bS_plot, RooFit::MarkerColor(kBlue));
  bS_fit.plotOn(bS_plot, RooFit::LineColor(kBlue));
  std::cout<<" === === "<<std::endl;
  std::cout<<"chi^2/ndof of bC = "<<bC_plot->chiSquare()<<std::endl;
  std::cout<<"chi^2/ndof of bS = "<<bS_plot->chiSquare()<<std::endl;
  std::cout<<" === === "<<std::endl;
  TCanvas *c_bC=comparePlots2(bC_plot, bS_plot, h_mMMMMb_3Tag_SR, h_mMMMMb_3Tag_CR24, "Kinematic Extrapolation in "+tags+" Validation Region of Data; m_{X} GeV");
  double x_mean_bC=bC_p0.getVal();
  double x_kHi_bC=bC_p0.getVal()+bC_p2.getVal()*bC_p1.getVal();
  double x_kLo_bC=bC_p4.getVal();
  TLine *l_mean_bC=new TLine(x_mean_bC, 0, x_mean_bC, h_mMMMMb_3Tag_CR24->GetMaximum()); l_mean_bC->SetLineColor(kRed); l_mean_bC->Draw();
  TLine *l_kHi_bC=new TLine(x_kHi_bC, 0, x_kHi_bC, h_mMMMMb_3Tag_CR24->GetMaximum()); l_kHi_bC->SetLineColor(kRed); l_kHi_bC->SetLineStyle(9); l_kHi_bC->Draw();
  TLine *l_kLo_bC=new TLine(x_kLo_bC, 0, x_kLo_bC, h_mMMMMb_3Tag_CR24->GetMaximum()); l_kLo_bC->SetLineColor(kRed); l_kLo_bC->SetLineStyle(9); l_kLo_bC->Draw();
  double x_mean_bS=bS_p0.getVal();
  double x_kHi_bS=bS_p0.getVal()+bS_p2.getVal()*bS_p1.getVal();
  double x_kLo_bS=bS_p4.getVal();
  TLine *l_mean_bS=new TLine(x_mean_bS, 0, x_mean_bS, h_mMMMMb_3Tag_SR->GetMaximum()); l_mean_bS->SetLineColor(kBlue); l_mean_bS->Draw();
  TLine *l_kHi_bS=new TLine(x_kHi_bS, 0, x_kHi_bS, h_mMMMMb_3Tag_SR->GetMaximum()); l_kHi_bS->SetLineColor(kBlue); l_kHi_bS->SetLineStyle(9); l_kHi_bS->Draw();
  TLine *l_kLo_bS=new TLine(x_kLo_bS, 0, x_kLo_bS, h_mMMMMb_3Tag_SR->GetMaximum()); l_kLo_bS->SetLineColor(kBlue); l_kLo_bS->SetLineStyle(9); l_kLo_bS->Draw();
  c_bC->SaveAs(("c_compareData_"+tags+"_VR_RooFit_GaussExp.png").c_str());
  
   return;
  
  /*// Calculate Pi and DPi and dPi -- for shape systematics
  double PbC_0=bC_p0.getVal();
  double PbC_1=bC_p1.getVal();
  double PbC_2=bC_p2.getVal();
  double PbC_3=bC_p3.getVal();
  double PbC_4=bC_p4.getVal();
  double dPbC_0=bC_p0.getError();
  double dPbC_1=bC_p1.getError();
  double dPbC_2=bC_p2.getError();
  double dPbC_3=bC_p3.getError();
  double dPbC_4=bC_p4.getError();
  double PbS_0=bS_p0.getVal();
  double PbS_1=bS_p1.getVal();
  double PbS_2=bS_p2.getVal();
  double PbS_3=bS_p3.getVal();
  double PbS_4=bS_p4.getVal();
  double dPbS_0=bS_p0.getError();
  double dPbS_1=bS_p1.getError();
  double dPbS_2=bS_p2.getError();
  double dPbS_3=bS_p3.getError();
  double dPbS_4=bS_p4.getError();
  
  std::cout<<" = MMMM Background Prediction ==== "<<std::endl;
  TH1F *h_mMMMMa_3Tag_CR2=(TH1F*)f_MMMM_a->Get("h_mX_CR2");
  TH1F *h_mMMMMa_3Tag_CR4=(TH1F*)f_MMMM_a->Get("h_mX_CR4");
  TH1F *h_mMMMMa_3Tag_SR;
  h_mMMMMa_3Tag_SR=(TH1F*)f_MMMM_a->Get("h_mX_SR");
  TH1F *h_mMMMMa_3Tag_CR24=(TH1F*)h_mMMMMa_3Tag_CR2->Clone("h_mX_CR24");
  h_mMMMMa_3Tag_CR24->Add(h_mMMMMa_3Tag_CR4);
  h_mMMMMa_3Tag_CR24->Rebin(rebin);
  h_mMMMMa_3Tag_CR24->SetLineColor(kBlack);
  h_mMMMMa_3Tag_SR->Rebin(rebin);
  h_mMMMMa_3Tag_SR->SetLineColor(kBlue);
  TH1F *h_mMMMMa_3Tag_SR_Prediction=(TH1F*)h_mMMMMa_3Tag_CR24->Clone("h_mMMMMa_3Tag_SR_Prediction");
  double aC=h_mMMMMa_3Tag_CR24->Integral(h_mMMMMa_3Tag_CR24->FindBin(VR_lo), h_mMMMMa_3Tag_CR24->FindBin(VR_hi));
  // Get the scale of the prediction right
  // std::cout<<"bS/bC = "<<bS/bC<<std::endl;
  double ratioAt125 = h_mMMMMa_3Tag_SR->GetSumOfWeights() / h_mMMMMa_3Tag_CR24->GetSumOfWeights();

  std::cout<<"ratioAt125 = "<<ratioAt125<<", +- "<<errorAt125<<" (fract unc.) = "<<1.+errorAt125/ratioAt125<<std::endl;
  std::cout<<"Number of aC events in 18.6 /fb within the mass window = "<<h_mMMMMa_3Tag_SR_Prediction->Integral(h_mMMMMa_3Tag_SR_Prediction->FindBin(SR_lo), h_mMMMMa_3Tag_SR_Prediction->FindBin(SR_hi))<<std::endl;
  std::cout<<"Number of predicted events in 18.6 /fb within the mass window = "<<h_mMMMMa_3Tag_SR_Prediction->Integral(h_mMMMMa_3Tag_SR_Prediction->FindBin(SR_lo), h_mMMMMa_3Tag_SR_Prediction->FindBin(SR_hi))*ratioAt125<<std::endl;
  RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);
  RooRealVar bg_p0("bg_p0", "bg_p0", 250., 500.);
  RooRealVar bg_p1("bg_p1", "bg_p1", 30., 100.);
  RooRealVar bg_p2("bg_p2", "bg_p2", 0.0, 10.1);
  RooRealVar bg_p3("bg_p3", "bg_p3", 0., 1.);
  RooRealVar bg_p4("bg_p4", "bg_p4", 250., 410.);
  LowMassBackgroundPDF bg("bg", "bg", x, bg_p0, bg_p1, bg_p2, bg_p3, bg_p4);
  RooDataHist pred("pred", "Prediction from SB", RooArgList(x), h_mMMMMa_3Tag_SR_Prediction);
  bg.fitTo(pred, RooFit::Range(SR_lo, SR_hi));
  RooPlot *aC_plot=x.frame();
  pred.plotOn(aC_plot, RooFit::LineColor(kRed), RooFit::MarkerColor(kRed));
  bg.plotOn(aC_plot, RooFit::LineColor(kRed));
  TCanvas *c_rooFit=new TCanvas("c_rooFit", "c_rooFit", 1000, 700);
  if (tags!="MMMM") h_mMMMMa_3Tag_SR->Draw("Ep9 SAME");
  aC_plot->Draw();
  double x_mean_aC=bg_p0.getVal();
  double x_kHi_aC=bg_p0.getVal()+bg_p2.getVal()*bg_p1.getVal();
  double x_kLo_aC=bg_p4.getVal();
  TLine *l_mean_aC=new TLine(x_mean_aC, 0, x_mean_aC, h_mMMMMa_3Tag_SR_Prediction->GetMaximum()); l_mean_aC->SetLineColor(kRed); l_mean_aC->Draw();
  TLine *l_kHi_aC=new TLine(x_kHi_aC, 0, x_kHi_aC, h_mMMMMb_3Tag_SR->GetMaximum()); l_kHi_aC->SetLineColor(kRed); l_kHi_aC->SetLineStyle(9); l_kHi_aC->Draw();
  TLine *l_kLo_aC=new TLine(x_kLo_aC, 0, x_kLo_aC, h_mMMMMb_3Tag_SR->GetMaximum()); l_kLo_aC->SetLineColor(kRed); l_kLo_aC->SetLineStyle(9); l_kLo_aC->Draw();
  
  std::cout<<"chi^2/ndof of aC = "<<aC_plot->chiSquare()<<std::endl;
  //c_rooFit->SaveAs(("c_compareData_"+tags+"_SR_RooFit_GaussExp.png").c_str());
  
  // Prediction Curve with Shape Systematics
  double PaC_0=bg_p0.getVal();
  double PaC_1=bg_p1.getVal();
  double PaC_2=bg_p2.getVal();
  double PaC_3=bg_p3.getVal();
  double PaC_4=bg_p4.getVal();
  double dPaC_0=bg_p0.getError();
  double dPaC_1=bg_p1.getError();
  double dPaC_2=bg_p2.getError();
  double dPaC_3=bg_p3.getError();
  double dPaC_4=bg_p4.getError();
  double PaS_0=PaC_0*PbS_0/PbC_0;
  double PaS_1=PaC_1*PbS_1/PbC_1;
  double PaS_2=PaC_2*PbS_2/PbC_2;
  double PaS_3=PaC_3*PbS_3/PbC_3;
  double PaS_4=PaC_4*PbS_4/PbC_4;
  double dPaS_0=PaS_0*quad((dPaC_0/PaC_0), (dPbS_0/PbS_0), (dPbC_0/PbC_0));
  double dPaS_1=PaS_1*quad((dPaC_1/PaC_1), (dPbS_1/PbS_1), (dPbC_1/PbC_1));
  double dPaS_2=PaS_2*quad((dPaC_2/PaC_2), (dPbS_2/PbS_2), (dPbC_2/PbC_2));
  double dPaS_3=PaS_3*quad((dPaC_3/PaC_3), (dPbS_3/PbS_3), (dPbC_3/PbC_3));
  double dPaS_4=PaS_4*quad((dPaC_4/PaC_4), (dPbS_4/PbS_4), (dPbC_4/PbC_4));
  std::cout<<"(dPaC_0/PaC_0) = ("<<dPaC_0<<"/"<<PaC_0<<") = "<<(dPaC_0/PaC_0)<<"; (dPbS_0/PbS_0) = ("<<dPbS_0<<"/"<<PbS_0<<") = "<<(dPbS_0/PbS_0)<<"; (dPbC_0/PbC_0) = ("<<dPbC_0<<"/"<<PbC_0<<") = "<<(dPbC_0/PbC_0)<<std::endl;
  std::cout<<"(dPaC_1/PaC_1) = ("<<dPaC_1<<"/"<<PaC_1<<") = "<<(dPaC_1/PaC_1)<<"; (dPbS_1/PbS_1) = ("<<dPbS_1<<"/"<<PbS_1<<") = "<<(dPbS_1/PbS_1)<<"; (dPbC_1/PbC_1) = ("<<dPbC_1<<"/"<<PbC_1<<") = "<<(dPbC_1/PbC_1)<<std::endl; 
  std::cout<<"(dPaC_2/PaC_2) = ("<<dPaC_2<<"/"<<PaC_2<<") = "<<(dPaC_2/PaC_2)<<"; (dPbS_2/PbS_2) = ("<<dPbS_2<<"/"<<PbS_2<<") = "<<(dPbS_2/PbS_2)<<"; (dPbC_2/PbC_2) = ("<<dPbC_2<<"/"<<PbC_2<<") = "<<(dPbC_2/PbC_2)<<std::endl; 
  std::cout<<"(dPaC_3/PaC_3) = ("<<dPaC_3<<"/"<<PaC_3<<") = "<<(dPaC_3/PaC_3)<<"; (dPbS_3/PbS_3) = ("<<dPbS_3<<"/"<<PbS_3<<") = "<<(dPbS_3/PbS_3)<<"; (dPbC_3/PbC_3) = ("<<dPbC_3<<"/"<<PbC_3<<") = "<<(dPbC_3/PbC_3)<<std::endl;
  std::cout<<"(dPaC_4/PaC_4) = ("<<dPaC_4<<"/"<<PaC_4<<") = "<<(dPaC_4/PaC_4)<<"; (dPbS_4/PbS_4) = ("<<dPbS_4<<"/"<<PbS_4<<") = "<<(dPbS_4/PbS_4)<<"; (dPbC_4/PbC_4) = ("<<dPbC_4<<"/"<<PbC_4<<") = "<<(dPbC_4/PbC_4)<<std::endl;
  std::cout<<" Predicted PaS_0 = "<<PaS_0<<" +- "<<dPaS_0<<std::endl;
  std::cout<<" Predicted PaS_1 = "<<PaS_1<<" +- "<<dPaS_1<<std::endl;
  std::cout<<" Predicted PaS_2 = "<<PaS_2<<" +- "<<dPaS_2<<std::endl;
  std::cout<<" Predicted PaS_3 = "<<PaS_3<<" +- "<<dPaS_3<<std::endl;
  std::cout<<" Predicted PaS_4 = "<<PaS_4<<" +- "<<dPaS_4<<std::endl;
  RooRealVar bg_pred0;
  RooRealVar bg_pred1;
  RooRealVar bg_pred2;
  RooRealVar bg_pred3;
  RooRealVar bg_pred4;
   //RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);
   LowMassBackgroundPDF bg_pred("bg_pred", "Background Predicted for Signal Region", x, bg_pred0, bg_pred1, bg_pred2, bg_pred3, bg_pred4);

  if (tags!="MMMM")
  {
    bg_pred0=new RooRealVar("bg_pred0", "bg_pred0", PaS_0-dPaS_0/2., PaS_0+dPaS_0/2.);
    bg_pred1=new RooRealVar("bg_pred1", "bg_pred1", PaS_1-dPaS_1/2., PaS_1+dPaS_1/2.);
    bg_pred2=new RooRealVar("bg_pred2", "bg_pred2", PaS_2-dPaS_2/2., PaS_2+dPaS_2/2.);
    bg_pred3=new RooRealVar("bg_pred3", "bg_pred3", PaS_3-dPaS_3/2., PaS_3+dPaS_3/2.);
    bg_pred4=new RooRealVar("bg_pred4", "bg_pred4", PaS_4-dPaS_4/2., PaS_4+dPaS_4/2.);
  }
  else
  {
    bg_pred0=new RooRealVar("bg_pred0", "bg_pred0", PaS_0);  bg_pred0.setError(dPaS_0);
    bg_pred1=new RooRealVar("bg_pred1", "bg_pred1", PaS_1);  bg_pred1.setError(dPaS_1);
    bg_pred2=new RooRealVar("bg_pred2", "bg_pred2", PaS_2);  bg_pred2.setError(dPaS_2);
    bg_pred3=new RooRealVar("bg_pred3", "bg_pred3", PaS_3);  bg_pred3.setError(dPaS_3);
    bg_pred4=new RooRealVar("bg_pred4", "bg_pred4", PaS_4);  bg_pred4.setError(dPaS_4);
  }
  //RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);
 // LowMassBackgroundPDF bg_pred("bg_pred", "Background Predicted for Signal Region", x, bg_pred0, bg_pred1, bg_pred2, bg_pred3, bg_pred4);
  RooPlot *aS_plot=x.frame();
  if (tags!="MMMM")
  {
    RooDataHist unblind("unblind", "Signal Region", RooArgList(x), h_mMMMMa_3Tag_SR);
    unblind.plotOn(aS_plot, RooFit::LineColor(kBlue), RooFit::MarkerColor(kBlue));
    bg_pred.fitTo(unblind, RooFit::Range(SR_lo, SR_hi));
    bg_pred.plotOn(aS_plot, RooFit::LineColor(kBlue));
  //   std::cout<<" asd "<<std::endl;
    aS_plot->Draw("same");
  }
  else
  {
    bg_pred.plotOn(aC_plot, RooFit::LineColor(kBlue), RooFit::Range(SR_lo, SR_hi));
    std::cout<<" asd "<<std::endl;	
    aC_plot->Draw("same");
     std::cout<<" asd 1"<<std::endl; 
   //generating pseudo-data

   int num_toys = 0;
   

  TFile *data_8TeVData2012=new TFile("MMMM/a/Histograms_BJetPlusX_Run2012BCD_Skim.root");
  TH1F *h_mX_CR2=(TH1F*)data_8TeVData2012->Get("h_mX_CR2");
  TH1F *h_mX_CR4=(TH1F*)data_8TeVData2012->Get("h_mX_CR4");
  TH1F *h_mX_CR24=(TH1F*)h_mX_CR2->Clone("h_mX_CR24");
  h_mX_CR24->Add(h_mX_CR4);
  TH1F *h_mX_SR=(TH1F*)data_8TeVData2012->Get("h_mX_SR");
// ***** fit fake data with same function ***** /

  TH1F *h_mX_SR_fakeData=(TH1F*)h_mX_CR24->Clone("h_mX_SR_fakeData");
  h_mX_SR_fakeData->Reset();
  h_mX_SR_fakeData->SetLineColor(kBlack);
  std::cout<<"Background CR->SR Scale, 277.542/h_mX_CR24->GetSumOfWeights() = "<<h_mX_SR->GetSumOfWeights()/h_mX_CR24->GetSumOfWeights()<<"  num " << h_mX_SR->GetSumOfWeights()<<"  den "<<h_mX_CR24->GetSumOfWeights()<< std::endl;
  // fillUnweightedTH1(h_mX_CR24, h_mX_SR_fakeData, (277.542/h_mX_CR24->GetSumOfWeights()));

  fillByThrowTH1(h_mX_CR24, h_mX_SR_fakeData, (h_mX_SR->GetSumOfWeights()/h_mX_CR24->GetSumOfWeights()));



  std::string filename="MMMM_nominal/a_KinFit/Histograms_DiJetPt_glugluToX400ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root";
  int rebin=1;
  std::string filename_old = "/gpfs/ddn/cms/user/cvernier/H4b_step2/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/DiJetPt_glugluToX400ToHHTobbbb_8TeV_width1MeV_GEN_FASTSIM_HLT.root";
  TFile *signal=new TFile(filename.c_str());
  TH1F *h_mX_SR_signal=(TH1F*)signal->Get("h_mX_SR"); 
 
  TFile *signal_old=new TFile(filename_old.c_str());
  // Calculate nSignal events given production cross section, branching fractions and efficiency
  TH1F *h_CountWithPU=(TH1F*)signal_old->Get("CountWithPU");
  double nSignal_init=h_CountWithPU->GetBinContent(1);
  double nSignal_now=h_mX_SR_signal->GetSumOfWeights();
  double eff_signal=nSignal_now/nSignal_init;
  std::cout<<"Signal efficiency = "<<eff_signal<<std::endl;
  // double totalLumi=13241.968; // /pb
  double totalLumi=18600.0; // /pb
  double prodXsec_1=1.0; // pb
//h_mX_SR_signal->Draw();
  // fillUnweightedTH1(h_mX_SR_signal, h_mX_SR_fakeData, prodXsec_1*totalLumi/nSignal_init);
  //fillByThrowTH1(h_mX_SR_signal, h_mX_SR_fakeData, prodXsec_1*totalLumi/nSignal_init);
//h_mX_SR_fakeData->Draw();
  RooDataHist unblind("unblind", "Signal Region", RooArgList(x), h_mX_SR_fakeData);
  unblind.plotOn(aS_plot, RooFit::LineColor(kGreen), RooFit::MarkerColor(kGreen));
  bg.fitTo(unblind, RooFit::Range(SR_lo, SR_hi));
  bg.plotOn(aS_plot, RooFit::LineColor(kGreen));
  std::cout<<" ehi "<<std::endl;
  aS_plot->Draw("same");




  num_toys = h_mX_SR_fakeData->Integral(h_mX_SR_fakeData->FindBin(SR_lo), h_mX_SR_fakeData->FindBin(SR_hi)-1); //GetSumOfWeights();///h_mX_CR24->GetSumOfWeights();

  //RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);
  //LowMassBackgroundPDF bg_pred3("bg_pred3", "Background Predicted for Signal Region", x, bg_pred0, bg_pred1, bg_pred2, bg_pred3, bg_pred4);
  h_mX_SR_signal-> Scale(prodXsec_1*totalLumi/nSignal_init); 
  std::cout<<h_mX_SR_signal->GetIntegral()<<"  asd  "<<std::endl;
   RooDataHist data_sig("data_sig", "Data", RooArgList(x), h_mX_SR_signal);
   


   RooDataHist* data_ = bg.generateBinned(x,num_toys, RooAbsData::Poisson);
   data_->add(data_sig);
   RooDataHist data_obs(*data_, "data_obs");
   std::cout<<" asd 2 "<<std::endl;
   TH1 * h_fake = data_obs.createHistogram("h_fake",x);
   //h_fake->Draw();
  TH1F *h_fakeData=(TH1F*)h_fake->Clone("h_fakeData");
  h_fakeData->Reset();
    
  fillByThrowTH1(h_fake, h_fakeData, 1);


  RooDataHist obs("obs", "Data", RooArgList(x), h_fakeData);
  RooPlot *plot=x.frame();
  obs.plotOn(plot);
  TCanvas *c_data=new TCanvas("c_data", "c_data", 500, 500);
  plot->Draw();
  RooWorkspace *w1=new RooWorkspace("HbbHbb");
  w1->import(obs);
  std::string point="400";
  std::cout<<"mass = "<<point<<std::endl;
  w1->SaveAs(("w_data_"+point+".root").c_str());
 



  }
  double x_mean_aS=bg_pred0.getVal();
  double x_k_aS=bg_pred0.getVal()+bg_pred2.getVal()*bg_pred1.getVal();
  TLine *l_mean_aS=new TLine(x_mean_aS, 0, x_mean_aS, h_mMMMMa_3Tag_SR_Prediction->GetMaximum()); l_mean_aS->SetLineColor(kBlue); l_mean_aS->Draw();
  TLine *l_k_aS=new TLine(x_k_aS, 0, x_k_aS, h_mMMMMa_3Tag_SR_Prediction->GetMaximum()); l_k_aS->SetLineColor(kBlue); l_k_aS->SetLineStyle(9); l_k_aS->Draw();
  
  std::cout<<" === === "<<std::endl;
  std::cout<<"chi^2/ndof of bC = "<<bC_plot->chiSquare()<<std::endl;
  std::cout<<"chi^2/ndof of bS = "<<bS_plot->chiSquare()<<std::endl;
  std::cout<<"chi^2/ndof of aC = "<<aC_plot->chiSquare()<<std::endl;
 // std::cout<<"chi^2/ndof of aS = "<<aS_plot->chiSquare()<<std::endl;
  std::cout<<" === === "<<std::endl;
   
  c_rooFit->SaveAs(("c_compareData_"+tags+"_SR_RooFit_GaussExp.png").c_str());

  RooWorkspace *w=new RooWorkspace("HbbHbb");
  //std::cout<<"  asd "<<std::endl;
  //bg_pred.fitTo(pred, RooFit::Range(SR_lo, SR_hi));
  w->import(bg);
  w->SaveAs("w_background_LowMass.root"); 
  */   
}




  
  
  

