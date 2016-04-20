
#include <iostream>
#include<iomanip>
#include <fstream>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"

#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "Input/binning_util.h"
#include "../mBinning_HT.h"
using namespace std;

const Int_t nTrg=2;
const Double_t deta=1.4;

void Draw_Efficiency();
void Draw_Pt_spectra(TH1F *a1[],TH1F *a2[],TH1F *a3[],TH1F *a4[],TH1F *a5[],TH1F *a6[],TH1F *a7[],TH1F *a8[]);
void Efficiency_Correction_sts(TH1F *HT0_Raw_NPE_sts,TH1F *HT2_Raw_NPE_sts,TH1F *&HT0_NPE_sts,TH1F *&HT2_NPE_sts);
void Efficiency_Correction_sys(TH1F *HT0_Raw_NPE_sys,TH1F *HT2_Raw_NPE_sys,TH1F *&HT0_NPE_sys,TH1F *&HT2_NPE_sys);

void Cross_Section_sts(TH1F *HT0_NPE_sts,TH1F *HT2_NPE_sts,TH1F *&HT0_NPE_spectra_sts,TH1F *&HT2_NPE_spectra_sts,Double_t ,Double_t );
void Cross_Section_sys(TH1F *HT0_NPE_sys,TH1F *HT2_NPE_sys,TH1F *&HT0_NPE_spectra_sys,TH1F *&HT2_NPE_spectra_sys,Double_t ,Double_t );


void Draw_Cross_Section(TH1F *HT0_NPE_sts,TH1F *HT0_NPE_sys,TH1F *HT2_NPE_sts,TH1F *HT2_NPE_sys,TH1F *HT_NPE_sts,TH1F *HT_NPE_sys,TH1F *&,TH1F *&);


void Draw_NPE_PHE_ratio(TH1F *,TH1F *,TH1F *,TH1F *,TH1F *&,TH1F *&,TH1F *&,TH1F *&);  
void Raw_NPE_spectra_sts(TH1F *,TH1F *,TH1F *,TH1F *,TH1F *&,TH1F *&);  
void Raw_NPE_spectra_sys(TH1F *,TH1F *,TH1F *,TH1F *,TH1F *&,TH1F *&);  


TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
void setpad(TPad *,float , float , float , float );

//TH1F * Correc=new TH1F("Correc","",Nbins_HT,Pt_bin_HT);

TFile *infile_data =new TFile("Input/hist_4_13.root","read");
TFile *infile_BEMC =new TFile("Input/BEMC_efficiency.root","read");
TFile *infile_dEdx_cut =new TFile("Input/nSigma_Cut_efficiency.root","read");
TFile *infile_Trigger =new TFile("Input/TrigEfficiency_HT.root","read");
TFile *infile_Tracking =new TFile("Input/Tracking_efficiency_HT.root","read");
TFile *infile_PHE_re =new TFile("Input/Photonic_re_Efficiency.root","read");
TFile *infile_purity =new TFile("Input/purity_HT.root","read");

TH1F *eff_BEMC=(TH1F *) infile_BEMC->Get("efficiency_bemc");
TH1F *eff_dedx=(TH1F *) infile_dEdx_cut->Get("nsigmaE_HT");
TH1F *eff_Trigger_HT0=(TH1F *) infile_Trigger->Get("TrigEfficiency_HT0");
TH1F *eff_Trigger_HT2=(TH1F *) infile_Trigger->Get("TrigEfficiency_HT2");
TH1F *eff_Tracking=(TH1F *) infile_Tracking->Get("Tracking_efficiency_HT");
TH1F *eff_PHE_re=(TH1F *) infile_PHE_re->Get("PHE_re_efficiency");
TH1F *eff_purity=(TH1F *) infile_purity->Get("purity_HT");

TH1F *mh1MB_Nevents_psTrg0=(TH1F *) infile_data->Get("mh1MB_Nevents_psTrg0");
TH1F *mh1MB_Nevents_psTrg1=(TH1F *) infile_data->Get("mh1MB_Nevents_psTrg1");

TFile *File_fonll=new TFile("Input/fonll200gev.root","READ");
TFile *File_Jpsi=new TFile("Input/Jpsi.root","READ");

void NPE_Cross_section()
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  char buf[1024];
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1,"Y");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13); 

  Double_t  HT0_NmbEvents=1.95584425792000000e+11;// mh1MB_Nevents_psTrg0->GetBinContent(2);
  Double_t  HT2_NmbEvents=1.99543894176000000e+11;// mh1MB_Nevents_psTrg1->GetBinContent(2);

  // Double_t  HT0_NmbEvents= mh1MB_Nevents_psTrg0->GetBinContent(2);
  // Double_t  HT2_NmbEvents= mh1MB_Nevents_psTrg1->GetBinContent(2);


  TH1F *Inclusive[nTrg];
  TH1F *Inclusive_ps[nTrg];

  TH2F *Photonic_unlike[nTrg];
  TH2F *Photonic_like[nTrg];
  TH2F *Photonic_unlike_like[nTrg];
  
  TH2F *Photonic_ps_unlike[nTrg];
  TH2F *Photonic_ps_like[nTrg];
  TH2F *Photonic_ps_unlike_like[nTrg];
  

  TH1F *Photonic_unlike_pt[nTrg];
  TH1F *Photonic_like_pt[nTrg];
  TH1F *Photonic_unlike_like_pt[nTrg];

  TH1F *Photonic_unlike_pt_ps[nTrg];
  TH1F *Photonic_like_pt_ps[nTrg];
  TH1F *Photonic_unlike_like_pt_ps[nTrg];



  
  
  TString Photonic_unlikename[nTrg]={"mh2InvMassUnlikeTrg0","mh2InvMassUnlikeTrg1"};
  TString Photonic_likename[nTrg]={"mh2InvMasslikeTrg0","mh2InvMasslikeTrg1"};
  TString Photonic_unlikename_ps[nTrg]={"mh2InvMassUnlike_psTrg0","mh2InvMassUnlike_psTrg1"};
  TString Photonic_likename_ps[nTrg]={"mh2InvMasslike_psTrg0","mh2InvMasslike_psTrg1"};

  TString Inclusive_name[nTrg]={"mh1electronPtTrg0","mh1electronPtTrg1"};
  TString Inclusive_name_ps[nTrg]={"mh1electronPt_psTrg0","mh1electronPt_psTrg1"};
  

  for(Int_t iTrg=0;iTrg<nTrg;iTrg++)
    {
      
      // incusve 
      Inclusive[iTrg]=(TH1F *) infile_data->Get(Inclusive_name[iTrg]);
      Inclusive_ps[iTrg]=(TH1F *) infile_data->Get(Inclusive_name_ps[iTrg]);


      //photonic 

      Photonic_unlike[iTrg]=(TH2F *) infile_data->Get(Photonic_unlikename[iTrg]);
      Photonic_like[iTrg]=(TH2F *) infile_data->Get(Photonic_likename[iTrg]);
      sprintf(buf,"unlike_%i",iTrg);
      Photonic_unlike_pt[iTrg]=(TH1F *) Photonic_unlike[iTrg]->ProjectionY(buf);
      sprintf(buf,"like_%i",iTrg);
      Photonic_like_pt[iTrg]=(TH1F *) Photonic_like[iTrg]->ProjectionY(buf);

      sprintf(buf,"unlike_like%i",iTrg);
      Photonic_unlike_like_pt[iTrg]=(TH1F *) Photonic_unlike_pt[iTrg]->Clone(buf); 
      Photonic_unlike_like_pt[iTrg]->Sumw2();
      Photonic_unlike_like_pt[iTrg]->Add(Photonic_like_pt[iTrg],-1);

           
      Photonic_ps_unlike[iTrg]=(TH2F *) infile_data->Get(Photonic_unlikename_ps[iTrg]);
      Photonic_ps_like[iTrg]=(TH2F *) infile_data->Get(Photonic_likename_ps[iTrg]);

      sprintf(buf,"unlike_ps%i",iTrg);
      Photonic_unlike_pt_ps[iTrg]=(TH1F *) Photonic_ps_unlike[iTrg]->ProjectionY(buf);
      sprintf(buf,"like_ps%i",iTrg);
      Photonic_like_pt_ps[iTrg]=(TH1F *) Photonic_ps_like[iTrg]->ProjectionY(buf);

      sprintf(buf,"unlike_like_ps%i",iTrg);
      Photonic_unlike_like_pt_ps[iTrg]=(TH1F *) Photonic_unlike_pt_ps[iTrg]->Clone(buf); 
      Photonic_unlike_like_pt_ps[iTrg]->Sumw2();
      Photonic_unlike_like_pt_ps[iTrg]->Add(Photonic_like_pt_ps[iTrg],-1);

    }
  

  Draw_Efficiency();
  
  Draw_Pt_spectra(Inclusive,Photonic_unlike_pt,Photonic_like_pt,Photonic_unlike_like_pt,Inclusive_ps,Photonic_unlike_pt_ps,Photonic_like_pt_ps,Photonic_unlike_like_pt_ps);
  
  TH1F *HT0_inclusive=(TH1F *) Inclusive_ps[0]->Rebin(NpT_bins_run12_HT,"HT0_inclusive",pt_run12_HT);
  TH1F *HT2_inclusive=(TH1F *) Inclusive_ps[1]->Rebin(NpT_bins_run12_HT,"HT2_inclusive",pt_run12_HT);

  TH1F *HT0_phe=(TH1F *) Photonic_unlike_like_pt_ps[0]->Rebin(NpT_bins_run12_HT,"HT0_phe",pt_run12_HT);
  TH1F *HT2_phe=(TH1F *) Photonic_unlike_like_pt_ps[1]->Rebin(NpT_bins_run12_HT,"HT2_phe",pt_run12_HT);
  
  TH1F *HT0_Raw_NPE_sts;
  TH1F *HT0_Raw_NPE_sys;
  TH1F *HT2_Raw_NPE_sts;
  TH1F *HT2_Raw_NPE_sys;

  // Get Raw NPE yeild 
   Raw_NPE_spectra_sts(HT0_inclusive,HT2_inclusive,HT0_phe,HT2_phe,HT0_Raw_NPE_sts,HT2_Raw_NPE_sts);
   Raw_NPE_spectra_sys(HT0_inclusive,HT2_inclusive,HT0_phe,HT2_phe,HT0_Raw_NPE_sys,HT2_Raw_NPE_sys);

   TH1F *HT0_NPE_sts;
   TH1F *HT0_NPE_sys;
   TH1F *HT2_NPE_sts;
   TH1F *HT2_NPE_sys;
   
   // Get  NPE yeild corrected by the efficiency 
   Efficiency_Correction_sts(HT0_Raw_NPE_sts,HT2_Raw_NPE_sts,HT0_NPE_sts,HT2_NPE_sts);
   Efficiency_Correction_sys(HT0_Raw_NPE_sys,HT2_Raw_NPE_sys,HT0_NPE_sys,HT2_NPE_sys);
   
  TH1F *HT0_NPE_spectra_sts;
  TH1F *HT0_NPE_spectra_sys;
  
  TH1F *HT2_NPE_spectra_sts;
  TH1F *HT2_NPE_spectra_sys;

  //   Caculate the cross section
  Cross_Section_sts(HT0_NPE_sts,HT2_NPE_sts,HT0_NPE_spectra_sts,HT2_NPE_spectra_sts,HT0_NmbEvents,HT2_NmbEvents);
  Cross_Section_sys(HT0_NPE_sys,HT2_NPE_sys,HT0_NPE_spectra_sys,HT2_NPE_spectra_sys,HT0_NmbEvents,HT2_NmbEvents);
  
  
  TH1F *run12_HT_NPE_sts=new TH1F("run12_HT_NPE_sts","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *run12_HT_NPE_sys=new TH1F("run12_HT_NPE_sys","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *run12_HT_NPE_FONLL_sts=new TH1F("run12_HT_NPE_FONLL_sts","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *run12_HT_NPE_FONLL_sys=new TH1F("run12_HT_NPE_FONLL_sys","",NpT_bins_run12_HT,pt_run12_HT);

  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      if(ipt<3)
	{
      run12_HT_NPE_sts->SetBinContent(ipt+1,HT0_NPE_spectra_sts->GetBinContent(ipt+1));
      run12_HT_NPE_sts->SetBinError(ipt+1,HT0_NPE_spectra_sts->GetBinError(ipt+1));

      run12_HT_NPE_sys->SetBinContent(ipt+1,HT0_NPE_spectra_sys->GetBinContent(ipt+1));
      run12_HT_NPE_sys->SetBinError(ipt+1,HT0_NPE_spectra_sys->GetBinError(ipt+1));
      	}
      else
	{
	  run12_HT_NPE_sts->SetBinContent(ipt+1,HT2_NPE_spectra_sts->GetBinContent(ipt+1));
	  run12_HT_NPE_sts->SetBinError(ipt+1,HT2_NPE_spectra_sts->GetBinError(ipt+1));

	  run12_HT_NPE_sys->SetBinContent(ipt+1,HT2_NPE_spectra_sys->GetBinContent(ipt+1));
	  run12_HT_NPE_sys->SetBinError(ipt+1,HT2_NPE_spectra_sys->GetBinError(ipt+1));
	}
    }

  // run12_HT_NPE_sts->Draw();
  // run12_HT_NPE_sys->SetLineColor(2);
  // run12_HT_NPE_sys->Draw("same");
  
  Draw_Cross_Section(HT0_NPE_spectra_sts,HT2_NPE_spectra_sys,HT2_NPE_spectra_sts,HT2_NPE_spectra_sys,run12_HT_NPE_sts,run12_HT_NPE_sys,run12_HT_NPE_FONLL_sts,run12_HT_NPE_FONLL_sys);
}
void Draw_Cross_Section(TH1F *HT0_NPE_sts,TH1F *HT0_NPE_sys,TH1F *HT2_NPE_sts,TH1F *HT2_NPE_sys,TH1F *HT_NPE_sts,TH1F *HT_NPE_sys,TH1F *&run12_HT_NPE_FONLL_sts,TH1F *&run12_HT_NPE_FONLL_sys)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *FONLL=new TH1F("FONLL","",NpT_bins_run12_HT,pt_run12_HT);


  TGraphErrors *gFONLLu=( TGraphErrors *) File_fonll->Get("gFONLLu");
  TGraphErrors *gFONLLc=( TGraphErrors *) File_fonll->Get("gFONLLc");
  TGraphErrors *gFONLLl=( TGraphErrors *) File_fonll->Get("gFONLLl");

  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      FONLL->SetBinContent(ipt+1,gFONLLc->Eval(FONLL->GetBinCenter(ipt+1)));
      FONLL->SetBinError(ipt+1,0);
    }
  
  TH1F *Jpsi_cross_section=(TH1F *) File_Jpsi->Get("JPsi");  

  TGraphErrors *gR_FONLLu= (TGraphErrors *) File_fonll->Get("gR_FONLLu");
  TGraphErrors *gR_FONLLc= (TGraphErrors *) File_fonll->Get("gR_FONLLc");
  TGraphErrors *gR_FONLLl= (TGraphErrors *) File_fonll->Get("gR_FONLLl");
  

  TCanvas *c6=new TCanvas("c6","",800,600);
  gPad->SetLogy();

  TH2F *h6=new TH2F("h6","",100,2,10,100,1e-13,1e-2);
  h6->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h6->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");

  h6->Draw();
  // HT_NPE_sys->Scale(0.4);

  gFONLLc->Draw("same");
  gFONLLu->Draw("same");
  gFONLLl->Draw("same");

  TF1* f1 = new TF1("f1","[0]/(pow(x,[1])+[2])",2,10);
  f1->SetParameter(0,10);
  f1->SetParameter(1,8);
  f1->SetParameter(2,0.7); 
  HT_NPE_sts->Fit("f1","R","same" ,2,10);
  HT_NPE_sts->Fit("f1","R","same" ,2,10);

 
  const int nbins = NpT_bins_run12_HT+1;
  


  
  TH1D* hCorrection = get_bin_shift((unsigned long int)1e9,f1,nbins,pt_run12_HT);
  HT_NPE_sts->Divide(hCorrection);
  HT_NPE_sys->Divide(hCorrection);

  HT_NPE_sts->Add(Jpsi_cross_section,-1);
  HT_NPE_sys->Add(Jpsi_cross_section,-1);


  TGraphErrors  *gr_run12_HT_NPE_sts=new TGraphErrors(HT_NPE_sts);
  TGraphErrors  *gr_run12_HT_NPE_sys=new TGraphErrors(HT_NPE_sys);

  run12_HT_NPE_FONLL_sts->Divide(HT_NPE_sts,FONLL,1,1);
  run12_HT_NPE_FONLL_sys->Divide(HT_NPE_sys,FONLL,1,1);
  
  TGraphErrors  *gr_run12_HT_NPE_FONLL_sts=new TGraphErrors(run12_HT_NPE_FONLL_sts);
  TGraphErrors  *gr_run12_HT_NPE_FONLL_sys=new TGraphErrors(run12_HT_NPE_FONLL_sys);
  
  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      gr_run12_HT_NPE_sts->SetPointError(ipt,0,gr_run12_HT_NPE_sts->GetErrorY(ipt));
      gr_run12_HT_NPE_sys->SetPointError(ipt,0.05,gr_run12_HT_NPE_sys->GetErrorY(ipt));
      
      gr_run12_HT_NPE_FONLL_sts->SetPointError(ipt,0,gr_run12_HT_NPE_FONLL_sts->GetErrorY(ipt));
      gr_run12_HT_NPE_FONLL_sys->SetPointError(ipt,0.05,gr_run12_HT_NPE_FONLL_sys->GetErrorY(ipt));
    }  


  gr_run12_HT_NPE_sts->SetMarkerStyle(20);
  gr_run12_HT_NPE_sys->SetMarkerColor(2);
  gr_run12_HT_NPE_sts->SetLineColor(2);
  gr_run12_HT_NPE_sys->SetLineColor(2);


  gr_run12_HT_NPE_sys->Draw("sameE1p[]");
  gr_run12_HT_NPE_sts->Draw("samePE");
  c6->SaveAs("Cross_section.pdf");
  
  TCanvas *c7=new TCanvas("c7","",800,600);
  //  gPad->SetLogy();
  TH2F *h7=new TH2F("h7","",100,2,10,100,0,3);
  h7->Draw();
  h7->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h7->GetYaxis()->SetTitle("Data/FONLL");
  

  //run12_HT_NPE_FONLL_sys->Draw("same");

  gr_run12_HT_NPE_FONLL_sts->SetMarkerStyle(20);
  gr_run12_HT_NPE_FONLL_sts->SetMarkerColor(2);
  gr_run12_HT_NPE_FONLL_sys->SetMarkerColor(2);

  gr_run12_HT_NPE_FONLL_sts->SetLineColor(2);
  gr_run12_HT_NPE_FONLL_sys->SetLineColor(2);




  gr_run12_HT_NPE_FONLL_sys->Draw("samepe1[]");
  gr_run12_HT_NPE_FONLL_sts->Draw("samePE");
  
  gR_FONLLc->Draw("same");
  gR_FONLLu->Draw("same");
  gR_FONLLl->Draw("same");
 c7->SaveAs("Cross_section_fonll.pdf");



 TFile *file=new TFile("run12_Npe_HT.root","RECREATE");

 HT_NPE_sts->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 HT_NPE_sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 run12_HT_NPE_FONLL_sts->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 run12_HT_NPE_FONLL_sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");

 HT_NPE_sts->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");
 HT_NPE_sys->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");
 run12_HT_NPE_FONLL_sts->GetYaxis()->SetTitle("Data/FONLL");
 run12_HT_NPE_FONLL_sys->GetYaxis()->SetTitle("Data/FONLL");

 HT_NPE_sts->Write();
 HT_NPE_sys->Write();
 run12_HT_NPE_FONLL_sts->Write();
 run12_HT_NPE_FONLL_sys->Write();

 file->Close();

}
void Cross_Section_sts(TH1F *HT0_NPE_sts,TH1F *HT2_NPE_sts,TH1F *&HT0_NPE_spectra_sts,TH1F *&HT2_NPE_spectra_sts,Double_t HT0_NmbEvents,Double_t HT2_NmbEvents)
{

  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *HT0_Npe_sts=(TH1F *) HT0_NPE_sts->Clone("HT0_Npe_sts");
  TH1F *HT2_Npe_sts=(TH1F *) HT2_NPE_sts->Clone("HT2_Npe_sts");
  //2Pi
  HT0_Npe_sts->Scale(1./(2*TMath::Pi()));
  HT2_Npe_sts->Scale(1./(2*TMath::Pi()));
  // delta eta
  HT0_Npe_sts->Scale(1./deta);
  HT2_Npe_sts->Scale(1./deta);

  // delta pT
  HT0_Npe_sts->Scale(1.,"width");
  HT2_Npe_sts->Scale(1.,"width");

  //pT
  TH1F *HT0_Norm_sts=new TH1F("HT0_Norm_sts","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *HT2_Norm_sts=new TH1F("HT2_Norm_sts","",NpT_bins_run12_HT,pt_run12_HT);

  HT0_Npe_sts->Scale(1./HT0_NmbEvents);
  HT2_Npe_sts->Scale(1./HT2_NmbEvents);

  //Pt*NSD*BBCTrigger_efficiency*2(eplus and eminus)
  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      HT0_Norm_sts->SetBinContent(ipt+1,HT0_Npe_sts->GetBinCenter(ipt+1)*2/30.);
      HT2_Norm_sts->SetBinContent(ipt+1,HT2_Npe_sts->GetBinCenter(ipt+1)*2/30.);
      
      HT0_Norm_sts->SetBinError(ipt+1,0);
      HT2_Norm_sts->SetBinError(ipt+1,0);
    }
  
  HT0_Npe_sts->Divide(HT0_Norm_sts);
  HT2_Npe_sts->Divide(HT2_Norm_sts);

  HT0_NPE_spectra_sts=(TH1F *) HT0_Npe_sts->Clone("HT0_NPE_spectra_sts");
  HT2_NPE_spectra_sts=(TH1F *) HT2_Npe_sts->Clone("HT2_NPE_spectra_sts");
  

  // HT0_Npe_sts->Draw();
  // HT2_Npe_sts->SetMarkerColor(2);
  // HT2_Npe_sts->Draw("same");

}

void Cross_Section_sys(TH1F *HT0_NPE_sys,TH1F *HT2_NPE_sys,TH1F *&HT0_NPE_spectra_sys,TH1F *&HT2_NPE_spectra_sys,Double_t HT0_NmbEvents,Double_t HT2_NmbEvents)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *HT0_Npe_sys=(TH1F *) HT0_NPE_sys->Clone("HT0_Npe_sys");
  TH1F *HT2_Npe_sys=(TH1F *) HT2_NPE_sys->Clone("HT2_Npe_sys");
  //2Pi
  HT0_Npe_sys->Scale(1./(2*TMath::Pi()));
  HT2_Npe_sys->Scale(1./(2*TMath::Pi()));
  // delta eta
  HT0_Npe_sys->Scale(1./deta);
  HT2_Npe_sys->Scale(1./deta);

  // delta pT
  HT0_Npe_sys->Scale(1.,"width");
  HT2_Npe_sys->Scale(1.,"width");

  HT0_Npe_sys->Scale(1./HT0_NmbEvents);
  HT2_Npe_sys->Scale(1./HT2_NmbEvents);

  //pT
  TH1F *HT0_pt_sys=new TH1F("HT0_pt_sys","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *HT2_pt_sys=new TH1F("HT2_pt_sys","",NpT_bins_run12_HT,pt_run12_HT);

  TH1F *NSD_sys=new TH1F("NSD_sys","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *BBC_trigger=new TH1F("BBC_Trigger","",NpT_bins_run12_HT,pt_run12_HT);


  //Pt*2(eplus and eminus)

  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      HT0_pt_sys->SetBinContent(ipt+1,HT0_Npe_sys->GetBinCenter(ipt+1)*2);  //*0.866*2/30.);
      HT2_pt_sys->SetBinContent(ipt+1,HT2_Npe_sys->GetBinCenter(ipt+1)*2) ; //*0.866*2/30.);
      
      HT0_pt_sys->SetBinError(ipt+1,0);
      HT2_pt_sys->SetBinError(ipt+1,0);

      NSD_sys->SetBinContent(ipt+1,30);
      NSD_sys->SetBinError(ipt+1,2.4);
      BBC_trigger->SetBinContent(ipt+1,0.866);
      BBC_trigger->SetBinError(ipt+1,0.08);
    }
  
  HT0_Npe_sys->Divide(HT0_pt_sys);
  HT2_Npe_sys->Divide(HT2_pt_sys);


  //  HT0_Npe_sys->Divide(BBC_trigger);
  HT0_Npe_sys->Multiply(NSD_sys);
  
  HT2_Npe_sys->Multiply(NSD_sys);
  // HT2_Npe_sys->Divide(BBC_trigger);



  HT0_NPE_spectra_sys=(TH1F *) HT0_Npe_sys->Clone("HT0_NPE_spectra_sys");
  HT2_NPE_spectra_sys=(TH1F *) HT2_Npe_sys->Clone("HT2_NPE_spectra_sys");
  

  // HT0_Npe_sts->Draw();
  // HT2_Npe_sts->SetMarkerColor(2);
  // HT2_Npe_sts->Draw("same");

}

// sts
void Efficiency_Correction_sts(TH1F *HT0_Raw_NPE_sts,TH1F *HT2_Raw_NPE_sts,TH1F *&HT0_NPE_sts,TH1F *&HT2_NPE_sts)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *HT0_Raw_Npe_sts=(TH1F *) HT0_Raw_NPE_sts->Clone("HT0_Raw_Npe_sts");
  TH1F *HT2_Raw_Npe_sts=(TH1F *) HT2_Raw_NPE_sts->Clone("HT2_Raw_Npe_sts");

  TH1F *Eff_BEMC_sts=(TH1F *)eff_BEMC->Clone("Eff_BEMC_sts");
  TH1F *Eff_dedx_sts=(TH1F *)eff_dedx->Clone("Eff_dedx_sts");
  TH1F *Eff_Trigger_HT0_sts=(TH1F *)eff_Trigger_HT0->Clone("Eff_Trigger_HT0_sts");
  TH1F *Eff_Trigger_HT2_sts=(TH1F *)eff_Trigger_HT2->Clone("Eff_Trigger_HT2_sts");
  TH1F *Eff_Tracking_sts=(TH1F *)eff_Tracking->Clone("Eff_Tracking_sts");

  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      Eff_BEMC_sts->SetBinError(ipt+1,0);
      Eff_dedx_sts->SetBinError(ipt+1,0);
      Eff_Trigger_HT0_sts->SetBinError(ipt+1,0);
      Eff_Trigger_HT2_sts->SetBinError(ipt+1,0);
      Eff_Tracking_sts->SetBinError(ipt+1,0);
    }
  
  HT0_Raw_Npe_sts->Divide(Eff_dedx_sts);
  HT0_Raw_Npe_sts->Divide(Eff_BEMC_sts);
  HT0_Raw_Npe_sts->Divide(Eff_Trigger_HT0_sts);
  HT0_Raw_Npe_sts->Divide(Eff_Tracking_sts);

  HT2_Raw_Npe_sts->Divide(Eff_dedx_sts);
  HT2_Raw_Npe_sts->Divide(Eff_BEMC_sts);
  HT2_Raw_Npe_sts->Divide(Eff_Trigger_HT2_sts);
  HT2_Raw_Npe_sts->Divide(Eff_Tracking_sts);

  HT0_NPE_sts=(TH1F *) HT0_Raw_Npe_sts->Clone("HT0_Npe_sts") ;
  HT2_NPE_sts=(TH1F *) HT2_Raw_Npe_sts->Clone("HT2_Npe_sts") ;

}
// sys
void Efficiency_Correction_sys(TH1F *HT0_Raw_NPE_sys,TH1F *HT2_Raw_NPE_sys,TH1F *&HT0_NPE_sys,TH1F *&HT2_NPE_sys)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *HT0_Raw_Npe_sys=(TH1F *) HT0_Raw_NPE_sys->Clone("HT0_Raw_Npe_sys");
  TH1F *HT2_Raw_Npe_sys=(TH1F *) HT2_Raw_NPE_sys->Clone("HT2_Raw_Npe_sys");

  TH1F *Eff_BEMC_sys=(TH1F *)eff_BEMC->Clone("Eff_BEMC_sys");
  TH1F *Eff_dedx_sys=(TH1F *)eff_dedx->Clone("Eff_dedx_sys");
  TH1F *Eff_Trigger_HT0_sys=(TH1F *)eff_Trigger_HT0->Clone("Eff_Trigger_HT0_sys");
  TH1F *Eff_Trigger_HT2_sys=(TH1F *)eff_Trigger_HT2->Clone("Eff_Trigger_HT2_sys");
  TH1F *Eff_Tracking_sys=(TH1F *)eff_Tracking->Clone("Eff_Tracking_sys");

  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      HT0_Raw_Npe_sys->SetBinError(ipt+1,0);
      HT2_Raw_Npe_sys->SetBinError(ipt+1,0);
    }
  
  HT0_Raw_Npe_sys->Divide(Eff_dedx_sys);
  HT0_Raw_Npe_sys->Divide(Eff_BEMC_sys);
  HT0_Raw_Npe_sys->Divide(Eff_Trigger_HT0_sys);
  HT0_Raw_Npe_sys->Divide(Eff_Tracking_sys);

  HT2_Raw_Npe_sys->Divide(Eff_dedx_sys);
  HT2_Raw_Npe_sys->Divide(Eff_BEMC_sys);
  HT2_Raw_Npe_sys->Divide(Eff_Trigger_HT2_sys);
  HT2_Raw_Npe_sys->Divide(Eff_Tracking_sys);

  HT0_NPE_sys=(TH1F *) HT0_Raw_Npe_sys->Clone("HT0_Npe_sys") ;
  HT2_NPE_sys=(TH1F *) HT2_Raw_Npe_sys->Clone("HT2_Npe_sys") ;

}
void Raw_NPE_spectra_sts(TH1F *HT0_inclusive,TH1F *HT2_inclusive,TH1F *HT0_phe,TH1F *HT2_phe,TH1F *&HT0_Raw_NPE_sts,TH1F *&HT2_Raw_NPE_sts) 
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *HT0_inclusive_sts=(TH1F *) HT0_inclusive->Clone("HT0_inclusive_sts");
  TH1F *HT2_inclusive_sts=(TH1F *) HT2_inclusive->Clone("HT2_inclusive_sts");

  TH1F *HT0_phe_sts=(TH1F *) HT0_phe->Clone("HT0_phe_sts");
  TH1F *HT2_phe_sts=(TH1F *) HT2_phe->Clone("HT2_phe_sts");

  TH1F *purity_sts= (TH1F *)eff_purity->Clone("purity_sts");
  TH1F *phe_re_sts= (TH1F *)eff_PHE_re->Clone("phe_re_sts");
  
  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      purity_sts->SetBinError(ipt+1,0);
      phe_re_sts->SetBinError(ipt+1,0);
    }
  
  HT0_inclusive_sts->Multiply(purity_sts);
  HT0_phe_sts->Divide(phe_re_sts);
  HT0_Raw_NPE_sts=(TH1F *) HT0_inclusive_sts->Clone("HT0_Raw_NPE_sts");
  HT0_Raw_NPE_sts->Add(HT0_phe_sts,-1);
  //  TH1F *HT0_S_B_sts=(TH1F *)HT0_Raw_NPE_sts->Clone("HT0_S_B_sts");

  HT2_inclusive_sts->Multiply(purity_sts);
  HT2_phe_sts->Divide(phe_re_sts);
  HT2_Raw_NPE_sts=(TH1F *) HT2_inclusive_sts->Clone("HT2_Raw_NPE_sts");
  HT2_Raw_NPE_sts->Add(HT2_phe_sts,-1);
  //  TH1F *HT0_S_B_sts=(TH1F *)HT0_Raw_NPE_sts->Clone("HT0_S_B_sts");



}

void Raw_NPE_spectra_sys(TH1F *HT0_inclusive,TH1F *HT2_inclusive,TH1F *HT0_phe,TH1F *HT2_phe,TH1F *&HT0_Raw_NPE_sys,TH1F *&HT2_Raw_NPE_sys) 
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *HT0_inclusive_sys=(TH1F *) HT0_inclusive->Clone("HT0_inclusive_sys");
  TH1F *HT2_inclusive_sys=(TH1F *) HT2_inclusive->Clone("HT2_inclusive_sys");

  TH1F *HT0_phe_sys=(TH1F *) HT0_phe->Clone("HT0_phe_sys");
  TH1F *HT2_phe_sys=(TH1F *) HT2_phe->Clone("HT2_phe_sys");

  TH1F *purity_sys= (TH1F *)eff_purity->Clone("purity_sys");
  TH1F *phe_re_sys= (TH1F *)eff_PHE_re->Clone("phe_re_sys");
  
  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      HT0_inclusive_sys->SetBinError(ipt+1,0);
      HT2_inclusive_sys->SetBinError(ipt+1,0);
      HT0_phe_sys->SetBinError(ipt+1,0);
      HT2_phe_sys->SetBinError(ipt+1,0);
    }
  
  HT0_inclusive_sys->Multiply(purity_sys);
  HT0_phe_sys->Divide(phe_re_sys);
  HT0_Raw_NPE_sys=(TH1F *) HT0_inclusive_sys->Clone("HT0_Raw_NPE_sys");
  HT0_Raw_NPE_sys->Add(HT0_phe_sys,-1);
  //  TH1F *HT0_S_B_sys=(TH1F *)HT0_Raw_NPE_sys->Clone("HT0_S_B_sys");

  HT2_inclusive_sys->Multiply(purity_sys);
  HT2_phe_sys->Divide(phe_re_sys);
  HT2_Raw_NPE_sys=(TH1F *) HT2_inclusive_sys->Clone("HT2_Raw_NPE_sys");
  HT2_Raw_NPE_sys->Add(HT2_phe_sys,-1);
  //  TH1F *HT0_S_B_sys=(TH1F *)HT0_Raw_NPE_sys->Clone("HT0_S_B_sys");

}

void Draw_NPE_PHE_ratio(TH1F *HT0_inclusive,TH1F *HT2_inclusive,TH1F *HT0_phe,TH1F *HT2_phe,TH1F *&HT0_Raw_NPE_sts,TH1F *&HT0_Raw_NPE_sys,TH1F *&HT2_Raw_NPE_sts,TH1F *&HT2_Raw_NPE_sys)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  
  
  TH1F *HT0_inclusive_sts=(TH1F *) HT0_inclusive->Clone("HT0_inclusive_sts");
  TH1F *HT0_inclusive_sys=(TH1F *) HT0_inclusive->Clone("HT0_inclusive_sys");

  TH1F *HT2_inclusive_sts=(TH1F *) HT2_inclusive->Clone("HT2_inclusive_sts");
  TH1F *HT2_inclusive_sys=(TH1F *) HT2_inclusive->Clone("HT2_inclusive_sys");


  TH1F *HT0_phe_sts=(TH1F *) HT0_phe->Clone("HT0_phe_sts");
  TH1F *HT0_phe_sys=(TH1F *) HT0_phe->Clone("HT0_phe_sys");

  TH1F *HT2_phe_sts=(TH1F *) HT2_phe->Clone("HT2_phe_sts");
  TH1F *HT2_phe_sys=(TH1F *) HT2_phe->Clone("HT2_phe_sys");

  TH1F *purity_sts= (TH1F *)eff_purity->Clone("purity_sts");
  TH1F *purity_sys= (TH1F *)eff_purity->Clone("purity_sys");


  TH1F *phe_re_sts= (TH1F *)eff_PHE_re->Clone("phe_re_sts");
  TH1F *phe_re_sys= (TH1F *)eff_PHE_re->Clone("phe_re_sys");

  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      HT0_inclusive_sys->SetBinError(ipt+1,0);
      HT2_inclusive_sys->SetBinError(ipt+1,0);
      HT0_phe_sys->SetBinError(ipt+1,0);
      HT2_phe_sys->SetBinError(ipt+1,0);

      purity_sts->SetBinError(ipt+1,0);
      phe_re_sts->SetBinError(ipt+1,0);
    }
  //--------------------------------------sts-----------------------------------
  //HT0
  HT0_inclusive_sts->Multiply(purity_sts);
  HT0_phe_sts->Divide(phe_re_sts);
  HT0_Raw_NPE_sts=(TH1F *) HT0_inclusive_sts->Clone("HT0_Raw_NPE_sts");
  HT0_Raw_NPE_sts->Add(HT0_phe_sts,-1);
  TH1F *HT0_S_B_sts=(TH1F *)HT0_Raw_NPE_sts->Clone("HT0_S_B_sts");
  HT0_S_B_sts->Divide(HT0_phe_sts);  
  //HT2
  HT2_inclusive_sts->Multiply(purity_sts);
  HT2_phe_sts->Divide(phe_re_sts);
  HT2_Raw_NPE_sts=(TH1F *) HT2_inclusive_sts->Clone("HT2_Raw_NPE_sts");
  HT2_Raw_NPE_sts->Add(HT2_phe_sts,-1);
  TH1F *HT2_S_B_sts=(TH1F *)HT2_Raw_NPE_sts->Clone("HT2_S_B_sts");
  HT2_S_B_sts->Divide(HT2_phe_sts);  

  //--------------------------------------sys-----------------------------------

  HT0_inclusive_sys->Multiply(purity_sys);
  HT0_phe_sys->Divide(phe_re_sys);
  HT0_Raw_NPE_sys=(TH1F *) HT0_inclusive_sys->Clone("HT0_Raw_NPE_sys");
  HT0_Raw_NPE_sys->Add(HT0_phe_sys,-1);
  TH1F *HT0_S_B_sys=(TH1F *)HT0_Raw_NPE_sys->Clone("HT0_S_B_sys");
  HT0_S_B_sys->Divide(HT0_phe_sys);  
  //HT2
  HT2_inclusive_sys->Multiply(purity_sys);
  HT2_phe_sys->Divide(phe_re_sys);
  HT2_Raw_NPE_sys=(TH1F *) HT2_inclusive_sys->Clone("HT2_Raw_NPE_sys");
  HT2_Raw_NPE_sys->Add(HT2_phe_sys,-1);
  TH1F *HT2_S_B_sys=(TH1F *)HT2_Raw_NPE_sys->Clone("HT2_S_B_sys");
  HT2_S_B_sys->Divide(HT2_phe_sys);  
  
    
  TCanvas *c5=new TCanvas("c5","",1000,800);
  TH2F *h5=new TH2F("h5","h5",100,2,10,100,0,2);
  h5->Draw();
  
  TGraphErrors  *gr_HT0_S_B_sts= new TGraphErrors(HT0_S_B_sts);
  TGraphErrors  *gr_HT2_S_B_sts= new TGraphErrors(HT2_S_B_sts);

  TGraphErrors  *gr_HT0_S_B_sys= new TGraphErrors(HT0_S_B_sys);
  TGraphErrors  *gr_HT2_S_B_sys= new TGraphErrors(HT2_S_B_sys);
 
  // gr_HT0_S_B_sys->Draw();
  // //  return; 

  gr_HT0_S_B_sys->SetMarkerStyle(20);
  gr_HT0_S_B_sts->SetMarkerStyle(20);

  gr_HT0_S_B_sys->SetMarkerColor(4);
  gr_HT0_S_B_sts->SetMarkerStyle(4);

  gr_HT0_S_B_sys->SetLineColor(2);
  gr_HT0_S_B_sys->SetLineWidth(3);

  gr_HT0_S_B_sts->SetMarkerStyle(4);

  // gr_HT0_S_B_sts->Draw("samePE");
  gr_HT0_S_B_sys->Draw("samepE[]");


  
 
  
  // HT0_S_B_sys->Draw("same");
  // HT2_S_B_sys->Draw("same");


 

  



 





}
void Draw_Pt_spectra(TH1F * Inclusive[nTrg],TH1F * Photonic_unlike_pt[nTrg],TH1F * Photonic_like_pt[nTrg],TH1F * Photonic_unlike_like_pt[nTrg],TH1F * Inclusive_ps[nTrg],TH1F * Photonic_unlike_pt_ps[nTrg],TH1F * Photonic_like_pt_ps[nTrg],TH1F * Photonic_unlike_like_pt_ps[nTrg])
{

  gStyle->SetOptStat(00000);

  //  Inclusive[0]->Draw();
  Int_t bMarkerstyle=20;
  
  for(Int_t iTrg=0;iTrg<nTrg;iTrg++)
    {
      Inclusive[iTrg]->SetLineColor(1);
      Inclusive_ps[iTrg]->SetLineColor(1);
      
      Photonic_unlike_like_pt[iTrg]->SetLineColor(2);
      Photonic_unlike_like_pt_ps[iTrg]->SetLineColor(2);

      Inclusive[iTrg]->SetMarkerColor(1);
      Inclusive_ps[iTrg]->SetMarkerColor(1);
      
      Photonic_unlike_like_pt[iTrg]->SetMarkerColor(2);
      Photonic_unlike_like_pt_ps[iTrg]->SetMarkerColor(2);


      Inclusive[iTrg]->SetMarkerStyle(bMarkerstyle+iTrg);
      Inclusive_ps[iTrg]->SetMarkerStyle(bMarkerstyle+iTrg);

      Photonic_unlike_like_pt[iTrg]->SetMarkerStyle(bMarkerstyle+4+iTrg);
      Photonic_unlike_like_pt_ps[iTrg]->SetMarkerStyle(bMarkerstyle+4+iTrg);
    }
  TLegend *legend = new TLegend(0.7,0.65,0.85,0.85);
  legend->AddEntry(Inclusive[0],"HT0 inclusive e","lp");
  legend->AddEntry(Photonic_unlike_like_pt[0],"HT0 photonic e","lp");
  legend->AddEntry(Inclusive[1],"HT2 inclusive e","lp");
  legend->AddEntry(Photonic_unlike_like_pt[1],"HT2 photonic e","lp");
  legend->SetTextFont(62);
  legend->SetTextSize(0.035);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);   
  TCanvas *c4=new TCanvas("c4","",1200,1000);
  gPad->SetLogy();
  TH2F *h4=new TH2F("h4","",100,2,10,100,1,3e4);
  h4->Draw();
  h4->GetXaxis()->SetTitle("p_{T} GeV/c");
  h4->GetYaxis()->SetTitle("Counts");
  Inclusive[0]->Draw("same");
  Photonic_unlike_like_pt[0]->Draw("same");
  Inclusive[1]->Draw("same");
  Photonic_unlike_like_pt[1]->Draw("same");
  legend->Draw();
  c4->SaveAs("pT_raw_spectra_after_prescale.pdf");

  TCanvas *c3=new TCanvas("c2","",1200,1000);
  gPad->SetLogy();
  TH2F *h3=new TH2F("h3","",100,2,10,100,1,3e6);
  h3->Draw();
  h3->GetXaxis()->SetTitle("p_{T} GeV/c");
  h3->GetYaxis()->SetTitle("Counts");

  Inclusive_ps[0]->Draw("same");
  Photonic_unlike_like_pt_ps[0]->Draw("same");
  Inclusive_ps[1]->Draw("same");
  Photonic_unlike_like_pt_ps[1]->Draw("same");
  legend->Draw();
  c3->SaveAs("pT_raw_spectra_before_prescale.pdf");
}
void Draw_Efficiency()
{
  TH1F *Eff_BEMC=(TH1F *)eff_BEMC->Clone("Eff_BEMC");
  TH1F *Eff_dedx=(TH1F *)eff_dedx->Clone("Eff_dedx");
  TH1F *Eff_Trigger_HT0=(TH1F *)eff_Trigger_HT0->Clone("Eff_Trigger_HT0");
  TH1F *Eff_Trigger_HT2=(TH1F *)eff_Trigger_HT2->Clone("Eff_Trigger_HT2");
  TH1F *Eff_Tracking=(TH1F *)eff_Tracking->Clone("Eff_Tracking");
  TH1F *Eff_PHE_re=(TH1F *)eff_PHE_re->Clone("Eff_PHE_re");
  TH1F *Eff_purity=(TH1F *)eff_purity->Clone("Eff_purity");

  Eff_BEMC->SetMarkerStyle(20);
  Eff_dedx->SetMarkerStyle(20);
  Eff_Trigger_HT0->SetMarkerStyle(20);
  Eff_Trigger_HT2->SetMarkerStyle(20);
  Eff_Tracking->SetMarkerStyle(20);
  Eff_PHE_re->SetMarkerStyle(20);
  Eff_purity->SetMarkerStyle(20);

  Eff_BEMC->SetMarkerColor(1);
  Eff_dedx->SetMarkerColor(2);
  Eff_Trigger_HT0->SetMarkerColor(3);
  Eff_Trigger_HT2->SetMarkerColor(4);
  Eff_Tracking->SetMarkerColor(5);
  Eff_PHE_re->SetMarkerColor(6);
  Eff_purity->SetMarkerColor(9);

  Eff_BEMC->SetLineColor(1);
  Eff_dedx->SetLineColor(2);
  Eff_Trigger_HT0->SetLineColor(3);
  Eff_Trigger_HT2->SetLineColor(4);
  Eff_Tracking->SetLineColor(5);
  Eff_PHE_re->SetLineColor(6);
  Eff_purity->SetLineColor(8);


  gStyle->SetOptStat(0);
  TCanvas *c2=new TCanvas("c2","",1200,800);
  c2->cd();

  TH2F *h2=new TH2F("h2","",100,2,10,100,-0.1,1.8);
  h2->Draw();
  h2->GetYaxis()->SetTitle("efficiency");
  h2->GetXaxis()->SetTitle("p_{T} (Gev/c)");

  Eff_BEMC->Draw("samePE");
  Eff_dedx->Draw("samePE");
  Eff_Trigger_HT0->Draw("samePE");
  Eff_Trigger_HT2->Draw("samePE");
  Eff_Tracking->Draw("samePE");
  Eff_PHE_re->Draw("samePE");
  Eff_purity->Draw("samePE");
  
  TLegend *legend = new TLegend(0.2,0.65,0.35,0.85);
  //  TLegend *legend = new TLegend(0.15,0.65,0.55,0.85);
  
  legend->AddEntry(Eff_BEMC,"BEMC","lp");
  legend->AddEntry(Eff_dedx,"dEdx","lp");
  legend->AddEntry(Eff_Trigger_HT0,"HT0 Trigger","lp");
  legend->AddEntry(Eff_Trigger_HT2,"HT2 Trigger","lp");
  legend->AddEntry(Eff_Tracking,"Tracking","lp");
  legend->AddEntry(Eff_PHE_re,"PHE_re","lp");
  legend->AddEntry(Eff_purity,"purity","lp");

  legend->SetTextFont(62);
  legend->SetTextSize(0.035);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);   
  
  legend->Draw();
  c2->SaveAs("Efficiency.pdf");
}


void setpad(TPad *pad,float left, float right, float top, float bottom)
{
  pad->SetFillColor(10);
  pad->SetBorderMode(0);
  pad->SetBorderSize(10);
  pad->SetFrameFillColor(10);
  pad->SetFrameBorderMode(0);
  pad->SetFrameBorderSize(0);
  pad->SetLeftMargin(left);
  pad->SetRightMargin(right);
  pad->SetTopMargin(top);
  pad->SetRightMargin(right);                                                  
  
  pad->SetBottomMargin(bottom);
}
TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
  TLatex *latex = new TLatex(x,y,text);
  latex->SetNDC();
  latex->SetTextFont(textFont);
  latex->SetTextSize(textSize);
  latex->SetTextColor(colorIndex);
  latex->Draw("same");
  return latex;
}
