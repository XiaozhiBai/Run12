
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
#include "Input/binning_shift.h"

#include "../mBinning_MB.h"
using namespace std;


const Double_t deta=1.4;

void Draw_Efficiency();
void Draw_Pt_spectra(TH1F *a1,TH1F *a2,TH1F *a3,TH1F *a4,TH1F *a5,TH1F *a6,TH1F *a7,TH1F *a8);
 void Efficiency_Correction_sts(TH1F *Raw_NPE_sts,TH1F *&NPE_sts);
 void Efficiency_Correction_sys(TH1F *Raw_NPE_sys,TH1F *&NPE_sys);

void Cross_Section_sts(TH1F *NPE_sts,TH1F *&NPE_spectra_sts,Double_t );
void Cross_Section_sys(TH1F *NPE_sys,TH1F *&NPE_spectra_sys,Double_t );



void Draw_Cross_Section(TH1F *NPE_sts,TH1F *NPE_sys,TH1F *HT_NPE_sts,TH1F *HT_NPE_sys,TH1F *&,TH1F *&);

//void Omega_Phi_Correction(TH1F *MB_NPE_sts,TH1F *MB_NPE_sys,TH1F *&,TH1F *&);

// void Draw_NPE_PHE_ratio(TH1F *,TH1F *,TH1F *,TH1F *,TH1F *&,TH1F *&,TH1F *&,TH1F *&);  
 void Raw_NPE_spectra_sts(TH1F *,TH1F *,TH1F *&);  
 void Raw_NPE_spectra_sys(TH1F *,TH1F *,TH1F *&);  


TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
void setpad(TPad *,float , float , float , float );

//TH1F * Correc=new TH1F("Correc","",Nbins_HT,Pt_bin_HT);

TFile *infile_data =new TFile("Input/hist_5_1.root","read");
TFile *infile_Tof_match =new TFile("Input/TOF_Match_efficiency.root","read");
TFile *infile_dEdx_cut =new TFile("Input/nSigma_Cut_efficiency.root","read");
TFile *infile_Tof_cut =new TFile("Input/Tof_Beta_cuts_efficiency.root","read");
TFile *infile_poe =new TFile("Input/Poe_efficiency.root","read");

TFile *infile_Tracking =new TFile("Input/Tracking_efficiency_MB.root","read");
TFile *infile_PHE_re =new TFile("Input/Photonic_re_Efficiency.root","read");
TFile *infile_purity =new TFile("Input/purity_MB_com.root","read");

TH1F *eff_Tof_match=(TH1F *) infile_Tof_match->Get("efficiency_tof_match");
TH1F *eff_dedx=(TH1F *) infile_dEdx_cut->Get("nsigmaE_MB");

TH1F *eff_Tof_cuts=(TH1F *) infile_Tof_cut->Get("Tof_beta_cut_MB");
TH1F *eff_poe=(TH1F *) infile_poe->Get("Poe_efficiency");


TH1F *eff_Tracking=(TH1F *) infile_Tracking->Get("Tracking_efficiency_MB");
TH1F *eff_PHE_re=(TH1F *) infile_PHE_re->Get("PHE_re_efficiency");
TH1F *eff_purity=(TH1F *) infile_purity->Get("purity_MB");

TH1F *mh1MB_Nevents_ps=(TH1F *) infile_data->Get("mh1MB_Nevents_psTrg2");


TFile *File_fonll=new TFile("Input/fonll200gev.root","READ");
TFile *File_Jpsi=new TFile("Input/Jpsi.root","READ");


TFile *File_Omega_phi=new TFile("Input/Omega_phi.root","READ");


void NPE_Cross_section()
{
  // eff_Tof_macth->Draw();
  // eff_dedx->Draw();
  // eff_Tof_cuts->Draw();
  // eff_Tracking->Draw();
  // eff_PHE_re->Draw();
  // eff_purity->Draw();
  // return;
  
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


   Double_t  NmbEvents= mh1MB_Nevents_ps->GetBinContent(2);
   cout<< NmbEvents<<endl;

  TH1F *Inclusive;
  TH1F *Inclusive_ps;

  TH2F *Photonic_unlike;
  TH2F *Photonic_like;
  TH2F *Photonic_unlike_like;
  
  TH2F *Photonic_ps_unlike;
  TH2F *Photonic_ps_like;
  TH2F *Photonic_ps_unlike_like;
  

  TH1F *Photonic_unlike_pt;
  TH1F *Photonic_like_pt;
  TH1F *Photonic_unlike_like_pt;

  TH1F *Photonic_unlike_pt_ps;
  TH1F *Photonic_like_pt_ps;
  TH1F *Photonic_unlike_like_pt_ps;



  
  
  TString Photonic_unlikename={"mh2InvMassUnlikeTrg2"};
  TString Photonic_likename={"mh2InvMasslikeTrg2"};
  TString Photonic_unlikename_ps={"mh2InvMassUnlike_psTrg2"};
  TString Photonic_likename_ps={"mh2InvMasslike_psTrg2"};

  TString Inclusive_name={"mh1electronPtTrg2"};
  TString Inclusive_name_ps={"mh1electronPt_psTrg2"};
  

      
      // incusve 
      Inclusive=(TH1F *) infile_data->Get(Inclusive_name);
      Inclusive_ps=(TH1F *) infile_data->Get(Inclusive_name_ps);


      //photonic 

      Photonic_unlike=(TH2F *) infile_data->Get(Photonic_unlikename);
      Photonic_like=(TH2F *) infile_data->Get(Photonic_likename);


      Photonic_unlike_pt=(TH1F *) Photonic_unlike->ProjectionY("unlike");
      Photonic_like_pt=(TH1F *) Photonic_like->ProjectionY("like");
      Photonic_unlike_like_pt=(TH1F *) Photonic_unlike_pt->Clone("unlike_like"); 
      Photonic_unlike_like_pt->Sumw2();
      Photonic_unlike_like_pt->Add(Photonic_like_pt,-1);

           
      Photonic_ps_unlike=(TH2F *) infile_data->Get(Photonic_unlikename_ps);
      Photonic_ps_like=(TH2F *) infile_data->Get(Photonic_likename_ps);


      Photonic_unlike_pt_ps=(TH1F *) Photonic_ps_unlike->ProjectionY("unlike_ps");
      Photonic_like_pt_ps=(TH1F *) Photonic_ps_like->ProjectionY("like_ps");
      Photonic_unlike_like_pt_ps=(TH1F *) Photonic_unlike_pt_ps->Clone("unlike_like_ps"); 
      Photonic_unlike_like_pt_ps->Sumw2();
      Photonic_unlike_like_pt_ps->Add(Photonic_like_pt_ps,-1);


  

   Draw_Efficiency();


   
   Draw_Pt_spectra(Inclusive,Photonic_unlike_pt,Photonic_like_pt,Photonic_unlike_like_pt,Inclusive_ps,Photonic_unlike_pt_ps,Photonic_like_pt_ps,Photonic_unlike_like_pt_ps);


  TH1F *inclusive=(TH1F *) Inclusive_ps->Rebin(NpT_bins_run12_MB,"inclusive",pt_run12_MB);
  TH1F *phe=(TH1F *) Photonic_unlike_like_pt_ps->Rebin(NpT_bins_run12_MB,"phe",pt_run12_MB);
  
  TH1F *Raw_NPE_sts;
  TH1F *Raw_NPE_sys;

  // Get Raw NPE yeild 
   Raw_NPE_spectra_sts(inclusive,phe,Raw_NPE_sts);
   Raw_NPE_spectra_sys(inclusive,phe,Raw_NPE_sys);

   // Raw_NPE_sts->Scale(1,"width");
   // Raw_NPE_sts->Draw();
   
   // return;
   
   TH1F *NPE_sts;
   TH1F *NPE_sys;
   TH1F *HT2_NPE_sts;
   TH1F *HT2_NPE_sys;
   
   // Get  NPE yeild corrected by the efficiency 
   Efficiency_Correction_sts(Raw_NPE_sts,NPE_sts);
   Efficiency_Correction_sys(Raw_NPE_sys,NPE_sys);


   


  TH1F *NPE_spectra_sts;
  TH1F *NPE_spectra_sys;
  

  //   Caculate the cross section

  Cross_Section_sts(NPE_sts,NPE_spectra_sts,NmbEvents);

  Cross_Section_sys(NPE_sys,NPE_spectra_sys,NmbEvents);

  NPE_spectra_sts->Draw();
  

  
  
  TH1F *run12_MB_NPE_sts=new TH1F("run12_MB_NPE_sts","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F *run12_MB_NPE_sys=new TH1F("run12_MB_NPE_sys","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F *run12_MB_NPE_FONLL_sts=new TH1F("run12_MB_NPE_FONLL_sts","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F *run12_MB_NPE_FONLL_sys=new TH1F("run12_MB_NPE_FONLL_sys","",NpT_bins_run12_MB,pt_run12_MB);

  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {


      run12_MB_NPE_sts->SetBinContent(ipt+1,NPE_spectra_sts->GetBinContent(ipt+1));
      run12_MB_NPE_sts->SetBinError(ipt+1,NPE_spectra_sts->GetBinError(ipt+1));

      run12_MB_NPE_sys->SetBinContent(ipt+1,NPE_spectra_sys->GetBinContent(ipt+1));
      run12_MB_NPE_sys->SetBinError(ipt+1,NPE_spectra_sys->GetBinError(ipt+1));


    }



  
  // run12_MB_NPE_sts->Draw();
  // run12_MB_NPE_sys->SetLineColor(2);
  // run12_MB_NPE_sys->Draw("same");
  
  Draw_Cross_Section(NPE_spectra_sts,NPE_spectra_sys,run12_MB_NPE_sts,run12_MB_NPE_sys,run12_MB_NPE_FONLL_sts,run12_MB_NPE_FONLL_sys);
}
void Draw_Cross_Section(TH1F *NPE_sts,TH1F *NPE_sys,TH1F *MB_NPE_sts,TH1F *MB_NPE_sys,TH1F *&run12_MB_NPE_FONLL_sts,TH1F *&run12_MB_NPE_FONLL_sys)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *FONLL=new TH1F("FONLL","",NpT_bins_run12_MB,pt_run12_MB);


  TGraphErrors *gFONLLu=( TGraphErrors *) File_fonll->Get("gFONLLu");
  TGraphErrors *gFONLLc=( TGraphErrors *) File_fonll->Get("gFONLLc");
  TGraphErrors *gFONLLl=( TGraphErrors *) File_fonll->Get("gFONLLl");

  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
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

  TH2F *h6=new TH2F("h6","",100,0.2,5,100,1e-13,1e-2);
  h6->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h6->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");

  h6->Draw();
  // MB_NPE_sys->Scale(0.4);

  gFONLLc->Draw("same");
  gFONLLu->Draw("same");
  gFONLLl->Draw("same");

  TF1* f1 = new TF1("f1","[0]/(pow(x,[1])+[2])",0.2,4);
  f1->SetParameter(0,10);
  f1->SetParameter(1,8);
  f1->SetParameter(2,0.7); 

  TH1F *MB_NPE_sts_fit=(TH1F *)  MB_NPE_sts->Clone("MB_NPE_sts_fit"); 
  
  MB_NPE_sts_fit->Fit("f1","R","same" ,0.2,4);
  MB_NPE_sts_fit->Fit("f1","R","same" ,0.2,4);
  
  const int nbins = NpT_bins_run12_MB+1;

  TH1D* hCorrection = get_bin_shift((unsigned long int)1e8,f1,nbins,pt_run12_MB);
   MB_NPE_sts->Divide(hCorrection);
   MB_NPE_sys->Divide(hCorrection);

   TH1F *Omega_sts=(TH1F *) File_Omega_phi->Get("Omega_sts");
   TH1F *Omega_sys=(TH1F *) File_Omega_phi->Get("Omega_sys");

   TH1F *Phi_sts=(TH1F *) File_Omega_phi->Get("Phi_sts");
   TH1F *Phi_sys=(TH1F *) File_Omega_phi->Get("Phi_sys");


   MB_NPE_sts->Add(Omega_sts,-1);
   MB_NPE_sys->Add(Omega_sys,-1);

   MB_NPE_sts->Add(Phi_sts,-1);
   MB_NPE_sys->Add(Phi_sys,-1);


  TGraphErrors  *gr_run12_MB_NPE_sts=new TGraphErrors(MB_NPE_sts);
  TGraphErrors  *gr_run12_MB_NPE_sys=new TGraphErrors(MB_NPE_sys);

  run12_MB_NPE_FONLL_sts->Divide(MB_NPE_sts,FONLL,1,1);
  run12_MB_NPE_FONLL_sys->Divide(MB_NPE_sys,FONLL,1,1);
  
  TGraphErrors  *gr_run12_MB_NPE_FONLL_sts=new TGraphErrors(run12_MB_NPE_FONLL_sts);
  TGraphErrors  *gr_run12_MB_NPE_FONLL_sys=new TGraphErrors(run12_MB_NPE_FONLL_sys);
  
  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      gr_run12_MB_NPE_sts->SetPointError(ipt,0,gr_run12_MB_NPE_sts->GetErrorY(ipt));
      gr_run12_MB_NPE_sys->SetPointError(ipt,0.05,gr_run12_MB_NPE_sys->GetErrorY(ipt));
      
      gr_run12_MB_NPE_FONLL_sts->SetPointError(ipt,0,gr_run12_MB_NPE_FONLL_sts->GetErrorY(ipt));
      gr_run12_MB_NPE_FONLL_sys->SetPointError(ipt,0.05,gr_run12_MB_NPE_FONLL_sys->GetErrorY(ipt));
    }  


  gr_run12_MB_NPE_sts->SetMarkerStyle(20);
  gr_run12_MB_NPE_sys->SetMarkerColor(2);
  gr_run12_MB_NPE_sts->SetLineColor(2);
  gr_run12_MB_NPE_sys->SetLineColor(2);


  gr_run12_MB_NPE_sys->Draw("sameE1p[]");
  gr_run12_MB_NPE_sts->Draw("samePE");
  c6->SaveAs("Cross_section.pdf");
  
  TCanvas *c7=new TCanvas("c7","",800,600);
  //  gPad->SetLogy();
  TH2F *h7=new TH2F("h7","",100,0.2,5,100,0,3);
  h7->Draw();
  h7->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h7->GetYaxis()->SetTitle("Data/FONLL");
  

  //run12_MB_NPE_FONLL_sys->Draw("same");

  gr_run12_MB_NPE_FONLL_sts->SetMarkerStyle(20);
  gr_run12_MB_NPE_FONLL_sts->SetMarkerColor(2);
  gr_run12_MB_NPE_FONLL_sys->SetMarkerColor(2);

  gr_run12_MB_NPE_FONLL_sts->SetLineColor(2);
  gr_run12_MB_NPE_FONLL_sys->SetLineColor(2);




  gr_run12_MB_NPE_FONLL_sys->Draw("samepe1[]");
  gr_run12_MB_NPE_FONLL_sts->Draw("samePE");
  
  gR_FONLLc->Draw("same");
  gR_FONLLu->Draw("same");
  gR_FONLLl->Draw("same");
 c7->SaveAs("Cross_section_fonll.pdf");



 TFile *file=new TFile("run12_Npe_MB.root","RECREATE");

 MB_NPE_sts->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 MB_NPE_sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 run12_MB_NPE_FONLL_sts->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 run12_MB_NPE_FONLL_sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");

 MB_NPE_sts->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");
 MB_NPE_sys->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");
 run12_MB_NPE_FONLL_sts->GetYaxis()->SetTitle("Data/FONLL");
 run12_MB_NPE_FONLL_sys->GetYaxis()->SetTitle("Data/FONLL");

 MB_NPE_sts->Write();
 MB_NPE_sys->Write();
 run12_MB_NPE_FONLL_sts->Write();
 run12_MB_NPE_FONLL_sys->Write();

 file->Close();

}
  
void Cross_Section_sts(TH1F *NPE_sts,TH1F *&NPE_spectra_sts,Double_t NmbEvents)
{

  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *Npe_sts=(TH1F *) NPE_sts->Clone("Npe_sts");

  //2Pi
  Npe_sts->Scale(1./(2*TMath::Pi()));
  
  // delta eta
  Npe_sts->Scale(1./deta);
  
  // delta pT
  Npe_sts->Scale(1.,"width");
  

  //pT
  TH1F *Norm_sts=new TH1F("Norm_sts","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F *Trigger_bias_sts=new TH1F("Trigger_bias_sts","",NpT_bins_run12_MB,pt_run12_MB);

  Npe_sts->Scale(1./NmbEvents);
  

  //Pt*NSD*BBCTrigger_efficiency*2(eplus and eminus)

  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      Norm_sts->SetBinContent(ipt+1,Npe_sts->GetBinCenter(ipt+1)*2/30.);
      Norm_sts->SetBinError(ipt+1,0);
      Trigger_bias_sts->SetBinContent(ipt+1,0.64);
      Trigger_bias_sts->SetBinError(ipt+1,0);
      
      
    }
  
  Npe_sts->Divide(Norm_sts);
  Npe_sts->Multiply(Trigger_bias_sts);


  NPE_spectra_sts=(TH1F *) Npe_sts->Clone("NPE_spectra_sts");

  

  // Npe_sts->Draw();
  // HT2_Npe_sts->SetMarkerColor(2);
  // HT2_Npe_sts->Draw("same");

}

void Cross_Section_sys(TH1F *NPE_sys,TH1F *&NPE_spectra_sys,Double_t NmbEvents)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *Npe_sys=(TH1F *) NPE_sys->Clone("Npe_sys");

  //2Pi
  Npe_sys->Scale(1./(2*TMath::Pi()));

  // delta eta
  Npe_sys->Scale(1./deta);


  // delta pT
  Npe_sys->Scale(1.,"width");


  Npe_sys->Scale(1./NmbEvents);

  //pT
  TH1F *pt_sys=new TH1F("pt_sys","",NpT_bins_run12_MB,pt_run12_MB);


  TH1F *NSD_sys=new TH1F("NSD_sys","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F *BBC_trigger=new TH1F("BBC_Trigger","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F *Trigger_bias_sys=new TH1F("Trigger_bias_sys","",NpT_bins_run12_MB,pt_run12_MB);

  //Pt*2(eplus and eminus)

  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      pt_sys->SetBinContent(ipt+1,Npe_sys->GetBinCenter(ipt+1)*2);  //*0.866*2/30.);
      pt_sys->SetBinError(ipt+1,0);


      NSD_sys->SetBinContent(ipt+1,30);
      NSD_sys->SetBinError(ipt+1,2.4);
      BBC_trigger->SetBinContent(ipt+1,0.866);
      BBC_trigger->SetBinError(ipt+1,0.08);

      Trigger_bias_sys->SetBinContent(ipt+1,0.64);
      Trigger_bias_sys->SetBinError(ipt+1,0.08);

    }
  
  Npe_sys->Divide(pt_sys);
  //  Npe_sys->Divide(pt_sys);



  //  Npe_sys->Divide(BBC_trigger);
  Npe_sys->Multiply(NSD_sys);
  Npe_sys->Multiply(Trigger_bias_sys);
  

  // HT2_Npe_sys->Divide(BBC_trigger);



  NPE_spectra_sys=(TH1F *) Npe_sys->Clone("NPE_spectra_sys");

  

  // Npe_sts->Draw();
  // HT2_Npe_sts->SetMarkerColor(2);
  // HT2_Npe_sts->Draw("same");

}

// sts
void Efficiency_Correction_sts(TH1F *Raw_NPE_sts,TH1F *&NPE_sts)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

 
  
  TH1F *Raw_Npe_sts=(TH1F *) Raw_NPE_sts->Clone("Raw_Npe_sts");

  TH1F *Eff_tof_match_sts=(TH1F *)eff_Tof_match->Clone("Eff_tof_match_sts");
  TH1F *Eff_poe_sts=(TH1F *)eff_poe->Clone("Eff_poe_sts");

  TH1F *Eff_dedx_sts=(TH1F *)eff_dedx->Clone("Eff_dedx_sts");
  TH1F *Eff_tof_cut_sts=(TH1F *)eff_Tof_cuts->Clone("Eff_tof_cut_sts");
  TH1F *Eff_Tracking_sts=(TH1F *)eff_Tracking->Clone("Eff_Tracking_sts");

  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      Eff_tof_match_sts->SetBinError(ipt+1,0);
      Eff_dedx_sts->SetBinError(ipt+1,0);
      Eff_tof_cut_sts->SetBinError(ipt+1,0);
      Eff_poe_sts->SetBinError(ipt+1,0);
      
      Eff_Tracking_sts->SetBinError(ipt+1,0);
    }
  
  Raw_Npe_sts->Divide(Eff_dedx_sts);
  Raw_Npe_sts->Divide(Eff_tof_match_sts);
   Raw_Npe_sts->Divide(Eff_tof_cut_sts);
   Raw_Npe_sts->Divide(Eff_Tracking_sts);
   Raw_Npe_sts->Divide(Eff_poe_sts);

  NPE_sts=(TH1F *) Raw_Npe_sts->Clone("Npe_sts") ;


}
// sys
void Efficiency_Correction_sys(TH1F *Raw_NPE_sys,TH1F *&NPE_sys)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *Raw_Npe_sys=(TH1F *) Raw_NPE_sys->Clone("Raw_Npe_sys");


  TH1F *Eff_tof_match_sys=(TH1F *)eff_Tof_match->Clone("Eff_tof_match_sys");
  TH1F *Eff_dedx_sys=(TH1F *)eff_dedx->Clone("Eff_dedx_sys");
  TH1F *Eff_tof_cut_sys=(TH1F *)eff_Tof_cuts->Clone("Eff_tof_cut_sys");
  TH1F *Eff_poe_sys=(TH1F *)eff_poe->Clone("Eff_poe_sys");
  TH1F *Eff_Tracking_sys=(TH1F *)eff_Tracking->Clone("Eff_Tracking_sys");

  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      Raw_Npe_sys->SetBinError(ipt+1,0);

    }
  
  Raw_Npe_sys->Divide(Eff_dedx_sys);
  Raw_Npe_sys->Divide(Eff_tof_match_sys);
  Raw_Npe_sys->Divide(Eff_tof_cut_sys);
  Raw_Npe_sys->Divide(Eff_Tracking_sys);
  Raw_Npe_sys->Divide(Eff_poe_sys);
  

  NPE_sys=(TH1F *) Raw_Npe_sys->Clone("Npe_sys") ;


}

void Raw_NPE_spectra_sts(TH1F *inclusive,TH1F *phe,TH1F *&Raw_NPE_sts) 
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *inclusive_sts=(TH1F *) inclusive->Clone("inclusive_sts");


  TH1F *phe_sts=(TH1F *) phe->Clone("phe_sts");


  TH1F *purity_sts= (TH1F *)eff_purity->Clone("purity_sts");
  TH1F *phe_re_sts= (TH1F *)eff_PHE_re->Clone("phe_re_sts");
  
  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      purity_sts->SetBinError(ipt+1,0);
      phe_re_sts->SetBinError(ipt+1,0);
    }
  
  inclusive_sts->Multiply(purity_sts);
  phe_sts->Divide(phe_re_sts);
  Raw_NPE_sts=(TH1F *) inclusive_sts->Clone("Raw_NPE_sts");
  Raw_NPE_sts->Add(phe_sts,-1);
  //  TH1F *S_B_sts=(TH1F *)Raw_NPE_sts->Clone("S_B_sts");




}

void Raw_NPE_spectra_sys(TH1F *inclusive,TH1F *phe,TH1F *&Raw_NPE_sys) 
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

  TH1F *inclusive_sys=(TH1F *) inclusive->Clone("inclusive_sys");
  TH1F *phe_sys=(TH1F *) phe->Clone("phe_sys");

  TH1F *purity_sys= (TH1F *)eff_purity->Clone("purity_sys");
  TH1F *phe_re_sys= (TH1F *)eff_PHE_re->Clone("phe_re_sys");
  
  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      inclusive_sys->SetBinError(ipt+1,0);
      phe_sys->SetBinError(ipt+1,0);
    }
  
  inclusive_sys->Multiply(purity_sys);
  phe_sys->Divide(phe_re_sys);
  Raw_NPE_sys=(TH1F *) inclusive_sys->Clone("Raw_NPE_sys");
  Raw_NPE_sys->Add(phe_sys,-1);
  //  TH1F *S_B_sys=(TH1F *)Raw_NPE_sys->Clone("S_B_sys");


}
/*
void Draw_NPE_PHE_ratio(TH1F *inclusive,TH1F *HT2_inclusive,TH1F *phe,TH1F *HT2_phe,TH1F *&Raw_NPE_sts,TH1F *&Raw_NPE_sys,TH1F *&HT2_Raw_NPE_sts,TH1F *&HT2_Raw_NPE_sys)
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  
  
  TH1F *inclusive_sts=(TH1F *) inclusive->Clone("inclusive_sts");
  TH1F *inclusive_sys=(TH1F *) inclusive->Clone("inclusive_sys");

  TH1F *HT2_inclusive_sts=(TH1F *) HT2_inclusive->Clone("HT2_inclusive_sts");
  TH1F *HT2_inclusive_sys=(TH1F *) HT2_inclusive->Clone("HT2_inclusive_sys");


  TH1F *phe_sts=(TH1F *) phe->Clone("phe_sts");
  TH1F *phe_sys=(TH1F *) phe->Clone("phe_sys");

  TH1F *HT2_phe_sts=(TH1F *) HT2_phe->Clone("HT2_phe_sts");
  TH1F *HT2_phe_sys=(TH1F *) HT2_phe->Clone("HT2_phe_sys");

  TH1F *purity_sts= (TH1F *)eff_purity->Clone("purity_sts");
  TH1F *purity_sys= (TH1F *)eff_purity->Clone("purity_sys");


  TH1F *phe_re_sts= (TH1F *)eff_PHE_re->Clone("phe_re_sts");
  TH1F *phe_re_sys= (TH1F *)eff_PHE_re->Clone("phe_re_sys");

  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      inclusive_sys->SetBinError(ipt+1,0);
      HT2_inclusive_sys->SetBinError(ipt+1,0);
      phe_sys->SetBinError(ipt+1,0);
      HT2_phe_sys->SetBinError(ipt+1,0);

      purity_sts->SetBinError(ipt+1,0);
      phe_re_sts->SetBinError(ipt+1,0);
    }
  //--------------------------------------sts-----------------------------------
  //HT0
  inclusive_sts->Multiply(purity_sts);
  phe_sts->Divide(phe_re_sts);
  Raw_NPE_sts=(TH1F *) inclusive_sts->Clone("Raw_NPE_sts");
  Raw_NPE_sts->Add(phe_sts,-1);
  TH1F *S_B_sts=(TH1F *)Raw_NPE_sts->Clone("S_B_sts");
  S_B_sts->Divide(phe_sts);  
  //HT2
  HT2_inclusive_sts->Multiply(purity_sts);
  HT2_phe_sts->Divide(phe_re_sts);
  HT2_Raw_NPE_sts=(TH1F *) HT2_inclusive_sts->Clone("HT2_Raw_NPE_sts");
  HT2_Raw_NPE_sts->Add(HT2_phe_sts,-1);
  TH1F *HT2_S_B_sts=(TH1F *)HT2_Raw_NPE_sts->Clone("HT2_S_B_sts");
  HT2_S_B_sts->Divide(HT2_phe_sts);  

  //--------------------------------------sys-----------------------------------

  inclusive_sys->Multiply(purity_sys);
  phe_sys->Divide(phe_re_sys);
  Raw_NPE_sys=(TH1F *) inclusive_sys->Clone("Raw_NPE_sys");
  Raw_NPE_sys->Add(phe_sys,-1);
  TH1F *S_B_sys=(TH1F *)Raw_NPE_sys->Clone("S_B_sys");
  S_B_sys->Divide(phe_sys);  
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
  
  TGraphErrors  *gr_S_B_sts= new TGraphErrors(S_B_sts);
  TGraphErrors  *gr_MB2_S_B_sts= new TGraphErrors(HT2_S_B_sts);

  TGraphErrors  *gr_S_B_sys= new TGraphErrors(S_B_sys);
  TGraphErrors  *gr_MB2_S_B_sys= new TGraphErrors(HT2_S_B_sys);
 
  // gr_S_B_sys->Draw();
  // //  return; 

  gr_S_B_sys->SetMarkerStyle(20);
  gr_S_B_sts->SetMarkerStyle(20);

  gr_S_B_sys->SetMarkerColor(4);
  gr_S_B_sts->SetMarkerStyle(4);

  gr_S_B_sys->SetLineColor(2);
  gr_S_B_sys->SetLineWidth(3);

  gr_S_B_sts->SetMarkerStyle(4);

  // gr_S_B_sts->Draw("samePE");
  gr_S_B_sys->Draw("samepE[]");


  
 
  
  // S_B_sys->Draw("same");
  // HT2_S_B_sys->Draw("same");


 

  



 





}
  */
void Draw_Pt_spectra(TH1F * Inclusive,TH1F * Photonic_unlike_pt,TH1F * Photonic_like_pt,TH1F * Photonic_unlike_like_pt,TH1F * Inclusive_ps,TH1F * Photonic_unlike_pt_ps,TH1F * Photonic_like_pt_ps,TH1F * Photonic_unlike_like_pt_ps)
{

  gStyle->SetOptStat(00000);

  //  Inclusive[0]->Draw();
  Int_t bMarkerstyle=20;
  

      Inclusive->SetLineColor(1);
      Inclusive_ps->SetLineColor(1);
      
      Photonic_unlike_like_pt->SetLineColor(2);
      Photonic_unlike_like_pt_ps->SetLineColor(2);

      Inclusive->SetMarkerColor(1);
      Inclusive_ps->SetMarkerColor(1);
      
      Photonic_unlike_like_pt->SetMarkerColor(2);
      Photonic_unlike_like_pt_ps->SetMarkerColor(2);


      Inclusive->SetMarkerStyle(bMarkerstyle);
      Inclusive_ps->SetMarkerStyle(bMarkerstyle);

      Photonic_unlike_like_pt->SetMarkerStyle(bMarkerstyle+4);
      Photonic_unlike_like_pt_ps->SetMarkerStyle(bMarkerstyle+4);

      
  TLegend *legend = new TLegend(0.7,0.65,0.85,0.85);
  legend->AddEntry(Inclusive,"MB inclusive e","lp");
  legend->AddEntry(Photonic_unlike_like_pt,"MB photonic e","lp");
  legend->SetTextFont(62);
  legend->SetTextSize(0.035);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);   
  TCanvas *c4=new TCanvas("c4","",1200,1000);
  gPad->SetLogy();
  TH2F *h4=new TH2F("h4","",100,1,5,100,0.1,1e5);
  h4->Draw();
  h4->GetXaxis()->SetTitle("p_{T} GeV/c");
  h4->GetYaxis()->SetTitle("Counts");
  Inclusive->Draw("same");
  Photonic_unlike_like_pt->Draw("same");
  Inclusive->Draw("same");
  Photonic_unlike_like_pt->Draw("same");
   legend->Draw();
  c4->SaveAs("pT_raw_spectra_after_prescale.pdf");

  TCanvas *c3=new TCanvas("c2","",1200,1000);
  gPad->SetLogy();
  TH2F *h3=new TH2F("h3","",100,1,5,100,1,3e7);
  h3->Draw();
  h3->GetXaxis()->SetTitle("p_{T} GeV/c");
  h3->GetYaxis()->SetTitle("Counts");

  Inclusive_ps->Draw("same");
  Photonic_unlike_like_pt_ps->Draw("same");
  Inclusive_ps->Draw("same");
  Photonic_unlike_like_pt_ps->Draw("same");
  legend->Draw();
  c3->SaveAs("pT_raw_spectra_before_prescale.pdf");
}
  
void Draw_Efficiency()
{
  // eff_Tof_match->Draw();
  // eff_dedx->Draw();
  // eff_Tof_cuts->Draw();
  // eff_Tracking->Draw();
  // eff_PHE_re->Draw();
  // eff_purity->Draw();

  TH1F *Eff_Tof_match=(TH1F *)eff_Tof_match->Clone("Eff_Tof_match");
  TH1F *Eff_dedx=(TH1F *)eff_dedx->Clone("Eff_dedx");
  TH1F *Eff_Tof_cut=(TH1F *)eff_Tof_cuts->Clone("Eff_Tof_cut");
  TH1F *Eff_Tracking=(TH1F *)eff_Tracking->Clone("Eff_Tracking");
  TH1F *Eff_PHE_re=(TH1F *)eff_PHE_re->Clone("Eff_PHE_re");
  TH1F *Eff_purity=(TH1F *)eff_purity->Clone("Eff_purity");

  Eff_Tof_match->SetMarkerStyle(20);
  Eff_dedx->SetMarkerStyle(20);
  Eff_Tof_cut->SetMarkerStyle(20);

  Eff_Tracking->SetMarkerStyle(20);
  Eff_PHE_re->SetMarkerStyle(20);
  Eff_purity->SetMarkerStyle(20);

  Eff_Tof_match->SetMarkerColor(1);
  Eff_dedx->SetMarkerColor(2);
  Eff_Tof_cut->SetMarkerColor(3);

  Eff_Tracking->SetMarkerColor(5);
  Eff_PHE_re->SetMarkerColor(6);
  Eff_purity->SetMarkerColor(9);

  Eff_Tof_match->SetLineColor(1);
  Eff_dedx->SetLineColor(2);
  Eff_Tof_cut->SetLineColor(3);

  Eff_Tracking->SetLineColor(5);
  Eff_PHE_re->SetLineColor(6);
  Eff_purity->SetLineColor(8);


  gStyle->SetOptStat(0);
  TCanvas *c2=new TCanvas("c2","",1200,800);
  c2->cd();

  TH2F *h2=new TH2F("h2","",100,0.2,4,100,-0.1,1.8);
  h2->Draw();
  h2->GetYaxis()->SetTitle("efficiency");
  h2->GetXaxis()->SetTitle("p_{T} (Gev/c)");

  Eff_Tof_match->Draw("samePE");
  Eff_dedx->Draw("samePE");

  Eff_Tof_cut->Draw("samePE");
  Eff_Tracking->Draw("samePE");
  Eff_PHE_re->Draw("samePE");
  Eff_purity->Draw("samePE");
  
  TLegend *legend = new TLegend(0.2,0.65,0.35,0.85);
  //  TLegend *legend = new TLegend(0.15,0.65,0.55,0.85);
  
  legend->AddEntry(Eff_Tof_match,"TOF macth","lp");
  legend->AddEntry(Eff_dedx,"dEdx","lp");
  legend->AddEntry(Eff_Tof_cut,"Tof Cuts","lp");

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
