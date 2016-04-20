



#include <iostream>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "../mBinning_HT.h"
using namespace std;

void Draw_efficiency(TH1F *);
void Tracking_efficiency()
{

  TH1F * mh1Track_noCut;
  TH1F * mh1Track_Cut20;
  TH1F * mh1Track_Cut25;
  TH1F * Tracking_efficiency_HT=new TH1F("Tracking_efficiency_HT","",NpT_bins_run12_HT,pt_run12_HT);
  
  TFile * inFile= new TFile("Tracking_eff.root","READ");
  TH1F * mh1Track_MC= (TH1F *) inFile->Get("mNTrack_nocuts");
  TH1F * mh1Track_MC20= (TH1F *) inFile->Get("mNTrack_cut");
  TH1F * mh1Track_MC25= (TH1F *) inFile->Get("mNTrack_cut_25");

  TH1F * mh1_Track=(TH1F *) mh1Track_MC->Rebin(NpT_bins_run12_HT,"mh1_Track",pt_run12_HT);
  TH1F * mh1_Track_20Cut=(TH1F *) mh1Track_MC20->Rebin(NpT_bins_run12_HT,"mh1_Track_20Cut",pt_run12_HT);
  TH1F * mh1_Track_25Cut=(TH1F *) mh1Track_MC25->Rebin(NpT_bins_run12_HT,"mh1_Track_25Cut",pt_run12_HT);

  TGraphAsymmErrors *Tracking_Ef=new TGraphAsymmErrors(mh1_Track_20Cut,mh1_Track,"N");
  TGraphAsymmErrors *Tracking_Ef_sys=new TGraphAsymmErrors(mh1_Track_25Cut,mh1_Track,"N");
 
  Tracking_Ef->Draw();
  Tracking_Ef_sys->Draw("same");

  for(Int_t i=0;i<NpT_bins_run12_HT;i++)
    {
      Double_t x_20=0,y_20=0,y_err_20=0;
      Double_t x_25=0,y_25=0,y_err_25=0;

      Tracking_Ef->GetPoint(i,x_20,y_20);
      y_err_20=Tracking_Ef->GetErrorY(i);
      y_err_20=Tracking_Ef_sys->GetErrorY(i);

      Tracking_Ef_sys->GetPoint(i,x_25,y_25);
      y_err_25=Tracking_Ef_sys->GetErrorY(i);
      y_err_25=Tracking_Ef_sys->GetErrorY(i);
      

      Tracking_efficiency_HT->SetBinContent(i+1,y_20);
      Tracking_efficiency_HT->SetBinError(i+1,y_err_20+fabs(y_20-y_25));
    }
   Draw_efficiency(Tracking_efficiency_HT);  
}
void Draw_efficiency(TH1F * efficiency)
{
  gStyle->SetOptStat(0);  
  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd();
  TH2F * hh=new TH2F("hh","",10,2,10,10,0,1);
  hh->GetXaxis()->SetTitle("P_{T} GeV/c");
  hh->GetYaxis()->SetTitle("Tracking efficiency");
  
  hh->Draw();

  efficiency->SetMarkerStyle(20);
  efficiency->SetMarkerColor(4);
  efficiency->Draw("samePE1");
  

  // TLegend *legend = new TLegend(0.3,0.2,0.8,0.4); 
  // legend->AddEntry(Tracking_high,"Tracking efficiency from High Tower Trigger ","lpe"); 
  // legend->AddEntry(efficiency,"Tracking efficiency from VPDMB Trigger ","lpe"); 
  // legend->SetTextSize(0.03); 
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);   
  // legend->Draw("same");
  
  c2->SaveAs("Tracking_efficiency_HT.pdf");
  
  TFile *file=new TFile("Tracking_efficiency_HT.root","RECREATE");
  efficiency->Write();
  file->Close();
}



