/*****************************************************
* this macro for the trigger Efficiency 
*
*   4/7/2014 by Xiaozhi
*
*
*****************************************************
*/

#include <iostream>
#include <fstream> 
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

void TrgEfficiency()
{
  TH2F::SetDefaultSumw2();
  TH1F::SetDefaultSumw2();

  gStyle->SetOptStat(false);
  gStyle->SetOptTitle(true);
  
  TFile * file=new TFile("Trigger_efficiency.root","READ");

  TH1F * TrigEfficiency_HT0=new TH1F("TrigEfficiency_HT0","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F * TrigEfficiency_HT2=new TH1F("TrigEfficiency_HT2","",NpT_bins_run12_HT,pt_run12_HT);
  
  TH1F *  m_DSMadc=(TH1F*)file->Get("m_DSMadc");  
  TH1F *  m_DSMadc11_cut=(TH1F* )file->Get("m_DSMadc11_cut");
  TH1F *  m_DSMadc18_cut=(TH1F*)file->Get("m_DSMadc18_cut");
  
  TH1F *  hh=(TH1F *) m_DSMadc->Rebin(NpT_bins_run12_HT,"hh",pt_run12_HT);
  TH1F *  hh_11=(TH1F *) m_DSMadc11_cut->Rebin(NpT_bins_run12_HT,"hh_11",pt_run12_HT);
  TH1F *  hh_18=(TH1F *) m_DSMadc18_cut->Rebin(NpT_bins_run12_HT,"hh_18",pt_run12_HT);
  
  TGraphAsymmErrors *h11=new TGraphAsymmErrors(hh_11,hh,"N");
  TGraphAsymmErrors *h18=new TGraphAsymmErrors(hh_18,hh,"N");
  

  for(Int_t i=0;i<NpT_bins_run12_HT;i++)
    {
      Double_t x=0,y=0,x_err=0,y_err=0;
      h11->GetPoint(i,x,y);
      y_err=h11->GetErrorY(i);
      TrigEfficiency_HT0->SetBinContent(i+1,y);
      TrigEfficiency_HT0 ->SetBinError(i+1,y_err);

      h18->GetPoint(i,x,y);
      y_err=h18->GetErrorY(i);
      
      TrigEfficiency_HT2->SetBinContent(i+1,y);
      TrigEfficiency_HT2->SetBinError(i+1,y_err);
    }
  
 
   TCanvas *c2=new TCanvas("c2","",1200,1000);

   gPad->SetLogy();  
   TH2F * h2=new TH2F("h2","",10,2,10,10,1e-4,10);
   h2->GetXaxis()->SetTitle("P_{T} GeV/c");
   h2->GetYaxis()->SetTitle("Trigger efficiency");
   h2->Draw();
  
 
   TrigEfficiency_HT0->SetMarkerStyle(20);
   TrigEfficiency_HT0->SetMarkerColor(2);

   TrigEfficiency_HT2->SetMarkerStyle(20);
   TrigEfficiency_HT2->SetMarkerColor(4);

   TrigEfficiency_HT0->Draw("SAMEPE1");
   TrigEfficiency_HT2->Draw("SAMEPE1");

   

   
   TLegend *legend = new TLegend(0.6,0.4,0.85,0.6);
   legend->AddEntry( TrigEfficiency_HT0," HT0 Trigger efficiency","lp");
   legend->AddEntry( TrigEfficiency_HT2," HT2 Trigger efficiency","lp");
   
   legend->SetBorderSize(0);
   legend->SetTextSize(0.035);
   legend->SetFillStyle(0);
   legend->SetTextFont(62);
   legend->Draw();   

     c2->SaveAs("Trigger_efficiency.pdf");   
   TFile *File=new TFile("TrigEfficiency_HT.root","recreate");
   TrigEfficiency_HT0->Write();
   TrigEfficiency_HT2->Write();
   File->Close();
   

}
