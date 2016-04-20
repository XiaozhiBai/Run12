
/************************************************

  This macro is for the photonic electron reconstruction efficiency
  xiaozhi bai
  xiaozhi@uic.edu
Thu Jul  9 11:47:12 CDT 2015
 ************************************************/
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
#include "../mBinning_MB.h"

using namespace std;

TH2F * InvMass_Pt[3];

TH1F * InvMass[3];
TH1F * mh1Pt[3];

TString histname[3]={"Gamma_conversion","Pi0Dalitz","EtaDalitzDecay"};
void Draw_contribution(TH1F *,TH1F *,TH1F * );
//TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
//void setpad(TPad *,float , float , float , float );
void Photonic_RE_relative_contribution(TString Gamma_conversion_MB="Root_File_PHE/Gamma_MB.root",TString Pi0Dalitz_MB="Root_File_PHE/Pi0_dalitz_MB.root",TString EtaDalitz_MB="Root_File_PHE/Eta_dalitz_MB.root") 
{
  char buf[1024];
  gStyle->SetOptStat(0);                                                    
  gStyle->SetOptFit(kTRUE);                                                    
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1,"Y");
  //  gStyle->SetLabelSize(0.08,"xy");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13); 
  
  TH1F::SetDefaultSumw2();  

  TFile * inFile_Gamma_MB= new TFile(Gamma_conversion_MB); //gammma  
  TFile * inFile_Pi0Dalitz_MB= new TFile(Pi0Dalitz_MB); //  Pi0
  TFile * inFile_EtaDalitz_MB= new TFile(EtaDalitz_MB);  //eta
  
  InvMass_Pt[0]=(TH2F *) inFile_Gamma_MB->Get("mInv_massVspt");
  InvMass_Pt[1]=(TH2F *) inFile_Pi0Dalitz_MB->Get("mInv_massVspt");
  InvMass_Pt[2]=(TH2F *) inFile_EtaDalitz_MB->Get("mInv_massVspt");



  TH1F * mh1_Gamma=(TH1F *) InvMass_Pt[0]->Rebin(NpT_bins_run12_MB,"mh1_Gamma",pt_run12_MB); 
  TH1F * mh1_pi0=(TH1F *) InvMass_Pt[1]->Rebin(NpT_bins_run12_MB,"mh1_pi0",pt_run12_MB); 
  TH1F * mh1_Eta=(TH1F *) InvMass_Pt[2]->Rebin(NpT_bins_run12_MB,"mh1_Eta",pt_run12_MB); 


  mh1_Gamma->Scale(1./(291530*100));
  mh1_pi0->Scale(1./((452725*5)));
  mh1_Eta->Scale(1./(242066*5));

  mh1_Gamma->Scale(1);
  mh1_pi0->Scale(0.01174);
  mh1_Eta->Scale(0.0069);

  
  mh1_Gamma->Scale(2);
  mh1_pi0->Scale(3);
  mh1_Eta->Scale(3);





  TH1F *total=new TH1F("total","",NpT_bins_run12_MB,pt_run12_MB);
  total->Sumw2();
 
  total->Add(mh1_Gamma);
  total->Add(mh1_pi0);
  total->Add(mh1_Eta);
  

  mh1_Gamma->Divide(total);
  mh1_pi0->Divide(total);
  mh1_Eta->Divide(total);
  

  Draw_contribution(mh1_Gamma,mh1_pi0, mh1_Eta);
}

void Draw_contribution(TH1F * mh1_Gamma,TH1F *mh1_pi0,TH1F * mh1_Eta)
{
  
  char buf[1024];
  gStyle->SetOptStat(00000);    
  TCanvas *c2 = new TCanvas("c2","",0,0,800,600);
  c2->cd(); 

  TH2F * hh=new TH2F("hh","",100,1,4,100,0,1);
  
  hh->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hh->GetYaxis()->SetTitle("Relative Contribution");
  hh->GetYaxis()->CenterTitle(true);
  hh->GetXaxis()->CenterTitle(true);
  hh->Draw();

  mh1_Gamma->SetLineColor(1);
  mh1_pi0->SetLineColor(2);
  mh1_Eta->SetLineColor(3);

  mh1_Gamma->SetMarkerColor(1);
  mh1_pi0->SetMarkerColor(2);
  mh1_Eta->SetMarkerColor(3);
 
  mh1_Gamma->SetMarkerStyle(20);
  mh1_pi0->SetMarkerStyle(20);
  mh1_Eta->SetMarkerStyle(20);
  

  mh1_Gamma->Draw("samePE");
  mh1_pi0->Draw("samePL");
  mh1_Eta->Draw("samePL");
  
  TLegend *legend_MB = new TLegend(0.2,0.75,0.4,0.85);
  legend_MB->AddEntry(mh1_Gamma, " HT  #gamma #rightarrow e^{+} e^{-}","lp"); 
  legend_MB->AddEntry(mh1_pi0, " HT  #pi #rightarrow #gamma  e^{+} e^{-}","lp"); 
  legend_MB->AddEntry(mh1_Eta, " HT  #eta #rightarrow #gamma  e^{+} e^{-}","lp"); 
  

  legend_MB->SetTextSize(0.03); 
  legend_MB->SetBorderSize(0);
  legend_MB->SetFillStyle(0);
  legend_MB->SetTextFont(62);
  legend_MB->Draw();

  c2->SaveAs("ralative_phe.pdf");

  
  TFile *file=new TFile("PHE_Contribution_MB.root","RECREATE");
  
  mh1_Gamma->Write();
  mh1_pi0->Write();
  mh1_Eta->Write();
  
  file->Close();

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
