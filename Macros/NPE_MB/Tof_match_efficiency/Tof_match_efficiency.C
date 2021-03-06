
// By xiaozhi
#include <iostream>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include<vector>
#include "../mBinning_MB.h"

using namespace std;

const int mh2Nhist= 1;

TH2F *mh2hist_unlike_pass[mh2Nhist];
TH2F *mh2hist_like_pass[mh2Nhist];

TH1F *mh1hist_pT_unlike_pass[mh2Nhist];
TH1F *mh1hist_pT_like_pass[mh2Nhist];
TH1F *mh1hist_pT_unlike_like_pass[mh2Nhist];

TH2F *mh2hist_unlike_total[mh2Nhist];
TH2F *mh2hist_like_total[mh2Nhist];
TH1F *mh1hist_pT_unlike_total[mh2Nhist];
TH1F *mh1hist_pT_like_total[mh2Nhist];
TH1F *mh1hist_pT_unlike_like_total[mh2Nhist];

TH3F *mh3hist_unlike_partnernSigmaE_pt[mh2Nhist];
TH3F *mh3hist_like_partnernSigmaE_pt[mh2Nhist];

TH1F *mh1hist_unlike[mh2Nhist][NpT_bins_run12_MB];
TH1F *mh1hist_like[mh2Nhist][NpT_bins_run12_MB];
TH1F *mh1hist_unlike_like[mh2Nhist][NpT_bins_run12_MB];



TString mh2HistName_Photonic_unlike_total[mh2Nhist]={"mh2Part_Ele_MassVspT_noTofMatchcut_unlikeTrg2"};
TString mh2HistName_Photonic_like_total[mh2Nhist]={"mh2Part_Ele_MassVspT_noTofMatchcut_likeTrg2"};

TString mh2HistName_Photonic_unlike_pass[mh2Nhist]={"mh2Part_Ele_MassVspT_TofMatchcut_unlikeTrg2"};
TString mh2HistName_Photonic_like_pass[mh2Nhist]={"mh2Part_Ele_MassVspT_TofMatchcut_likeTrg2"};

TString mh1_TitleY[mh2Nhist]={"Counts"};
TString mh1_TitleX[mh2Nhist]={"Mass ee (GeV/c^2)"};



char buf[1024];
void Tof_match_efficiency()
{
  cout<< sizeof(pt_run12_MB)/sizeof(pt_run12_MB[0])<<endl;

  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  

   gStyle->SetOptFit(1111);
   gStyle->SetTitleSize(0.05,"XY");
  // gStyle->SetTitleFontSize(0.06);
   gStyle->SetTitleOffset(1,"X");
   gStyle->SetTitleOffset(1.4,"Y");
  // gStyle->SetPadTopMargin(0.13);
  // gStyle->SetPadRightMargin(0.02);
  // gStyle->SetPadBottomMargin(0.13);
   gStyle->SetPadLeftMargin(0.15); 
  // gStyle->SetEndErrorSize(4);
  

     TFile *File=new TFile("../RootFile/Root_File_5_2/hist_5_2.root","READ");
  
  for(Int_t ihist=0;ihist<mh2Nhist;ihist++)
    {
      mh2hist_unlike_pass[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_unlike_pass[ihist]);
      mh2hist_like_pass[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_like_pass[ihist]);
      mh2hist_unlike_total[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_unlike_total[ihist]);
      mh2hist_like_total[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_like_total[ihist]);
      
      sprintf(buf,"unlike_passTrg%i",ihist);
      mh1hist_pT_unlike_pass[ihist]=(TH1F *)mh2hist_unlike_pass[ihist]->ProjectionX(buf,1,10,"");
      sprintf(buf,"like_passTrg%i",ihist);
      mh1hist_pT_like_pass[ihist]=(TH1F *)mh2hist_like_pass[ihist]->ProjectionX(buf,1,10,"");
      sprintf(buf,"unlike_like_passTrg%i",ihist);
      mh1hist_pT_unlike_like_pass[ihist]=(TH1F *)mh1hist_pT_unlike_pass[ihist]->Clone(buf);
      mh1hist_pT_unlike_like_pass[ihist]->Add(mh1hist_pT_like_pass[ihist],-1);

      sprintf(buf,"unlike_totalTrg%i",ihist);
      mh1hist_pT_unlike_total[ihist]=(TH1F *)mh2hist_unlike_total[ihist]->ProjectionX(buf,1,10,"");
      sprintf(buf,"like_totalTrg%i",ihist);
      mh1hist_pT_like_total[ihist]=(TH1F *)mh2hist_like_total[ihist]->ProjectionX(buf,1,10,"");
      sprintf(buf,"unlike_like_totalTrg%i",ihist);
      mh1hist_pT_unlike_like_total[ihist]=(TH1F *)mh1hist_pT_unlike_total[ihist]->Clone(buf);
      mh1hist_pT_unlike_like_total[ihist]->Add(mh1hist_pT_like_total[ihist],-1);
      
    }

  // mh1hist_pT_unlike_like_total[0]->Draw();
  // mh1hist_pT_unlike_like_pass[0]->Draw("same");
  // mh1hist_pT_unlike_like_pass[0]->SetLineColor(2);

  // return;

  TH1F *Pt_pass_rebin=(TH1F *)  mh1hist_pT_unlike_like_pass[0]->Rebin(NpT_bins_run12_MB,"Pt_pass_rebin",pt_run12_MB);
  TH1F *Pt_total_rebin=(TH1F *) mh1hist_pT_unlike_like_total[0]->Rebin(NpT_bins_run12_MB,"Pt_total_rebin",pt_run12_MB);

  // for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
  //   {
  //     cout<< " psss " <<Pt_pass_rebin->GetBinContent(ipt+1)<< " total "<<Pt_total_rebin->GetBinContent(ipt+1)<< " ratio"<< Pt_pass_rebin->GetBinContent(ipt+1)/Pt_total_rebin->GetBinContent(ipt+1)<<endl;;
  //   }
  
  TGraphAsymmErrors *efficiency=new TGraphAsymmErrors(Pt_pass_rebin,Pt_total_rebin,"N");
  
  TH1F *efficiency_tof_match=new TH1F("efficiency_tof_match","",NpT_bins_run12_MB,pt_run12_MB);
  
  for(Int_t i=0;i<NpT_bins_run12_MB;i++)
    {
      Double_t x=0,y=0,x_err=0,y_err=0;
      efficiency->GetPoint(i,x,y);
      y_err=efficiency->GetErrorY(i);
      efficiency_tof_match->SetBinContent(efficiency_tof_match->FindBin(x),y);
      efficiency_tof_match->SetBinError(i+1,y_err);

    }
  // return;
  efficiency_tof_match->SetMarkerStyle(20);
  efficiency_tof_match->SetMarkerSize(1);
  efficiency_tof_match->SetMarkerColor(2);
  
  TH2F *h2=new TH2F("h2","",100,0.2,4,100,0,1);
  h2->GetXaxis()->SetTitle("P_{T} GeV/c");
  h2->GetYaxis()->SetTitle("Tof Match efficiency"); 
  
  gStyle->SetOptStat(0);     
  TCanvas *c2= new TCanvas("c2","",800,600);
  efficiency_tof_match->Draw("P");
  h2->Draw();
  efficiency_tof_match->Draw("same");
  c2->SaveAs("Tof_macth_efficiency.pdf");
  
  gStyle->SetOptStat(1111);     
  TFile *file_tof_match=new TFile("TOF_Match_efficiency.root","RECREATE");
  efficiency_tof_match->Write();
  file_tof_match->Close();
  


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




