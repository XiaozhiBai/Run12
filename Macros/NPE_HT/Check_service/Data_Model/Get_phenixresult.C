#include <iostream>
#include<iomanip>
#include <fstream>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

// #include "../mBinning_HT.h"
// #include "../mBinning_MB.h"

TH1F *Phex=new TH1F();

const Int_t Npt=28;
Double_t ptbin[Npt+1]={0.325,0.375,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,6,7,8,9};
void Get_phenixresult()
{
  TH1F *Phenix_sts=new TH1F("Phenix_sts","",Npt,ptbin);
  TH1F *Phenix_sys=new TH1F("Phenix_sys","",Npt,ptbin);
  
  TH1F *Phenix_fonll=new TH1F("Phenix_fonll","",Npt,ptbin);


  
  
  ifstream file;
  file.open("phenixresult.dat",ios::in);
  Double_t a=0,b=0,c=0,d=0;
  while(!file.eof())
    { 
      file>>a>>b>>c>>d;
      cout<<a<<",";// "  "<<b<<"  "<<c<<endl;
      Phenix_sts->SetBinContent(Phenix_sts->FindBin(a),b);
      Phenix_sts->SetBinError(Phenix_sts->FindBin(a),c);
      
      Phenix_sys->SetBinContent(Phenix_sys->FindBin(a),b);
      Phenix_sys->SetBinError(Phenix_sys->FindBin(a),d);
    }

  // Phenix_sts->Draw();
  
  // return;
  
  TFile *file_data_fonll=new TFile("fonll200gev.root","READ");


  TGraphErrors *gFONLLu=( TGraphErrors *) file_data_fonll->Get("gFONLLu");
  TGraphErrors *gFONLLc=( TGraphErrors *) file_data_fonll->Get("gFONLLc");
  TGraphErrors *gFONLLl=( TGraphErrors *) file_data_fonll->Get("gFONLLl");

  TGraphErrors *gR_FONLLu= (TGraphErrors *) file_data_fonll->Get("gR_FONLLu");
  TGraphErrors *gR_FONLLc= (TGraphErrors *) file_data_fonll->Get("gR_FONLLc");
  TGraphErrors *gR_FONLLl= (TGraphErrors *) file_data_fonll->Get("gR_FONLLl");
  // gFONLLc->Draw();
  // return;
  
  for(int ipt=0;ipt<Npt;ipt++)
    {
      Phenix_fonll->SetBinContent(ipt+1,gFONLLc->Eval(Phenix_fonll->GetBinCenter(ipt+1)));
      Phenix_fonll->SetBinError(ipt+1,0);
      
    }
  // Phenix_fonll->Draw();
  
  //  return;

  TH1F *  Phenix_fonll_ratio_sts=(TH1F *) Phenix_sts->Clone("Phenix_fonll_ratio_sts");
  TH1F *  Phenix_fonll_ratio_sys=(TH1F *) Phenix_sys->Clone("Phenix_fonll_ratio_sys"); 
  

    Phenix_fonll->Draw();
   
    gFONLLu->Draw("same");
    gFONLLc->Draw("same");
    gFONLLl->Draw("same");

  
   Phenix_fonll_ratio_sts->Divide(Phenix_fonll);
   Phenix_fonll_ratio_sys->Divide(Phenix_fonll);

   Phenix_fonll_ratio_sts->Draw();
   gR_FONLLu->Draw("same");
   gR_FONLLc->Draw("same");
   gR_FONLLl->Draw("same");
   
  
  
  // Phenix_sts->Draw();

   TFile *file_Pheix=new TFile("Npe_Phenix.root","RECREATE");
   Phenix_sts->Write();
      
   Phenix_sys->Write();

   file_Pheix->Close();

}
