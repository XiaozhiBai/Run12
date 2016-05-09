#include <iostream>
#include<iomanip>
#include <fstream>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TLine.h"

#include "TMath.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "../mBinning_MB.h"


void Get_Omiga_phi()
{


  TH1F  *Omiga_NPE_CC_sts=new TH1F("Omiga_NPE_CC_sts","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F  *Phi_NPE_CC_sts=new TH1F("Phi_NPE_CC_sts","",NpT_bins_run12_MB,pt_run12_MB);

  TH1F  *Omiga_NPE_CC_sys=new TH1F("Omiga_NPE_CC_sys","",NpT_bins_run12_MB,pt_run12_MB);
  TH1F  *Phi_NPE_CC_sys=new TH1F("Phi_NPE_CC_sys","",NpT_bins_run12_MB,pt_run12_MB);


ifstream file_2;
file_2.open("omegaphiccratio.dat",ios::in);

 Double_t a1=0,b1=0,c1=0,d1=0,e1=0,f1=0;
 while(!file_2.eof())
   { 
     file_2>>a1>>b1>>c1>>d1>>e1>>f1;
     
     cout<<a1<<"  "<<b1<<"  "<<c1<<"  "<<d1 <<e1<<"  " <<f1 <<endl;
     Omiga_NPE_CC_sys->SetBinContent(Omiga_NPE_CC_sys->FindBin(a1),c1);
     Omiga_NPE_CC_sys->SetBinError(Omiga_NPE_CC_sys->FindBin(a1),d1);
     
     Phi_NPE_CC_sys->SetBinContent(Phi_NPE_CC_sys->FindBin(a1),e1);
     Phi_NPE_CC_sys->SetBinError(Phi_NPE_CC_sys->FindBin(a1),f1);


     Omiga_NPE_CC_sts->SetBinContent(Omiga_NPE_CC_sts->FindBin(a1),c1);
     Omiga_NPE_CC_sts->SetBinError(Omiga_NPE_CC_sts->FindBin(a1),0);
     
     Phi_NPE_CC_sts->SetBinContent(Phi_NPE_CC_sts->FindBin(a1),e1);
     Phi_NPE_CC_sts->SetBinError(Phi_NPE_CC_sts->FindBin(a1),0);

     
   }
 Phi_NPE_CC_sts->Draw();
 Phi_NPE_CC_sys->Draw("same");



 TFile *file_npe=new TFile("run12_Npe_MB.root","READ");

 TH1F *run12_MB_NPE_sts=(TH1F *) file_npe->Get("run12_MB_NPE_sts");
 TH1F *run12_MB_NPE_sys=(TH1F *) file_npe->Get("run12_MB_NPE_sys");
 
 
 /*
 gStyle->SetOptStat(00);

 TCanvas *c4=new TCanvas("c4","",800,600);
 TH2F *hhh=new TH2F("hhh","",10,0,4,10,0,0.15);
 hhh->GetXaxis()->SetTitle("p_{T} GeV/c");
 hhh->GetYaxis()->SetTitle("Ratio");

 hhh->Draw();
 Omiga_NPE_CC->SetMarkerStyle(20);
 Phi_NPE_CC->SetMarkerStyle(20);
 Omiga_NPE_CC->SetMarkerColor(2);
 Phi_NPE_CC->SetMarkerColor(3);
 Omiga_NPE_CC->SetLineColor(2);
 Phi_NPE_CC->SetLineColor(3);

 Omiga_NPE_CC->Draw("same");
 Phi_NPE_CC->Draw("same");

 TLegend *legend  = new TLegend(0.4,0.65,0.8,0.85);
  legend ->AddEntry(Omiga_NPE_CC,"#omega #rightarrow e/ c #rightarrow e","lpe");
  legend ->AddEntry(Phi_NPE_CC,"#phi #rightarrow e/ c #rightarrow e","lpe");
  legend ->SetBorderSize(0);
  legend ->SetTextSize(0.04);
  legend ->SetFillStyle(0);
  legend ->SetTextFont(62);
  legend ->Draw("same");

  c4->SaveAs("Omega_phi.pdf");
  // return;
  */

 // run12_MB_NPE_sys->Draw();
 // return;
 
 TH1F*  Omega_sts=(TH1F *) run12_MB_NPE_sts->Clone("Omega_sts");
 TH1F*  Omega_sys=(TH1F *) run12_MB_NPE_sys->Clone("Omega_sys");
 
 TH1F*  Phi_sts=(TH1F *) run12_MB_NPE_sts->Clone("Phi_sts");
 TH1F*  Phi_sys=(TH1F *) run12_MB_NPE_sys->Clone("Phi_sys");

 
 Omega_sts->Multiply(Omiga_NPE_CC_sts);
 Omega_sys->Multiply(Omiga_NPE_CC_sys);
 
 Phi_sts->Multiply(Omiga_NPE_CC_sts);
 Phi_sys->Multiply(Omiga_NPE_CC_sys);
 


 

 // TCanvas *c4=new TCanvas("c4","",800,600);
 //  TH2F *hh=new TH2F("hh","",10,0,4,10,0,0.15);
 // hh->GetXaxis()->SetTitle("p_{T} GeV/c");
 // hh->GetYaxis()->SetTitle("Ratio");

 //->Draw();
 // Phi->Draw();
 // Omiga->Draw("same");
 // Omiga_NPE_CC->SetMarkerStyle(20);
 // Phi_NPE_CC->SetMarkerStyle(20);
 // Omiga_NPE_CC->SetMarkerColor(2);
 // Phi_NPE_CC->SetMarkerColor(3);
 // Omiga_NPE_CC->SetLineColor(2);
 // Phi_NPE_CC->SetLineColor(3);
 
 // Omiga_NPE_CC->Draw("same");
 // Phi_NPE_CC->Draw("same");

 // TLegend *legend  = new TLegend(0.4,0.65,0.8,0.85);
 //  legend ->AddEntry(Omiga_NPE_CC,"#omega #rightarrow e/ c #rightarrow e","lpe");
 //  legend ->AddEntry(Phi_NPE_CC,"#phi #rightarrow e/ c #rightarrow e","lpe");
 //  legend ->SetBorderSize(0);
 //  legend ->SetTextSize(0.04);
 //  legend ->SetFillStyle(0);
 //  legend ->SetTextFont(62);
 //  legend ->Draw("same");


 

 
 TFile * file_3=new TFile("Omega_phi.root","RECREATE");


 Omega_sts->Write();
 Phi_sts->Write();

 Omega_sys->Write();
 Phi_sys->Write();
 
 // Omiga_NPE_CC->Write();
 // Phi_NPE_CC->Write();


 file_3->Close(); 

 

 return;



// TH2F  *hh=new TH2F ("hh","",100,0,5,100,0,0.2);
// hh->Draw();

// Omiga_NPE_CC->Draw("same");
// //Phi_NPE_CC->Draw("same");





}
