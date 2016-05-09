
#include <fstream>
#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

Int_t const NHIST=2;

TString Histname_constrain[NHIST]={"Mean_Fit","Sigma_Fit"};
TString Histname_fit[NHIST]={"mh1Electron_Mean","mh1Electron_Sigma"};

  TH1F *mh1Constrain[NHIST];
TH1F *mh1Fit[NHIST];




char buf[1024];

void Draw_Mean_Sigma()
{
  gStyle->SetOptStat(0);
  
  TFile *File_constrain=new TFile("nSigmaE_calibration_primary.root","READ");
  TFile *File_Fit=new TFile("Mean_Sigma_Hadron.root","READ");

  for(Int_t ihist=0;ihist<NHIST;ihist++)
    {

      mh1Constrain[ihist]=(TH1F *)File_constrain->Get(Histname_constrain[ihist]);
      mh1Fit[ihist]=(TH1F *) File_Fit->Get(Histname_fit[ihist]);

      mh1Constrain[ihist]->SetMarkerStyle(24);
      mh1Fit[ihist]->SetMarkerStyle(20);
      

      mh1Constrain[ihist]->SetMarkerColor(ihist+1);
      mh1Fit[ihist]->SetMarkerColor(ihist+1);

      mh1Constrain[ihist]->SetLineColor(ihist+1);
      mh1Fit[ihist]->SetLineColor(ihist+1);
      

      
    }

  mh1Fit_Constant=(TH1F *) File_Fit->Get("mh1Electron_Constant");



  TCanvas *c5=new TCanvas("c5","",800,600);
  TH2F *hh_mean=new TH2F("hh_mean","",20,1.2,4,20,-1,2);
  hh_mean->GetXaxis()->SetTitle("p_{T} (GeV)");
  hh_mean->GetYaxis()->SetTitle("Nsigma E Mean and sigma ");

  hh_mean->Draw();
  for(Int_t ihist=0;ihist<NHIST;ihist++)
    {
      mh1Constrain[ihist]->Draw("samePE");
      mh1Fit[ihist]->Draw("samePE");
    }


  
    TLegend *legend  = new TLegend(0.3,0.75,0.5,0.85);
    legend ->AddEntry(mh1Constrain[0],"Photonic  electron  #sigma_{e} mean","lp");
    legend ->AddEntry(mh1Constrain[1],"Photonic  electron  #sigma_{e} sigma","lp");

    legend ->SetBorderSize(0);
    legend ->SetTextSize(0.03);
    legend ->SetFillStyle(0);


    TLegend *legend_fit  = new TLegend(0.3,0.65,0.5,0.75);

    legend_fit ->AddEntry(mh1Fit[0]," #sigma_{e} mean from fit","lp");
    legend_fit ->AddEntry(mh1Fit[1]," #sigma_{e} sigma from fit","lp");

    legend_fit ->SetBorderSize(0);
    legend_fit ->SetTextSize(0.03);
    legend_fit ->SetFillStyle(0);

    legend ->Draw("same");
    legend_fit ->Draw("same");

    c5->SaveAs("mean_sigma.pdf");




    TCanvas *c6=new TCanvas("c6","",800,600);

    c6->cd();
    gPad->SetLogy(1);
    TH2F *hh_C=new TH2F("hh_C","",20,1.2,4,20,10,40000);
    hh_C->GetXaxis()->SetTitle("p_{T} (GeV)");
    hh_C->GetYaxis()->SetTitle("Electron Fit Constant ");
    
    mh1Fit_Constant->Scale(1,"width");

    hh_C->Draw();
    mh1Fit_Constant->SetMarkerStyle(20);
    mh1Fit_Constant->SetMarkerColor(4);

    mh1Fit_Constant->Draw("SameP");
    
  c6->SaveAs("Constant.pdf");
  
  

  
}
/*
  return;

  TH1F * Pion_mean_Inv= (TH1F *)File_Re_InvM->Get("Pion_mean");
  TH1F * Pion_sigma_Inv=(TH1F *)File_Re_InvM->Get("Pion_sigma");
  TH1F * Proton_mean_Inv=(TH1F *)File_Re_InvM->Get("Proton_mean");;
  TH1F *  Proton_sigma_Inv=(TH1F *)File_Re_InvM->Get("Proton_sigma");;




  TCanvas *c5=new TCanvas("c5","",800,600);
  TH2F *hh_mean=new TH2F("hh_mean","",20,0,4,20,-10,25);
  hh_mean->GetXaxis()->SetTitle("p_{T} (GeV)");
  hh_mean->GetYaxis()->SetTitle("Nsigma E Mean");

  Pion_mean_Inv->SetMarkerStyle(20);
  Pion_sigma_Inv->SetMarkerStyle(24);
  Proton_mean_Inv->SetMarkerStyle(20);
  Proton_sigma_Inv->SetMarkerStyle(24);

  Pion_mean_Inv->SetMarkerColor(1);
  Pion_sigma_Inv->SetMarkerColor(1);
  Proton_mean_Inv->SetMarkerColor(2);
  Proton_sigma_Inv->SetMarkerColor(2);
  
  Pion_mean_Inv->SetLineColor(1);
  Pion_sigma_Inv->SetLineColor(1);
  Proton_mean_Inv->SetLineColor(2);
  Proton_sigma_Inv->SetLineColor(2);

  // Fit

  TFile *File_Fit=new TFile("Mean_Sigma_Hadron.root","READ");


  TH1F * Pion_mean_Fit= (TH1F *)File_Fit->Get("mh1Pion_Mean");
  TH1F * Pion_sigma_Fit=(TH1F *)File_Fit->Get("mh1Pion_Sigma");
  TH1F * Proton_mean_Fit=(TH1F *)File_Fit->Get("mh1K_P_Mean");;
  TH1F *  Proton_sigma_Fit=(TH1F *)File_Fit->Get("mh1K_P_Sigma");;

  TH1F * Kaon_mean_Fit=(TH1F *)File_Fit->Get("mh1K_K_Mean");;
  TH1F *  Kaon_sigma_Fit=(TH1F *)File_Fit->Get("mh1K_K_Sigma");;

  TH1F * Electron_mean_Fit=(TH1F *)File_Fit->Get("mh1Electron_Mean");
  TH1F *  Electron_sigma_Fit=(TH1F *)File_Fit->Get("mh1Electron_Sigma");

  //---------------------------------
  TH1F * Pion_Constant_Fit= (TH1F *)File_Fit->Get("mh1Pion_Constant");
  TH1F * Proton_Constant_Fit=(TH1F *)File_Fit->Get("mh1K_P_Constant");;
  TH1F * Kaon_Constarnt_Fit=(TH1F *)File_Fit->Get("mh1K_K_Constant");;
  TH1F * Electron_Constant_Fit=(TH1F *)File_Fit->Get("mh1Electron_Constant");

  for(Int_t ibin=0;ibin<9;ibin++)
    {
      if(ibin<9)
	{
	  // Proton_mean_Fit->SetBinContent(ibin+1,-999);
	  // Proton_sigma_Fit->SetBinContent(ibin+1,-999);
	  //Kaon_mean_Fit->SetBinContent(ibin+1,-999);
	  // Kaon_sigma_Fit->SetBinContent(ibin+1,-999);
	  cout<<" 22"<<endl;
	}
    }  


  Electron_mean_Fit->SetMarkerStyle(23);
  Electron_sigma_Fit->SetMarkerStyle(25);  

  Pion_mean_Fit->SetMarkerStyle(23);
  Pion_sigma_Fit->SetMarkerStyle(24);
  Proton_mean_Fit->SetMarkerStyle(23);
  Proton_sigma_Fit->SetMarkerStyle(24);

  Kaon_mean_Fit->SetMarkerStyle(23);
  Kaon_sigma_Fit->SetMarkerStyle(24);
  

  Pion_mean_Fit->SetMarkerColor(3);
  Pion_sigma_Fit->SetMarkerColor(3);
  Proton_mean_Fit->SetMarkerColor(4);
  Proton_sigma_Fit->SetMarkerColor(4);
  Kaon_mean_Fit->SetMarkerColor(15);
  Kaon_sigma_Fit->SetMarkerColor(15);
   


  Pion_mean_Fit->SetLineColor(3);
  Pion_sigma_Fit->SetLineColor(3);
  Proton_mean_Fit->SetLineColor(4);
  Proton_sigma_Fit->SetLineColor(4);
  Kaon_mean_Fit->SetLineColor(15);
  Kaon_sigma_Fit->SetLineColor(15);



  Pion_Constant_Fit->SetMarkerStyle(20);
  Proton_Constant_Fit->SetMarkerStyle(21);
  Kaon_Constarnt_Fit->SetMarkerStyle(22);
  Electron_Constant_Fit->SetMarkerStyle(23);

  Pion_Constant_Fit->SetMarkerColor(1);
  Proton_Constant_Fit->SetMarkerColor(2);
  Kaon_Constarnt_Fit->SetMarkerColor(4);
  Electron_Constant_Fit->SetMarkerColor(5);



  // Tof only
 // TFile *File_Tof=new TFile("file_Tof_only.root","READ");


 //  TH1F * Pion_mean_Tof= (TH1F *)File_Tof->Get("mh1_Pion_mean_Tof");
 //  TH1F * Pion_sigma_Tof=(TH1F *)File_Tof->Get("mh1_Pion_sigma_Tof");
 //  TH1F * Proton_mean_Tof=(TH1F *)File_Tof->Get("mh1_Proton_mean_Tof");;
 //  TH1F *  Proton_sigma_Tof=(TH1F *)File_Tof->Get("mh1_Proton_sigma_Tof");

  TFile *File_Tof=new TFile("TOF_ONlY_DATA.root","READ");
  
  
  TH1F * Pion_mean_Tof= (TH1F *)File_Tof->Get("mh1_Pion_mean");
  TH1F * Pion_sigma_Tof=(TH1F *)File_Tof->Get("mh1_Pion_sigma");
  TH1F * Proton_mean_Tof=(TH1F *)File_Tof->Get("mh1_Proton_mean");
  TH1F *  Proton_sigma_Tof=(TH1F *)File_Tof->Get("mh1_Proton_sigma");
  TH1F * Kaon_mean_Tof=(TH1F *)File_Tof->Get("mh1_Kaon_mean");
  TH1F *  Kaon_sigma_Tof=(TH1F *)File_Tof->Get("mh1_Kaon_sigma");
  
  
  Pion_mean_Tof->SetMarkerStyle(20);
  Pion_sigma_Tof->SetMarkerStyle(24);
  Proton_mean_Tof->SetMarkerStyle(20);
  Proton_sigma_Tof->SetMarkerStyle(24);
  Kaon_mean_Tof->SetMarkerStyle(20);
  Kaon_sigma_Tof->SetMarkerStyle(24);
  
  Pion_mean_Tof->SetMarkerColor(5);
  Pion_sigma_Tof->SetMarkerColor(5);
  Proton_mean_Tof->SetMarkerColor(6);
  Proton_sigma_Tof->SetMarkerColor(6);
  Kaon_mean_Tof->SetMarkerColor(7);
  Kaon_sigma_Tof->SetMarkerColor(7);
  
  Pion_mean_Tof->SetLineColor(5);
  Pion_sigma_Tof->SetLineColor(5);
  Proton_mean_Tof->SetLineColor(6);
  Proton_sigma_Tof->SetLineColor(6);

  Kaon_mean_Tof->SetLineColor(7);
  Kaon_sigma_Tof->SetLineColor(7);
  
  // ------------------------------
  // shift nsigma e
 TFile *File_Electron_mean_sigma=new TFile("Electron_Mean_igma.root ","READ");


  TH1F * Electron_mean= (TH1F *)File_Electron_mean_sigma->Get("mh1mean");
  TH1F * Electron_sigma=(TH1F *)File_Electron_mean_sigma->Get("mh1sigma");
  //TH1F * Proton_mean_Shift=(TH1F *)File_Shift->Get("mh1_Proton_mean_Shift");;
  // TH1F *  Proton_sigma_Shift=(TH1F *)File_Shift->Get("mh1_Proton_sigma_Shift");




  Electron_mean->SetMarkerStyle(25);
  Electron_sigma ->SetMarkerStyle(22);
  // Proton_mean_Shift->SetMarkerStyle(20);
  // Proton_sigma_Shift->SetMarkerStyle(24);

  Electron_mean->SetMarkerColor(11);
  Electron_sigma->SetMarkerColor(11);
  // Proton_mean_Shift->SetMarkerColor(8);
  // Proton_sigma_Shift->SetMarkerColor(8);
  
  Electron_mean->SetLineColor(7);
  Electron_sigma->SetLineColor(7);
  // Proton_mean_Shift->SetLineColor(8);
  // Proton_sigma_Shift->SetLineColor(8);
  //------------------------------------------


  
  c5->cd();
  hh_mean->Draw();

    Pion_mean_Inv->Draw("sameP");
   Proton_mean_Inv->Draw("sameP");
  Pion_mean_Fit->Draw("sameP");
  
  Proton_mean_Fit->Draw("sameP");
  Kaon_mean_Fit->Draw("sameP");
  Electron_mean_Fit->Draw("sameP");
  Pion_mean_Tof->Draw("sameP");
  Proton_mean_Tof->Draw("sameP");
  Kaon_mean_Tof->Draw("sameP");
    
    Electron_mean->Draw("sameP");
    // Electron_sigma->Draw("same");
  // Proton_mean_Shift->Draw("sameP");
   //Proton_sigma_Shift->Draw("sameP");
  

    TLegend *legend  = new TLegend(0.6,0.4,0.85,0.8);
    legend ->AddEntry(Pion_mean_Inv," Pion From KShort n  #sigma_{e} mean","lpe");

    legend ->AddEntry(Proton_mean_Inv," Proton From #Lambda #sigma_{e} mean","lpe");


 


    legend ->AddEntry(Electron_sigma," Photonic electron","lpe");
    legend ->AddEntry(Electron_mean_Fit,"Electron 4 Gauss Fit  #sigma_{e} mean","lpe");
 


  legend ->AddEntry(Pion_mean_Tof," Pion From Tof #sigma_{e} mean","lpe");
  legend ->AddEntry(Pion_mean_Fit," Pion From 4 Gauss Fit  #sigma_{e} mean","lpe");
  // legend ->AddEntry(Pion_sigma_Tof," Pion From Tof  sigma","lpe");
  legend ->AddEntry(Proton_mean_Tof," Proton Tof #sigma_{e} mean","lpe");
  legend ->AddEntry(Proton_mean_Fit," Proton 4 Gauss Fit  #sigma_{e} mean","lpe");

  legend ->AddEntry(Kaon_mean_Tof," Kaon Tof #sigma_{e} mean","lpe");
  legend ->AddEntry(Kaon_mean_Fit," Kaon 4 Gauss Fit  #sigma_{e} mean","lpe");
  
  //  legend ->AddEntry(Proton_sigma_Tof," Proton From Tof sigma ","lpe");

  // legend ->AddEntry(Electron_mean," Pion From Shift   #sigma_{e} mean","lpe");
  
  // legend ->AddEntry(Proton_mean_Shift," Proton Shift  #sigma_{e} mean","lpe");
  // legend ->AddEntry(Proton_sigma_Shift," Proton Shift   sigma ","lpe");
  
  legend ->SetBorderSize(0);
  legend ->SetTextSize(0.03);
  legend ->SetFillStyle(0);
  legend ->Draw("same");
  

  c5->SaveAs("mean.pdf");
  
  TCanvas *c6=new TCanvas("c6","",800,600);
 
  TH2F *hh_sigma=new TH2F("hh_sigma","",20,0,4,20,0.001,3);
  hh_sigma->GetXaxis()->SetTitle("p_{T} (GeV)");
  hh_sigma->GetYaxis()->SetTitle("Nsigma E Sigma");
  
  hh_sigma->Draw();


  // Pion_sigma_Inv->Draw("sameP");
  // Proton_sigma_Inv->Draw("sameP");
  Pion_sigma_Fit->Draw("sameP");
  Kaon_sigma_Fit->Draw("sameP");
  Proton_sigma_Fit->Draw("sameP");
  
  Pion_sigma_Tof->Draw("sameP");
  Proton_sigma_Tof->Draw("sameP");
  Kaon_sigma_Tof->Draw("sameP");
  Electron_sigma_Fit->Draw("sameP");  
  

  TLegend *legend_sigma  = new TLegend(0.4,0.6,0.65,0.8);
  
  // legend_sigma ->AddEntry(Pion_sigma_Inv," Pion From KShort n  #sigma_{e} sigma","lpe");
  // legend_sigma ->AddEntry(Proton_sigma_Inv," Proton From #Lambda sigma ","lpe");
  legend_sigma ->AddEntry(Pion_sigma_Fit," Pion From 4 Gauss  Fit sigma","lpe");
  legend_sigma ->AddEntry(Proton_sigma_Fit," Proton From 4 Gauss Fit  sigma ","lpe");
  legend_sigma ->AddEntry(Kaon_sigma_Fit," Kaon From 4 Gauss Fit  sigma ","lpe");
  legend_sigma ->AddEntry(Electron_sigma_Fit," Electron From 4 Gauss Fit  sigma ","lpe");
  
  legend_sigma ->AddEntry(Pion_sigma_Tof," Pion From Tof  sigma","lpe");
  legend_sigma ->AddEntry(Proton_sigma_Tof," Proton From Tof  sigma","lpe");
  legend_sigma ->AddEntry(Kaon_sigma_Tof," Kaon From Tof  sigma","lpe");
  
  
  
  //  legend ->AddEntry(Proton_sigma_Tof," Proton From Tof sigma ","lpe");

  // legend ->AddEntry(Electron_mean," Pion From Shift   #sigma_{e} mean","lpe");
  // legend ->AddEntry(Electron_sigma," Pion From Shift  sigma","lpe");
  // legend ->AddEntry(Proton_mean_Shift," Proton Shift  #sigma_{e} mean","lpe");
  // legend ->AddEntry(Proton_sigma_Shift," Proton Shift   sigma ","lpe");
  
  legend_sigma ->SetBorderSize(0);
  legend_sigma ->SetTextSize(0.03);
  legend_sigma ->SetFillStyle(0);
  legend_sigma ->Draw("same");
  

  c6->SaveAs("sigma.pdf");
  

  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd();
  gPad->SetLogy();

  TH2F *h2=new TH2F("h2","",10,0,4,1000,1,1e8);
  h2->Draw();
  h2->SetTitle(" Multi Gauss Fit Normalization Constant");
  h2->GetXaxis()->SetTitle("p_{T} GeV/c");
  h2->GetYaxis()->SetTitle("Fit Constant");

  Pion_Constant_Fit->Draw("SAMEP");
  Proton_Constant_Fit->Draw("SAMEPE1");
  Kaon_Constarnt_Fit->Draw("SAMEPE");
  Electron_Constant_Fit->Draw("SAMEPE1");


  TLegend *legend_Constant  = new TLegend(0.25,0.1,0.5,0.3);
  
  // legend_Constant ->AddEntry(Pion_Constant_Inv," Pion From KShort n  #Constant_{e} Constant","lpe");
  // legend_Constant ->AddEntry(Proton_Constant_Inv," Proton From #Lambda Constant ","lpe");
  legend_Constant ->AddEntry(Pion_Constant_Fit," Pion From 4 Gauss  Fit Constant","lpe");
  legend_Constant ->AddEntry(Proton_Constant_Fit," Proton From 4 Gauss Fit  Constant ","lpe");
  legend_Constant ->AddEntry(Kaon_Constarnt_Fit," Kaon From 4 Gauss Fit  Constant ","lpe");
  legend_Constant ->AddEntry(Electron_Constant_Fit," Electron From 4 Gauss Fit  Constant ","lpe");
  
  
  
  
  //  legend ->AddEntry(Proton_sigma_Tof," Proton From Tof sigma ","lpe");

  // legend ->AddEntry(Electron_mean," Pion From Shift   #sigma_{e} mean","lpe");
  // legend ->AddEntry(Electron_sigma," Pion From Shift  sigma","lpe");
  // legend ->AddEntry(Proton_mean_Shift," Proton Shift  #sigma_{e} mean","lpe");
  // legend ->AddEntry(Proton_sigma_Shift," Proton Shift   sigma ","lpe");
  
  legend_Constant ->SetBorderSize(0);
  legend_Constant ->SetTextSize(0.03);
  legend_Constant ->SetFillStyle(0);
  legend_Constant ->SetTextFont(62);
  
  legend_Constant ->Draw("same");

  c2->SaveAs("Constant.pdf");



}
*/
