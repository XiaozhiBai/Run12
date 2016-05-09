  
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
#include"../mBinning_MB.h"
//const int NPtBin= 28;


using namespace std;

TH1F * mh1Sinle_Track_MB[3];
TH1F * mh1Pair_MB[3];




TH1F * PHE_RE_MB[3];


TString histname[3]={"Gamma_conversion","Pi0Dalitz","EtaDalitzDecay"};
//void Draw_effciency(TH1F *a[],TH1F *b[]);
void Draw_effciency_Total( TH1F * const  PHE_RE_MB[]);
TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
void setpad(TPad *,float , float , float , float );


void Photonic_RE_efficiency(TString Gamma_conversion_MB="Root_File_PHE/Gamma_MB.root",TString Pi0Dalitz_MB="Root_File_PHE/Pi0_dalitz_MB.root",TString EtaDalitz_MB="Root_File_PHE/Eta_dalitz_MB.root")
{

  char buf[1024];
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
  
  // Get single and pair for gamma, Pi0, Eta 
  mh1Sinle_Track_MB[0]=(TH1F *) inFile_Gamma_MB->Get("mSingle_Track")->Clone("single_gamma_MB");
  mh1Pair_MB[0]=(TH1F *)inFile_Gamma_MB->Get("mPair")->Clone("Pair_gamma_MB");

  mh1Sinle_Track_MB[1]=(TH1F *) inFile_Pi0Dalitz_MB->Get("mSingle_Track")->Clone("single_pi0_MB");
  mh1Pair_MB[1]=(TH1F *)inFile_Pi0Dalitz_MB->Get("mPair")->Clone("Pair_pi0_MB");
  
  mh1Sinle_Track_MB[2]=(TH1F *) inFile_EtaDalitz_MB->Get("mSingle_Track")->Clone("single_Eta_MB");
  mh1Pair_MB[2]=(TH1F *)inFile_EtaDalitz_MB->Get("mPair")->Clone("Pair_Eta_MB");


  //  mh1Pair_MB[0]->Draw();
  // return ;

  TString histname_single_MB[3]={"Gamma_conversion_single_MB","Pi0Dalitz_single_MB","EtaDalitzDecay_single_MB"};
  TString histname_pair_MB[3]={"Gamma_conversion_pair_MB","Pi0Dalitz_pair_MB","EtaDalitzDecay_pair_MB"};
  

  TH1F *  mh1Sinle_Track_rebin_MB[3];
  TH1F * mh1Pair_rebin_MB[3];
     
  


  TGraphAsymmErrors  * Photonic_re_eff_MB[3];   


  
  for(Int_t i=0;i<3;i++)
    {
      mh1Sinle_Track_rebin_MB[i]=(TH1F *) mh1Sinle_Track_MB[i]->Rebin(NpT_bins_run12_MB,histname_single_MB[i],pt_run12_MB);
      mh1Pair_rebin_MB[i]=(TH1F * ) mh1Pair_MB[i]->Rebin(NpT_bins_run12_MB,histname_pair_MB[i],pt_run12_MB);
      
      
      Photonic_re_eff_MB[i]=new TGraphAsymmErrors(mh1Pair_rebin_MB[i],mh1Sinle_Track_rebin_MB[i],"N");
 
    
      
      PHE_RE_MB[i]=new TH1F(histname[i]+TString("MB"),"Photonic electron reconstruction efficiency",NpT_bins_run12_MB,pt_run12_MB);
      PHE_RE_MB[i]=new TH1F(histname[i]+TString("HT"),"Photonic electron reconstruction efficiency",NpT_bins_run12_MB,pt_run12_MB);
      for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
	{
	  Double_t x=0,y=0,x_err=0,y_err=0;
	  Photonic_re_eff_MB[i]  ->GetPoint(ipt,x,y);
	  y_err=Photonic_re_eff_MB[i]->GetErrorY(ipt);
	  
	  PHE_RE_MB[i]->SetBinContent(ipt+1,y);
	  PHE_RE_MB[i]->SetBinError(ipt+1,y_err);
	}
    }
 



  Draw_effciency_Total(PHE_RE_MB);
}
void Draw_effciency_Total(TH1F * const PHE_RE_MB[] )
{
  gStyle->SetOptStat(0);
  cout<<"  Draw _efficiency"<<endl;
  TFile * file_PHE=TFile::Open("PHE_Contribution_MB.root");

  TH1F*	mh1_Gamma=(TH1F *)file_PHE->Get("mh1_Gamma");	
  TH1F*	mh1_pi0=(TH1F *)file_PHE->Get("mh1_pi0");	
  TH1F*	mh1_Eta=(TH1F *)file_PHE->Get("mh1_Eta");
  

  TH1F *mh1_MB=new TH1F("mh1_MB","",NpT_bins_run12_MB,pt_run12_MB);

  
  
  TH1F *Gamma_MB=(TH1F *) PHE_RE_MB[0]->Clone("Gamma_MB");
  TH1F *pi0_MB=(TH1F *) PHE_RE_MB[1]->Clone("pi0_MB");
  TH1F *Eta_MB=(TH1F *) PHE_RE_MB[2]->Clone("Eta_MB");

  Gamma_MB->Multiply(mh1_Gamma);
  pi0_MB->Multiply(mh1_pi0);
  Eta_MB->Multiply(mh1_Eta);

  mh1_MB->Add(pi0_MB);
  mh1_MB->Add(Gamma_MB);
  mh1_MB->Add(Eta_MB);

  cout<<"Total "<<endl;
  
  mh1_MB->Print("all");

  TCanvas *c2 = new TCanvas("c2"," Rhotonic_Re_eff",0,0,800,600);
  c2->cd(); 
  
  TH2F * hh=new TH2F("hh","",100,0,4,100,0,1);
  
  hh->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hh->GetYaxis()->SetTitle(" #epsilon Photonic");
  hh->GetYaxis()->CenterTitle(true);
  hh->GetXaxis()->CenterTitle(true);

  hh->Draw();

  //  mh1_MB->SetMarkerStyle(20);
  mh1_MB->SetMarkerStyle(20);
  TH1F *PHE_re_efficiency=(TH1F *) mh1_MB->Clone("PHE_re_efficiency");
  //  mh1_MB->SetMarkerStyle(20);

  TF1 *f_1 = new TF1("f_1","[0]*(1.0 -TMath::Exp(-(x+[1])/[2]))",0,20);
  f_1->SetParameters(0.8,9.37659e-02,4.40904e-01);
  f_1->SetLineColor(2);

  mh1_MB->Fit("f_1","R+","same",0.2,14);
  
  // TFile *file=new TFile("Photonic_re_Efficiency.root","READ");
  // TH1F *mh1_MB_test=(TH1F *) file->Get("mh1_MB");
  // mh1_MB_test->SetMarkerColor(2);
  
  //  mh1_MB_test->Draw("sameP");
   
  c2->SaveAs("PHE_Efficiency.pdf");
  
  TFile *file=new TFile("Photonic_re_Efficiency.root","RECREATE");
  PHE_re_efficiency->Write();
  file->Close();
  
  // TF1 *f_2 = new TF1("f_2","[0]/(TMath::Exp(-(x-[1])/[2]) + 1)",0,20);
  // f_2->FixParameter(0,0.8);
  // f_2->SetParameter(1,5.0);
  // f_2->SetParameter(2,1.5);

  // f_2->SetLineColor(4);

  //  mh1_MB->Fit("f_2","R+","same",0.2,14);


  //  mh1_MB->Draw("sameCL");
}

void Draw_effciency( TH1F *  PHE_RE_MB[])
{
  
  
  // TH1F::SetDefaultSumw2();
  char buf[1024];
  gStyle->SetOptStat(00000);    
  TCanvas *c2 = new TCanvas("c2"," Rhotonic_Re_eff",0,0,800,600);
  c2->cd(); 

  TH2F * hh=new TH2F("hh","Photonic Electron Rec. Efficiency",100,0,14,100,0,1);
  
  hh->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hh->GetYaxis()->SetTitle(" #epsilon Photonic");
  hh->GetYaxis()->CenterTitle(true);
  hh->GetXaxis()->CenterTitle(true);
  //Gamma
  PHE_RE_MB[0]->SetMarkerStyle(20);
  PHE_RE_MB[0]->SetMarkerColor(1);
  PHE_RE_MB[0]->SetLineColor(1);
  PHE_RE_MB[0]->SetMarkerStyle(24);
  PHE_RE_MB[0]->SetMarkerColor(1);
  PHE_RE_MB[0]->SetLineColor(1);
  //Pi0
  PHE_RE_MB[1]->SetMarkerStyle(20);
  PHE_RE_MB[1]->SetMarkerColor(2);
  PHE_RE_MB[1]->SetLineColor(2);
  PHE_RE_MB[1]->SetMarkerStyle(24);
  PHE_RE_MB[1]->SetMarkerColor(2);
  PHE_RE_MB[1]->SetLineColor(2);
  //Eta
  
  PHE_RE_MB[2]->SetMarkerStyle(20);
  PHE_RE_MB[2]->SetMarkerColor(4);
  PHE_RE_MB[2]->SetLineColor(4);
  PHE_RE_MB[2]->SetMarkerStyle(24);
  PHE_RE_MB[2]->SetMarkerColor(4);
  PHE_RE_MB[2]->SetLineColor(4);
  
  hh->Draw();
  PHE_RE_MB[0]->Draw("samePE1");
  PHE_RE_MB[1]->Draw("samePE1");
  PHE_RE_MB[2]->Draw("samePE1");

 
  PHE_RE_MB[0]->Draw("samePE1");
  PHE_RE_MB[1]->Draw("samePE1");
  PHE_RE_MB[2]->Draw("samePE1");
  
  


  // TLegend *legend_MB = new TLegend(0.65,0.35,0.85,0.5);   
  // legend_MB->AddEntry(PHE_RE_MB[0], "HT  #gamma conversion","lp"); 
  // legend_MB->AddEntry(PHE_RE_MB[1], "HT  #pi0 #rightarrow->#gamma e^{+} e^{-}","lp"); 
  // legend_MB->AddEntry(PHE_RE_MB[2], "HT  #eta #rightarrow->#gamma e^{+} e^{-}","lp"); 

  TLegend *legend_MB = new TLegend(0.5,0.35,0.65,0.5);

  legend_MB->AddEntry(PHE_RE_MB[0], " MB  #gamma conversion","lp"); 
  legend_MB->AddEntry(PHE_RE_MB[1], " MB  #pi #rightarrow->#gamma e^{+} e^{-}","lp"); 
  legend_MB->AddEntry(PHE_RE_MB[2], " MB  #eta #rightarrow->#gamma e^{+} e^{-}","lp"); 
  
  legend_MB->SetTextSize(0.03); 
  legend_MB->SetBorderSize(0);
  legend_MB->SetFillStyle(0);   
  legend_MB->Draw();

  legend_MB->SetTextSize(0.03); 
  legend_MB->SetBorderSize(0);
  legend_MB->SetFillStyle(0);   
  legend_MB->Draw();


  c2->SaveAs("PHE_Plots/Photonic_re_efficiency_La0080.pdf");    
  
  // fit 

  TCanvas *c3=new TCanvas("c3","",800,600);
  c3->cd();

  hh->Draw();
  TF1 *f0 = new TF1("f0","[0]*(1.0 -TMath::Exp(-(x+[1])/[2]))",0,20);
  f0->SetParameters(0.8,9.37659e-02,4.40904e-01);

  f0->SetLineColor(1);
  TF1 *f1 = new TF1("f1","[0]*(1.0 -TMath::Exp(-(x+[1])/[2]))",0,20);
  f1->SetParameters(0.8,9.37659e-02,4.40904e-01);
  f1->SetLineColor(2);

  TF1 *f2 = new TF1("f2","[0]*(1.0 -TMath::Exp(-(x+[1])/[2]))",0,20);
  f2->SetParameters(0.8,9.37659e-02,4.40904e-01);

  f2->SetLineColor(4);

  PHE_RE_MB[0]->Fit("f0","R","same",0.2,14);
  PHE_RE_MB[1]->Fit("f1","R","same",0.2,14);
  PHE_RE_MB[2]->Fit("f2","R","same",0.2,14);
  legend_MB->Draw();


  c3->SaveAs("PHE_Plots/Photonic_MB.pdf");


  TCanvas *c4=new TCanvas("c4","",800,600);
  c4->cd();

  hh->Draw();
  PHE_RE_MB[0]->Fit("f0","R","same",0.2,14);
  PHE_RE_MB[1]->Fit("f1","R","same",0.2,14);
  PHE_RE_MB[2]->Fit("f2","R","same",0.2,14);
  legend_MB->Draw();
  c4->SaveAs("PHE_Plots/Photonic_MB.pdf");



  
  // TFile *file_2=new TFile("Photonic_ReconStructionEfficiency_MB.root","RECREATE");
  
  // PHE_RE_MB[0]->Write();
  // PHE_RE_MB[1]->Write();
  // PHE_RE_MB[2]->Write();
  
  // file_2->Close();

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




