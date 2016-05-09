#include<fstream>
#include <iostream>
#include<iomanip>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "../mBinning_MB.h"

using namespace std;
void Draw_purity(TH1D * );
void Draw_purity_sys(TH1D *,TH1D *,TH1D * );

void Draw_purity()
{

  
  
  
  TH1D *purity_MB=new TH1D("purity_MB","",NpT_bins_run12_MB,pt_run12_MB);
  TH1D *purity_oneSigma=new TH1D("purity_oneSigma","",NpT_bins_run12_MB,pt_run12_MB);
  TH1D *purity_twoSigma=new TH1D("purity_twoSigma","",NpT_bins_run12_MB,pt_run12_MB);
  TH1D *purity_threeSigma=new TH1D("purity_threeSigma","",NpT_bins_run12_MB,pt_run12_MB);

  TFile *file_lowpt=new TFile("Combine_Effiicency/Inclusive_purity_MB_EMC.root","READ");
  TFile *file_highpt=new TFile("purity_MB_highpT.root","READ");

  TH1D *purity_MB_l=(TH1D *) file_lowpt->Get("mh1purity_44")->Clone("purity_MB_l");
  // TH1D *purity_oneSigma_l=(TH1D *) file_lowpt->Get("purity_oneSigma")->Clone("purity_oneSigma_l");
  // TH1D *purity_twoSigma_l=(TH1D *) file_lowpt->Get("purity_twoSigma")->Clone("purity_twoSigma_l");
  // TH1D *purity_threeSigma_l=(TH1D *) file_lowpt->Get("purity_threeSigma")->Clone("purity_threeSigma_l");

  TH1D *purity_MB_h=(TH1D *) file_highpt->Get("purity_MB")->Clone("purity_MB_h");
  // TH1D *purity_oneSigma_h=(TH1D *) file_highpt->Get("purity_oneSigma")->Clone("purity_oneSigma_h");
  // TH1D *purity_twoSigma_h=(TH1D *) file_highpt->Get("purity_twoSigma")->Clone("purity_twoSigma_h");
  // TH1D *purity_threeSigma_h=(TH1D *) file_highpt->Get("purity_threeSigma")->Clone("purity_threeSigma_h");


  for(int ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      if(ipt<NpT_bins_run12_MB-8)
      	{
      	  purity_MB->SetBinContent(ipt+1,purity_MB_l->GetBinContent(ipt+1));
      	  purity_MB->SetBinError(ipt+1,purity_MB_l->GetBinError(ipt+1));
	  
      	  // purity_oneSigma->SetBinContent(ipt+1,purity_oneSigma_l->GetBinContent(ipt+1));
      	  // purity_oneSigma->SetBinError(ipt+1,purity_oneSigma_l->GetBinError(ipt+1));

      	  // purity_twoSigma->SetBinContent(ipt+1,purity_twoSigma_l->GetBinContent(ipt+1));
      	  // purity_twoSigma->SetBinError(ipt+1,purity_twoSigma_l->GetBinError(ipt+1));

      	  // purity_threeSigma->SetBinContent(ipt+1,purity_threeSigma_l->GetBinContent(ipt+1));
      	  // purity_threeSigma->SetBinError(ipt+1,purity_threeSigma_l->GetBinError(ipt+1));


	  
      	}
      else
      	{

	  
  	  purity_MB->SetBinContent(ipt+1,purity_MB_h->GetBinContent(ipt+1));
	  purity_MB->SetBinError(ipt+1,purity_MB_h->GetBinError(ipt+1));

	  // purity_oneSigma->SetBinContent(ipt+1,purity_oneSigma_h->GetBinContent(ipt+1));
	  // purity_oneSigma->SetBinError(ipt+1,purity_oneSigma_h->GetBinError(ipt+1));

	  // purity_twoSigma->SetBinContent(ipt+1,purity_twoSigma_h->GetBinContent(ipt+1));
	  // purity_twoSigma->SetBinError(ipt+1,purity_twoSigma_h->GetBinError(ipt+1));

	  // purity_threeSigma->SetBinContent(ipt+1,purity_threeSigma_h->GetBinContent(ipt+1));
	  // purity_threeSigma->SetBinError(ipt+1,purity_threeSigma_h->GetBinError(ipt+1));

	  
	     	}
    }
  
  //  Draw_purity_sys(purity_oneSigma,purity_twoSigma,purity_threeSigma);
   Draw_purity(purity_MB);

   TFile *file=new TFile("purity_MB_com.root","RECREATE");
   purity_MB->Write();
   file->Close();
   
   
}
void Draw_purity(TH1D * purity)
{
  gStyle->SetOptStat(00000);
  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd();
  TH2F *h2=new TH2F("h2","",10,0,4,10,0,1.1);
  h2->Draw();
  h2->GetXaxis()->SetTitle("p_{T}");
  h2->GetYaxis()->SetTitle("purity");
  purity->SetMarkerStyle(20);
  purity->SetMarkerColor(2);
  purity->Draw("samePE1");
  c2->SaveAs("purity_MB.pdf");
}

void Draw_purity_sys(TH1D * purity_oneSigma,TH1D *purity_twoSigma,TH1D *purity_threeSigma)
{
  gStyle->SetOptStat(00000);
  TCanvas *c3=new TCanvas("c3","",800,600);
  c3->cd();
  TH2F *h3=new TH2F("h3","",10,1,4,10,0.6,1.1);
  h3->Draw();
  h3->GetXaxis()->SetTitle("p_{T}");
  h3->GetYaxis()->SetTitle("purity");
  purity_oneSigma->SetMarkerStyle(20);
  purity_twoSigma->SetMarkerStyle(20);
  purity_threeSigma->SetMarkerStyle(20);

  purity_oneSigma->SetMarkerColor(1);
  purity_twoSigma->SetMarkerColor(2);
  purity_threeSigma->SetMarkerColor(3);


  purity_oneSigma->Draw("samePE1");
  purity_twoSigma->Draw("samePE1");
  purity_threeSigma->Draw("samePE1");
  TLegend *legend = new TLegend(0.25,0.35,0.45,0.6);
  legend->AddEntry(purity_oneSigma,"one sigma standard deviation","lp");
  legend->AddEntry(purity_twoSigma,"two sigma standard deviation","lp");
  legend->AddEntry(purity_threeSigma,"three sigma standard deviation","lp");
  
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03);
  legend ->SetTextFont(62);       
  legend->Draw("same");

  c3->SaveAs("purity_MB_sys.pdf");
}
