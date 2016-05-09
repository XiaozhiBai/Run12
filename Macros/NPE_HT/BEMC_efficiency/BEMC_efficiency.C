
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
#include "../mBinning_HT.h"

using namespace std;

const int mh2Nhist= 2;

TH2F *mh2hist_unlike_pass[mh2Nhist];
TH2F *mh2hist_like_pass[mh2Nhist];
TH1F *mh1hist_pT_unlike_pass[mh2Nhist];
TH1F *mh1hist_pT_like_pass[mh2Nhist];
TH1F *mh1hist_pT_unlike_like_pass[mh2Nhist];


TH1F *mh1hist_pT_unlike_pass_sys[mh2Nhist];
TH1F *mh1hist_pT_like_pass_sys[mh2Nhist];
TH1F *mh1hist_pT_unlike_like_pass_sys[mh2Nhist];


TH2F *mh2hist_unlike_total[mh2Nhist];
TH2F *mh2hist_like_total[mh2Nhist];

TH1F *mh1hist_pT_unlike_total[mh2Nhist];
TH1F *mh1hist_pT_like_total[mh2Nhist];
TH1F *mh1hist_pT_unlike_like_total[mh2Nhist];

TH1F *mh1hist_pT_unlike_total_sys[mh2Nhist];
TH1F *mh1hist_pT_like_total_sys[mh2Nhist];
TH1F *mh1hist_pT_unlike_like_total_sys[mh2Nhist];



TH3F *mh3hist_unlike_partnernSigmaE_pt[mh2Nhist];
TH3F *mh3hist_like_partnernSigmaE_pt[mh2Nhist];

TH1F *mh1hist_unlike[mh2Nhist][NpT_bins_run12_HT];
TH1F *mh1hist_like[mh2Nhist][NpT_bins_run12_HT];
TH1F *mh1hist_unlike_like[mh2Nhist][NpT_bins_run12_HT];



TString mh2HistName_Photonic_unlike_total[mh2Nhist]={"mh2Prim_Ele_MassVspT_noBEMCcut_unlikeTrg0","mh2Prim_Ele_MassVspT_noBEMCcut_unlikeTrg1"};
TString mh2HistName_Photonic_like_total[mh2Nhist]={"mh2Prim_Ele_MassVspT_noBEMCcut_likeTrg0","mh2Prim_Ele_MassVspT_BEMCcut_likeTrg1"};

TString mh2HistName_Photonic_unlike_pass[mh2Nhist]={"mh2Prim_Ele_MassVspT_BEMCcut_unlikeTrg0","mh2Prim_Ele_MassVspT_BEMCcut_unlikeTrg1"};
TString mh2HistName_Photonic_like_pass[mh2Nhist]={"mh2Prim_Ele_MassVspT_BEMCcut_likeTrg0","mh2Prim_Ele_MassVspT_BEMCcut_likeTrg1"};


TString mh3HistName_Photonic_unlike_total[mh2Nhist]={"mh3nSigmaEPart_pT_Mass_unlike_totalTrg0","mh3nSigmaEPart_pT_Mass_unlike_totalTrg1"};
TString mh3HistName_Photonic_like_total[mh2Nhist]={"mh3nSigmaEPart_pT_Mass_like_totalTrg0","mh3nSigmaEPart_pT_Mass_like_totalTrg1"};



TString mh1_TitleY[mh2Nhist]={"Counts","Counts"};
TString mh1_TitleX[mh2Nhist]={"Mass ee (GeV/c^2)","partner Nsigmae"};



char buf[1024];
void BEMC_efficiency()
{
  cout<< sizeof(pt_run12_HT)/sizeof(pt_run12_HT[0])<<endl;

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
  
   //    TFile *File=new TFile("../RootFile/Root_File_2_26/hist_2_26.root","READ");
   //dev/   TFile *File=new TFile("../RootFile/Root_File_3_30/hist_3_30.root","READ");
     TFile *File=new TFile("../RootFile/Root_File_5_1/hist_5_1.root","READ");
  
  for(Int_t ihist=0;ihist<mh2Nhist;ihist++)
    {
      mh2hist_unlike_pass[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_unlike_pass[ihist]);
      mh2hist_like_pass[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_like_pass[ihist]);
      mh2hist_unlike_total[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_unlike_total[ihist]);
      mh2hist_like_total[ihist]=(TH2F *) File->Get(mh2HistName_Photonic_like_total[ihist]);
      
      sprintf(buf,"unlike_passTrg%i",ihist);
      mh1hist_pT_unlike_pass[ihist]=(TH1F *)mh2hist_unlike_pass[ihist]->ProjectionX(buf,1,10);
      sprintf(buf,"like_passTrg%i",ihist);
      mh1hist_pT_like_pass[ihist]=(TH1F *)mh2hist_like_pass[ihist]->ProjectionX(buf,1,10);
      sprintf(buf,"unlike_like_passTrg%i",ihist);
      mh1hist_pT_unlike_like_pass[ihist]=(TH1F *)mh1hist_pT_unlike_pass[ihist]->Clone(buf);
      mh1hist_pT_unlike_like_pass[ihist]->Add(mh1hist_pT_like_pass[ihist],-1);

      // sys
      sprintf(buf,"unlike_passTrg_sys%i",ihist);
      mh1hist_pT_unlike_pass_sys[ihist]=(TH1F *)mh2hist_unlike_pass[ihist]->ProjectionX(buf,1,20);
      sprintf(buf,"like_passTrg_sys%i",ihist);
      mh1hist_pT_like_pass_sys[ihist]=(TH1F *)mh2hist_like_pass[ihist]->ProjectionX(buf,1,20);
      sprintf(buf,"unlike_like_passTrg_sys%i",ihist);
      mh1hist_pT_unlike_like_pass_sys[ihist]=(TH1F *)mh1hist_pT_unlike_pass_sys[ihist]->Clone(buf);
      mh1hist_pT_unlike_like_pass_sys[ihist]->Add(mh1hist_pT_like_pass_sys[ihist],-1);




      
      sprintf(buf,"unlike_totalTrg%i",ihist);
      mh1hist_pT_unlike_total[ihist]=(TH1F *)mh2hist_unlike_total[ihist]->ProjectionX(buf,1,10);
      sprintf(buf,"like_totalTrg%i",ihist);
      mh1hist_pT_like_total[ihist]=(TH1F *)mh2hist_like_total[ihist]->ProjectionX(buf,1,10);
      sprintf(buf,"unlike_like_totalTrg%i",ihist);
      mh1hist_pT_unlike_like_total[ihist]=(TH1F *)mh1hist_pT_unlike_total[ihist]->Clone(buf);
      mh1hist_pT_unlike_like_total[ihist]->Add(mh1hist_pT_like_total[ihist],-1);

      //sys
      sprintf(buf,"unlike_totalTrg_sys%i",ihist);
      mh1hist_pT_unlike_total_sys[ihist]=(TH1F *)mh2hist_unlike_total[ihist]->ProjectionX(buf,1,20);
      sprintf(buf,"like_totalTrg_sys%i",ihist);
      mh1hist_pT_like_total_sys[ihist]=(TH1F *)mh2hist_like_total[ihist]->ProjectionX(buf,1,20);
      sprintf(buf,"unlike_like_totalTrg_sys%i",ihist);
      mh1hist_pT_unlike_like_total_sys[ihist]=(TH1F *)mh1hist_pT_unlike_total_sys[ihist]->Clone(buf);
      mh1hist_pT_unlike_like_total_sys[ihist]->Add(mh1hist_pT_like_total_sys[ihist],-1);

      
      mh3hist_unlike_partnernSigmaE_pt[ihist]=(TH3F *) File->Get(mh3HistName_Photonic_unlike_total[ihist]);
      mh3hist_like_partnernSigmaE_pt[ihist]=(TH3F *) File->Get(mh3HistName_Photonic_like_total[ihist]);

      // mh3hist_unlike_partnernSigmaE_pt[ihist]->RebinX(4);
      // mh3hist_like_partnernSigmaE_pt[ihist]->RebinX(4);

    }
  
  //  mh3hist_unlike_partnernSigmaE_pt[0]->Draw();
  // mh3hist_unlike_mass_pt[0]->Draw();
  // mh3hist_like_mass_pt[0]->Draw();
  //  mh3hist_unlike_partnernSigmaE_pt[0]->ProjectionX()->Draw();
  
  // return; 

  mh1hist_pT_unlike_like_pass[0]->Add(mh1hist_pT_unlike_like_pass[1],1);
  mh1hist_pT_unlike_like_total[0]->Add(mh1hist_pT_unlike_like_total[1],1);
  
  mh1hist_pT_unlike_like_pass_sys[0]->Add(mh1hist_pT_unlike_like_pass_sys[1],1);
  mh1hist_pT_unlike_like_total_sys[0]->Add(mh1hist_pT_unlike_like_total_sys[1],1);



  TFile *file=new TFile("Tracking_eff_HT.root","READ");

  TH1F *mNBEMC_nocuts=(TH1F *) file->Get("mNBEMC_nocuts");
  TH1F *mNBEMC_cut=(TH1F *) file->Get("mNBEMC_cut");

  TH1F *mNBEMC_nocuts_rebin=(TH1F *)  mNBEMC_nocuts->Rebin(NpT_bins_run12_HT,"mNBEMC_nocuts_rebin",pt_run12_HT);
  TH1F *mNBEMC_cut_rebin=(TH1F *) mNBEMC_cut->Rebin(NpT_bins_run12_HT,"mNBEMC_cut_rebin",pt_run12_HT);
  

  TH1F *Pt_pass_rebin=(TH1F *)  mh1hist_pT_unlike_like_pass[0]->Rebin(NpT_bins_run12_HT,"Pt_pass_rebin",pt_run12_HT);
  TH1F *Pt_total_rebin=(TH1F *) mh1hist_pT_unlike_like_total[0]->Rebin(NpT_bins_run12_HT,"Pt_total_rebin",pt_run12_HT);

  TH1F *Pt_pass_rebin_sys=(TH1F *)  mh1hist_pT_unlike_like_pass_sys[0]->Rebin(NpT_bins_run12_HT,"Pt_pass_rebin_sys",pt_run12_HT);
  TH1F *Pt_total_rebin_sys=(TH1F *) mh1hist_pT_unlike_like_total_sys[0]->Rebin(NpT_bins_run12_HT,"Pt_total_rebin_sys",pt_run12_HT);


  
  // for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
  //   {
  //     cout<< " psss " <<Pt_pass_rebin->GetBinContent(ipt+1)<< " total "<<Pt_total_rebin->GetBinContent(ipt+1)<< " ratio"<< Pt_pass_rebin->GetBinContent(ipt+1)/Pt_total_rebin->GetBinContent(ipt+1)<<endl;;
  //   }
  //TGraphAsymmErrors *efficiency_mc=new TGraphAsymmErrors(mNBEMC_cut_rebin,mNBEMC_nocuts_rebin,"N");
  TGraphAsymmErrors *efficiency=new TGraphAsymmErrors(Pt_pass_rebin,Pt_total_rebin,"N");
  TGraphAsymmErrors *efficiency_sys=new TGraphAsymmErrors(Pt_pass_rebin_sys,Pt_total_rebin_sys,"N");

  // efficiency->GetYaxis()->SetRangeUser(0,1);
  // efficiency->Draw();
  // efficiency_mc->Draw("same");
  // return ;
  
  
  TH1F *efficiency_bemc=new TH1F("efficiency_bemc","",NpT_bins_run12_HT,pt_run12_HT);

  
  for(Int_t i=0;i<NpT_bins_run12_HT;i++)
    {
      Double_t x=0,y=0,x_err=0,y_err=0;
      Double_t x_sys=0,y_sys=0,x_err_sys=0,y_err_sys=0;

      efficiency->GetPoint(i,x,y);
      y_err=efficiency->GetErrorY(i);

      efficiency_sys->GetPoint(i,x_sys,y_sys);
      y_err_sys=efficiency_sys->GetErrorY(i);

      Double_t bemc_err=sqrt(y_err*y_err+(y-y_sys)*(y-y_sys));
      efficiency_bemc->SetBinContent(efficiency_bemc->FindBin(x),y);
      efficiency_bemc->SetBinError(i+1,bemc_err);

    }
  // return;
  efficiency_bemc->SetMarkerStyle(20);
  efficiency_bemc->SetMarkerSize(1);
  efficiency_bemc->SetMarkerColor(2);
  
  TH2F *h2=new TH2F("h2","",100,2,14,100,0,1);
  h2->GetXaxis()->SetTitle("P_{T} GeV/c");
  h2->GetYaxis()->SetTitle("BEMC efficiency"); 
  
  gStyle->SetOptStat(0);     
  TCanvas *c2= new TCanvas("c2","",800,600);
  efficiency_bemc->Draw("P");
  h2->Draw();
  efficiency_bemc->Draw("same");
  //  efficiency_sys->Draw("same");

  c2->SaveAs("Bemc_efficiency.pdf");
  
  gStyle->SetOptStat(1111);     
  TFile *file_bemc=new TFile("BEMC_efficiency.root","RECREATE");
  efficiency_bemc->Write();
  file_bemc->Close();
  




  // mh3hist_unlike_mass_pt[0]->Add(mh3hist_unlike_mass_pt[1],1);
  // mh3hist_like_mass_pt[0]->Add(mh3hist_like_mass_pt[1],1);


  mh3hist_unlike_partnernSigmaE_pt[0]->Sumw2();
  mh3hist_like_partnernSigmaE_pt[0]->Sumw2();


  mh3hist_unlike_partnernSigmaE_pt[0]->Add(mh3hist_unlike_partnernSigmaE_pt[1],1);
  mh3hist_like_partnernSigmaE_pt[0]->Add(mh3hist_like_partnernSigmaE_pt[1],1);

  mh3hist_unlike_partnernSigmaE_pt[0]->RebinX(4);
  mh3hist_like_partnernSigmaE_pt[0]->RebinX(4);
  
  mh2hist_unlike_pass[0]->Add(mh2hist_unlike_pass[1],1);
  mh2hist_like_pass[0]->Add(mh2hist_like_pass[1],1);

  
  for(Int_t ihist=0;ihist<mh2Nhist;ihist++)
    {
      for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
	{
	  // inv Mass
	  if(ihist==0){
	    sprintf(buf,"mass_unlike%ihist%ipt",ihist,ipt);
	    mh1hist_unlike[ihist][ipt]=(TH1F *)mh2hist_unlike_pass[0]->ProjectionY(buf,ptBinX_low_HT[ipt],ptBinX_high_HT[ipt],"eo");
	    sprintf(buf,"mass_like%ihist%ipt",ihist,ipt);
	    mh1hist_like[ihist][ipt]=(TH1F *)mh2hist_like_pass[0]->ProjectionY(buf,ptBinX_low_HT[ipt],ptBinX_high_HT[ipt],"eo");
 	    sprintf(buf,"mass_unlike_like%ihist%ipt",ihist,ipt);
	    mh1hist_unlike_like[ihist][ipt]=(TH1F *)mh1hist_unlike[ihist][ipt]->Clone(buf);
	    mh1hist_unlike_like[ihist][ipt]->Add(mh1hist_like[ihist][ipt],-1);

	    mh1hist_unlike[ihist][ipt]->GetXaxis()->SetTitle(mh1_TitleX[ihist]);
	    mh1hist_unlike[ihist][ipt]->GetYaxis()->SetTitle(mh1_TitleY[ihist]);
	  }
	  
	  if(ihist==1){
	    sprintf(buf,"partnerNsigmaE_unlike%ihist%ipt",ihist,ipt);
	    mh1hist_unlike[ihist][ipt]=(TH1F *) mh3hist_unlike_partnernSigmaE_pt[0]->ProjectionX(buf,ptBinX_low_HT[ipt],ptBinX_high_HT[ipt],1,20,"eo");
	    sprintf(buf,"partnerNsigmaE_like%ihist%ipt",ihist,ipt);
	    mh1hist_like[ihist][ipt]=(TH1F *) mh3hist_like_partnernSigmaE_pt[0]->ProjectionX(buf,ptBinX_low_HT[ipt],ptBinX_high_HT[ipt],1,20,"eo");
	    // mh1hist_like[ihist][ipt]=(TH1F *)mh2hist_like_pass[0]->ProjectionY(buf,ptBinX_low_HT[ipt],ptBinX_high_HT[ipt],"eo");
 	    // sprintf(buf,"mass_unlike_like%ihist%ipt",ihist,ipt);
	    sprintf(buf,"partnerNsigmaE_unlike_like%ihist%ipt",ihist,ipt);
	    mh1hist_unlike_like[ihist][ipt]=(TH1F *)mh1hist_unlike[ihist][ipt]->Clone(buf);
	    mh1hist_unlike_like[ihist][ipt]->Add(mh1hist_like[ihist][ipt],-1);
	    
	    mh1hist_unlike[ihist][ipt]->GetXaxis()->SetTitle(mh1_TitleX[ihist]);
	    mh1hist_unlike[ihist][ipt]->GetYaxis()->SetTitle(mh1_TitleY[ihist]);

	    // mh1hist_unlike[ihist][ipt]->RebinX(10);
	    // mh1hist_unlike_like[ihist][ipt]->RebinX(10);
	    // mh1hist_like[ihist][ipt]->RebinX(10);

	      mh1hist_unlike[ihist][ipt]->GetXaxis()->SetRangeUser(-3,3);

	  }
	  
	}
    }
  TLegend *legend  = new TLegend(0.25,0.4,0.75,0.65);
  legend ->AddEntry( mh1hist_unlike[0][0],"Unlike","lpe");
  legend ->AddEntry( mh1hist_like[0][0],"like","lpe");
  legend ->AddEntry( mh1hist_unlike_like[0][0],"Unlike - Like","lpe");
  legend ->SetBorderSize(0);
  legend ->SetTextSize(0.07);
  legend ->SetFillStyle(0);
  legend->SetTextFont(62);
  
  for(Int_t ihist=0;ihist<mh2Nhist;ihist++)
    {  
      TCanvas *c5=new TCanvas("c5","",1200,800);
      c5->Divide(4,3);
      int ipad=1;
      for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++){
	gPad->SetLogy(0);
	c5->cd(ipt+1);
	//gPad->Clear();	
	
	mh1hist_unlike[ihist][ipt]->SetLineColor(1);
	mh1hist_like[ihist][ipt]->SetLineColor(2);
	mh1hist_unlike_like[ihist][ipt]->SetLineColor(3);	
	mh1hist_unlike[ihist][ipt]->SetTitle(mh1_pT_Title_HT[ipt]);
	mh1hist_unlike[ihist][ipt]->Draw();
	mh1hist_like[ihist][ipt]->Draw("same");
	mh1hist_unlike_like[ihist][ipt]->Draw("same");
      }

      c5->cd(NpT_bins_run12_HT+1);
      legend ->Draw();
      sprintf(buf,mh2HistName_Photonic_unlike_pass[ihist]+TString(".pdf"));
      c5->SaveAs(buf);
    }
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




