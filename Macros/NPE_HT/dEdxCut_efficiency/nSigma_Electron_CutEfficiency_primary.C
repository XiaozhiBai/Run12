/*-------------------------------------------------------------- 

  run14 HT Npe analyia 

  This macro is for the photo electron nsigma electron calibration  for the nsigma Electran calibration in different Pt bin
  ---------------------------------------------------------------

  Mon Feb 15 00:40:29 CST 2016

  by Xiaozhi

*/
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
#include "TProfile.h"
#include "../mBinning_HT.h"
#include "TVirtualFitter.h"
void Fit_electron_nsigma_Mean(TH1F * a[], TH1F *b[],TH1F *c[]);
TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);

using namespace std;

const int nTrg=2;

char buf[1024];

TFile *file=new TFile("../RootFile/Root_File_5_1/hist_5_1.root","READ");

TH2F *nsigmaE_pi[nTrg];
TH2F *nsigmaE_proton[nTrg];
TH2F *nsigmaE_kaon[nTrg];
ofstream outdata("nsigemeanAndsigmaeff.dat");
void nSigma_Electron_CutEfficiency_primary()
{
  TH3F::SetDefaultSumw2();  
  TH2F::SetDefaultSumw2();  
  TH1F::SetDefaultSumw2();  
 
  // gStyle->SetOptStat(00000);
  // gStyle->SetTitleSize(0.05,"XY");
  // gStyle->SetTitleFontSize(0.06);
  // gStyle->SetTitleOffset(1,"X");
  // gStyle->SetTitleOffset(1,"Y");
  // gStyle->SetPadTopMargin(0.13);
  // gStyle->SetPadRightMargin(0.02);
  // gStyle->SetPadBottomMargin(0.13);
  // gStyle->SetPadLeftMargin(0.1); 
  
  TH3F *mh3nSigmaE_pT_Mass_unlike[nTrg];
  TH3F *mh3nSigmaE_pT_Mass_like[nTrg];
  TH3F *mh3nSigmaE_pT_Mass_unlike_like[nTrg];
  
  
// TString mh3nSigmaE_pT_Mass_unlike_name[nTrg]={"mh3nSigmaEPart_pT_Mass_unlike_passTrg0","mh3nSigmaEPart_pT_Mass_unlike_passTrg1"};
//   TString mh3nSigmaE_pT_Mass_like_name[nTrg]={"mh3nSigmaEPart_pT_Mass_like_passTrg0","mh3nSigmaEPart_pT_Mass_like_passTrg1"};
  TString mh3nSigmaE_pT_Mass_unlike_name[nTrg]={"mh3nSigmaE_pT_Mass_unlikeTrg0","mh3nSigmaE_pT_Mass_unlikeTrg1"};
  TString mh3nSigmaE_pT_Mass_like_name[nTrg]={"mh3nSigmaE_pT_Mass_likeTrg0","mh3nSigmaE_pT_Mass_likeTrg1"};

  // TH2F *mh2NisgmaE_unlike;
  // TH2F *mh2NisgmaE_like;

TString nsigmaE_pi_name[nTrg]={"mh2Pion_nSigmaElecTrg0","mh2Pion_nSigmaElecTrg1"};
TString nsigmaE_proton_name[nTrg]={"mh2Proton_nSigmaElecTrg0","mh2Proton_nSigmaElecTrg1"};
TString nsigmaE_kaon_name[nTrg]={"mh2Kaon_nSigmaElecTrg0","mh2Kaon_nSigmaElecTrg1"};
  
  
  TH1F *mh1nSigmaE_unlike[NpT_bins_run12_HT];
  TH1F *mh1NisgmaE_like[NpT_bins_run12_HT];
  TH1F *mh1nSigmaE_unlike_like[NpT_bins_run12_HT];

 

  
  // mh2NisgmaE_unlike=(TH2F *) file->Get("mh2partnerNsigmaeVspt_unlikeTrg1");
  // mh2NisgmaE_like=(TH2F *) file->Get("mh2partnerNsigmaeVspt_likeTrg1");

  for(Int_t iTrg=0;iTrg<nTrg;iTrg++)
    {
      mh3nSigmaE_pT_Mass_unlike[iTrg]=(TH3F *) file->Get(mh3nSigmaE_pT_Mass_unlike_name[iTrg]);
      mh3nSigmaE_pT_Mass_like[iTrg]=(TH3F *) file->Get(mh3nSigmaE_pT_Mass_like_name[iTrg]);
      sprintf(buf,"mh3nSigmaE_pT_Mass_unlike_likeTrg%i",iTrg);
      mh3nSigmaE_pT_Mass_unlike_like[iTrg]=(TH3F *)mh3nSigmaE_pT_Mass_unlike[iTrg]->Clone(buf);
      mh3nSigmaE_pT_Mass_unlike_like[iTrg]->Add(mh3nSigmaE_pT_Mass_like[iTrg],-1);
      nsigmaE_pi[iTrg]=(TH2F *) file->Get(nsigmaE_pi_name[iTrg]);
      nsigmaE_proton[iTrg]=(TH2F *) file->Get(nsigmaE_proton_name[iTrg]);
      nsigmaE_kaon[iTrg]=(TH2F *) file->Get(nsigmaE_kaon_name[iTrg]);
    }
  
  mh3nSigmaE_pT_Mass_unlike[0]->Sumw2();
  mh3nSigmaE_pT_Mass_like[0]->Sumw2();
  mh3nSigmaE_pT_Mass_unlike_like[0]->Sumw2();
  nsigmaE_pi[0]->Sumw2();
  nsigmaE_proton[0]->Sumw2();
  nsigmaE_kaon[0]->Sumw2();


  mh3nSigmaE_pT_Mass_unlike[0]->Add( mh3nSigmaE_pT_Mass_unlike[1],1);
  mh3nSigmaE_pT_Mass_like[0]->Add( mh3nSigmaE_pT_Mass_like[1],1);
  mh3nSigmaE_pT_Mass_unlike_like[0]->Add( mh3nSigmaE_pT_Mass_unlike_like[1],1);
  nsigmaE_pi[0]->Add(nsigmaE_pi[1],1);
  nsigmaE_proton[0]->Add(nsigmaE_proton[1],1);
  nsigmaE_kaon[0]->Add(nsigmaE_kaon[1],1);



  mh3nSigmaE_pT_Mass_unlike[0]->RebinX(2);
  mh3nSigmaE_pT_Mass_like[0]->RebinX(2);
  mh3nSigmaE_pT_Mass_unlike_like[0]->RebinX(2);
  
  
  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      sprintf(buf,"nSigmaE_unlike_pt%i",ipt);
      mh1nSigmaE_unlike[ipt]=(TH1F *) mh3nSigmaE_pT_Mass_unlike[0]->ProjectionX(buf,ptBinX_low_HT[ipt],ptBinX_high_HT[ipt],1,20);
      sprintf(buf,"nSigmaE_like_pt%i",ipt);
      mh1NisgmaE_like[ipt]=(TH1F *) mh3nSigmaE_pT_Mass_like[0]->ProjectionX(buf,ptBinX_low_HT[ipt],ptBinX_high_HT[ipt],1,20);

      sprintf(buf,"nSigmaE_unlike_like_pt%i",ipt);
      mh1nSigmaE_unlike_like[ipt]=(TH1F *)   mh1nSigmaE_unlike[ipt]->Clone(buf);
      mh1nSigmaE_unlike_like[ipt]->Add(mh1NisgmaE_like[ipt],-1);
    }
  Fit_electron_nsigma_Mean(mh1nSigmaE_unlike,mh1NisgmaE_like,mh1nSigmaE_unlike_like);
}
void Fit_electron_nsigma_Mean(TH1F *mh1NisgmaE_unlike[NpT_bins_run12_HT],TH1F *mh1NisgmaE_like[NpT_bins_run12_HT],TH1F *mh1NisgmaE_unlike_like[NpT_bins_run12_HT])
{
  gStyle->SetOptFit(1111);
  
  TH1F *Mean=new TH1F("Mean","",    NpT_bins_run12_HT,pt_run12_HT);
  TH1F *Mean_u=new TH1F("Mean_u","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *Mean_d=new TH1F("Mean_d","",NpT_bins_run12_HT,pt_run12_HT);
  
  TH1F *Sigma=new TH1F("Sigma","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *Sigma_u=new TH1F("Sigma_u","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *Sigma_d=new TH1F("Sigma_d","",NpT_bins_run12_HT,pt_run12_HT);
  
  
  TCanvas *c2=new TCanvas("c2","",1200,1000);
  TCanvas *c3=new TCanvas("c3","",1200,1000);
  
  c2->Divide(3,2,0.001,0.001);
  c3->Divide(3,2,0.001,0.001);
  
  int Npad=1;
  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      TF1 *f1 = new TF1(TString("f1"),"gaus",-5,5);      
      if(ipt<6)
        {
          c2->cd(Npad++);
          gPad->SetLogy(0);
        }
      else 
        {
          c3->cd(Npad++);
          gPad->SetLogy(0);
        }
      if(Npad==7) Npad=1;

      mh1NisgmaE_unlike_like[ipt]->SetTitle(mh1_pT_Title_HT[ipt]);
      mh1NisgmaE_unlike_like[ipt]->GetXaxis()->SetTitle("nSigmaE");
      mh1NisgmaE_unlike_like[ipt]->GetYaxis()->SetTitle("Counts");

      mh1NisgmaE_unlike_like[ipt]->GetYaxis()->SetRangeUser(1,1.2*mh1NisgmaE_unlike[ipt]->GetMaximum());
      mh1NisgmaE_unlike_like[ipt]->GetXaxis()->SetRangeUser(-3,3);      


      mh1NisgmaE_unlike[ipt]->SetLineColor(1);
      mh1NisgmaE_like[ipt]->SetLineColor(3);
      mh1NisgmaE_unlike_like[ipt]->SetLineColor(4);


      mh1NisgmaE_unlike_like[ipt]->Fit(f1,"R","same",-3,3);
      mh1NisgmaE_unlike[ipt]->Draw("same");
      mh1NisgmaE_unlike_like[ipt]->Draw("same");
      mh1NisgmaE_like[ipt]->Draw("same");

      TLegend *legend = new TLegend(0.15,0.65,0.4,0.8);
      legend->AddEntry(mh1NisgmaE_unlike[ipt],"Unlike ","lpe");
      legend->AddEntry(mh1NisgmaE_like[ipt],"Like ","lpe");
      legend->AddEntry(mh1NisgmaE_unlike_like[ipt],"Unlike - Like ","lpe");
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.035);
      legend->Draw("same");

      TVirtualFitter * fitter = TVirtualFitter::GetFitter();
      assert(fitter != 0);
      double * cov =fitter->GetCovarianceMatrix();
      
      cout<<cov[0]<<" "<<sqrt(cov[4])<<" "<<sqrt(cov[8])<<" "<<cov[7]<<" meanerr="<<f1->GetParError(1)<<"sigmaerr "<<f1->GetParError(2)<<endl;
      outdata << (pt_run12_HT[ipt]+pt_run12_HT[ipt+1])/2       << "  " << 0.5*(pt_run12_HT[ipt+1]-pt_run12_HT[ipt])    << "  "
	      << f1->GetParameter(1)     << "  " << cov[4]    << "  " 
	      << f1->GetParameter(2)     << "  " << cov[8]    << "  "<<cov[7]<< endl; 
      

      // if(Npad==7)
      //   Npad=1;     
      // cout<< " mean: "<< f1->GetParameter(1)<<" mean uncertainty:  "<< f1->GetParError(1)<<endl;
      


    }



  // Mean->SetMarkerColor(1);
  // Mean->SetLineColor(1);
  // Mean->SetMarkerStyle(20);
  
  // Sigma->SetMarkerColor(1);
  // Sigma->SetLineColor(1);
  // Sigma->SetMarkerStyle(24);

  c2->SaveAs("nsigmaE_c2_partner.pdf");
  c3->SaveAs("nsigmaE_c3_partner.pdf");

 
  
}
TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
  TLatex *latex = new TLatex(x,y,text);                                        
  latex->SetNDC();  
  // latex->SetTextFont(textFont);                                                
  latex->SetTextSize(textSize);
  latex->SetTextColor(colorIndex);
  latex->Draw("same");   
  return latex; 
}




