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

void Fit_electron_nsigma_Mean(TH1F * a[], TH1F *b[],TH1F *c[]);
TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);

using namespace std;

const int nTrg=2;

char buf[1024];

//TFile *file=new TFile("../RootFile/Root_File_2_26/hist_2_26.root","READ");
//TFile *file=new TFile("../RootFile/Root_File_3_30/hist_3_30.root","READ");
TFile *file=new TFile("../RootFile/Root_File_5_2/hist_5_2.root","READ");
TH2F *nsigmaE_pi[nTrg];
TH2F *nsigmaE_proton[nTrg];
TH2F *nsigmaE_kaon[nTrg];

void nSigma_Electron_Calibration_partner()
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
  
  TString mh3nSigmaE_pT_Mass_unlike_name[nTrg]={"mh3nSigmaEPart_pT_Mass_unlike_passTrg0","mh3nSigmaEPart_pT_Mass_unlike_passTrg1"};
  TString mh3nSigmaE_pT_Mass_like_name[nTrg]={"mh3nSigmaEPart_pT_Mass_like_passTrg0","mh3nSigmaEPart_pT_Mass_like_passTrg1"};
  
  
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
  TH1F *Mean_Fit=new TH1F("Mean_Fit","",NpT_bins_run12_HT,pt_run12_HT);

  TH1F *Sigma=new TH1F("Sigma","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *Sigma_u=new TH1F("Sigma_u","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *Sigma_d=new TH1F("Sigma_d","",NpT_bins_run12_HT,pt_run12_HT);
  TH1F *Sigma_Fit=new TH1F("Sigma_Fit","",NpT_bins_run12_HT,pt_run12_HT);
  
  
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


      
      // if(Npad==7)
      //   Npad=1;     
      // cout<< " mean: "<< f1->GetParameter(1)<<" mean uncertainty:  "<< f1->GetParError(1)<<endl;
      
      Mean->SetBinContent(ipt+1,f1->GetParameter(1));
      Mean->SetBinError(ipt+1,f1->GetParError(1));

      Mean_u->SetBinContent(ipt+1,f1->GetParameter(1)+f1->GetParError(1));
      Mean_u->SetBinError(ipt+1,f1->GetParError(1));
      
      Mean_d->SetBinContent(ipt+1,f1->GetParameter(1)-f1->GetParError(1));
      Mean_d->SetBinError(ipt+1,f1->GetParError(1));

      Sigma->SetBinContent(ipt+1,f1->GetParameter(2));
      Sigma->SetBinError(ipt+1,f1->GetParError(2));

      Sigma_u->SetBinContent(ipt+1,f1->GetParameter(2)+f1->GetParError(2));
      Sigma_u->SetBinError(ipt+1,f1->GetParError(2));
      
      Sigma_d->SetBinContent(ipt+1,f1->GetParameter(2)-f1->GetParError(2));
      Sigma_d->SetBinError(ipt+1,f1->GetParError(2));

    }


  Mean->SetMarkerColor(1);
  Mean->SetLineColor(1);
  Mean->SetMarkerStyle(20);
  
  Sigma->SetMarkerColor(1);
  Sigma->SetLineColor(1);
  Sigma->SetMarkerStyle(24);

  c2->SaveAs("nsigmaE_c2_partner.pdf");
  c3->SaveAs("nsigmaE_c3_partner.pdf");


  gStyle->SetOptStat(00000);
  
  TCanvas *c5=new TCanvas("c5","",1200,1000);
  TH2F *h5_a=new TH2F("h5_a","",10,1.5,11,100,-1,1.5);
  h5_a->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h5_a->GetYaxis()->SetTitle("Electron nsigma Mean");
    
  
  TF1 *fmean=new TF1(TString("fmean"),"[0]",2,10);
  TF1 *fmean_u=new TF1(TString("fmean_u"),"[0]",2,10);
  TF1 *fmean_d=new TF1(TString("fmean_d"),"[0]",2,10);
  TF1 *fsigma=new TF1(TString("fsigma"),"[0]",2,10);
  TF1 *fsigma_u=new TF1(TString("fsigma_u"),"[0]",2,10);
  TF1 *fsigma_d=new TF1(TString("fsigma_d"),"[0]",2,10);  
  
  fmean_u->SetLineStyle(7);
  fmean_d->SetLineStyle(7);

  fsigma_u->SetLineStyle(7);
  fsigma_d->SetLineStyle(7);

  h5_a->Draw();
  Mean->Draw("same");
  Sigma->Draw("same");
  

  Mean->Fit(fmean,"R","same",2,10); 
  Mean_u->Fit(fmean_u,"R0","same",2,10); 
  Mean_d->Fit(fmean_d,"R0","same",2,10); 
  
  Sigma->Fit(fsigma,"R","same",2,10);
  Sigma_u->Fit(fsigma_u,"R0","same",2,10);
  Sigma_d->Fit(fsigma_d,"R0","same",2,10);
  
  fmean_u->SetParameter(0,fmean->GetParameter(0)+0.04);
  fmean_d->SetParameter(0,fmean->GetParameter(0)-0.04);

  fsigma_u->SetParameter(0,fsigma->GetParameter(0)+0.035);
  fsigma_d->SetParameter(0,fsigma->GetParameter(0)-0.035);
  

  for(Int_t ipt=0;ipt<NpT_bins_run12_HT;ipt++)
    {
      Mean_Fit->SetBinContent(ipt+1,fmean->GetParameter(0));
      Mean_Fit->SetBinError(ipt+1,0.04);
      
      Sigma_Fit->SetBinContent(ipt+1,fsigma->GetParameter(0));
      Sigma_Fit->SetBinError(ipt+1,0.035);
    }

  fmean_u->Draw("same");
  fmean_d->Draw("same");
  fsigma_u->Draw("same");
  fsigma_d->Draw("same");
  

  // sprintf(buf,"pol0:  %2.2f +- %4.5f",fmean->GetParameter(0),fmean->GetParameter(0));
  // drawLatex(0.45,0.75,buf,70,0.035,2);
  // sprintf(buf,"pol0:  %2.2f +- %4.5f",fsigma->GetParameter(0),fsigma->GetParameter(0));
  // drawLatex(0.45,0.55,buf,70,0.035,2);
  
  sprintf(buf,"Electron Mean: pol0:  %2.2f +- %4.3f",fmean->GetParameter(0),0.04);
  drawLatex(0.15,0.35,buf,70,0.035,2);
  sprintf(buf,"Electron Sigma: pol0:  %2.2f +- %4.3f",fsigma->GetParameter(0),0.035);
  drawLatex(0.15,0.75,buf,70,0.035,2);

  c5->SaveAs("Nsigma_mean_sigma_partner.pdf");


   // nsigmaE_pi[0]->RebinX(2);
   // nsigmaE_proton[0]->RebinX(2);
   // nsigmaE_kaon[0]->RebinX(2);
   
  TCanvas *c6=new TCanvas("c6","",1000,800);
  TH2F *h6=new TH2F("h6","",10,2,10,100,-7,1);
  h6->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h6->GetYaxis()->SetTitle("");
  
  
  TProfile *nsigma_E_pi=  (TProfile *) nsigmaE_pi[0]->ProfileX("nsigma_E_pi",1,-1,"");
  TProfile *nsigma_E_proton=  (TProfile *) nsigmaE_proton[0]->ProfileX("nsigma_E_proton",1,-1,"");
  TProfile *nsigma_E_kaon=  (TProfile *) nsigmaE_kaon[0]->ProfileX("nsigma_E_kaon",1,-1,"");
  
  nsigma_E_pi->SetLineColor(1);
  nsigma_E_proton->SetLineColor(3);
  nsigma_E_kaon->SetLineColor(4);
  
  
  // nsigma_E_pi->Draw("same");
  // nsigma_E_proton->Draw("same");
  // nsigma_E_kaon->Draw("same");
  // Mean->Draw("same");
  
  TH1F * nsigma_e_pi=(TH1F *) nsigma_E_pi->Rebin(NpT_bins_run12_HT,"nsigma_e_pi",pt_run12_HT); 
  TH1F * nsigma_e_proton=(TH1F *) nsigma_E_proton->Rebin(NpT_bins_run12_HT,"nsigma_E_proton",pt_run12_HT); 
  TH1F * nsigma_e_kaon=(TH1F *) nsigma_E_kaon->Rebin(NpT_bins_run12_HT,"nsigma_E_kaon",pt_run12_HT); 
  
   nsigma_e_pi->SetMarkerStyle(24);
   nsigma_e_pi->SetMarkerColor(6);
   nsigma_e_pi->SetLineColor(6);
   
   nsigma_e_proton->SetMarkerStyle(24);
   nsigma_e_proton->SetMarkerColor(4);
   nsigma_e_proton->SetLineColor(4);
   
   nsigma_e_kaon->SetMarkerStyle(24);
   nsigma_e_kaon->SetMarkerColor(3);
   nsigma_e_kaon->SetLineColor(3);
   
   h6->Draw("");   
   Mean->Draw("same");
   nsigma_e_pi->Draw("same");
   nsigma_e_proton->Draw("same");
   nsigma_e_kaon->Draw("same");
   
   
   TLegend *legend_6  = new TLegend(0.15,0.55,0.4,0.75);
   legend_6 ->AddEntry(Mean,"n#sigma_{e}","lp");
   legend_6 ->AddEntry(nsigma_e_pi,"n#sigma_{e}-n#sigma_{#pi}","lp");
   legend_6 ->AddEntry(nsigma_e_proton,"n#sigma_{e}-n#sigma_{proton}","lp");
   legend_6 ->AddEntry(nsigma_e_kaon,"n#sigma_{e}-n#sigma_{kaon}","lp");
   
   legend_6 ->SetBorderSize(0);
   legend_6 ->SetTextSize(0.035);
   legend_6 ->SetFillStyle(0);
   legend_6 ->SetTextFont(62);
   legend_6->Draw();
   c6->SaveAs("NsigmaE_hadron_partner.pdf");

   
   TH1F *pi_nsigmaE= (TH1F *) Mean->Clone("pi_nsigmaE");
   TH1F *proton_nsigmaE= (TH1F *) Mean->Clone("proton_nsigmaE");
   TH1F *kaon_nsigmaE= (TH1F *) Mean->Clone("kaon_nsigmaE");

   TCanvas *c7=new TCanvas("c7","",1000,800);
   TH2F *h7=new TH2F("h7","",10,2,10,100,-7,1);
   h7->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h7->GetYaxis()->SetTitle("");
   
   pi_nsigmaE->Add(nsigma_e_pi,1);
   proton_nsigmaE->Add(nsigma_e_proton,1);
   kaon_nsigmaE->Add(nsigma_e_kaon,1);


   pi_nsigmaE->SetMarkerStyle(24);
   pi_nsigmaE->SetMarkerColor(6);
   pi_nsigmaE->SetLineColor(6);
  
   proton_nsigmaE->SetMarkerStyle(24);
   proton_nsigmaE->SetMarkerColor(4);
   proton_nsigmaE->SetLineColor(4);
   
   kaon_nsigmaE->SetMarkerStyle(24);
   kaon_nsigmaE->SetMarkerColor(3);
   kaon_nsigmaE->SetLineColor(3);


   h7->Draw();
   Mean->Draw("same");
   pi_nsigmaE->Draw("same");
   proton_nsigmaE->Draw("same");
   kaon_nsigmaE->Draw("same");

   
   TLegend *legend_7  = new TLegend(0.15,0.45,0.4,0.65);
   legend_7 ->AddEntry(Mean,"n#sigma_{e}","lp");
   legend_7 ->AddEntry(pi_nsigmaE,"n#sigma_{#pi}","lp");
   legend_7 ->AddEntry(proton_nsigmaE,"n#sigma_{proton}","lp");
   legend_7 ->AddEntry(kaon_nsigmaE,"n#sigma_{kaon}","lp");

   legend_7 ->SetBorderSize(0);
   legend_7 ->SetTextSize(0.035);
   legend_7 ->SetFillStyle(0);
   legend_7 ->SetTextFont(62);
   legend_7->Draw();
   
   c7->SaveAs("hadron_nsigmaE_calibration_partner.pdf");

   TFile *File=new TFile("nSigmaE_calibration_partner.root","RECREATE");
   
   Mean->Write();
   Sigma->Write();
   Mean_Fit->Write();
   Sigma_Fit->Write();

   Mean_u->Write();
   Sigma_u->Write();

   Mean_d->Write();
   Sigma_d->Write();

   fmean->Write();
   fsigma->Write();

   fmean_u->Write();
   fsigma_u->Write();

   fmean_d->Write();
   fsigma_d->Write();

   pi_nsigmaE->Write();
   proton_nsigmaE->Write();
   kaon_nsigmaE->Write();

   File->Close();
  
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




