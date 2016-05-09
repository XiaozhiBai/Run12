/*--------------------Xiaozhi------------------------------------ 

---------------------------------------------------------------
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
#include "TRandom3.h"
#include "../mBinning_MB.h"

using namespace std;



void Draw_purity(TH1D * );
void Draw_purity_sys(TH1D *,TH1D *,TH1D * );

Double_t Get_sts_Error(TH1F *const inclusive,Int_t, TF1 *, TF1 *,TF1 *,TF1 *,TF1 *,Int_t,Int_t,Double_t ,Double_t);
//        mPurity_sts= Get_sts_Error(nsigmaE_inclusive[ipt],ipt,total_3,g1,g2,g3,g4,3,0,purity_mean,mPurity_sys);
void Fit_electron_purity(TH1F * a[]);

TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
Double_t Fit_purity(TH1F * const inclusive,Int_t ipt,Float_t E_mean_err,Float_t E_sigma_err, TF1 *,TF1 *, TF1 *,TF1 *,TF1*, int ,Int_t);

//TFile *File=new TFile("nSigmaE_calibration.root","READ");
TFile *File=new TFile("nSigmaE_calibration_primary.root","READ");
//TFile *File=new TFile("nSigmaE_calibration_partner.root","READ");


TH1F *pi_Mean=(TH1F *) File->Get("pi_nsigmaE");
TH1F *proton_Mean=(TH1F *) File->Get("proton_nsigmaE");
TH1F *kaon_Mean=(TH1F *) File->Get("kaon_nsigmaE");

// TH1F *electron_Sigma=(TH1F *) File->Get("Sigma");
// TH1F *electron_Mean=(TH1F *) File->Get("Mean");  

TH1F *electron_Sigma=(TH1F *) File->Get("Sigma_Fit");
TH1F *electron_Mean=(TH1F *) File->Get("Mean_Fit");  


//
TFile *File_Tof=new TFile("TOF_ONlY_DATA.root","READ");

TH1F *pi_mean_tof=(TH1F *) File_Tof->Get("mh1_Pion_mean");
TH1F *proton_mean_tof=(TH1F *) File_Tof->Get("mh1_Proton_mean");
TH1F *kaon_mean_tof=(TH1F *) File_Tof->Get("mh1_Kaon_mean");

TH1F *pi_sigma_tof=(TH1F *) File_Tof->Get("mh1_Pion_sigma");
TH1F *proton_sigma_tof=(TH1F *) File_Tof->Get("mh1_Proton_sigma");
TH1F *kaon_sigma_tof=(TH1F *) File_Tof->Get("mh1_Kaon_sigma");

char buf[1024];


TRandom3 *gRnd= new TRandom3(0);
int purity_MB_lowpT(){

TH3F::SetDefaultSumw2();  
TH2F::SetDefaultSumw2();  
TH1F::SetDefaultSumw2();  

//  gStyle->SetOptStat(00000);
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1,"Y");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.1); 

  //  electron_Mean->Draw();

  TH1F *nsigmaE_inclusive[NpT_bins_run12_MB];

  TFile * inFile=new TFile("../RootFile/Root_File_5_1/hist_5_1.root","READ");
  TH2F * nsigmaE_inclusive_Trig2=(TH2F *) inFile->Get("mh2nSigmaElecTrg2"); 



  cout<< "  GetFile"<<endl;
  for(Int_t ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      sprintf(buf,"nSigmaE_pt%i",ipt);
      nsigmaE_inclusive[ipt]=(TH1F *) nsigmaE_inclusive_Trig2->ProjectionX(buf,ptBinX_low[ipt],ptBinX_high[ipt]);
      nsigmaE_inclusive[ipt]->Rebin(3);
      nsigmaE_inclusive[ipt]->SetMarkerStyle(24);
      nsigmaE_inclusive[ipt]->SetMarkerSize(0.5);


    }
  
   
    Fit_electron_purity(nsigmaE_inclusive);
  return 0;
  

}
void Fit_electron_purity(TH1F * nsigmaE_inclusive[])
{
  TH1D *purity_MB=new TH1D("purity_MB","",NpT_bins_run12_MB,pt_run12_MB);
  TH1D *purity_oneSigma=new TH1D("purity_oneSigma","",NpT_bins_run12_MB,pt_run12_MB);
  TH1D *purity_twoSigma=new TH1D("purity_twoSigma","",NpT_bins_run12_MB,pt_run12_MB);
  TH1D *purity_threeSigma=new TH1D("purity_threeSigma","",NpT_bins_run12_MB,pt_run12_MB);

   gStyle->SetOptFit(1111);
  TCanvas *c2=new TCanvas("c2","",1200,800);
  TCanvas *c3=new TCanvas("c3","",1200,800);
  TCanvas *c22=new TCanvas("c22","",1200,800);
  TCanvas *c33=new TCanvas("c33","",1200,800);

    TCanvas *c333=new TCanvas("c333","",1200,800);
  
  c2->Divide(3,2);
  c3->Divide(3,2);
  
  c22->Divide(3,2);
  c33->Divide(3,2);

  
  int Npad=1;

  for(Int_t ipt=20;ipt<NpT_bins_run12_MB-5;ipt++)
     {
       sprintf(buf,"%ipt",ipt);
       TF1 *total_3 =new TF1(buf+TString("total_3"),"[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)+[9]*TMath::Gaus(x,[10],[11],1)",-15,15);// fit the total 
      TF1 *  g1=new TF1(buf+TString("g1"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15);   
      TF1 *  g2=new TF1(buf+TString("g2"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15);   
      TF1 *  g3=new TF1(buf+TString("g3"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15);
      TF1 *  g4=new TF1(buf+TString("g4"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15); 
      //       total_3->SetParNames("#pi C", "#pi #mu", "#pi #sigma", "proton C", "proton #mu", "proton #sigma","kaon C","kaon #mu","kaon sigma");
      total_3->SetLineColor(kBlue);
      

      g1->SetLineColor(3);
      g2->SetLineColor(6);
      g3->SetLineColor(7);
      g4->SetLineColor(2);

      g1->SetLineStyle(7);
      g2->SetLineStyle(7);
      g3->SetLineStyle(7);


      
      
      // if(ipt<23)
      // 	{
      // 	  c2->cd(Npad++);
      // 	  gPad->SetLogy(1);
      // 	}
      // else (ipt<29) 
      //   {
      //     c22->cd(Npad++);
      //     gPad->SetLogy(1);
      //   }

      // else if(ipt<18) 
      //   {
      //     c22->cd(Npad++);
      //     gPad->SetLogy(1);
      //   }

      if(ipt==18)
	c333->cd();
      else continue;
      

      //      if(Npad==7) Npad=1;
      nsigmaE_inclusive[ipt]->SetTitle(mh1_pT_Title[ipt]);
      nsigmaE_inclusive[ipt]->GetXaxis()->SetTitle("nSigmaE");
      nsigmaE_inclusive[ipt]->GetYaxis()->SetTitle("Counts");
      
      
      nsigmaE_inclusive[ipt]->GetXaxis()->SetRangeUser(-15,15);  









      

      // if(ipt<21)
      // 	{
      // 	  total_3->FixParameter(1,pi_mean_tof->GetBinContent(ipt+1));
      
      // 	  total_3->SetParLimits(2,0.7,0.95);

      // 	  if(ipt==17)
      // 	    total_3->SetParLimits(3,50,4e3);
      // 	  else
      // 	    total_3->SetParLimits(3,50,2e8);
      
      // 	  total_3->SetParLimits(4,proton_mean_tof->GetBinContent(ipt+1)-0.5,proton_mean_tof->GetBinContent(ipt+1)+0.5);
      
      // 	  total_3->SetParLimits(5,0.7,1);
      // 	  total_3->SetParLimits(7,kaon_mean_tof->GetBinContent(ipt+1)-0.5,kaon_mean_tof->GetBinContent(ipt+1)+0.5);
      // 	  total_3->SetParLimits(8,0.7,1);
      // 	}
      // else {

      // 	total_3->FixParameter(1,pi_mean_tof->GetBinContent(ipt+1));
      // 	total_3->FixParameter(2,pi_sigma_tof->GetBinContent(ipt+1));


      // 	total_3->SetParLimits(3,100,1e7);
      // 	total_3->FixParameter(4,proton_mean_tof->GetBinContent(ipt+1));
      // 	//   if(ip==)
      // 	total_3->SetParLimits(5,0.8,1);
      // 	// else
      // 	//   total_3->FixParameter(5,proton_sigma_tof->GetBinContent(ipt+1));
  
      // 	total_3->SetParLimits(6,100,1e6);

      // 	total_3->FixParameter(7,kaon_mean_tof->GetBinContent(ipt+1));
      // 	total_3->FixParameter(8,kaon_sigma_tof->GetBinContent(ipt+1));
      // }



      // //electron mean
      // total_3->SetParameter(10,electron_Mean->GetBinContent(ipt+1));
      // total_3->SetParLimits(10,electron_Mean->GetBinContent(ipt+1)-electron_Mean->GetBinError(ipt+1),electron_Mean->GetBinContent(ipt+1)+electron_Mean->GetBinError(ipt+1));
      // //electron sigma
      // total_3->SetParameter(11,electron_Sigma->GetBinContent(ipt+1));
      // total_3->SetParLimits(11,electron_Sigma->GetBinContent(ipt+1)-electron_Sigma->GetBinError(ipt+1),electron_Sigma->GetBinContent(ipt+1)+electron_Sigma->GetBinError(ipt+1));
	  


      // if(ipt<21)
      // 	{
	  total_3->FixParameter(1,pi_mean_tof->GetBinContent(ipt+1));
	  total_3->SetParameter(2,0.85);
	   
	   total_3->SetParLimits(2,0.7,0.95);


	  // if(ipt==17)
	  //   total_3->SetParLimits(3,50,4e3);
	  //	  else
	   //   total_3->SetParLimits(3,50,2e8);

	   total_3->SetParameter(4,proton_mean_tof->GetBinContent(ipt+1));
		  
	  total_3->SetParLimits(4,proton_mean_tof->GetBinContent(ipt+1)-0.5,proton_mean_tof->GetBinContent(ipt+1)+0.5);
	  total_3->SetParameter(5,0.85);
	  total_3->SetParLimits(5,0.7,1);

	  total_3->SetParameter(7,kaon_mean_tof->GetBinContent(ipt+1));
	  
	  total_3->SetParLimits(7,kaon_mean_tof->GetBinContent(ipt+1)-0.5,kaon_mean_tof->GetBinContent(ipt+1)+0.5);
	  total_3->SetParLimits(8,0.7,1);
	  //	}
      // else {
	
      // 	total_3->FixParameter(1,pi_mean_tof->GetBinContent(ipt+1));
      // 	total_3->FixParameter(2,pi_sigma_tof->GetBinContent(ipt+1));
	
	
      // 	total_3->SetParLimits(3,100,1e7);
      // 	total_3->FixParameter(4,proton_mean_tof->GetBinContent(ipt+1));
      // 	//   if(ip==)
      // 	total_3->SetParLimits(5,0.8,1);
      // 	// else
      // 	//   total_3->FixParameter(5,proton_sigma_tof->GetBinContent(ipt+1));
	
      // 	total_3->SetParLimits(6,100,1e6);
	
      // 	total_3->FixParameter(7,kaon_mean_tof->GetBinContent(ipt+1));
      // 	total_3->FixParameter(8,kaon_sigma_tof->GetBinContent(ipt+1));
      // }
      

      
  //electron mean
      total_3->SetParameter(10,electron_Mean->GetBinContent(ipt+1));
      total_3->SetParLimits(10,electron_Mean->GetBinContent(ipt+1)-electron_Mean->GetBinError(ipt+1),electron_Mean->GetBinContent(ipt+1)+electron_Mean->GetBinError(ipt+1));
      //electron sigma
      total_3->SetParameter(11,electron_Sigma->GetBinContent(ipt+1));
      total_3->SetParLimits(11,electron_Sigma->GetBinContent(ipt+1)-electron_Sigma->GetBinError(ipt+1),electron_Sigma->GetBinContent(ipt+1)+electron_Sigma->GetBinError(ipt+1));
      

      
      nsigmaE_inclusive[ipt]->Fit(total_3,"R","",-10,4);

      
	        continue;            
      g1->SetParameter(0,total_3->GetParameter(0));
      g1->SetParameter(1,total_3->GetParameter(1));
      g1->SetParameter(2,total_3->GetParameter(2));
      
      g2->SetParameter(0,total_3->GetParameter(3));
      g2->SetParameter(1,total_3->GetParameter(4));
      g2->SetParameter(2,total_3->GetParameter(5));
      
      g3->SetParameter(0,total_3->GetParameter(6));
      g3->SetParameter(1,total_3->GetParameter(7));
      g3->SetParameter(2,total_3->GetParameter(8));

      g4->SetParameter(0,total_3->GetParameter(9));
      g4->SetParameter(1,total_3->GetParameter(10));
      g4->SetParameter(2,total_3->GetParameter(11));
      
       
      g1->Draw("same");
      g2->Draw("same");
      g3->Draw("same");
      g4->Draw("same");
      total_3->Draw("same");      
      
      TLegend *legend = new TLegend(0.15,0.65,0.35,0.8);
      legend->AddEntry(g1,"#pi ","lp");
      legend->AddEntry(g2,"proton ","lp");
      legend->AddEntry(g3,"kaon","lp");
      legend->AddEntry(g4,"e","lp");

      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.055);
      legend ->SetTextFont(62);       
      legend->Draw("same");

      c2->SaveAs("purity_fit_c2.pdf");
           c333->SaveAs("purity_fit_c333.pdf");
      c22->SaveAs("purity_fit_c22.pdf");


           continue;
      
            
      Double_t mPurity=g4->Integral(-1,3)/total_3->Integral(-1,3);
      purity_MB->SetBinContent(ipt+1,mPurity);
          
      
      
      Double_t mPurity_temp_oneSigma= Fit_purity(nsigmaE_inclusive[ipt],ipt,electron_Mean->GetBinError(ipt+1),electron_Sigma->GetBinError(ipt+1),total_3,g1,g2,g3,g4,1,1); 
      Double_t mPurity_temp_twoSigma= Fit_purity(nsigmaE_inclusive[ipt],ipt,electron_Mean->GetBinError(ipt+1),electron_Sigma->GetBinError(ipt+1),total_3,g1,g2,g3,g4,2,1); 
      Double_t mPurity_temp_threeSigma= Fit_purity(nsigmaE_inclusive[ipt],ipt,electron_Mean->GetBinError(ipt+1),electron_Sigma->GetBinError(ipt+1),total_3,g1,g2,g3,g4,3,1);
            
      purity_oneSigma->SetBinContent(ipt+1,mPurity_temp_oneSigma);
      purity_twoSigma->SetBinContent(ipt+1,mPurity_temp_twoSigma);
      purity_threeSigma->SetBinContent(ipt+1,mPurity_temp_threeSigma);
      
      
      
      Double_t purity_mean=(mPurity_temp_oneSigma+mPurity_temp_twoSigma+mPurity_temp_threeSigma)/3.;

      
      Double_t mPurity_sys=0;
      if(abs(purity_mean-mPurity_temp_oneSigma)>mPurity_sys)
      	mPurity_sys=abs(purity_mean-mPurity_temp_oneSigma);
      
      if(abs(purity_mean-mPurity_temp_twoSigma)>mPurity_sys)
      	mPurity_sys=abs(purity_mean-mPurity_temp_twoSigma);
      
      if(abs(purity_mean-mPurity_temp_threeSigma)>mPurity_sys)
      	mPurity_sys=abs(purity_mean-mPurity_temp_threeSigma);
      



      Double_t mPurity_sts=0.0;

      // mPurity_sts= Get_sts_Error(nsigmaE_inclusive[ipt],ipt,total_3,g1,g2,g3,g4,3,0,purity_mean,mPurity_sys);

      purity_MB->SetBinContent(ipt+1,purity_mean);
      purity_MB->SetBinError(ipt+1,sqrt(mPurity_sys*mPurity_sys+mPurity_sts*mPurity_sts));

      Fit_purity(nsigmaE_inclusive[ipt],ipt,electron_Mean->GetBinError(ipt+1),electron_Sigma->GetBinError(ipt+1),total_3,g1,g2,g3,g4,1,1); 
      Fit_purity(nsigmaE_inclusive[ipt],ipt,electron_Mean->GetBinError(ipt+1),electron_Sigma->GetBinError(ipt+1),total_3,g1,g2,g3,g4,2,1); 
      Fit_purity(nsigmaE_inclusive[ipt],ipt,electron_Mean->GetBinError(ipt+1),electron_Sigma->GetBinError(ipt+1),total_3,g1,g2,g3,g4,3,1); 
      
       TFile *file_purity=new TFile("purity_MB_lowpT.root","RECREATE");
       purity_MB->Write();
       purity_oneSigma->Write();
       purity_twoSigma->Write();
       purity_threeSigma->Write();
       file_purity->Close();
    }

  //  Draw_purity_sys(purity_oneSigma,purity_twoSigma,purity_threeSigma);

   Draw_purity(purity_MB);
}
Double_t Fit_purity(TH1F * const inclusive,Int_t ipt,Float_t E_mean_err,Float_t E_sigma_err,TF1 *total_3,TF1 * g1,TF1 *g2,TF1 *g3, TF1 *g4,Int_t Sigma_flag,Int_t Draw_flag)
{
  TH1F *Inclusive= (TH1F *)  inclusive->Clone("Inclusive");
  TCanvas * c4=new TCanvas("c4","",1200,800);
  gPad->SetLogy();
  total_3->SetParNames("#pi C", "#pi #mu", "#pi #sigma", "proton C", "proton #mu", "proton #sigma","kaon C","kaon #mu","kaon sigma");  
  total_3->SetLineColor(kBlue);

  g1->SetLineColor(3);
  g2->SetLineColor(6);
  g3->SetLineColor(7);
   
  g4->SetLineColor(2);
   
  g1->SetLineStyle(7);
  g2->SetLineStyle(7);
  g3->SetLineStyle(7);
   
  
  
  Inclusive->SetTitle(mh1_pT_Title[ipt]);
  Inclusive->GetXaxis()->SetTitle("nSigmaE");
  Inclusive->GetYaxis()->SetTitle("Counts");
      
 
  // total_3->FixParameter(1,pi_mean_tof->GetBinContent(ipt+1));
  // total_3->FixParameter(2,pi_sigma_tof->GetBinContent(ipt+1));
  // total_3->SetParLimits(3,100,1e6);


  // total_3->FixParameter(4,proton_mean_tof->GetBinContent(ipt+1));

  // if(ipt<17)
  //   total_3->SetParLimits(5,0.8,1);
  // else
  //   total_3->FixParameter(5,proton_sigma_tof->GetBinContent(ipt+1));

  // total_3->SetParLimits(6,100,1e6);
  // total_3->FixParameter(7,kaon_mean_tof->GetBinContent(ipt+1));
  // total_3->FixParameter(8,kaon_sigma_tof->GetBinContent(ipt+1));





  if(ipt<21)
    {
      total_3->FixParameter(1,pi_mean_tof->GetBinContent(ipt+1));
      
      total_3->SetParLimits(2,0.7,0.95);

      if(ipt==17)
	total_3->SetParLimits(3,50,4e3);
      else
	total_3->SetParLimits(3,50,2e8);
      
      total_3->SetParLimits(4,proton_mean_tof->GetBinContent(ipt+1)-0.5,proton_mean_tof->GetBinContent(ipt+1)+0.5);
      
     total_3->SetParLimits(5,0.7,1);
     total_3->SetParLimits(7,kaon_mean_tof->GetBinContent(ipt+1)-0.5,kaon_mean_tof->GetBinContent(ipt+1)+0.5);
     total_3->SetParLimits(8,0.7,1);
    }
  else {

    total_3->FixParameter(1,pi_mean_tof->GetBinContent(ipt+1));
    total_3->FixParameter(2,pi_sigma_tof->GetBinContent(ipt+1));


    total_3->SetParLimits(3,100,1e7);
    total_3->FixParameter(4,proton_mean_tof->GetBinContent(ipt+1));
  //   if(ip==)
    total_3->SetParLimits(5,0.8,1);
  // else
    //   total_3->FixParameter(5,proton_sigma_tof->GetBinContent(ipt+1));
  
    total_3->SetParLimits(6,100,1e6);

    total_3->FixParameter(7,kaon_mean_tof->GetBinContent(ipt+1));
    total_3->FixParameter(8,kaon_sigma_tof->GetBinContent(ipt+1));
  }
  


  //electron mean
  total_3->SetParameter(10,electron_Mean->GetBinContent(ipt+1));
  total_3->SetParLimits(10,electron_Mean->GetBinContent(ipt+1)-Sigma_flag*electron_Mean->GetBinError(ipt+1),electron_Mean->GetBinContent(ipt+1)+Sigma_flag*electron_Mean->GetBinError(ipt+1));
  //electron sigma
  total_3->SetParameter(11,electron_Sigma->GetBinContent(ipt+1));
  total_3->SetParLimits(11,electron_Sigma->GetBinContent(ipt+1)-Sigma_flag*electron_Sigma->GetBinError(ipt+1),electron_Sigma->GetBinContent(ipt+1)+Sigma_flag*electron_Sigma->GetBinError(ipt+1));

  //    total_3->FixParameter(10,electron_Mean->GetBinContent(ipt+1));
  //    //    total_3->SetParLimits(10,electron_Mean->GetBinContent(ipt+1)-Sigma_flag*electron_Mean->GetBinError(ipt+1),electron_Mean->GetBinContent(ipt+1)+Sigma_flag*electron_Mean->GetBinError(ipt+1));
  // //electron sigma
  // total_3->FixParameter(11,electron_Sigma->GetBinContent(ipt+1));
  // total_3->SetParLimits(11,electron_Sigma->GetBinContent(ipt+1)-Sigma_flag*electron_Sigma->GetBinError(ipt+1),electron_Sigma->GetBinContent(ipt+1)+Sigma_flag*electron_Sigma->GetBinError(ipt+1));

	
  
  if(Draw_flag==0)
    { 
      Inclusive->Fit(total_3,"R0","",-10,4);
    }  
  else
    Inclusive->Fit(total_3,"R","same",-10,4);
  
  g1->SetParameter(0,total_3->GetParameter(0));
  g1->SetParameter(1,total_3->GetParameter(1));
  g1->SetParameter(2,total_3->GetParameter(2));
      
  g2->SetParameter(0,total_3->GetParameter(3));
  g2->SetParameter(1,total_3->GetParameter(4));
  g2->SetParameter(2,total_3->GetParameter(5));
      
  g3->SetParameter(0,total_3->GetParameter(6));
  g3->SetParameter(1,total_3->GetParameter(7));
  g3->SetParameter(2,total_3->GetParameter(8));

  g4->SetParameter(0,total_3->GetParameter(9));
  g4->SetParameter(1,total_3->GetParameter(10));
  g4->SetParameter(2,total_3->GetParameter(11));
    
  if(Draw_flag!=0){
    TLegend *legend = new TLegend(0.15,0.6,0.35,0.8);
    legend->AddEntry(g1,"#pi ","lp");
    legend->AddEntry(g2,"proton ","lp");
    legend->AddEntry(g3,"kaon","lp");
    legend->AddEntry(g4,"e","lp");
	
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.055);
    legend ->SetTextFont(62);       
    legend->Draw("same");
    
    g1->Draw("same");
    g2->Draw("same");
    g3->Draw("same");
    g4->Draw("same");
    
    sprintf(buf,"purity_fit_%ipt_%i_sigma_constrain.pdf",ipt,Sigma_flag);
    c4->SaveAs(buf);
        
  }
  return g4->Integral(-1,3)/total_3->Integral(-1,3);

}
Double_t Get_sts_Error(TH1F * const inclusive,Int_t ipt,TF1 * total_3, TF1 *g1,TF1 *g2,TF1 *g3,TF1 *g4,Int_t Sigma_flag,Int_t draw_flag,Double_t purity_mean,Double_t purity_err)
{
  return 0.;
char buf[1024];
  TH1F *Inclusive= (TH1F *)  inclusive->Clone("Inclusive");
  TRandom3 *gRnd= new TRandom3(0);
  sprintf(buf,"purity_fit_%ipt",ipt);
  
  Double_t bin_low=purity_mean-50*purity_err;
  Double_t bin_high=purity_mean+50*purity_err;


  if(ipt==1)
    {
      bin_low=purity_mean-100*purity_err;
      bin_high=purity_mean+100*purity_err;
    }
    if(ipt==5)
    {
      bin_low=purity_mean-100000*purity_err;
      bin_high=purity_mean+100000*purity_err;
    }


  TH1F *mh1purity_sts=new TH1F(buf,"",800,bin_low,bin_high);
  TF1 *f_Fit = new TF1(TString("f_Fit"),"gaus",-1,1);
  f_Fit->SetParameter(1,purity_mean);
  TH1F *nsigmae_inclusive_shift=(TH1F *) Inclusive->Clone("nsigmae_inclusive_shift");

  for(Int_t i=1;i<=1000;i++)
    {
      for(Int_t bin=1; bin<=Inclusive->GetNbinsX();bin++)
        {
          if(Inclusive->GetBinContent(bin) && Inclusive->GetBinError(bin))
            
            nsigmae_inclusive_shift->SetBinContent(bin,gRnd->Gaus(Inclusive->GetBinContent(bin),Inclusive->GetBinError(bin)));
          
        } 
      Double_t purity =Fit_purity(nsigmae_inclusive_shift,ipt,electron_Mean->GetBinError(ipt+1), electron_Sigma->GetBinError(ipt+1),total_3,g1,g2,g3,g4,Sigma_flag,draw_flag);
      mh1purity_sts->Fill(purity);
      cout<< " STS!!!!!"<<purity<<"  ipt="<<ipt<<endl;
   
    }

  TCanvas *c5=new TCanvas("c5","",600,800);
  f_Fit->SetParameter(1, mh1purity_sts->GetMean());
  f_Fit->SetParLimits(1, mh1purity_sts->GetMean()- 50*mh1purity_sts->GetMeanError(), mh1purity_sts->GetMean()+50*mh1purity_sts->GetMeanError());

  mh1purity_sts->Fit(f_Fit,"R","");//0.5*bin_low,0.5*bin_high);
  
  mh1purity_sts->GetXaxis()->SetTitle("purity");
  mh1purity_sts->GetYaxis()->SetTitle("Counts");
  c5->SaveAs(buf+TString("purity_sts.pdf"));
  
  return f_Fit->GetParameter(2);
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

  TLegend *legend = new TLegend(0.15,0.35,0.35,0.6);
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

//Double_t Get_purity_sys(TH1F * const in)
void Draw_purity(TH1D * purity)
{
  gStyle->SetOptStat(00000);
  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd();
  TH2F *h2=new TH2F("h2","",10,0.2,4,10,0,1.1);
  h2->Draw();
  h2->GetXaxis()->SetTitle("p_{T}");
  h2->GetYaxis()->SetTitle("purity");
  purity->SetMarkerStyle(20);
  purity->SetMarkerColor(2);
  purity->Draw("samePE1");
  c2->SaveAs("purity_MB.pdf");
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
