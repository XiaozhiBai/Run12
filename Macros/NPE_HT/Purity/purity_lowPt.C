/*--------------------Xiaozhi------------------------------------ 

This macro is for the photo electron nsigma in different Pt bin
for the nsigma Electran calibration 
low Pt electron purity

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
#include "/Users/xiaozhi/NPE_work/NPE_low_pt/Binning/Bin.h"

using namespace std;

//const int Nbins_MB_ab=26;

//Double_t pt_bin_MB_46[]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,2,2.5,3,3.5,4};

//Double_t Pt_xx_low_46[]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,2,2.5,3,3.5,4};

//Double_t Pt_xx_high_46[]={0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,2,2.5,3,3.5,4};

//Int_t  Ptbin_low_46[]={4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,40,50,60,70};
// Int_t  Ptbin_high_46[]={5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,40,50,60,70,80};

//Int_t Ptbin_low_46[]={2,4,6,12,20,26};
//Int_t Ptbin_high_46[]={4,6,12,20,26,40};

void Get_electron_nsigma_Mean( TH3F* ,TH3F* ,TH2F *);
void Draw_mean_sigma(TH1F *a[],Float_t *,Float_t *,Float_t*,Float_t *);
void Get_purity(TH1F * const a[],Float_t *,Float_t *,Float_t*,Float_t *,Float_t ,Float_t ,Float_t ,Float_t );
void Get_sys_Error( TH1F   * const a[], Float_t,Float_t ,Float_t ,Float_t ,Float_t b[],Float_t c[] ,Float_t d[],TH1D *&,Int_t);
void Get_sts_Error(TH1F * const a[]  , Float_t,Float_t ,Float_t ,Float_t , TH1D *& ,TH1D *, Int_t);

Float_t Fit_purity(TH1F * const, Float_t,Float_t ,Float_t ,Float_t ,Int_t ,Int_t);

void Draw_purity(TH1D * );

void nSigma_cut_efficiency( TH1F *,TH1F *,TH1F *,TH1F *,TH1F *,TH1F *,Float_t,Float_t,Float_t,Float_t);
TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);


int purity_lowPt(){

  gStyle->SetOptStat(00000);
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1,"Y");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.1); 

  //  gStyle->SetEndErrorSize(4);



  
   TFile * inFile= new TFile("/Users/xiaozhi/NPE_work/NPE_low_pt/Root_File/Root_File_5_31/hist_5_31.root","READ");

  TH3F* nsigmaEUnlike=(TH3F *) inFile->Get("mh3nSigmaEUnlike_VPD");
  TH3F* nsigmaElike=(TH3F *) inFile->Get("mh3nSigmaElike_VPD");
  TH2F * nsigmaE_inclusive=(TH2F *) inFile->Get("mh2nSigmaElec_VPD"); 

   nsigmaE_inclusive->Draw();


  //return 0;



  //nsigmaEUnlike->Draw();

     Get_electron_nsigma_Mean(nsigmaEUnlike ,nsigmaElike,nsigmaE_inclusive);
return 0;


}
void Get_electron_nsigma_Mean( TH3F* mh3nsigmaEUnlike ,TH3F* mh3nsigmaElike,TH2F *mh2nsigmaE_inclusive)
{

  //  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);

  char buf[1024];
  TH1F * nsigmaEUnlike[Nbins_MB_46];
  TH1F * nsigmaElike[Nbins_MB_46];
  TH1F * nsigmaEUnlike_like[Nbins_MB_46];
  TH1F* nsigmaE_inclusive[Nbins_MB_46];

  mh3nsigmaEUnlike->GetXaxis()->SetRangeUser(-3.5,3.5); // unlike sign
  mh3nsigmaElike->GetXaxis()->SetRangeUser(-3.5,3.5); // like sign
   mh2nsigmaE_inclusive->GetXaxis()->SetRangeUser(-15,20);  // total

   for(Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
     {
      sprintf(buf,"nsigmaE_Pt%2.2f",Ptbin_low_46[ipt]/10.);
      cout<< buf<<endl;     
      nsigmaEUnlike[ipt]=(TH1F *) mh3nsigmaEUnlike->ProjectionX(buf+TString("Unlike"),Ptbin_low_46[ipt],Ptbin_high_46[ipt],0,1);
      nsigmaElike[ipt]=(TH1F *) mh3nsigmaElike->ProjectionX(buf+TString("like"),Ptbin_low_46[ipt],Ptbin_high_46[ipt],0,1);
      
      
      nsigmaEUnlike_like[ipt]=(TH1F *) nsigmaEUnlike[ipt]->Clone(buf+TString("Unlike-like"));
      nsigmaEUnlike_like[ipt]->Add(nsigmaElike[ipt],-1);

      sprintf(buf,"nsigmaE_inclusive%2.2f",Ptbin_low_46[ipt]/10.);
      nsigmaE_inclusive[ipt]=(TH1F *) mh2nsigmaE_inclusive->ProjectionX(buf,Ptbin_low_46[ipt],Ptbin_high_46[ipt]);    
      sprintf(buf,"Pt from %2.2fGeV - %2.2fGeV",Pt_xx_low_46[ipt],Pt_xx_high_46[ipt]);
      nsigmaE_inclusive[ipt]->SetTitle(buf);
      
    }
  
  //  nsigmaE_inclusive[20]->Draw();
  
  // return;
  
  
  // nsigmaEUnlike[2]->Draw();
  // nsigmaEUnlike_like[2]->SetLineColor(2);
  
  // nsigmaEUnlike_like[2]->Draw("same");
  // nsigmaElike[3]->Draw();
  // return;
  
  TCanvas *c1=new TCanvas("c1","",800,600);
  TCanvas *c2=new TCanvas("c2","",800,600);
  TCanvas *c3=new TCanvas("c3","",800,600);
  TCanvas *c4=new TCanvas("c4","",800,600);
  TCanvas *c5=new TCanvas("c5","",800,600);

  c1->Divide(3,3);
  c2->Divide(3,3);
  c3->Divide(3,3);
  c4->Divide(3,3);
  c5->Divide(3,3);
  
  
  //  c1->Draw();
  // return;
  int Npad=1;
  TF1 *f1 = new TF1(TString("f1"),"gaus",-5,5);
  for(Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
    {
      nsigmaEUnlike[ipt]->SetLineColor(1);
      nsigmaElike[ipt]->SetLineColor(6);
      nsigmaEUnlike_like[ipt]->SetLineColor(3);
      nsigmaEUnlike_like[ipt]->SetMarkerStyle(20);
      nsigmaEUnlike_like[ipt]->SetMarkerColor(3);
      
      sprintf(buf,"nsigma Electron P_{T} %2.2fGeV - %2.2fGev",Pt_xx_low_46[ipt],Pt_xx_high_46[ipt]);
      
      nsigmaEUnlike[ipt]->SetTitle(buf); 
      nsigmaEUnlike[ipt]->GetXaxis()->SetTitle("NSigmaE");
      nsigmaEUnlike[ipt]->GetYaxis()->SetTitle("Counts");
      
      if(ipt<9)
	{
	  c1->cd(Npad++);
	  gPad->SetLogy(0);
	}
      else if(ipt<18)
	{
	  c2->cd(Npad++);
	  gPad->SetLogy(0);
	}
      else  if(ipt<27)
	{
	  c3->cd(Npad++);
	  gPad->SetLogy(0);
	}

      else  if(ipt<36)
	{
	  c4->cd(Npad++);
	  gPad->SetLogy(0);
	}
      else  if(ipt<45)
	{
	  c5->cd(Npad++);
	  gPad->SetLogy(0);
	}

      
      if(Npad==10)
	Npad=1;    
      // gPad->SetLogy();

        nsigmaEUnlike_like[ipt]->Fit(f1,"RLL","",-3,3);
      nsigmaEUnlike[ipt]->Draw("sameE1");
      nsigmaEUnlike_like[ipt]->Draw("sameE1");
      //nsigmaEUnlike[ipt]->Draw("sameE1");
      nsigmaElike[ipt]->Draw("sameE1");

      TLegend *legend  = new TLegend(0.15,0.65,0.35,0.85);
      legend ->AddEntry(nsigmaEUnlike[ipt]," Unlike sign","lpe");
      legend ->AddEntry(nsigmaElike[ipt]," Like sign","lpe");
      legend ->AddEntry(nsigmaEUnlike_like[ipt]," Unlike-like","lpe");
      legend ->SetBorderSize(0);
      legend ->SetTextSize(0.045);
      legend ->SetFillStyle(0);
      legend ->Draw("same");


    }

     c1->SaveAs("nsigmaE_c1.pdf");
     c2->SaveAs("nsigmaE_c2.pdf");
     c3->SaveAs("nsigmaE_c3.pdf");
     c4->SaveAs("nsigmaE_c4.pdf");
     c5->SaveAs("nsigmaE_c5.pdf");

     // return ;

     TF1 *f = new TF1(TString("f"),"gaus",-5,5);

     Float_t mean[Nbins_MB_46];
     Float_t sigma[Nbins_MB_46];
     Float_t meanError[Nbins_MB_46];
     Float_t sigmaErr[Nbins_MB_46];

 Float_t purity_sts_Error=0;
     
     for( Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
       {
	 nsigmaEUnlike_like[ipt]->Fit(f,"R0q","",-2,3);
	 mean[ipt]=f->GetParameter(1);
	 sigma[ipt]=f->GetParameter(2);
	 meanError[ipt]=f->GetParError(1);
	 sigmaErr[ipt]=f->GetParError(2);	
	 
		 // cout<< "mean"<<"  "<< mean[ipt]<< " sigma "<< sigma[ipt]<<" meanErroor "<< meanError[ipt]<< " sigmaError "<< sigmaErr[ipt] << endl;
       }
     
     
     Draw_mean_sigma(nsigmaE_inclusive,mean,meanError,sigma,sigmaErr);
    
     //    nSigma_cut_efficiency(nsigmaEUnlike_like);
}




void Get_purity(TH1F * const nsigmaE_inclusive[],Float_t Mean[],Float_t MeanErr[],Float_t Sigma[],Float_t SigmaErr[],Float_t mean_central,Float_t sigma_central ,Float_t Mean_err ,Float_t Sigma_err )
{

 
  
  // gStyle->SetOptFit(0); 
 Double_t Nsigma=5;  
  TCanvas *c1=new TCanvas("c1","",800,600);
  TCanvas *c2=new TCanvas("c2","",800,600);
  TCanvas *c3=new TCanvas("c3","",800,600);
  TCanvas *c4=new TCanvas("c4","",800,600);
  TCanvas *c5=new TCanvas("c5","",800,600);
  TCanvas *c6=new TCanvas("c6","",800,600);
 
  c1->Divide(3,3);
  c2->Divide(3,3);
  c3->Divide(3,3);
  c4->Divide(3,3);
  c5->Divide(3,3);
  c6->Divide(3,3);

 
  // nsigmaE_inclusive[2]->Draw();
  //return; 
  TH1D * mh1purity_sys=new TH1D("mh1purity_sys","Inclusive electron purity",Nbins_MB_46,pt_bin_MB_46);
  TH1D * mh1purity_sts=new TH1D("mh1purity_sts","Inclusive electron purity",Nbins_MB_46,pt_bin_MB_46);
  TH1D * mh1purity_all=new TH1D("mh1purity_all","Inclusive electron purity",Nbins_MB_46,pt_bin_MB_46);
  
  
  TH1F * mh1purity_one_sigma=new TH1F("mh1purity_one_sigma","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1purity_two_sigma=new TH1F("mh1purity_two_sigma","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1purity_three_sigma=new TH1F("mh1purity_three_sigma","",Nbins_MB_46,pt_bin_MB_46);
  

 

  //
  TH1F * mh1K_P_Constant=new TH1F("mh1K_P_Constant","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1K_P_Mean=new TH1F("mh1K_P_Mean","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1K_P_Sigma=new TH1F("mh1K_P_Sigma","",Nbins_MB_46,pt_bin_MB_46);

  TH1F * mh1K_K_Constant=new TH1F("mh1K_K_Constant","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1K_K_Mean=new TH1F("mh1K_K_Mean","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1K_K_Sigma=new TH1F("mh1K_K_Sigma","",Nbins_MB_46,pt_bin_MB_46);

  TH1F * mh1Pion_Constant=new TH1F("mh1Pion_Constant","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1Pion_Mean=new TH1F("mh1Pion_Mean","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1Pion_Sigma=new TH1F("mh1Pion_Sigma","",Nbins_MB_46,pt_bin_MB_46);
  
  
  TH1F * mh1Electron_Constant=new TH1F("mh1Electron_Constant","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1Electron_Mean=new TH1F("mh1Electron_Mean","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1Electron_Sigma=new TH1F("mh1Electron_Sigma","",Nbins_MB_46,pt_bin_MB_46);
 
  //


  
 Float_t purity=0;
 Float_t purity_sys_Error=0;
 Float_t purity_one_sigma_constrain[Nbins_MB_46]={0};
 Float_t purity_two_sigma_constrain[Nbins_MB_46]={0};
 Float_t purity_three_sigma_constrain[Nbins_MB_46]={0};
  
  
 TFile *file_2=new TFile("TOF_ONlY_DATA.root","read");
 
 TH1F *  mh1_Pion_mean = (TH1F *) file_2->Get("mh1_Pion_mean");
 TH1F *  mh1_Kaon_mean = (TH1F *) file_2->Get("mh1_Kaon_mean");
 TH1F *  mh1_Proton_mean = (TH1F *) file_2->Get("mh1_Proton_mean");
 TH1F *  mh1_Pion_sigma = (TH1F *) file_2->Get("mh1_Pion_sigma");
 TH1F *  mh1_Kaon_sigma = (TH1F *) file_2->Get("mh1_Kaon_sigma");
 TH1F *  mh1_Proton_sigma = (TH1F *) file_2->Get("mh1_Proton_sigma");

 TFile *file_4=new TFile("Constant_Constrain.root","READ");
   
   TH1F *mh1_Electron_Constant =(TH1F *) file_4->Get("mh1_Electron_Constant");
   TH1F *mh1_Pion_Constant =(TH1F *) file_4->Get("mh1_Pion_Constant");
   TH1F *mh1_Kaon_Constant =(TH1F *) file_4->Get("mh1_Kaon_Constant");
   TH1F *mh1_Proton_Constant =(TH1F *) file_4->Get("mh1_Proton_Constant");


   //     file_2->Close();
   
   Int_t  Npad=1;
   
   //dev/   for( Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
     for( Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
     {
       TF1 *g1=new TF1(nsigmaE_inclusive[ipt]->GetName()+TString("g1"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15);   //fit the ele
       TF1 *g2=new TF1(nsigmaE_inclusive[ipt]->GetName()+TString("g2"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15);  //fit the electron 
       TF1 *g3=new TF1(nsigmaE_inclusive[0]->GetName()+TString("g3"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15); 
       TF1 *g4=new TF1(nsigmaE_inclusive[0]->GetName()+TString("g4"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15); 
      
       TF1 *total_3 =new TF1(nsigmaE_inclusive[ipt]->GetName()+TString("total_3"),"[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)+[9]*TMath::Gaus(x,[10],[11],1)",-15,15);// fit the total 
       
       TF1 *total_2 =new TF1(nsigmaE_inclusive[ipt]->GetName()+TString("total_2"),"[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1)",-15,15);// fit the total 
       
       // total_3->SetParNames( "K C", "K #mu", "K #sigma", "#pi C", "#pi #mu", "#pi #sigma","e C","e #mu","e #sigma","P C","P #mu","P #sigma");
       total_2->SetParNames( "#pi C", "#pi #mu", "#pi #sigma", "e C", "e #mu", "e #sigma");
       
     g1->SetLineColor(kRed);
     g2->SetLineColor(kGreen+3);
     g3->SetLineColor(kMagenta);
     g4->SetLineColor(7);
 
    total_2->SetLineColor(kBlue);
     total_3->SetLineColor(kBlue);
     
     cout<< "Fit "<<endl;
     if( ipt<9)
       {
	 c1->cd(Npad++);	
	 gPad->SetLogy();
       }
     else if( ipt<18)
       {
	 c2->cd(Npad++);	
	 gPad->SetLogy();
       }
     else if( ipt<27)
       {
	 c3->cd(Npad++);	
	 gPad->SetLogy();
       }
     else if( ipt<36)
       {
	 c4->cd(Npad++);	
	 gPad->SetLogy();
       }

     else if( ipt<45)
       {
	 c5->cd(Npad++);	
	 gPad->SetLogy();
       }

     else if( ipt<46)
       {
	 c6->cd(Npad++);	
	 gPad->SetLogy();
       }
   
  
     if(Npad==10)
       Npad=1;

    
     /*  
	 Float_t purity=0;
	 Float_t purity_sys_Error=0;
	 Float_t purity_one_sigma_constrain[Nbins_MB_46];
	 Float_t purity_two_sigma_constrain[Nbins_MB_46];
	 Float_t purity_three_sigma_constrain[Nbins_MB_46];
	 
	 Float_t purity_sts_Error=0;
     */
     
     
     //     file_2->Close();
     

     if(ipt<0)
       {
	 
	 //pion constrain
	 total_2->SetParameter(1,mh1_Pion_mean->GetBinContent(ipt+1));
	 total_2->SetParLimits(1,mh1_Pion_mean->GetBinContent(ipt+1)-0.2,mh1_Pion_mean->GetBinContent(ipt+1)+0.2);
	 
	 total_2->SetParameter(2,mh1_Pion_sigma->GetBinContent(ipt+1));
	  total_2->SetParLimits(2,mh1_Pion_sigma->GetBinContent(ipt+1)-0.2,mh1_Pion_sigma->GetBinContent(ipt+1)+0.2);
	 
	 // total_2->SetParameter(2,1);
	 // total_2->SetParLimits(2,0.7,1.3);
	

 
	 
	 //electron constrain 	 
	 total_2->SetParameter(4,mean_central);
	 total_2->SetParameter(5,sigma_central);
	 total_2->SetParLimits(4, mean_central-Mean_err,mean_central+Mean_err);
	 total_2->SetParLimits(5, sigma_central-Sigma_err,sigma_central+Sigma_err);
	 



	 
	 nsigmaE_inclusive[ipt]->Fit(total_2,"R","same",-9,4);
	 
	 g1->SetParameter(0,total_2->GetParameter(3));
	 g1->SetParameter(1,total_2->GetParameter(4));
	 g1->SetParameter(2,total_2->GetParameter(5));
	 
	 g2->SetParameter(0,total_2->GetParameter(0));
	 g2->SetParameter(1,total_2->GetParameter(1));
	 g2->SetParameter(2,total_2->GetParameter(2));
	 
	 g1->Draw("same");
	 g2->Draw("same");
	 g3->Draw("same");
	 
	 
	 
	 
	 mh1Pion_Mean->SetBinContent(ipt+1,g2->GetParameter(1));
	 mh1Pion_Mean->SetBinError(ipt+1,g2->GetParError(1));
	 
	 
         
	 mh1Pion_Sigma->SetBinContent(ipt+1,g2->GetParameter(2));
	 mh1Pion_Sigma->SetBinError(ipt+1,g2->GetParError(2));




	 mh1Electron_Mean->SetBinContent(ipt+1,g1->GetParameter(1));
	 mh1Electron_Mean->SetBinError(ipt+1,g1->GetParError(1));

	 mh1Electron_Sigma->SetBinContent(ipt+1,g1->GetParameter(2));
	 mh1Electron_Sigma->SetBinError(ipt+1,g1->GetParError(2));


	 TLegend *legend  = new TLegend(0.15,0.65,0.35,0.85);
	 legend ->AddEntry(g1," e","lpe");
	 legend ->AddEntry(g2," #pi","lpe");
	 
	 if(16<=ipt)
	   legend ->AddEntry(g3," k+p","lpe");
	 
	 legend ->SetBorderSize(0);
	 legend ->SetTextSize(0.055);
	 legend ->SetFillStyle(0);
	 legend ->Draw("same");
	 purity=g1->Integral(-1,3)/total_2->Integral(-1,3);	  
	 
	 //	  purity_sys_Error=Get_sys_Error(nsigmaE_inclusive[ipt],mean_central,Mean_err,sigma_central,Sigma_err,purity_one_sigma_constrain,purity_two_sigma_constrain,purity_three_sigma_constrain,2);
	 
	 //	   purity_sts_Error=Get_sts_Error(nsigmaE_inclusive[ipt],mean_central,Mean_err,sigma_central,Sigma_err,ipt,2);	
	 
       }
          
     else
       {
	 //continue;
	
	 if(15<ipt && ipt<36)
	   {
	      total_3->SetParameter(3,mh1_Pion_Constant->GetBinContent(ipt+1));
	     // total_3->SetParLimits(3,mh1_Pion_Constant->GetBinContent(ipt+1)-1000,mh1_Pion_Constant->GetBinContent(ipt+1)+1000);
	   }
	 total_3->SetParameter(4,mh1_Pion_mean->GetBinContent(ipt+1));
	 total_3->SetParLimits(4,mh1_Pion_mean->GetBinContent(ipt+1)-Nsigma*mh1_Pion_mean->GetBinError(ipt+1)-0.3,mh1_Pion_mean->GetBinContent(ipt+1)+Nsigma*mh1_Pion_mean->GetBinError(ipt+1)+0.3);
	 
	

	 total_3->SetParameter(5,mh1_Pion_sigma->GetBinContent(ipt+1));
	 total_3->SetParLimits(5,mh1_Pion_sigma->GetBinContent(ipt+1)-Nsigma*mh1_Pion_sigma->GetBinError(ipt+1)-0.1,mh1_Pion_sigma->GetBinContent(ipt+1)+Nsigma*mh1_Pion_sigma->GetBinError(ipt+1)+0.1);
	 
     
     
	 // kaon --------------------------------------------
	 // 
	  if(4<ipt && ipt<8){
	      total_3->SetParameter(0,mh1_Kaon_Constant->GetBinContent(ipt+1));	  
	    total_3->SetParameter(0,500);	  
	    // total_3->SetParLimits(0,mh1_Kaon_Constant->GetBinContent(ipt+1)-3000,mh1_Kaon_Constant->GetBinContent(ipt+1)+3000);
	     total_3->SetParLimits(0,100,500);
	  }
	  
	  else if(10<ipt && ipt<30)
	  {
	    total_3->SetParameter(0,0.1*mh1_Pion_Constant->GetBinContent(ipt+1));	  	    
	    // total_3->SetParameter(0,1000);	  	    
	    //total_3->SetParLimits(0,mh1_Kaon_Constant->GetBinContent(35)-100,mh1_Kaon_Constant->GetBinContent(35)+100);
	  }
	
	  total_3->SetParameter(1,mh1_Kaon_mean->GetBinContent(ipt+1));
	  total_3->SetParLimits(1,mh1_Kaon_mean->GetBinContent(ipt+1)-Nsigma*mh1_Kaon_mean->GetBinError(ipt+1)-0.2,mh1_Kaon_mean->GetBinContent(ipt+1)+Nsigma*mh1_Kaon_mean->GetBinError(ipt+1)+0.2);
	  total_3->SetParameter(2,mh1_Kaon_sigma->GetBinContent(ipt+1));
	  total_3->SetParLimits(2,mh1_Kaon_sigma->GetBinContent(ipt+1)-Nsigma*mh1_Kaon_sigma->GetBinError(ipt+1)-0.1,mh1_Kaon_sigma->GetBinContent(ipt)+Nsigma*mh1_Kaon_sigma->GetBinError(ipt+1)+0.1);
	  // Proton
	  
	    if(10<ipt && ipt<23){
	   //   Double_t costant_proton=(mh1_Proton_Constant->GetBinContent(ipt+4)+mh1_Proton_Constant->GetBinContent(ipt-4))/2.0;
	      // total_3->SetParameter(9,mh1_Proton_Constant->GetBinContent(ipt+1));
	   total_3->SetParameter(9,1000);
	    
	   //  total_3->SetParLimits(9,mh1_Proton_Constant->GetBinContent(ipt+1)-1000,mh1_Proton_Constant->GetBinContent(ipt+1)+1000);
	   //   total_3->SetParLimits(9,costant_proton-400,costant_proton+400);
	    }
     
    
	  total_3->SetParameter(10,mh1_Proton_mean->GetBinContent(ipt+1));
	  total_3->SetParLimits(10,mh1_Proton_mean->GetBinContent(ipt+1)-Nsigma*mh1_Proton_mean->GetBinError(ipt+1)-0.2,mh1_Proton_mean->GetBinContent(ipt+1)+Nsigma*mh1_Proton_mean->GetBinError(ipt+1)+0.2);
	  total_3->SetParameter(11,mh1_Proton_sigma->GetBinContent(ipt+1));
	  total_3->SetParLimits(11,mh1_Proton_sigma->GetBinContent(ipt+1)-Nsigma*mh1_Proton_sigma->GetBinError(ipt+1)-0.1,mh1_Proton_sigma->GetBinContent(ipt)+Nsigma*mh1_Proton_sigma->GetBinError(ipt+1)+0.1);
     
     // electron ---------------------------------------
     if(5<ipt &&ipt<25){
       
       total_3->SetParameter(6,mh1_Electron_Constant->GetBinContent(ipt+1));
       
       total_3->SetParLimits(6,mh1_Electron_Constant->GetBinContent(ipt+1)-5,mh1_Electron_Constant->GetBinContent(ipt+1)+5);
     }
     total_3->SetParameter(7,mean_central);
     total_3->SetParameter(8,sigma_central);
     total_3->SetParLimits(7, mean_central-Mean_err,mean_central+Mean_err);
     total_3->SetParLimits(8, sigma_central-Sigma_err,sigma_central+Sigma_err);
     
     
     
     cout<< " aa"<< ipt<< endl;      
     
     nsigmaE_inclusive[ipt]->Fit(total_3,"R","same",-10,13);
     
     cout<< " bb"<< ipt<< endl;      
     
     g1->SetParameter(0,total_3->GetParameter(6));
     g1->SetParameter(1,total_3->GetParameter(7));
     g1->SetParameter(2,total_3->GetParameter(8));
     
     g2->SetParameter(0,total_3->GetParameter(3));
     g2->SetParameter(1,total_3->GetParameter(4));
     g2->SetParameter(2,total_3->GetParameter(5));
     
     g3->SetParameter(0,total_3->GetParameter(0));
     g3->SetParameter(1,total_3->GetParameter(1));
     g3->SetParameter(2,total_3->GetParameter(2));
     
     g4->SetParameter(0,total_3->GetParameter(9));
     g4->SetParameter(1,total_3->GetParameter(10));
     g4->SetParameter(2,total_3->GetParameter(11));
     
     // Set err 
     g1->SetParError(0,total_3->GetParError(6));
     g1->SetParError(1,total_3->GetParError(7));
     g1->SetParError(2,total_3->GetParError(8));
     
     g2->SetParError(0,total_3->GetParError(3));
     g2->SetParError(1,total_3->GetParError(4));
     g2->SetParError(2,total_3->GetParError(5));
     
     g3->SetParError(0,total_3->GetParError(0));
     g3->SetParError(1,total_3->GetParError(1));
     g3->SetParError(2,total_3->GetParError(2));
     
     g4->SetParError(0,total_3->GetParError(9));
     g4->SetParError(1,total_3->GetParError(10));
     g4->SetParError(2,total_3->GetParError(11));
     
     
     
      g1->Draw("same");
      g2->Draw("same");
      g3->Draw("same");
      g4->Draw("same");
     
	  
     mh1Pion_Constant->SetBinContent(ipt+1,g2->GetParameter(0));
     mh1Pion_Constant->SetBinError(ipt+1,g2->GetParError(0));
     
     mh1Pion_Mean->SetBinContent(ipt+1,g2->GetParameter(1));
     mh1Pion_Mean->SetBinError(ipt+1,g2->GetParError(1));
     
     mh1Pion_Sigma->SetBinContent(ipt+1,g2->GetParameter(2));
     mh1Pion_Sigma->SetBinError(ipt+1,g2->GetParError(2));
     
     mh1Electron_Constant->SetBinContent(ipt+1,g1->GetParameter(0));
     mh1Electron_Constant->SetBinError(ipt+1,g1->GetParError(0));
     
     mh1Electron_Mean->SetBinContent(ipt+1,g1->GetParameter(1));
     mh1Electron_Mean->SetBinError(ipt+1,g1->GetParError(1));
     
     mh1Electron_Sigma->SetBinContent(ipt+1,g1->GetParameter(2));
     mh1Electron_Sigma->SetBinError(ipt+1,g1->GetParError(2));
     
     mh1K_K_Constant->SetBinContent(ipt+1,g3->GetParameter(0));
     mh1K_K_Constant->SetBinError(ipt+1,g3->GetParError(0));
     
     mh1K_K_Mean->SetBinContent(ipt+1,g3->GetParameter(1));
     mh1K_K_Mean->SetBinError(ipt+1,g3->GetParError(1));
     
     mh1K_K_Sigma->SetBinContent(ipt+1,g3->GetParameter(2));
     mh1K_K_Sigma->SetBinError(ipt+1,g3->GetParError(2));
     
     mh1K_P_Constant->SetBinContent(ipt+1,g4->GetParameter(0));
     mh1K_P_Constant->SetBinError(ipt+1,g4->GetParError(0));
     
     mh1K_P_Mean->SetBinContent(ipt+1,g4->GetParameter(1));
     mh1K_P_Mean->SetBinError(ipt+1,g4->GetParError(1));
     
     mh1K_P_Sigma->SetBinContent(ipt+1,g4->GetParameter(2));
     mh1K_P_Sigma->SetBinError(ipt+1,g4->GetParError(2));
     
     
     
     
     
     TLegend *legend  = new TLegend(0.15,0.65,0.35,0.85);
     legend ->AddEntry(g1," e","lpe");
     legend ->AddEntry(g2," #pi","lpe");
     
     if(8<=ipt)
       legend ->AddEntry(g3," k","lpe");
     legend ->AddEntry(g4," p","lpe");
     
     legend ->SetBorderSize(0);
     legend ->SetTextSize(0.055);
     legend ->SetFillStyle(0);
     legend ->Draw("same");
     
     purity=g1->Integral(-1,3)/total_3->Integral(-1,3);	 
     
     // purity_sys_Error=Get_sys_Error(nsigmaE_inclusive[ipt],mean_central,Mean_err,sigma_central,Sigma_err,purity_one_sigma_constrain,purity_two_sigma_constrain,purity_three_sigma_constrain,mh1purity_sys,3);
     
     // purity_sts_Error=Get_sts_Error(nsigmaE_inclusive[ipt],mean_central,Mean_err,sigma_central,Sigma_err,ipt,3);	
      }
     
        if(Npad==10)
 	Npad=1;
 
       TLegend *legend  = new TLegend(0.15,0.65,0.35,0.85);
       legend ->AddEntry(g1," e","lpe");
       legend ->AddEntry(g2," #pi","lpe");

       if(9<=ipt)
	 legend ->AddEntry(g3," k+p","lpe");
       
      legend ->SetBorderSize(0);
      legend ->SetTextSize(0.055);
      legend ->SetFillStyle(0);
      // legend ->Draw("same");


      

      //  cout<<" purity <<<<<<<<<<<<<<<   "<<purity<<"  sts"<< purity_sys_Error<< endl;
     
      // mh1purity_all->SetBinContent(ipt+1,purity);
      // mh1purity->SetBinError(ipt+1,purity_sys_Error);
      // mh1purity_one_sigma->SetBinContent(ipt+1,purity_one_sigma_constrain);
      // mh1purity_one_sigma->SetBinError(ipt+1,0);
      // mh1purity_two_sigma->SetBinContent(ipt+1,purity_two_sigma_constrain);
      // mh1purity_two_sigma->SetBinError(ipt+1,0);
      // mh1purity_three_sigma->SetBinContent(ipt+1,purity_three_sigma_constrain);
      // mh1purity_three_sigma->SetBinError(ipt+1,0);

      

      //      fstream purity_out;
       fstream  purity_out("purity.dat",ios::trunc|ios::out);


      purity_out<< purity<< "  "<<purity_sys_Error<<endl;
   

     }


    
   c1->SaveAs("nsigmaE_inclusiveLog_c1Fit.pdf");
   c2->SaveAs("nsigmaE_inclusiveLog_c2Fit.pdf");
   c3->SaveAs("nsigmaE_inclusiveLog_c3Fit.pdf");
   c4->SaveAs("nsigmaE_inclusiveLog_c4Fit.pdf");
   c5->SaveAs("nsigmaE_inclusiveLog_c5Fit.pdf");
   c6->SaveAs("nsigmaE_inclusiveLog_c6Fit.pdf");
   

   // return;

  TFile *file_3=new TFile("Mean_Sigma_Hadron.root","RECREATE");



   mh1Pion_Mean->Write();
   mh1Pion_Sigma->Write();
   mh1Pion_Constant->Write();
   mh1Electron_Constant->Write();
   mh1Electron_Mean->Write();
   mh1Electron_Sigma->Write();

   mh1K_P_Constant->Write();
   mh1K_P_Mean->Write();
   mh1K_P_Sigma->Write();

   mh1K_K_Constant->Write();
   mh1K_K_Mean->Write();
   mh1K_K_Sigma->Write();
  
   // Get_sys_Error(nsigmaE_inclusive,mean_central,Mean_err,sigma_central,Sigma_err,purity_one_sigma_constrain,purity_two_sigma_constrain,purity_three_sigma_constrain,mh1purity_sys,3);

  // TH1D * mh1purity_sts=(TH1D *) mh1purity_sys->Clone("mh1purity_sts");
  //////////////////////888888888888888888888888
   // Get_sts_Error(nsigmaE_inclusive,mean_central,Mean_err,sigma_central,Sigma_err,mh1purity_sts,mh1purity_sys,3);
  ///////////////////////////////8888888888888888888888




    for(Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
   {
     Double_t purity_err=mh1purity_sys->GetBinError(ipt+1)*mh1purity_sys->GetBinError(ipt+1)+mh1purity_sts->GetBinError(ipt+1)*mh1purity_sts->GetBinError(ipt+1);
     
     mh1purity_all->SetBinContent(ipt+1,mh1purity_sys->GetBinContent(ipt+1));
     mh1purity_all->SetBinError(ipt+1,sqrt(purity_err));
     
     
   }
   
//   throw (-1);
   //  Draw_purity(mh1purity_all);
   Draw_purity(mh1purity_all);
 

   
}


void Draw_purity(TH1D * purity)
{
  
  TCanvas *c2=new TCanvas("c2","",800,600);
  
  c2->cd();

  TH2F *hh=new TH2F("hh","",10,0,4.5,10,0,1.2);

  hh->Draw();
  hh->GetXaxis()->SetTitle("p_{T}");
  hh->GetYaxis()->SetTitle("purity");

  purity->SetMarkerStyle(20);
  purity->SetMarkerColor(2);
 

  purity->Draw("samePE1");

  c2->SaveAs("purity.pdf");


  TH1D *purity_low_pt=new TH1D("purity_low_pt","",10,pt_bin_MB_46);

  for(Int_t ipt=0;ipt<10;ipt++)
    {
      purity_low_pt->SetBinContent(ipt+1,purity->GetBinContent(ipt+17));
      purity_low_pt->SetBinError(ipt+1,purity->GetBinError(ipt+17));

    }


  TFile *file=new TFile("Inclusive_purity.root","RECREATE");

  purity_low_pt->Write();

  file->Close();


}



Float_t Fit_purity( TH1F * const  hh, Float_t  mean_central,Float_t  Mean_err, Float_t sigma_central,Float_t  Sigma_err,Int_t  ipt ,Int_t Fit_flag)
{

  TCanvas *c1=new TCanvas("c1","",600,800);
  c1->cd();

  Float_t Purity=0;
 Double_t Nsigma=5;    

  TF1 *g1=new TF1(hh->GetName()+TString("g1"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15);   //fit the ele
  TF1 *g2=new TF1(hh->GetName()+TString("g2"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15);  //fit the electron 
  TF1 *g3=new TF1(hh->GetName()+TString("g3"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15); 
  TF1 *g4=new TF1(hh->GetName()+TString("g4"),"[0]*TMath::Gaus(x,[1],[2],1)",-15,15); 
      
  TF1 *total_3 =new TF1(hh->GetName()+TString("total_3"),"[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)+[9]*TMath::Gaus(x,[10],[11],1)",-15,15);// fit the total 



   
  TF1 *total_2 =new TF1(TString("total_2"),"[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1)",-15,15);// fit the total 
 
  // total->SetParNames( "e C", "e #mu", "e #sigma", "#pi C", "#pi #mu", "#pi #sigma","K C","K #mu","K #sigma");



     

  
  TFile *file_2=new TFile("TOF_ONlY_DATA.root","read");
 
  TH1F *  mh1_Pion_mean = (TH1F *) file_2->Get("mh1_Pion_mean");
  TH1F *  mh1_Kaon_mean = (TH1F *) file_2->Get("mh1_Kaon_mean");
  TH1F *  mh1_Proton_mean = (TH1F *) file_2->Get("mh1_Proton_mean");
  TH1F *  mh1_Pion_sigma = (TH1F *) file_2->Get("mh1_Pion_sigma");
  TH1F *  mh1_Kaon_sigma = (TH1F *) file_2->Get("mh1_Kaon_sigma");
  TH1F *  mh1_Proton_sigma = (TH1F *) file_2->Get("mh1_Proton_sigma");

  TFile *file_4=new TFile("Constant_Constrain.root","READ");
   
  TH1F *mh1_Electron_Constant =(TH1F *) file_4->Get("mh1_Electron_Constant");
  TH1F *mh1_Pion_Constant =(TH1F *) file_4->Get("mh1_Pion_Constant");
  TH1F *mh1_Kaon_Constant =(TH1F *) file_4->Get("mh1_Kaon_Constant");
  TH1F *mh1_Proton_Constant =(TH1F *) file_4->Get("mh1_Proton_Constant");

    

  if(Fit_flag==2)
    {	 
      //pion constrain
      total_2->SetParameter(1,mh1_Pion_mean->GetBinContent(ipt+1));
      total_2->SetParLimits(1,mh1_Pion_mean->GetBinContent(ipt+1)-0.3,mh1_Pion_mean->GetBinContent(ipt+1)+0.3);
	 
      //	 total_2->SetParameter(2,mh1_Pion_sigma->GetBinContent(ipt+1));
      //dev/ total_2->SetParLimits(2,mh1_Pion_sigma->GetBinContent(ipt+1)-0.2,mh1_Pion_sigma->GetBinContent(ipt+1)+0.2);

      total_2->SetParameter(2,1);
      total_2->SetParLimits(2,0.8,1.2);
	 


	 
      //electron constrain 	 
      total_2->SetParameter(4,mean_central);
      total_2->SetParameter(5,sigma_central);
      total_2->SetParLimits(4, mean_central-Mean_err,mean_central+Mean_err);
      total_2->SetParLimits(5, sigma_central-Sigma_err,sigma_central+Sigma_err);
	 
	 
      hh->Fit(total_2,"R","same",-10,4);
	 
	 
      g1->SetParameter(0,total_2->GetParameter(3));
      g1->SetParameter(1,total_2->GetParameter(4));
      g1->SetParameter(2,total_2->GetParameter(5));
	 
      g2->SetParameter(0,total_2->GetParameter(0));
      g2->SetParameter(1,total_2->GetParameter(1));
      g2->SetParameter(2,total_2->GetParameter(2));
	 
		 

      Purity=g1->Integral(-1,3)/total_2->Integral(-1,3);
      ////////////////////////////////////////////////////////


    }
  else   if(Fit_flag==3)
    {

       

      if(10<ipt && ipt<36)
	{
	  total_3->SetParameter(3,mh1_Pion_Constant->GetBinContent(ipt+1));
	  total_3->SetParLimits(3,mh1_Pion_Constant->GetBinContent(ipt+1)-500,mh1_Pion_Constant->GetBinContent(ipt+1)+500);
	}
       
      if(ipt<5)
	{
	  total_3->SetParameter(4,mh1_Pion_mean->GetBinContent(ipt+1));
	  total_3->SetParLimits(4,mh1_Pion_mean->GetBinContent(ipt+1)-0.3,mh1_Pion_mean->GetBinContent(ipt+1)+0.3);
	     
	}
      else{
	total_3->SetParameter(4,mh1_Pion_mean->GetBinContent(ipt+1));
	total_3->SetParLimits(4,mh1_Pion_mean->GetBinContent(ipt+1)-Nsigma*mh1_Pion_mean->GetBinError(ipt+1)+0.2,mh1_Pion_mean->GetBinContent(ipt+1)+Nsigma*mh1_Pion_mean->GetBinError(ipt+1)-0.2);
      }
	
      total_3->SetParameter(5,mh1_Pion_sigma->GetBinContent(ipt+1));
      total_3->SetParLimits(5,mh1_Pion_sigma->GetBinContent(ipt+1)-Nsigma*mh1_Pion_sigma->GetBinError(ipt+1)-0.1,mh1_Pion_sigma->GetBinContent(ipt+1)+Nsigma*mh1_Pion_sigma->GetBinError(ipt+1)+0.1);
	 
     
     
      // kaon --------------------------------------------
      total_3->SetParameter(0,mh1_Kaon_Constant->GetBinContent(ipt+1));
      if(5<ipt && ipt<36){
	total_3->SetParLimits(0,mh1_Kaon_Constant->GetBinContent(ipt+1)-1000,mh1_Kaon_Constant->GetBinContent(ipt+1)+1000);
      }
	  
      else if(36<=ipt)
	{
	  total_3->SetParLimits(0,mh1_Kaon_Constant->GetBinContent(35)-1000,mh1_Kaon_Constant->GetBinContent(35)+1000);
	}
	
      total_3->SetParameter(1,mh1_Kaon_mean->GetBinContent(ipt+1));
      total_3->SetParLimits(1,mh1_Kaon_mean->GetBinContent(ipt+1)-Nsigma*mh1_Kaon_mean->GetBinError(ipt+1)-0.2,mh1_Kaon_mean->GetBinContent(ipt+1)+Nsigma*mh1_Kaon_mean->GetBinError(ipt+1)+0.2);
      total_3->SetParameter(2,mh1_Kaon_sigma->GetBinContent(ipt+1));
      total_3->SetParLimits(2,mh1_Kaon_sigma->GetBinContent(ipt+1)-Nsigma*mh1_Kaon_sigma->GetBinError(ipt+1)-0.1,mh1_Kaon_sigma->GetBinContent(ipt)+Nsigma*mh1_Kaon_sigma->GetBinError(ipt+1)+0.1);
      // Proton
	  
      // if(18<ipt && ipt<23){
      //   Double_t costant_proton=(mh1_Proton_Constant->GetBinContent(ipt+4)+mh1_Proton_Constant->GetBinContent(ipt-4))/2.0;
      //   total_3->SetParameter(9,mh1_Proton_Constant->GetBinContent(ipt+1));
      //   total_3->SetParLimits(9,mh1_Proton_Constant->GetBinContent(ipt+1)-500,mh1_Proton_Constant->GetBinContent(ipt+1)+500);
      //   total_3->SetParLimits(9,costant_proton-400,costant_proton+400);
      // }
     
      total_3->SetParameter(10,mh1_Proton_mean->GetBinContent(ipt+1));
      total_3->SetParLimits(10,mh1_Proton_mean->GetBinContent(ipt+1)-Nsigma*mh1_Proton_mean->GetBinError(ipt+1)-0.2,mh1_Proton_mean->GetBinContent(ipt+1)+Nsigma*mh1_Proton_mean->GetBinError(ipt+1)-0.2);
      total_3->SetParameter(11,mh1_Proton_sigma->GetBinContent(ipt+1));
      total_3->SetParLimits(11,mh1_Proton_sigma->GetBinContent(ipt+1)-Nsigma*mh1_Proton_sigma->GetBinError(ipt+1)-0.1,mh1_Proton_sigma->GetBinContent(ipt)+Nsigma*mh1_Proton_sigma->GetBinError(ipt+1)+0.1);
     
      // electron ---------------------------------------
      if(5<ipt &&ipt<25){
       
	total_3->SetParameter(6,mh1_Electron_Constant->GetBinContent(ipt+1));
       
	total_3->SetParLimits(6,mh1_Electron_Constant->GetBinContent(ipt+1)-100,mh1_Electron_Constant->GetBinContent(ipt+1)+100);
      }
      total_3->SetParameter(7,mean_central);
      total_3->SetParameter(8,sigma_central);
      total_3->SetParLimits(7, mean_central-Mean_err,mean_central+Mean_err);
      total_3->SetParLimits(8, sigma_central-Sigma_err,sigma_central+Sigma_err);

     
	  
	  

      hh->Fit(total_3,"R","same",-10,15);
	  
     

    g1->SetParameter(0,total_3->GetParameter(6));
     g1->SetParameter(1,total_3->GetParameter(7));
     g1->SetParameter(2,total_3->GetParameter(8));
     
     g2->SetParameter(0,total_3->GetParameter(3));
     g2->SetParameter(1,total_3->GetParameter(4));
     g2->SetParameter(2,total_3->GetParameter(5));
     
     g3->SetParameter(0,total_3->GetParameter(0));
     g3->SetParameter(1,total_3->GetParameter(1));
     g3->SetParameter(2,total_3->GetParameter(2));
     
     g4->SetParameter(0,total_3->GetParameter(9));
     g4->SetParameter(1,total_3->GetParameter(10));
     g4->SetParameter(2,total_3->GetParameter(11));
     

      g1->SetParError(0,total_3->GetParError(6));
      g1->SetParError(1,total_3->GetParError(7));
      g1->SetParError(2,total_3->GetParError(8));
     
     g2->SetParError(0,total_3->GetParError(3));
     g2->SetParError(1,total_3->GetParError(4));
     g2->SetParError(2,total_3->GetParError(5));
     
     g3->SetParError(0,total_3->GetParError(0));
     g3->SetParError(1,total_3->GetParError(1));
     g3->SetParError(2,total_3->GetParError(2));
     
     g4->SetParError(0,total_3->GetParError(9));
     g4->SetParError(1,total_3->GetParError(10));
     g4->SetParError(2,total_3->GetParError(11));
 
     Purity=g1->Integral(-1,3)/total_3->Integral(-1,3);
     // cout<< " Pu"<<Purity<<endl;
     // g1->SetLineColor(4);
     // g1->Draw("same");

    
    }

  file_2->Close();

  return Purity;
     
  

}

void  Get_sts_Error(TH1F * const  mh1nsigmae_inclusive[], Float_t mean_central,Float_t Mean_err,Float_t sigma_central,Float_t Sigma_err,TH1D *& mh1purity_sts, TH1D * purity_sys ,Int_t Fit_Flag)
{
   
 
  return ;

  gStyle->SetOptFit(1111);
  cout<< purity_sys->GetBinContent(9)<< "  "<<  purity_sys->GetBinContent(10)<< endl;
  
  //  throw (-1);
  
  
  
  //    gStyle->SetOptStat(true);
  
  char buf[1024];
  TRandom3 *gRnd= new TRandom3(0);
  Int_t Nbins=400;
  
  // fstream purity_sys("purity_sys",)   
  
  TCanvas *c7=new TCanvas("c7","",800,600);
  TCanvas *c8=new TCanvas("c8","",800,600);
  TCanvas *c9=new TCanvas("c9","",800,600);
  c7->Divide(3,3);
  c8->Divide(3,3);
  c9->Divide(3,3);

  Int_t Npad=1;


  for(int ipt=0;ipt<Nbins_MB_46;ipt++)
    {
      
      Float_t BinCenter=  purity_sys->GetBinContent(ipt+1);
      
      Float_t Bin_next=purity_sys->GetBinContent(ipt+2);
      
      
      Double_t width=fabs(BinCenter-Bin_next);
      
      
      if(ipt==0)
	width=fabs(BinCenter-Bin_next)*2;
      
      if(ipt<0 && ipt <7)
	{
	  width=fabs(BinCenter-Bin_next)/6;
	  Nbins=1000;
	}
      if(ipt==7)
	
	width=fabs(BinCenter-Bin_next)/200.;
      
      if(ipt==8)
	
	width=fabs(BinCenter-Bin_next)/10.;
      
      if(ipt==9)
	
	
	width=fabs(BinCenter-Bin_next)/5.;
      
      if(ipt==10)
	
	width=fabs(BinCenter-Bin_next)/10.;
      
      if(ipt==12)
	width=fabs(BinCenter-Bin_next)*10.;
      if(ipt==13)
	width=fabs(BinCenter-Bin_next)*10.;
      
      if(ipt==14)
	width=fabs(BinCenter-Bin_next)*5.;
      if(ipt==15)
	width=fabs(BinCenter-Bin_next)*2.;
      
      
      if(ipt==16)
	width=fabs(BinCenter-Bin_next)*1.5;
      
      if(ipt==17)
	width=fabs(BinCenter-Bin_next)/2.;
      
      if(ipt==19)
	width=fabs(BinCenter-Bin_next)/2.;
      
      if(ipt==20)
	width=fabs(BinCenter-Bin_next)/4.;
      
      if(ipt==21)
	width=fabs(BinCenter-Bin_next);
      if(ipt==22)
	width=fabs(BinCenter-Bin_next)/4.;
      if(ipt==24)
	width=fabs(BinCenter-Bin_next);
      if(ipt==25)
	{
	  Bin_next=purity_sys->GetBinContent(24);
	}
      
      // Instead use some array width = widths[ipt], where widths has been filled in a preliminary step.
      
      
      
      // TH1D * purity_Fit=new TH1D("purity_Fit","purity_sts",100,BinCenter-width,BinCenter+width);
      TH1D * purity_Fit=new TH1D("purity_Fit","purity_sts",Nbins,BinCenter-width,BinCenter+width);
      purity_Fit->Sumw2();
      
      TF1 *f_Fit = new TF1(TString("f_Fit"),"gaus",BinCenter-width,BinCenter+width);
      
      
      
      if(ipt<8)
	Fit_Flag=2;
      else Fit_Flag=3;
      
      TH1F * nsigmae_inclusive_shift=(TH1F *) mh1nsigmae_inclusive[ipt]->Clone("nsigmae_inclusive_shift");
      
      for(Int_t i=1;i<1000;i++)
	// for(Int_t i=1;i<100;i++)
      	{
	  
	     
	
	  for(Int_t bin=1; bin<=mh1nsigmae_inclusive[ipt]->GetNbinsX();bin++)
	    {
	      if(mh1nsigmae_inclusive[ipt]->GetBinContent(bin) && mh1nsigmae_inclusive[ipt]->GetBinError(bin))
	    
		nsigmae_inclusive_shift->SetBinContent(bin,gRnd->Gaus(mh1nsigmae_inclusive[ipt]->GetBinContent(bin),mh1nsigmae_inclusive[ipt]->GetBinError(bin)));
	      
	    } 
	  
	  
	   Float_t purity =Fit_purity(nsigmae_inclusive_shift,  mean_central, Mean_err, sigma_central, Sigma_err,ipt,Fit_Flag);
	  
	  
	  //	  mh1mean->Fill(purity);
	  
	  //
	   purity_Fit->Fill(purity);
	}
      
    
      // TCanvas *c4=new TCanvas("c4","",800,600);
      
      // c4->cd();
      
      f_Fit->SetParameter(1,purity_sys->GetBinContent(ipt+1));
      if(ipt<9)
	c7->cd(Npad++);
      else if (ipt<18)
      c8->cd(Npad++);
      else if(18<=ipt)
	c9->cd(Npad++);
      if(Npad==10)
	Npad=1;
      purity_Fit->Fit(f_Fit,"RLL","",BinCenter-0.6*width,BinCenter+0.6*width);
      purity_Fit->Draw("same");
      
      mh1purity_sts->SetBinContent(ipt+1,purity_sys->GetBinContent(ipt+1));
      
      mh1purity_sts->SetBinError(ipt+1,f_Fit->GetParameter(2));
      //      f_Fit->Draw("same");
      
      //  mh1mean->Draw();
      
      // char buf[1024];
      //  sprintf(buf,"purity_Pt_%d",ipt+1);
      
      // sprintf(buf,"sts_purity%2.2f.gif",ipt);
      // c4->SaveAs(TString("Sts")+buf+TString(".gif"));
      
      //sprintf(buf,"purity_Pt_%d_to_%d_ThreeSigmaconstrain",int(pt_low*10),int(pt_high*10));
      
      // cout<<"--------------------------------------->"<<buf<<endl;
      // c3->SaveAs(TString(buf)+TString(".pdf"));    
    }
  
  c7->SaveAs("c1_test.pdf");
  c8->SaveAs("c2_test.pdf");
  c9->SaveAs("c3_test.pdf");
}	

void Get_sys_Error(TH1F * const  hh[], Float_t mean_central,Float_t Mean_err,Float_t sigma_central,Float_t Sigma_err,Float_t purity_one_sigma_constrain[],Float_t purity_two_sigma_constrain[],Float_t purity_three_sigma_constrain[],TH1D *&mh1purity,Int_t Fit_flag)
{

 
  TH1F * mh1purity_one_sigma=new TH1F("mh1purity_one_sigma","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1purity_two_sigma=new TH1F("mh1purity_two_sigma","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1purity_three_sigma=new TH1F("mh1purity_three_sigma","",Nbins_MB_46,pt_bin_MB_46);

  Float_t purity_average=0;
  Float_t purity_err;
  for(Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
    {

      //      cout<< mean_central<< "  "<< Mean_err<<"  "<< Sigma_err<<"  "<< "  "<< purity_two_sigma_constrain[2]<< endl;

      if(ipt<0)

	{      
	  Fit_flag=2;
	  purity_one_sigma_constrain[ipt]=Fit_purity(hh[ipt],  mean_central, Mean_err, sigma_central, Sigma_err, ipt,Fit_flag);
	  purity_two_sigma_constrain[ipt]=Fit_purity(hh[ipt],  mean_central, 2*Mean_err, sigma_central, 2*Sigma_err,ipt, Fit_flag);
	  purity_three_sigma_constrain[ipt]=Fit_purity(hh[ipt],  mean_central, 3*Mean_err, sigma_central, 3*Sigma_err, ipt,Fit_flag);
	  
	  cout << purity_one_sigma_constrain[ipt]<< "  "<<purity_two_sigma_constrain[ipt]<<" "<<purity_three_sigma_constrain[ipt]<<" sys  "<<endl;
	  mh1purity_one_sigma->SetBinContent(ipt+1,purity_one_sigma_constrain[ipt]);
	  mh1purity_one_sigma->SetBinError(ipt+1,0);

	  mh1purity_two_sigma->SetBinContent(ipt+1,purity_two_sigma_constrain[ipt]);

	  mh1purity_two_sigma->SetBinError(ipt+1,0);
	  mh1purity_three_sigma->SetBinContent(ipt+1,purity_three_sigma_constrain[ipt]);

	  mh1purity_three_sigma->SetBinError(ipt+1,0);
	  purity_average=(purity_one_sigma_constrain[ipt]+purity_two_sigma_constrain[ipt]+purity_three_sigma_constrain[ipt])/3.0;
	  mh1purity->SetBinContent(ipt+1,purity_average);
	  purity_err=fabs(purity_average-purity_one_sigma_constrain[ipt]);
	  
	  if(purity_err<fabs(purity_average-purity_two_sigma_constrain[ipt]))
	    purity_err=fabs(purity_average-purity_two_sigma_constrain[ipt]);
	  
	  if(purity_err<fabs(purity_average-purity_three_sigma_constrain[ipt]))
	    purity_err=fabs(purity_average-purity_three_sigma_constrain[ipt]);
	  
	 mh1purity->SetBinError(ipt+1,purity_err);

	 cout << purity_one_sigma_constrain[ipt]<< "  "<<purity_two_sigma_constrain[ipt]<<" "<<purity_three_sigma_constrain[ipt]<<" sys  "<< purity_average<<"  "<< purity_err<<endl;

	
	}
      else
	{   
	  Fit_flag=3;
	  purity_one_sigma_constrain[ipt]=Fit_purity(hh[ipt],  mean_central, Mean_err, sigma_central, Sigma_err, ipt ,Fit_flag);
	  purity_two_sigma_constrain[ipt]=Fit_purity(hh[ipt],  mean_central, 2*Mean_err, sigma_central, 2*Sigma_err, ipt ,Fit_flag);
	  purity_three_sigma_constrain[ipt]=Fit_purity(hh[ipt],  mean_central, 3*Mean_err, sigma_central, 3*Sigma_err, ipt, Fit_flag);

	  cout << purity_one_sigma_constrain[ipt]<< "  "<<purity_two_sigma_constrain[ipt]<<" "<<purity_three_sigma_constrain[ipt]<<" sys  "<<endl;
	 

	  mh1purity_one_sigma->SetBinContent(ipt+1,purity_one_sigma_constrain[ipt]);

	  mh1purity_one_sigma->SetBinError(ipt+1,0);
	  mh1purity_two_sigma->SetBinContent(ipt+1,purity_two_sigma_constrain[ipt]);
	  mh1purity_two_sigma->SetBinError(ipt+1,0);
	  mh1purity_three_sigma->SetBinContent(ipt+1,purity_three_sigma_constrain[ipt]);
	 
	  mh1purity_three_sigma->SetBinError(ipt+1,0);
	  purity_average=(purity_one_sigma_constrain[ipt]+purity_two_sigma_constrain[ipt]+purity_three_sigma_constrain[ipt])/3.0;
 mh1purity->SetBinContent(ipt+1,purity_average);	

 purity_err=fabs(purity_average-purity_one_sigma_constrain[ipt]);
	  
 if(purity_err<fabs(purity_average-purity_two_sigma_constrain[ipt]))
   purity_err=fabs(purity_average-purity_two_sigma_constrain[ipt]);
 
 if(purity_err<fabs(purity_average-purity_three_sigma_constrain[ipt]))
   purity_err=fabs(purity_average-purity_three_sigma_constrain[ipt]);
 
 mh1purity->SetBinError(ipt+1,purity_err);
	
}    
    }

  mh1purity_one_sigma->SetMarkerStyle(20);
  mh1purity_two_sigma->SetMarkerStyle(20);
  mh1purity_three_sigma->SetMarkerStyle(20);
  mh1purity->SetMarkerStyle(22);

  mh1purity_one_sigma->SetMarkerColor(1);
  mh1purity_two_sigma->SetMarkerColor(2);
  mh1purity_three_sigma->SetMarkerColor(3);
  mh1purity->SetMarkerColor(4);





  TCanvas *c1=new TCanvas("c1","",800,600);

  c1->cd();

  TH2F *H=new TH2F("H","H",20,0,4,10,0,1.2);
  H->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  H->GetYaxis()->SetTitle("purity");

  H->GetXaxis()->CenterTitle();
  H->GetYaxis()->CenterTitle();
  H->SetTitle("Inclusive Elelctron Purity differnet constrain");
  mh1purity->SetTitle("purity");
  H->Draw();
  mh1purity_one_sigma->Draw("samePX0");
  mh1purity_two_sigma->Draw("samePX0");
  mh1purity_three_sigma->Draw("samePX0");
  mh1purity->Draw("samePX0");
  
 TLegend *legend  = new TLegend(0.4,0.2,0.6,0.5);
 legend ->AddEntry(mh1purity_one_sigma,"One sigma devaration constrain","lpe");
 legend ->AddEntry(mh1purity_two_sigma,"Two sigma devaration constrain","lpe");
 legend ->AddEntry(mh1purity_three_sigma,"Three sigma devaration constrain","lpe");
 legend ->SetBorderSize(0);
  legend ->SetTextSize(0.045);
  legend ->SetFillStyle(0);
  legend ->Draw("same");


  c1->SaveAs("purity_Three_sigmaConstrain.pdf");



  


  /*
  
   purity_one_sigma_constrain=Fit_purity(hh,  mean_central, Mean_err, sigma_central, Sigma_err, Fit_flag);
   purity_two_sigma_constrain=Fit_purity(hh,  mean_central, Mean_err, sigma_central, 2*Sigma_err, Fit_flag);
   purity_three_sigma_constrain=Fit_purity(hh,  mean_central, Mean_err, sigma_central, 3*Sigma_err, Fit_flag);
  
  
  
  
  
  //  return   else   if(Fit_flag==3)
  

  Float_t purity_average=(purity_one_sigma_constrain+purity_two_sigma_constrain+purity_three_sigma_constrain)/3.0;

  
  Float_t sys=fabs(purity_one_sigma_constrain-purity_average);
  
  if(sys<fabs(purity_two_sigma_constrain-purity_average))
    sys=fabs(purity_two_sigma_constrain-purity_average);

  if(sys<fabs(purity_three_sigma_constrain-purity_average))
    sys=fabs(purity_three_sigma_constrain-purity_average);
  */
  // else 
  //  return 0.5;
  
}

void Draw_mean_sigma(TH1F * nsigmaE_inclusive[],Float_t Mean[],Float_t MeanErr[],Float_t Sigma[],Float_t SigmaErr[])
{

  char buf[1024]; 
  TH1F * mh1mean=new TH1F("mh1mean","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1meanUp=new TH1F("mh1meanUp","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1meanDown=new TH1F("mh1meanDown","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1sigma=new TH1F("mh1sigma","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1sigmaUp=new TH1F("mh1sigmaUp","",Nbins_MB_46,pt_bin_MB_46);
  TH1F * mh1sigmaDown=new TH1F("mh1sigmaDown","",Nbins_MB_46,pt_bin_MB_46);



  for(Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
    {
      mh1mean->SetBinContent(ipt+1,Mean[ipt]);
      mh1mean->SetBinError(ipt+1,MeanErr[ipt]);
      mh1meanUp->SetBinContent(ipt+1,Mean[ipt]+MeanErr[ipt]);
      mh1meanUp->SetBinError(ipt+1,MeanErr[ipt]);
      mh1meanDown->SetBinContent(ipt+1,Mean[ipt]-MeanErr[ipt]);
      mh1meanDown->SetBinError(ipt+1,MeanErr[ipt]);
      //
      mh1sigma->SetBinContent(ipt+1,Sigma[ipt]);
      mh1sigma->SetBinError(ipt+1,SigmaErr[ipt]);
      mh1sigmaUp->SetBinContent(ipt+1,Sigma[ipt]+SigmaErr[ipt]);
      mh1sigmaUp->SetBinError(ipt+1,SigmaErr[ipt]);
      mh1sigmaDown->SetBinContent(ipt+1,Sigma[ipt]-SigmaErr[ipt]);
      mh1sigmaDown->SetBinError(ipt+1,SigmaErr[ipt]);


    }  

  TF1 *pol_1=new TF1(TString("pol_1"),"[0]",0.2,3);
  TF1 *pol_2=new TF1(TString("pol_2"),"[0]",0.2,3);
  TF1 *pol_3=new TF1(TString("pol_3"),"[0]",0.2,3);
  TF1 *pol_4=new TF1(TString("pol_4"),"[0]",0.2,3);
  TF1 *pol_5=new TF1(TString("pol_5"),"[0]",0.2,3);
  TF1 *pol_6=new TF1(TString("pol_6"),"[0]",0.2,3);  
  
  pol_1->SetLineStyle(7);
  pol_3->SetLineStyle(7);
  pol_4->SetLineStyle(7);
  pol_6->SetLineStyle(7);
  
  pol_1->SetLineColor(1);
  pol_3->SetLineColor(1);
  pol_4->SetLineColor(2);
  pol_6->SetLineColor(2);
  
  TCanvas *c3=new TCanvas("c3","",800,600);
  c3->cd(); 
  TH2F *hh=new TH2F("hh","",10,0,4,100,-1,1.5);
  hh->SetTitle("Elecctron Mean and Sigma calibration");
  hh->GetXaxis()->SetTitle("p_{T}");


  hh->Draw();

  mh1mean->SetMarkerStyle(20);
  mh1sigma->SetMarkerStyle(4);
  
  
  mh1mean->Draw("same");
  mh1sigma->Draw("same");
 
  fstream Fit_Err("Electron_mean_sigma_Fit_err.dat",ios::trunc|ios::out);




  mh1mean->Fit(pol_2,"R","same",0.2,3);
  Fit_Err<<pol_2->GetParameter(0)<<"  "<<pol_2->GetParError(0)<<endl;
  mh1meanUp->Fit(pol_1,"R0","",0.2,3);
Fit_Err<<pol_1->GetParameter(0)<<"  "<<pol_1->GetParError(0)<<endl;
  mh1meanDown->Fit(pol_3,"R0","",0.2,3);
Fit_Err<<pol_3->GetParameter(0)<<"  "<<pol_3->GetParError(0)<<endl;

  mh1sigma->Fit(pol_5,"R","same",0.2,3);
Fit_Err<<pol_5->GetParameter(0)<<"  "<<pol_5->GetParError(0)<<endl;

  mh1sigmaUp->Fit(pol_4,"R0","",0.2,3);
Fit_Err<<pol_4->GetParameter(0)<<"  "<<pol_4->GetParError(0)<<endl;
  mh1sigmaDown->Fit(pol_6,"R0","",0.2,3);
Fit_Err<<pol_6->GetParameter(0)<<"  "<<pol_6->GetParError(0)<<endl;

  /*
  pol_1->SetLineStyle(7);
  pol_3->SetLineStyle(7);
  pol_4->SetLineStyle(7);
  pol_6->SetLineStyle(7);
  
  pol_1->SetLineColor(1);
  pol_3->SetLineColor(1);
  pol_4->SetLineColor(1);
  pol_6->SetLineColor(1);
  */

  pol_1->SetParameter(0,pol_1->GetParameter(0)+0.02);
  pol_3->SetParameter(0,pol_3->GetParameter(0)-0.02);

  pol_4->SetParameter(0,pol_4->GetParameter(0)+0.02);
  pol_6->SetParameter(0,pol_6->GetParameter(0)-0.02);


  pol_1->Draw("same");
  pol_3->Draw("same");
  pol_4->Draw("same");
  pol_6->Draw("same");

  Float_t  mean_central=pol_2->Eval(2);
  Float_t  mean_up=pol_1->Eval(2);
  Float_t  mean_down=pol_3->Eval(2);
  Float_t  sigma_central=pol_5->Eval(2);
  Float_t  sigma_up=pol_4->Eval(2);
  Float_t  sigma_down=pol_6->Eval(2);



  cout<< "\n"<< sigma_central<<"  "<< sigma_down<<endl;
  Float_t Mean_err=fabs(mean_central-mean_up)>fabs(mean_central-mean_down)?fabs(mean_central-mean_up):fabs(mean_central-mean_down);

  Float_t Sigma_err=fabs(sigma_central-sigma_up)>fabs(sigma_central-sigma_down)?fabs(sigma_central-sigma_up):fabs(sigma_central-sigma_down);

  sprintf(buf,"Mean:  %2.2f +- %4.5f",mean_central,Mean_err);
  drawLatex(0.6,0.8,buf,70,0.035,2);

sprintf(buf,"Sigma:  %2.2f +- %4.5f",sigma_central,Sigma_err);
  drawLatex(0.6,0.75,buf,70,0.035,2);


  cout<<"Mean_err  "<< Mean_err<< " sigma Error "<< Sigma_err<<endl;


 TFile *file=new TFile("Electron_Mean_igma.root","RECREATE");
 file->cd();
 pol_1->Write();
 pol_3->Write();
 pol_4->Write();
 pol_6->Write();

 mh1mean->Write();
 mh1sigma->Write();
 file->Close();

  /*  test the mean without nsigma e cut on partner

  TFile * file_1=new TFile("wiout_nsigmaE_cut.root","read");

  TF1 * mean=(TF1*) file_1->Get("pol_2");
 TF1 * sigma=(TF1*) file_1->Get("pol_6");

 mean->SetLineColor(4);
 sigma->SetLineColor(4);
 // mean ->Draw("same");
 // sigma->Draw("same");

 */

  TLegend *legend  = new TLegend(0.4,0.4,0.6,0.6);
  legend ->AddEntry(mh1mean,"Electron Mean","lpe");
  legend ->AddEntry(mh1sigma,"Electron Sigma","lpe");
  legend ->SetBorderSize(0);
  legend ->SetTextSize(0.045);
  legend ->SetFillStyle(0);
  legend ->Draw("same");
  c3->SaveAs("mean_sigma_electron.pdf");


  // return ;
  Get_purity(nsigmaE_inclusive,Mean,MeanErr,Sigma,SigmaErr,mean_central,sigma_central,Mean_err,Sigma_err); 



  // nSigma_cut_efficiency(mh1mean,mh1meanUp,mh1meanDown,mh1sigma,mh1sigmaUp,mh1sigmaDown,mean_central,sigma_central,Mean_err,Sigma_err);

}

void nSigma_cut_efficiency(TH1F * mh1mean,TH1F * mh1meanUp,TH1F * mh1meanDown,TH1F * mh1sigma,TH1F * mh1sigmaUp,TH1F * mh1sigmaDown,Float_t mean_central,Float_t sigma_central,Float_t Mean_err,Float_t Sigma_err)
{

  TH1F * nsigmaE_cut_efficiency=new TH1F("nsigmaE_cut_efficiency","",Nbins_MB_46,pt_bin_MB_46);
  // nsigma cut efficiency statistics uncertainty
  TF1 *f_nsigma=new TF1(TString("f_nsigma"),"gaus",-5,5);
  f_nsigma->SetParameter(0,1);
  f_nsigma->SetParError(0,0);
  f_nsigma->SetParameter(1,mean_central);
  f_nsigma->SetParError(1,Mean_err);
  f_nsigma->SetParameter(2,sigma_central);
  f_nsigma->SetParError(2,Sigma_err);

  //  f_nsigma->Draw();
  
  Float_t Nsigma_Cut_efficiency=f_nsigma->Integral(-1,3)/f_nsigma->Integral(-5,5);
  cout<< "\n "<< Nsigma_Cut_efficiency<<endl;
  
  f_nsigma->SetParameter(1,mean_central+Mean_err);
  Float_t Nsigma_Cut_efficiency_1=f_nsigma->Integral(-1,3)/f_nsigma->Integral(-5,5);

  f_nsigma->SetParameter(1,mean_central);
  f_nsigma->SetParameter(2,sigma_central+Sigma_err);
  Float_t Nsigma_Cut_efficiency_2=f_nsigma->Integral(-1,3)/f_nsigma->Integral(-5,5);
  cout<< "\n "<< Nsigma_Cut_efficiency<<"   "<< Nsigma_Cut_efficiency_1 << "  "<< Nsigma_Cut_efficiency_2<< endl;
  
  Float_t Nsigma_Cut_efficiency_sts=fabs(Nsigma_Cut_efficiency-Nsigma_Cut_efficiency_1)>fabs(Nsigma_Cut_efficiency-Nsigma_Cut_efficiency_2)?fabs(Nsigma_Cut_efficiency-Nsigma_Cut_efficiency_1):fabs(Nsigma_Cut_efficiency-Nsigma_Cut_efficiency_2);


  TF1 *pol_1=new TF1(TString("pol_1"),"[0]+[1]*x",0,4);
  TF1 *pol_2=new TF1(TString("pol_2"),"[0]+[1]*x",0,4);
  TF1 *pol_3=new TF1(TString("pol_3"),"[0]+[1]*x",0,4);
  TF1 *pol_4=new TF1(TString("pol_4"),"[0]+[1]*x",0,4);
  TF1 *pol_5=new TF1(TString("pol_5"),"[0]+[1]*x",0,4);
  TF1 *pol_6=new TF1(TString("pol_6"),"[0]+[1]*x",0,4); 

  TCanvas * c4=new TCanvas("c4","",800,600);
  c4->cd();
  TH2F *hh=new TH2F("hh","",10,0,4,100,-1,1.5);
  hh->Draw(); 
  mh1meanUp->Fit(pol_1,"R","same",0.2,3);
  mh1mean->Fit(pol_2,"R","same",0.2,3);
  mh1meanDown->Fit(pol_3,"R","same",0.2,3);
  mh1sigmaUp->Fit(pol_4,"R","same",0.2,3);
  mh1sigma->Fit(pol_5,"R","same",0.2,3);
 
  mh1sigmaDown->Fit(pol_6,"R","same",0.2,3);
 
  c4->SaveAs("pol_2.pdf");
  
  


     mh1meanUp->Draw("same");
  // nsigma cut systmatics uncertainty


  for(Int_t ipt=0;ipt<Nbins_MB_46;ipt++)
    {

      Float_t meanError=fabs(pol_1->Eval((pt_bin_MB_46[ipt]+pt_bin_MB_46[ipt+1])/2.0)-pol_3->Eval((pt_bin_MB_46[ipt]+pt_bin_MB_46[ipt+1])/2.0));
      Float_t SigmaError=fabs(pol_4->Eval((pt_bin_MB_46[ipt]+pt_bin_MB_46[ipt+1])/2.0)-pol_6->Eval((pt_bin_MB_46[ipt]+pt_bin_MB_46[ipt+1])/2.0));

      Float_t mean_central=pol_2->Eval((pt_bin_MB_46[ipt]+pt_bin_MB_46[ipt+1])/2.0);
      Float_t sigma_central=pol_5->Eval((pt_bin_MB_46[ipt]+pt_bin_MB_46[ipt+1])/2.0);
      f_nsigma->SetParameter(1,mean_central);
      
      f_nsigma->SetParError(1,meanError);
      f_nsigma->SetParameter(2,sigma_central);
      f_nsigma->SetParError(2,SigmaError);
      
      Float_t Nsigma_Cut_efficiency_pol_2=f_nsigma->Integral(-1,3)/f_nsigma->Integral(-5,5);
           
      Float_t  Nsigma_Cut_efficiency_sys=fabs(Nsigma_Cut_efficiency_pol_2-Nsigma_Cut_efficiency);
      //    cout<< meanError<< "  "<< SigmaError<<  "  /n"<< Nsigma_Cut_efficiency_pol_2<< "   "<<Nsigma_Cut_efficiency<<endl;
      
      // cout<< " mean err"<< meanError<< "   "<< SigmaError<<"  "<< Nsigma_Cut_efficiency_sys<< endl;
      Float_t       nsigmaE_err=sqrt(Nsigma_Cut_efficiency_sts*Nsigma_Cut_efficiency_sts+ Nsigma_Cut_efficiency_sys*Nsigma_Cut_efficiency_sys);
      cout<< " mean err"<< meanError<< "   "<< SigmaError<<"  "<< Nsigma_Cut_efficiency_sys<< "  "<< nsigmaE_err<< "  "<< endl;
      
      nsigmaE_cut_efficiency ->SetBinContent(ipt+1,Nsigma_Cut_efficiency);
      nsigmaE_cut_efficiency->SetBinError(ipt+1,nsigmaE_err);
      
    }
  
  
  
  
 cout<< " nsigma Cut"<<endl;
  
 TCanvas *c1=new TCanvas("c1","",800,600);
 c1->cd();
 TH2F *h=new TH2F("h","",10,0,4.5,20,0,1);
 h->Draw();
 h->SetTitle("nsigma Electron cut efficiency");
 h->GetXaxis()->SetTitle("p_{T}");
 h->GetYaxis()->SetTitle("efficiency");
 
 nsigmaE_cut_efficiency->SetMarkerStyle(20);
 nsigmaE_cut_efficiency->SetMarkerColor(2);
 nsigmaE_cut_efficiency->Draw("samePE1");
 
 c1->SaveAs("nsigma_cut_eff.pdf");
 
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
