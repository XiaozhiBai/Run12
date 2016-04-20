
#include <iostream>
#include<iomanip>
#include <fstream>
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
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "binning_util.h"
#include "/Users/xiaozhi/NPE_work/NPE_high_pt/Bin_highpt.h"
//#define Nbins_HT 13
using namespace std;

//Double_t xx[14]={2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.5,10,14};

//void Get_NoPhotonic_ele(const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *);
void Get_NoPhotonic_ele(const TH1F *, TH1F *, TH1F *,const TH1F *,TH1F *, TH1F*,const TH1F *,const TH1F *, TH1F *&, TH1F *&, TH1F *&);//,const TH1F *,const TH1F *);
void Draw_Inclusive_electron(const TH1F * ,const TH1F *);
void Draw_Photonic_electron(const TH1F * ,const TH1F * );
void GetSpectra_yeild(const TH1F *, TH1F *, TH1F *,const TH1F *,const TH1F *,const TH1F*,TH1F *,TH1F * );
void Draw_Inv_yield(TH1F* ,TH1F *,TH1F *,Double_t sys[],Double_t sts[],TH1F *,TH1F *);
void Draw_Efficiency(TH1F *,TH1F *,TH1F *,TH1F *);
void Get_sys_Error(const TH1F* ,const TH1F *,const TH1F * ,const TH1F *, const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F * ,TH1F* ,Double_t sys[]);
void Get_sts_Error(const TH1F* ,const TH1F *,const TH1F * ,const TH1F *, const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *,const TH1F *,TH1F *,Double_t sts[]);





TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
void setpad(TPad *,float , float , float , float );

 TH1F * Correc=new TH1F("Correc","",Nbins_HT,Pt_bin_HT);
void NPE_yeild()
{
  char buf[1024];


  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1,"Y");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13); 

  gStyle->SetEndErrorSize(4);

  //get the PHotonic electron 
  TFile * infile_Ph_eff =new TFile("Photonic_re_Efficiency_HT.root","read");
 
  TH1F * Ph_eff=(TH1F *)infile_Ph_eff->Get("mh1_MB");
    
    Ph_eff->Draw();

// return;
  TFile * infile_EMC_eff =new TFile("EMC_cut_efficiency_HT.root","read");
  TH1F * EMC_eff=(TH1F *)infile_EMC_eff->Get("EMC_cut_effciency");
 
    EMC_eff->Draw();
  // return;
  
  TFile *file_nsigma_cut =new TFile("Nsigma_cut_Efficiency_HT.root","read");
  TH1F * nsigma_Cut_efficiency=(TH1F *) file_nsigma_cut->Get("nsigma_cut_Efficiency");

    nsigma_Cut_efficiency->Draw();
  //  return;
 

  TFile *file_Tracking_eff =new TFile("Tracking_efficiency_HT.root","read");
  TH1F * Tracking_eff=(TH1F * )file_Tracking_eff->Get("Tracking_Efficiency");
  
   Tracking_eff->Draw();
// return;

  // this is the purity Photonic electron and Inclusive electron 
  TFile * infile_Ph_In=new TFile("Inclusive_Photo_HT2.root","read");
  
  // infile_Ph_In->Draw();
  // return;
 TH1F * mh1Inclusive_ele = (TH1F *) infile_Ph_In->Get("mh1Inclusive");
 TH1F * HT0_Inclusive_ele = (TH1F *) infile_Ph_In->Get("HT0_Inclusive");
 TH1F * HT2_Inclusive_ele = (TH1F *) infile_Ph_In->Get("HT2_Inclusive");


// mh1Inclusive_ele->Draw("pE");
// return; 
TH1F * mh1Photonic_ele = (TH1F *) infile_Ph_In->Get("mh1Photonic");
TH1F * HT0_Photonic_ele = (TH1F *) infile_Ph_In->Get("PhoHT0");
TH1F * HT2_Photonic_ele = (TH1F *) infile_Ph_In->Get("PhoHT2");
//mh1Photonic_ele->Draw("PE1");
//return;
 TFile * infile_purity=new TFile("Purity_HT.root","read");
 
 TH1F * mh1purity = (TH1F *) infile_purity->Get("purity");
 //mh1purity->Draw("p");
 


// return;
 TH1F *mh1NPE;
 TH1F *HT0_NPE;
 TH1F *HT2_NPE;
 Get_NoPhotonic_ele( mh1Inclusive_ele,HT0_Inclusive_ele,HT2_Inclusive_ele,mh1Photonic_ele,HT0_Photonic_ele,HT2_Photonic_ele,Ph_eff,mh1purity,mh1NPE,HT0_NPE,HT2_NPE);
 mh1NPE->Draw();

  HT0_NPE->Draw();

  HT0_Inclusive_ele->Draw();

  mh1NPE->Draw();

  //  return;
 
 GetSpectra_yeild( mh1NPE, HT0_NPE,HT2_NPE, nsigma_Cut_efficiency,EMC_eff,Tracking_eff,mh1purity,Ph_eff);
 //Tracking_eff  ->Draw();
 // cout<<mh1NPE->GetBinError(3)<<endl;
 // mh1NPE->Draw();
}
void Get_NoPhotonic_ele(const TH1F * mh1Inclusive_ele, TH1F * HT0_Inclusive,TH1F *HT2_Inclusive, const TH1F *mh1Photonic_ele, TH1F *HT0_Photonic_ele, TH1F * HT2_Photonic_ele, const TH1F *Ph_eff,const TH1F *mh1purity,TH1F * &NPE,TH1F * &HT0_NPE, TH1F * &HT2_NPE)
  
{

  

  
  TH1F::SetDefaultSumw2();
  
  gStyle->SetOptStat(00000);
  
  char buf[1024];
  
  TH1F * Inclusive_ele=(TH1F *) mh1Inclusive_ele->Clone("Inclusive_ele");
  TH1F * Photonic_ele=(TH1F *) mh1Photonic_ele->Clone("Photonic_ele");
  Inclusive_ele->Draw();
  
  // throw (-1);
  TH1F * PH_eff=(TH1F *) Ph_eff->Clone("PH_eff");
  TH1F * purity=(TH1F *) mh1purity->Clone("purity");
  Inclusive_ele->Multiply(purity);
  HT0_Inclusive->Multiply(purity);
  HT2_Inclusive->Multiply(purity);
  Photonic_ele->Divide( PH_eff);
  HT0_Photonic_ele->Divide(Ph_eff);
  HT2_Photonic_ele->Divide(Ph_eff);

  NPE=(TH1F*)  Inclusive_ele->Clone("NPE");
  HT0_NPE=(TH1F*)  HT0_Inclusive->Clone("HT0_NPE");
  HT2_NPE=(TH1F*)  HT2_Inclusive->Clone("HT2_NPE");
  NPE->Add(Photonic_ele,-1);
  HT0_NPE->Add(HT0_Photonic_ele,-1);
  HT2_NPE->Add(HT2_Photonic_ele,-1);

  TH1F * Inclusive_ele_sts=(TH1F *) mh1Inclusive_ele->Clone("Inclusive_ele_sts");
  TH1F * Photonic_ele_sts=(TH1F *) mh1Photonic_ele->Clone("Photonic_ele_sts");
  TH1F * PH_eff_sts=(TH1F *) Ph_eff->Clone("PH_eff_sts");
  TH1F * purity_sts=(TH1F *) mh1purity->Clone("purity_sts");

  TH1F * Inclusive_ele_sys=(TH1F *) mh1Inclusive_ele->Clone("Inclusive_ele_sys");
  TH1F * Photonic_ele_sys=(TH1F *) mh1Photonic_ele->Clone("Photonic_ele_sys");
  TH1F * PH_eff_sys=(TH1F *) Ph_eff->Clone("PH_eff_sys");
  TH1F * purity_sys=(TH1F *) mh1purity->Clone("purity_sys");
  



  for(Int_t i=1;i<=Nbins_HT;i++)
    {
      PH_eff_sts->SetBinError(i,0);
      purity_sts->SetBinError(i,0);
      Inclusive_ele_sys->SetBinError(i,0);
      Photonic_ele_sys->SetBinError(i,0);
    }

  Inclusive_ele_sts->Multiply(purity_sts);
  Photonic_ele_sts->Divide(PH_eff_sts);
  
TH1F *  NPE_sts= (TH1F *) Inclusive_ele_sts->Clone("NPE_sts");
  NPE_sts->Add(Photonic_ele_sts,-1);


    // systimatic uncertainty
    Inclusive_ele_sys->Multiply(purity_sys);
    Photonic_ele_sys->Divide(PH_eff_sys);

    TH1F *  NPE_sys= (TH1F *) Inclusive_ele_sts->Clone("NPE_sys");
    NPE_sys->Add(Photonic_ele_sys,-1);

    TH1F * Ration_NPE_PE=(TH1F*)NPE_sts->Clone("Ration_NPE_PE"); 
    Ration_NPE_PE->Divide( Photonic_ele_sts);

    NPE_sys->Divide(Photonic_ele_sys);

    fstream  outfile_Ratio12_sys("NPE_PE.txt",ios::trunc|ios::out);
    for(Int_t i=1;i<=Nbins_HT;i++)
      {
	outfile_Ratio12_sys<<NPE_sys->GetBinCenter(i)<<" "<<NPE_sys->GetBinContent(i)<<"  "<<Ration_NPE_PE->GetBinError(i)<<"  "<<NPE_sys->GetBinError(i)<<endl;
	
    }
    
    TGraphErrors* Npe_Pe_stat_err = new TGraphErrors("NPE_PE.txt","%lg %lg %lg");
    TGraphErrors* Npe_Pe_sys_err = new TGraphErrors("NPE_PE.txt","%lg %lg %*lg %lg");
    
    Npe_Pe_sys_err->SetLineColor(kBlue);
    
    TCanvas * c2=new TCanvas("c2","Inclusive Electron ",0,0,800,600);
    // setpad(c2,0.12,0.02,0.02,0.12);
    c2->cd();
    Ration_NPE_PE->SetMarkerSize(1);
  Ration_NPE_PE->SetMarkerStyle(20);
  Ration_NPE_PE->SetMarkerColor(2);
  Ration_NPE_PE->SetTitle("NPE and PHE Ratio  ");
  Ration_NPE_PE->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  Ration_NPE_PE->GetYaxis()->SetTitle("NPE/PHE");

  Ration_NPE_PE->GetXaxis()->CenterTitle(true);
  Ration_NPE_PE->GetYaxis()->CenterTitle(true);
  Ration_NPE_PE->GetXaxis()->SetTitleOffset(1);
  Ration_NPE_PE ->GetXaxis()->SetTitleOffset(1);
  Ration_NPE_PE->GetXaxis()->SetTitleSize(0.05);
  Ration_NPE_PE->GetYaxis()->SetTitleSize(0.05);


  Ration_NPE_PE->GetYaxis()->SetRangeUser(0,2);
  Ration_NPE_PE->Draw("PE1");
  Npe_Pe_sys_err->Draw("P[]:same");
  sprintf(buf,"NPE Over photonic");
  // drawLatex(0.35,0.85,buf,70,0.04,1);
  // c2->SetLogy(1);
  c2->SaveAs(TString("Ration_NPE_PHOE.gif"));

  //  throw (-1);  

}
// GetSpectra_yeild( mh1NPE, nsigma_Cut_efficiency,EMC_eff,Tracking_eff);
void GetSpectra_yeild(const TH1F  *NPE_h, TH1F * HT0_NPE, TH1F * HT2_NPE,const TH1F  * nsigmaE_Cut_efficiency_h,const TH1F  * EMC_Cut_efficiency_h,const TH1F *Track_recon_efficiency_h, TH1F * purity ,TH1F * PhEfficiency)
{
  TH1F:: SetDefaultSumw2();
    

  TCanvas *c4=new TCanvas("c4","c4",800,600);
  c4->cd();
  
  TH1F * N_npe=(TH1F *)NPE_h->Clone("N_npe");
  
  TH1F  * nsigmaE_Cut_efficiency=(TH1F *)nsigmaE_Cut_efficiency_h->Clone("nsigmaE_Cut_efficiency");
  
  TH1F  * EMC_Cut_efficiency=(TH1F  *) EMC_Cut_efficiency_h->Clone("EMC_Cut_efficiency");
  TH1F * Track_recon_efficiency=(TH1F*) Track_recon_efficiency_h->Clone("Track_recon_efficiency_h");
  // Draw_Efficiency(nsigmaE_Cut_efficiency,EMC_Cut_efficiency,Track_recon_efficiency);
    
  ///!!!!!!!!!!!!!!
  TFile *file=new TFile("Inclusive_Photo_HT2.root","read");


  
  // TH1F *HT0_Inclusive=(TH1F *) file->Get("HT0_Inclusive");
  // TH1F *HT2_Inclusive=(TH1F *) file->Get("HT2_Inclusive");
  TH1F *Inclusive=(TH1F *) file->Get("mh1Inclusive");
  TH1F *Photonic=(TH1F *) file->Get("mh1Photonic");
  
  // Photonic ->Draw();
  cout<<"11"<<endl;  
  // return;
  // HT0_Inclusive->Multiply(purity);
  // HT0_Inclusive->Divide(nsigmaE_Cut_efficiency); 
  
  // HT2_Inclusive->Multiply(purity);
  // HT2_Inclusive->Divide(nsigmaE_Cut_efficiency); 
  
  
  //  Float_t xx[9]={2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,16.5};
  
 
  TH1F * Correc=new TH1F("Correc","",Nbins_HT,Pt_bin_HT);
  TH1F * nor_BBC=new TH1F("nor_BBC","",Nbins_HT,Pt_bin_HT);
  TH1F * nor_BBCTOF=new TH1F("nor_BBCTOF","",Nbins_HT,Pt_bin_HT);
  TH1F * Combine_Nor=new TH1F("Combine_Nor","",Nbins_HT,Pt_bin_HT);     
  TH1F * BBC_trigger_eff=new TH1F("BBC_trigger_eff","",Nbins_HT,Pt_bin_HT);
  TH1F * Lumi=new TH1F("Lumi","",Nbins_HT,Pt_bin_HT);

  // return;

     TFile *Trg=new TFile("TrigEfficiency.root","read");
     TH1F *Trig_efficiency_HT0=(TH1F * ) Trg->Get("hh_11_0");
     TH1F *Trig_efficiency_HT2=(TH1F * ) Trg->Get("hh_18_0");
     TH1F * combile_Trg_efficiency=(TH1F *) Trig_efficiency_HT0->Clone("combile_Trg_efficiency");

     //     Trig_efficiency_HT2->Draw();

    
   
     TFile * Normalize=new TFile("hist_4_27.root","read");
     
     TH1F * mh1Vz_BBCMBPS =(TH1F *) Normalize->Get("mh1Vz_BBCMBPS");
     TH1F * mh1Vz_BBCMBTOF0PS =(TH1F *) Normalize->Get("mh1Vz_BBCMBTOF0PS");
     
     Double_t N_BBCMB=mh1Vz_BBCMBPS->Integral();
     Double_t N_BBCMBTOF=mh1Vz_BBCMBTOF0PS->Integral();


     
     Draw_Efficiency(nsigmaE_Cut_efficiency,EMC_Cut_efficiency,Track_recon_efficiency,combile_Trg_efficiency);
     
     for(Int_t ipt=1;ipt<=Nbins_HT;ipt++)
       {
	 nor_BBC->SetBinContent(ipt,N_BBCMB);
	 nor_BBC->SetBinError(ipt,0);
	 
	 nor_BBCTOF->SetBinContent(ipt,N_BBCMBTOF);
	 nor_BBCTOF->SetBinError(ipt,0);
	 
	 if( ipt<5)
	   {
	     Combine_Nor->SetBinContent(ipt,N_BBCMBTOF);
	     Combine_Nor->SetBinError(ipt,0);
	    }
	 else 
	   { 
	     Combine_Nor->SetBinContent(ipt,N_BBCMB);
	     Combine_Nor->SetBinError(ipt,0);
	     
	     combile_Trg_efficiency->SetBinContent(ipt,Trig_efficiency_HT2->GetBinContent(ipt));
	     combile_Trg_efficiency->SetBinError(ipt,Trig_efficiency_HT2->GetBinError(ipt));
	   }
	 Float_t pt=NPE_h->GetBinCenter(ipt);
	 Float_t correction_fector= 1/( 4*3.14159*1.4*pt*NPE_h->GetBinWidth(ipt));
	 Correc->SetBinContent(ipt,correction_fector);
	 Correc->SetBinError(ipt,0);
	 BBC_trigger_eff->SetBinContent(ipt,0.866);
	 BBC_trigger_eff->SetBinError(ipt,0.08);
	 Lumi->SetBinContent(ipt,30);
	 Lumi->SetBinError(ipt,2.4);
       }   
     
     TFile *File_jpsi=new TFile("Jpsi.root","read");
     
     TH1F * JPsi=(TH1F *) File_jpsi->Get("JPsi");
     TH1F * Before_Jpsi_Conrection =(TH1F *) N_npe->Clone("Before_Jpsi_Conrection");
     // N_npe->Add(JPsi,-1);  
     
     //  N_npe->Draw();
     // return;
     Double_t sys[Nbins_HT]={0,0,0,0,0,0,0,0,0,0,0};
     Double_t sts[Nbins_HT]={0,0,0,0,0,0,0,0,0,0,0};
     
      Get_sys_Error(Inclusive,Photonic,purity,PhEfficiency,nsigmaE_Cut_efficiency,EMC_Cut_efficiency,Track_recon_efficiency,combile_Trg_efficiency,Correc,Combine_Nor,Lumi,JPsi,sys);
         Get_sts_Error(Inclusive,Photonic,purity,PhEfficiency,nsigmaE_Cut_efficiency,EMC_Cut_efficiency,Track_recon_efficiency,combile_Trg_efficiency,Correc,Combine_Nor,Lumi,JPsi,sts);

	 cout<<      Inclusive->GetNbinsX()<<"  "<<Photonic->GetNbinsX() <<endl;

     cout<< " sys_error"<<sts[1]<<endl;
     cout<< " sts_error"<<sts[2]<<endl;


     // HT0_Inclusive->Multiply(purity);
     // HT2_Inclusive->Multiply(purity);
     
     // HT0_NPE->Divide(nsigmaE_Cut_efficiency); 
     // HT2_NPE->Divide(nsigmaE_Cut_efficiency);
     N_npe->Divide(nsigmaE_Cut_efficiency);
     
     // N_npe->Draw();
     // return;
     // HT0_NPE->Divide(EMC_Cut_efficiency); 
     // HT2_NPE->Divide(EMC_Cut_efficiency);
     N_npe->Divide(EMC_Cut_efficiency);     
     
     
     // HT0_NPE->Divide(Track_recon_efficiency); 
     // HT2_NPE->Divide(Track_recon_efficiency);
     N_npe->Divide(Track_recon_efficiency);
     
     // HT0_NPE->Divide(Trig_efficiency_HT0);
     //  HT2_NPE->Divide(Trig_efficiency_HT2);
        N_npe->Divide(combile_Trg_efficiency);
     
     
     // N_npe->Draw();
     // return;
     
     HT0_NPE->Multiply(Correc);
     HT2_NPE->Multiply(Correc);
     N_npe->Multiply(Correc);
     
    
     //  N_npe->Draw();
     // return;
      
     HT0_NPE->Divide(nor_BBCTOF);
     HT2_NPE->Divide(nor_BBC);
     N_npe->Divide(Combine_Nor);
     

     HT0_NPE->Multiply(Lumi);
     HT2_NPE->Multiply(Lumi);
     N_npe->Multiply(Lumi);

      N_npe->Add(JPsi,-1);  
     // HT0_NPE->Add(JPsi,-1);
     // HT2_NPE->Add(JPsi,-1);


     N_npe->Draw();
     TFile * File_npe=new TFile("NPE_spectra.root","recreate");
     N_npe->Write();
     File_npe->Close();
     

     Draw_Inv_yield(HT0_NPE,HT2_NPE,N_npe,sys,sts,Before_Jpsi_Conrection,JPsi);

}
void Get_sts_Error(const TH1F * mInclusive,const TH1F *mPhotonicE,const TH1F * mpurity,const TH1F * mPhotonicEReEff, const TH1F * mnsigmaE_Cut_efficiency,const TH1F * mEMC_Cut_efficiency,const TH1F * mTrack_recon_efficiency,const TH1F * mcombile_Trg_efficiency,const TH1F * mCorrec,const TH1F *mCombine_Nor,const TH1F * mLumi, TH1F * mJPsi, Double_t sts[])
{
    TH1::SetDefaultSumw2();
  
  TH1F * Inclusive=(TH1F *) mInclusive->Clone("Inclusive");
  TH1F * PhotonicE=(TH1F *) mPhotonicE->Clone("PhotonicE");
  TH1F * PhotonicEReEff=(TH1F *) mPhotonicEReEff->Clone("PhotonicEReEff,");
  TH1F * purity=(TH1F *) mpurity->Clone("purity");
  
  TH1F * nsigmaE_Cut_efficiency=(TH1F *) mnsigmaE_Cut_efficiency->Clone("nsigmaE_Cut_efficiency");
  TH1F * EMC_Cut_efficiency=(TH1F *) mEMC_Cut_efficiency->Clone("EMC_Cut_efficiency");
  TH1F * Track_recon_efficiency=(TH1F *) mTrack_recon_efficiency->Clone("Track_recon_efficiency");
  TH1F * combile_Trg_efficiency=(TH1F *) mcombile_Trg_efficiency->Clone("combile_Trg_efficiency");
  TH1F * Correc=(TH1F *) mCorrec->Clone("Correc");
  TH1F * Combine_Nor=(TH1F *) mCombine_Nor->Clone("Combine_Nor");
  TH1F * Lumi=(TH1F *) mLumi->Clone("mLumi");
  TH1F * JPsi=(TH1F *) mJPsi->Clone("JPsi");

  cout<< Inclusive->GetNbinsX()<<"  Nbinx"<<endl;
    
  for( Int_t i=1;i<=Inclusive->GetNbinsX();i++)
    {
      Inclusive->SetBinError(i,0);
      PhotonicE->SetBinError(i,0);
      
    }
  Inclusive->Sumw2();
  PhotonicE->Sumw2();
  Inclusive->Multiply(purity);
  PhotonicE->Divide(PhotonicEReEff);
  
  TH1F * NPE=(TH1F *)  Inclusive->Clone("NPE");
  NPE->Sumw2();
  NPE->Add(PhotonicE,-1);
  NPE->Divide(nsigmaE_Cut_efficiency);
  NPE->Divide(EMC_Cut_efficiency);
  NPE->Divide(Track_recon_efficiency);
  NPE->Divide(combile_Trg_efficiency);
  NPE->Multiply(Correc);
  NPE->Divide(Combine_Nor);
  NPE->Multiply(Lumi);
    NPE->Add(JPsi,-1); 

  /*
    TFile * file=new TFile("fonll200gev.root","read");
  TGraphErrors *gru=(TGraphErrors *) file->Get("gFONLLu");
  TGraphErrors *grc=(TGraphErrors *) file->Get("gFONLLc");
  TGraphErrors *grl=(TGraphErrors *) file->Get("gFONLLl");
  
  gru->Draw("same");
  grc->Draw("same");
  grl->Draw("same");

  //c1->SaveAs("HHH.gif");
  
  */

  for(Int_t i=1;i<=NPE->GetNbinsX();i++)
    
    {
      sts[i-1]=NPE->GetBinError(i);
      cout<<sts[i-1]<<endl;     
      
    }
  // return 0.1; 

} 

void Get_sys_Error(const TH1F * mInclusive,const TH1F *mPhotonicE,const TH1F * mpurity,const TH1F * mPhotonicEReEff, const TH1F * mnsigmaE_Cut_efficiency,const TH1F * mEMC_Cut_efficiency,const TH1F * mTrack_recon_efficiency,const TH1F * mcombile_Trg_efficiency,const TH1F * mCorrec,const TH1F *mCombine_Nor,const TH1F * mLumi,TH1F * mJPsi,Double_t sys[])

{
  //    return ;
  TH1F * Inclusive=(TH1F *) mInclusive->Clone("Inclusive");
  TH1F * PhotonicE=(TH1F *) mPhotonicE->Clone("PhotonicE");
 
  

  TH1F * PhotonicEReEff=(TH1F *) mPhotonicEReEff->Clone("PhotonicEReEff,");
  TH1F * purity=(TH1F *) mpurity->Clone("purity");
  
  TH1F * nsigmaE_Cut_efficiency=(TH1F *) mnsigmaE_Cut_efficiency->Clone("nsigmaE_Cut_efficiency");
  TH1F * EMC_Cut_efficiency=(TH1F *) mEMC_Cut_efficiency->Clone("EMC_Cut_efficiency");
  TH1F * Track_recon_efficiency=(TH1F *) mTrack_recon_efficiency->Clone("Track_recon_efficiency");
  TH1F * combile_Trg_efficiency=(TH1F *) mcombile_Trg_efficiency->Clone("combile_Trg_efficiency");
  TH1F * Correc=(TH1F *) mCorrec->Clone("Correc");
  TH1F * Combine_Nor=(TH1F *) mCombine_Nor->Clone("Combine_Nor");
  TH1F * Lumi=(TH1F *) mLumi->Clone("mLumi");
  TH1F * JPsi=(TH1F *) mJPsi->Clone("JPsi");

for( Int_t i=1;i<=Inclusive->GetNbinsX();i++)
    {
     PhotonicEReEff ->SetBinError(i,0);
     purity ->SetBinError(i,0);
     nsigmaE_Cut_efficiency ->SetBinError(i,0);
     EMC_Cut_efficiency ->SetBinError(i,0);
     Track_recon_efficiency ->SetBinError(i,0);
     combile_Trg_efficiency ->SetBinError(i,0);
     Correc ->SetBinError(i,0);
     Combine_Nor ->SetBinError(i,0);
     Lumi ->SetBinError(i,0);
     JPsi->SetBinError(i,0);
    }
  Inclusive->Sumw2();
  PhotonicE->Sumw2();
  Inclusive->Multiply(purity);
  PhotonicE->Divide(PhotonicEReEff);
 
  TH1F * NPE=(TH1F *)  Inclusive->Clone("NPE");
  NPE->Add(PhotonicE,-1);
  NPE->Divide(nsigmaE_Cut_efficiency);
  NPE->Divide(EMC_Cut_efficiency);
  NPE->Divide(Track_recon_efficiency);
  NPE->Divide(combile_Trg_efficiency);
  NPE->Multiply(Correc);
  NPE->Divide(Combine_Nor);
  NPE->Multiply(Lumi);
    NPE->Add(JPsi,-1);

for(Int_t i=1;i<=NPE->GetNbinsX();i++)

    {
      sys[i-1]=NPE->GetBinError(i);
      cout<<sys[i-1]<<endl;     
 
    }
  TCanvas *c1;
  NPE->Draw();

  // throw (-1);
  //  return 0.1;
} 

void Draw_Inv_yield(TH1F * HT0_NPE, TH1F * HT2_NPE,TH1F * N_npe,Double_t sys[],Double_t sts[],TH1F * before_Jpsi_Conrection,TH1F *jPsi)
{
  
  TH1F * Before_Jpsi_Conrection =( TH1F *) before_Jpsi_Conrection->Clone("Before_Jpsi_Conrection");
    TH1F * JPsi=(TH1F *) jPsi->Clone("JPsi");

  char buf[1024];
  TCanvas * c3=new TCanvas("c3","Inclusive Electron ",0,0,800,600);
 
  c3->SetLogy(1);  
 
  N_npe->SetMarkerSize(1);
  N_npe->SetMarkerStyle(20);
  N_npe->SetMarkerColor(2);
  
  
  //
  HT0_NPE->SetMarkerSize(1);
  HT0_NPE->SetMarkerStyle(20);
  HT0_NPE->SetMarkerColor(2);
   //
  HT2_NPE->SetMarkerSize(1);
  HT2_NPE->SetMarkerStyle(20);
  HT2_NPE->SetMarkerColor(3);

  N_npe->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  HT2_NPE  ->GetYaxis()->SetTitle("1/2#pi P_{t} dN^{2}/d_{y} d_{P_{t}}");
  HT2_NPE  ->GetXaxis()->SetTitle("P_{T}");
  
  N_npe->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");
  
  N_npe->SetTitle("Non-photonic Electron Cross Section");
  N_npe->GetYaxis()->CenterTitle(true);
  N_npe->GetXaxis()->CenterTitle(true);
  // N_npe->GetYaxis()->SetTitleOffset(1.35);
  N_npe->GetYaxis()->SetRangeUser(2.0e-12,5.0e-3);   
  N_npe->GetXaxis()->SetRangeUser(2+0.001,14-0.001);   
  TH2F *hh=new TH2F  ("hh","",20,1,16,2000,2.0e-12,5.0e-3);
  //N_npe->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");
  // hh->Draw();
  
  // N_npe->Print("all");
  
  TFile *file_2= new TFile("Bin_shift.root","recreate");
  TH1F * Before_Conrection= ( TH1F *)N_npe->Clone("Before_Conrection");
  Before_Conrection->Write();    

//----------------------------------------
  TF1* f1 = new TF1("f1","[0]/(pow(x,[1])+[2])",2,10);
  f1->SetParameter(0,10);
  f1->SetParameter(1,8);
  f1->SetParameter(2,0.7); 
  N_npe->Fit("f1","R0","" ,2,10);
  N_npe->Fit("f1","R0","" ,2,10);

  //c3->SaveAs("Fit.gif");


  N_npe->SetMarkerSize(1);
  N_npe->SetMarkerStyle(20);
  N_npe->SetMarkerColor(2);
 // return;
  const int nbins = Nbins_HT+1;
  // double xbins[nbins] = {2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,16.5};


  
     TH1D* hCorrection = get_bin_shift((unsigned long int)1e9,f1,nbins,Pt_bin_HT);
     N_npe->Divide(hCorrection);
  /**********/
   // binShift 
   // N_npe->Print("all");
   
 // return;
 // TH1F * After_Conrection= ( TH1F *)N_npe->Clone("After_Conrection");
 // After_Conrection->Write();
  file_2->Close();    
  //---------------------------------------
  
  TH1F * N_npe_sys=(TH1F *) N_npe->Clone("N_npe_sys");
  TH1F * N_npe_sts=(TH1F *) N_npe->Clone("N_npe_sts");

  
  fstream  outfile("NPE_cross_section.dat",ios::out|ios::trunc);  
  
  for(Int_t i=1;i<=N_npe->GetNbinsX();i++)
    {
      // N_npe->SetBinError(i,0);
      cout<<N_npe_sys->GetBinContent(i) <<endl;
      cout<<sys[i-1]<<endl;
      cout<<sts[i-1]<<endl;
      // N_npe_sys->SetBinError(i,0);
      // N_npe_sts->SetBinError(i,0);
      
      N_npe_sys->SetBinError(i,sys[i-1]);
      N_npe_sts->SetBinError(i,sts[i-1]);
      outfile<<N_npe->GetBinCenter(i)<<" "<<N_npe->GetBinContent(i)<<"  "<<
	sts[i-1]<<"  "<<sys[i-1]<<endl;      
    }
  
  TGraphErrors* data12_stat_err = new TGraphErrors("NPE_cross_section.dat","%lg %lg %lg");
  TGraphErrors* data12_sys_err = new TGraphErrors("NPE_cross_section.dat","%lg %lg %*lg %lg");

  data12_stat_err->SetLineColor(2);
  data12_sys_err->SetLineColor(4);


  data12_stat_err->SetMarkerStyle(20);
  data12_stat_err->SetMarkerColor(2);

  N_npe->GetXaxis()->SetTitleOffset(1);
  N_npe ->GetXaxis()->SetTitleOffset(1);
  N_npe->GetXaxis()->SetTitleSize(0.05);
  N_npe->GetYaxis()->SetTitleSize(0.05);  

  N_npe-> SetAxisRange(2.5+0.001,14-0.001,"x");   


  N_npe->Draw("PE1X0");
  data12_stat_err->Draw("SAMEP");
  data12_sys_err->Draw("P[]:same");
  
  
  // HT0_NPE->Draw("samePE1");
  // HT2_NPE->Draw("samePE1");
  
  

  TFile *Run_12=new TFile("Run_12.root","recreate");

  HT0_NPE->Write();
  HT2_NPE->Write();
  
  
  Run_12->Close();
  TGraphErrors* data09_stat_err = new TGraphErrors("NPERun09_corss_section.dat","%lg %lg %lg");
  TGraphErrors* data09_sys_err = new TGraphErrors("NPERun09_corss_section.dat","%lg %lg %*lg %lg");
  
  data09_stat_err->SetMarkerStyle(20);
  
  data09_stat_err->SetMarkerColor(4);
  data09_stat_err->SetLineColor(4);
  data09_sys_err->SetLineColor(4);
  
  data09_stat_err->Draw("SAMEP");
  data09_sys_err->Draw("P[]:same");
  
  
  Before_Jpsi_Conrection->SetMarkerStyle(20);
  Before_Jpsi_Conrection->SetMarkerColor(7);
  JPsi->SetMarkerStyle(20);
  JPsi->SetMarkerColor(5);

  // Before_Jpsi_Conrection->SetMarkerSize(0.8);
 
  // JPsi->SetMarkerSize(0.8);


  // Before_Jpsi_Conrection->Draw("samePE1X0");
  // JPsi->Draw("samePE1X0");



  
  HT2_NPE ->GetYaxis()->SetRangeUser(1.0e-13,1.0e-6); 
  // N_npe->Draw("PE1X0");
  N_npe_sys->SetLineColor(4);
  // N_npe_sys->SetFillColor(kRed);
  //      N_npe_sts->SetFillColor(kGreen);
  // N_npe_sys->Draw("samePEX0");
  // N_npe_sts->Draw("samePE1X0");  
  c3->SetLogy(1);
  // HT2_NPE->Draw();
  // HT0_NPE->Draw("same");
  // return ;

  
  
  TFile * file=new TFile("fonll200gev.root","read");
  TGraphErrors *gru=(TGraphErrors *) file->Get("gFONLLu");
  TGraphErrors *grc=(TGraphErrors *) file->Get("gFONLLc");
  TGraphErrors *grl=(TGraphErrors *) file->Get("gFONLLl");
  TGraphErrors *gRu=(TGraphErrors *) file->Get("gR_FONLLu");
  TGraphErrors *gRc=(TGraphErrors *) file->Get("gR_FONLLc");
  TGraphErrors *gRl=(TGraphErrors *) file->Get("gR_FONLLl");
  
  gru->Draw("same");
  grc->Draw("same");
  grl->Draw("same");

    //---------------------
  //  TGraphErrors* gr_stat_err = new TGraphErrors("star_published.dat","%lg %lg %lg");
  // TGraphErrors* gr_sys_err = new TGraphErrors("star_published.dat","%lg %lg %*lg %lg");
  TGraphErrors* gr_stat_err = new TGraphErrors("NPE_Run08_subtract_Jpsi.dat","%lg %lg %lg");
  TGraphErrors* gr_sys_err = new TGraphErrors("NPE_Run08_subtract_Jpsi.dat","%lg %lg %*lg %lg");
  
  gr_stat_err->SetMarkerStyle(20);
  
  // gr_stat_err->Draw();
  // throw (-1);

  gr_stat_err->Draw("SAMEP");
  gr_sys_err->Draw("P[]:same");
  
  //  HT2_NPE->Draw("samePE1");
  // HT0_NPE->Draw("samePE1");
  
  TLegend *legend = new TLegend(0.5,0.6,0.9,0.85);
  //legend->AddEntry( N_npe," NPE Run12 After Jpsi Correction","lpe");
  legend->AddEntry( N_npe," NPE Run12 ","lpe");
  legend->AddEntry( data09_stat_err," NPE Run09 ","lpe");
  legend->AddEntry(gr_stat_err ," NPE Run08^{[2]}","lpe");
  legend->AddEntry(  grc," FONLL","lpe");
  // legend->AddEntry(  Before_Jpsi_Conrection ," NPE Run12 Before Jpsi Correction","lpe");
  //      legend->AddEntry(  JPsi," Jpsi","lpe");
  
  legend->Draw("same");
  
  // legend->AddEntry( HT0_NPE," HT0 NPE Inv. Yield ","lpe");                
  //  legend->AddEntry( HT2_NPE,"HT2 NPE Inv. Yield","lpe");  
  // legend->SetTextFont(50);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);    
  legend->Draw("same");   
  
  //  TFile *File_jpsi=new TFile("Jpsi.root","read");
  
  // TH1F * JPsi=(TH1F *) File_jpsi->Get("JPsi");
  // JPsi->SetMarkerStyle(20);
  // JPsi->SetMarkerColor(3);
  // JPsi->Draw("samePE1");
  
  TH1F * FOll=new TH1F("FOll","",Nbins_HT,Pt_bin_HT);
  TH1F * Data_stat=new TH1F("Data_stat","",Nbins_HT,Pt_bin_HT);
  TH1F * Data_sys=new TH1F("Data_sys","",Nbins_HT,Pt_bin_HT);
  Data_sys->Sumw2();
  FOll->Sumw2();
  
  for(Int_t i=0;i<Nbins_HT;i++)
    {
      FOll->SetBinContent(i+1,grc->Eval((Pt_bin_HT[i]+Pt_bin_HT[i+1])/2.0));
      FOll->SetBinError(i+1,0);
      Data_stat->SetBinContent(i+1,N_npe_sts->GetBinContent(i+1));
      Data_stat->SetBinError(i+1,N_npe_sts->GetBinError(i+1));
      
      Data_sys->SetBinContent(i+1,N_npe_sys->GetBinContent(i+1));
      Data_sys->SetBinError(i+1,N_npe_sys->GetBinError(i+1));
      

      cout<< N_npe->GetBinContent(i+1)/FOll->GetBinContent(i+1)<<endl;


    }
  
  Data_stat->Divide(FOll);
  Data_sys->Divide(FOll);
  fstream  outfile_Ratio12("NPE_Ratio12.dat",ios::out|ios::trunc); 
  for(Int_t i=0;i<Nbins_HT;i++)
    {
      //outfile<<N_npe->GetBinCenter(i)<<" "<<N_npe->GetBinContent(i)<<"  "<<
      //	sts[i-1]<<"  "<<sys[i-1]<<endl;   
      
      outfile_Ratio12<<Data_sys->GetBinCenter(i+1)<<" "<<Data_sys->GetBinContent(i+1)<<"  "<<Data_stat->GetBinError(i+1)<<"  "<<Data_sys->GetBinError(i+1)<<endl;
      
    }
  
  

  c3->SaveAs(TString("HHT0_HT2_N_npe_yield.gif"));
  TCanvas * c5=new TCanvas("c5","",800,600);
  setpad(c5,0.12,0.03,0.12,0.12);
  //        N_npe->Draw("PE1");
  // FOll->Draw();
  //	TH1F * Data=(TH1F *) N_npe->Clone("Data");
    //	Data->Sumw2();
    // gStyle->SetTitleSize(0.05,"XY");
     //     gStyle->SetTitleFontSize(0.1);
     Data_stat->SetMarkerSize(1);
     Data_stat->SetMarkerStyle(20);
     Data_stat->SetMarkerColor(2);
     Data_stat->GetYaxis()->SetRangeUser(0,5);   
     Data_stat->GetXaxis()->CenterTitle(true);
     Data_stat->GetYaxis()->CenterTitle(true);
     Data_stat->GetXaxis()->SetTitleOffset(1);
     Data_stat->GetYaxis()->SetTitleOffset(1);
     //   Data->Setfont size
     //    Data->SetLabelSize(0.06);
    // Data->GetXaxis()->SetTitleSize(0.05);
    // Data->GetYaxis()->SetTitleSize(0.05);
     Data_stat->SetTitle(" Data and Fonll Ratio");
    // Data->SetTitleSize(0.3);



    Data_stat->GetYaxis()->SetTitle("DATA/FONLL");
    Data_stat->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	
 TGraphErrors* Ratio12_sys_err = new TGraphErrors("NPE_Ratio12.dat","%lg %lg %*lg %lg");
    Data_stat->SetLineColor(2);
    Ratio12_sys_err->SetLineColor(2);

    Data_stat->Draw("PE1X0");
    Ratio12_sys_err->Draw("P[]:same");
    
    TGraphErrors* Ratio08_stat_err = new TGraphErrors("Ratio_08.dat","%lg %lg %lg");
    TGraphErrors* Ratio08_sys_err = new TGraphErrors("Ratio_08.dat","%lg %lg %*lg %lg");
    
    Ratio08_stat_err->SetMarkerStyle(20);
    Ratio08_stat_err->Draw("samePE1");


    // Ratio08_sys_err->SetLineWidth(2);
    Ratio08_sys_err->Draw("P[]:same");

    TGraphErrors* Ratio09_stat_err = new TGraphErrors("Ratio09.dat","%lg %lg %lg");
    TGraphErrors* Ratio09_sys_err = new TGraphErrors("Ratio09.dat","%lg %lg %*lg %lg");

    Ratio09_stat_err->SetMarkerStyle(20);
    Ratio09_sys_err->SetLineColor(4);

    
    Ratio09_stat_err->SetMarkerColor(4);
    Ratio09_sys_err->SetMarkerColor(4);

    Ratio09_stat_err->Draw("samePE1");
    Ratio09_sys_err->Draw("P[]:same");
    
    TLegend *legend2 = new TLegend(0.15,0.65,0.55,0.85);
    legend2->AddEntry( Data_stat,"DATA/FONLL Run12","lpe");
    legend2->AddEntry(Ratio09_stat_err,"DATA/FONLL Run09","lpe");
    legend2->AddEntry(Ratio08_stat_err,"DATA/FONLL Run08^{[2]}","lpe");

    
    legend2->AddEntry( gRc," FONLL","lpe");
    

     
     // legend->AddEntry(gr_stat_err ," NPE 2008","lpe");
    //     legend2->SetTextFont(50);
     legend2->SetTextSize(0.04);
     legend2->SetBorderSize(0);
     legend2->SetFillStyle(0);    
     legend2->Draw("same");
     gRu->Draw("sameE");
     gRc->Draw("sameE");
     gRl->Draw("sameE");

   
	// c5->SetLogy(1);
    c5->SaveAs("NPE_Foll.gif");
    TCanvas * c6=new TCanvas("c6","",800,600);

    // c6->Divide(2,1);
    // c6->cd(1);
  



  
     TLegend *legend1 = new TLegend(0.45,0.7,0.85,0.6);
    
                  
      legend1->AddEntry( HT0_NPE," HT0 NPE Inv. Yield   ","lpe");                
     legend1->AddEntry( HT2_NPE,"HT2 NPE Inv. Yield ","lpe");  
     legend1->SetTextFont(50);
     legend1->SetTextSize(0.03);
     legend1->SetBorderSize(0);
     legend1->SetFillStyle(0);    
     legend1->Draw("same"); 
     c6->SetLogy(1);
     c6->SaveAs("NPEE.gif");    
 // throw (-1);
     TCanvas * c7=new TCanvas("c7","",800,600);
    
     TH1F *  Ration=(TH1F * ) HT2_NPE->Clone("Ration");
     Ration->GetYaxis()->SetRangeUser(0,2);
     Ration->SetTitle("Ratio of HT2 and HT0");
 
     Ration->Divide(HT0_NPE);
    
     Ration->Draw();
     Ration->GetYaxis()->SetTitle("HT2 NPE Inv. Yield over HT0 Inv.yield");

   c7->SaveAs("NPE_ratio.gif");

}
//void Draw_Efficiency(TH1F * ,TH1F *,TH1F *)

void  Draw_Efficiency(TH1F *nsigmaE_Cut_efficiency,TH1F *EMC_Cut_efficiency,TH1F * Track_recon_efficiency,TH1F *combile_Trg_efficiency)
{
  gStyle->SetOptStat(0);
  
  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->cd();
  nsigmaE_Cut_efficiency->SetMarkerSize(1);
  nsigmaE_Cut_efficiency->SetMarkerStyle(20);
  nsigmaE_Cut_efficiency ->SetMarkerColor(1);

  EMC_Cut_efficiency->SetMarkerSize(1);
  EMC_Cut_efficiency->SetMarkerStyle(20);
  EMC_Cut_efficiency ->SetMarkerColor(2);
  EMC_Cut_efficiency ->SetLineColor(2);


  Track_recon_efficiency ->SetMarkerSize(1);
  Track_recon_efficiency->SetMarkerStyle(20);
  Track_recon_efficiency ->SetMarkerColor(4);
  Track_recon_efficiency ->SetLineColor(4);

  combile_Trg_efficiency ->SetMarkerSize(1);
  combile_Trg_efficiency->SetMarkerStyle(20);
  combile_Trg_efficiency ->SetMarkerColor(7);
  combile_Trg_efficiency ->SetLineColor(7);



  

 TFile * infile_purity=new TFile("Purity.root","read");
 
 TH1F * mh1purity = (TH1F *) infile_purity->Get("purity");

 TFile * infile_Ph_eff =new TFile("Photonic_ReconStructionEfficiency.root","read");

 TH1F * mh1ReConEffienc = (TH1F *) infile_Ph_eff->Get("ReConEffienc");

 fstream  Error_Table("Error_table.txt",ios::trunc|ios::out);
    for(Int_t i=1;i<=Nbins_HT;i++)
    {
      //      outfile_Ratio12_sys<<NPE_sys->GetBinCenter(i)<<" "<<NPE_sys->GetBinContent(i)<<"  "<<Ration_NPE_PE->GetBinError(i)<<"  "<<NPE_sys->GetBinError(i)<<endl;

      Error_Table<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(4)<<mh1purity->GetBinCenter(i)<<"  "
	//		 <<mh1purity->GetBinContent(i)<<"  "<<mh1purity->GetBinError(i)<<"  "<<mh1purity->GetBinError(i)/mh1purity->GetBinContent(i)<<"  " 
		 <<mh1purity->GetBinError(i)/mh1purity->GetBinContent(i)<<"  " 
		 
    
			 <<mh1ReConEffienc->GetBinError(i)/mh1ReConEffienc->GetBinContent(i)<<" "

	//		 <<mh1ReConEffienc->GetBinError(i)/mh1ReConEffienc->GetBinContent(i)<<"  "
		
		 <<EMC_Cut_efficiency->GetBinError(i)/ EMC_Cut_efficiency->GetBinContent(i)<<"  "
		 <<Track_recon_efficiency->GetBinError(i)/Track_recon_efficiency->GetBinContent(i)<<"  "
		 <<combile_Trg_efficiency->GetBinError(i)/combile_Trg_efficiency->GetBinContent(i)<<"  "
		 <<nsigmaE_Cut_efficiency->GetBinError(i)/nsigmaE_Cut_efficiency->GetBinContent(i)<<"  "<<endl;


      
	
	}

    
 mh1purity ->SetMarkerSize(1);
 mh1purity->SetMarkerStyle(20);
 mh1purity ->SetMarkerColor(6);

 mh1ReConEffienc ->SetMarkerSize(1);
 mh1ReConEffienc->SetMarkerStyle(20);
 mh1ReConEffienc ->SetMarkerColor(5);


 TH1F *eID=(TH1F *) nsigmaE_Cut_efficiency->Clone("eID");
 eID->Sumw2();
 eID->Multiply(EMC_Cut_efficiency);
 eID->Multiply(Track_recon_efficiency);
 // eID->Multiply(combile_Trg_efficiency);
 
 eID ->SetLineColor(6);
 eID->SetMarkerStyle(20);
 eID ->SetMarkerColor(6);

 nsigmaE_Cut_efficiency->GetYaxis()->SetRangeUser(0,1.8);
 nsigmaE_Cut_efficiency->GetXaxis()->SetRangeUser(2.5,14);
 nsigmaE_Cut_efficiency->GetYaxis()->SetTitle("Electron selection efficiency");
 nsigmaE_Cut_efficiency->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 nsigmaE_Cut_efficiency->GetYaxis()->CenterTitle(true);
 nsigmaE_Cut_efficiency->GetXaxis()->CenterTitle(true);
 nsigmaE_Cut_efficiency->SetTitle(" Electron selection efficiency ");

 nsigmaE_Cut_efficiency->GetXaxis()->SetTitleOffset(1);
 nsigmaE_Cut_efficiency ->GetXaxis()->SetTitleOffset(1);
 nsigmaE_Cut_efficiency->GetXaxis()->SetTitleSize(0.05);
 nsigmaE_Cut_efficiency->GetYaxis()->SetTitleSize(0.05);  


   nsigmaE_Cut_efficiency->Draw("PE1");
   EMC_Cut_efficiency->Draw("samePE1");
   Track_recon_efficiency->Draw("samePE1");
   combile_Trg_efficiency->Draw("samePE1");

   // mh1purity->Draw("samePE1");
   // mh1ReConEffienc->Draw("samePE1");  
   eID->Draw("samePE1");
  
  TLegend *legend = new TLegend(0.45,0.6,0.75,0.85);
  legend->AddEntry(nsigmaE_Cut_efficiency,"#epsilon _{dEdx}   ","lpe");                 
  legend->AddEntry(EMC_Cut_efficiency,"#epsilon _{EMC}","lpe"); 
  legend->AddEntry(Track_recon_efficiency,"#epsilon _{Trk} ","lpe");   
   legend->AddEntry(combile_Trg_efficiency,"#epsilon _{Trg}","lpe");
   // legend->AddEntry(mh1purity,"#epsilon _{purity}  From Data","lpe");
   //  legend->AddEntry(mh1ReConEffienc,"#epsilon _{Pho_ReC}  From MC","lpe");
   
legend->AddEntry(eID,"#epsilon _{Total}","lpe");                
  
  
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(62);
  legend->Draw("same");     
 
  c1->SaveAs("cut_efficiency.pdf");

  TFile *file_efficiency=new TFile("Efficiency_HT.root","RECREATE");

  nsigmaE_Cut_efficiency->Write();
  EMC_Cut_efficiency->Write();
  Track_recon_efficiency->Write();
  combile_Trg_efficiency->Write();
  mh1purity->Write();
  mh1ReConEffienc->Write();
  eID->Write();

  file_efficiency->Close();


}
void Draw_Inclusive_electron(const TH1F *mh1Inclusive_ele_h,const TH1F * mh1purity_h)
{   
  gStyle->SetOptStat(00000);
  char buf[1024];
  TH1F * mh1Inclusive_ele=(TH1F *)mh1Inclusive_ele_h->Clone("mh1Inclusive_ele");
  TH1F * mh1purity=(TH1F *)mh1purity_h->Clone("mh1purity"); 
  TCanvas * c2=new TCanvas("c2","Inclusive Electron ",0,0,800,600);
  
 
  
  
  setpad(c2,0.12,0.02,0.02,0.12);
  c2->cd();
  mh1Inclusive_ele->SetMarkerSize(1);
  mh1Inclusive_ele->SetMarkerStyle(20);
  mh1Inclusive_ele->SetMarkerColor(1);
  mh1Inclusive_ele->GetXaxis()->SetTitle("Pt");
  mh1Inclusive_ele->GetYaxis()->SetTitle("dN/dPt");
  mh1Inclusive_ele->Draw("PE1");
  sprintf(buf,"Inclusive electron");
  drawLatex(0.55,0.85,buf,70,0.05,1);
  c2->SetLogy(1);
  c2->SaveAs(TString("Inclusive_electron.gif"));
  
  TCanvas * c3=new TCanvas("c3","Inclusive pure Electron ",0,0,800,600);
  
  setpad(c3,0.12,0.02,0.02,0.12);
  c3->cd();
  TH1F * hh_pur=(TH1F *) mh1Inclusive_ele->Clone("hh_pur");
  hh_pur-> Multiply(hh_pur,mh1purity);
  hh_pur->SetMarkerSize(1);
  hh_pur->SetMarkerStyle(20);  
  hh_pur->SetMarkerColor(2);
  hh_pur->GetXaxis()->SetTitle("Pt");
  hh_pur->GetYaxis()->SetTitle("dN/dPt");
  hh_pur->Draw("PE1");
  
  
  mh1Inclusive_ele->Draw("SAMEPE1");
  sprintf(buf,"Inclusive Pure electron");
  drawLatex(0.55,0.85,buf,70,0.05,1);
  
  TLegend *legend = new TLegend(0.5,0.9,0.8,0.6);
  legend->AddEntry(mh1Inclusive_ele," Inclusive  Electron ","lpe");                 
  legend->AddEntry(hh_pur,"After correction of the purity ","lpe");                
  
  legend->SetTextFont(50);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw("same");     
  c3->SetLogy(1);
  c3->SaveAs(TString("Inclusive_pure_electron.gif"));
  
  
}
void Draw_Photonic_electron(const TH1F * mh1Photonic_ele_h,const TH1F *Ph_eff_h)
{
  gStyle->SetOptStat(00000);
  char buf[1024];
  TH1F * mh1Photonic_ele=(TH1F *)mh1Photonic_ele_h->Clone("mh1Inclusive_ele");
  TH1F * mh1Photonic_eff=(TH1F *)Ph_eff_h->Clone("mh1Photonic_eff");
  TCanvas * c2=new TCanvas("c2","Photonic Electron ",0,0,800,600);
  
  setpad(c2,0.12,0.02,0.02,0.12);
  c2->cd();
  mh1Photonic_ele->SetMarkerSize(1);
  mh1Photonic_ele->SetMarkerStyle(20);
  mh1Photonic_ele->SetMarkerColor(1);
  mh1Photonic_ele->GetXaxis()->SetTitle("Pt");
  mh1Photonic_ele->GetYaxis()->SetTitle("dN/dPt");
  
  mh1Photonic_ele->Draw("PE1");
  sprintf(buf,"Photonic electron");
  drawLatex(0.55,0.85,buf,70,0.05,1);
  c2->SetLogy(1);
  c2->SaveAs(TString("Photonic_electron.gif"));
  
  TCanvas * c3=new TCanvas("c3","Photonic Recon Electron ",0,0,800,600);
  
  setpad(c3,0.12,0.02,0.02,0.12);
  c3->cd();
  TH1F * hh_pur=(TH1F *) mh1Photonic_ele->Clone("hh_pur");
  hh_pur-> Divide(hh_pur,mh1Photonic_eff);
  hh_pur->SetMarkerSize(1);
  hh_pur->SetMarkerStyle(20);
  hh_pur->SetMarkerColor(2);
  hh_pur->GetXaxis()->SetTitle("Pt");
  hh_pur->GetYaxis()->SetTitle("dN/dPt");
  
  hh_pur->Draw("PE1");
  mh1Photonic_ele->Draw("SAMEPE1");
  sprintf(buf,"Photonic Recon electron");
  drawLatex(0.55,0.85,buf,70,0.05,1);
  
  TLegend *legend = new TLegend(0.4,0.9,0.8,0.6);
  legend->AddEntry(mh1Photonic_ele,"  Photonic electron  ","lpe");
  legend->AddEntry(hh_pur,"After correction reconstruction efficiency ","lpe");
  
  legend->SetTextFont(50);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw("same");                                                                                                             
  c3->SetLogy(1);
  c3->SaveAs(TString("Photonic_Recon_electron.gif"));
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
