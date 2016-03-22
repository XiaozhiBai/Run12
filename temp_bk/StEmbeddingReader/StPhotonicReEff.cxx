/*
 *  
 *   3/17/2014/ xiaozhi bai
 *
 */




#include <iostream>
#include "StPhotonicReEff.h"
#include "StMcAccpTracks.h"
#include "StMcAccpEventCount.h"
#include "StPheRecoPair.h"
#include "StPheRecoTrack.h"
#include "TFile.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "cmath"
using namespace std;

ClassImp(StPhotonicReEff)
  
StPhotonicReEff::StPhotonicReEff(const char* outName)
{
  TH1F::SetDefaultSumw2();  
  mOutputFile=new TFile(outName,"RECREATE");
  //  weightFile=new TFile("pi0_gamma_fit.root","READ");
  weightFile=new TFile("/global/u1/x/xiao00/disk02/Weight.root","READ");
  
  
  Gamma_spectra=(TH1D *) weightFile->Get("Gamma_Total");
  pi0_spectra= (TH1D *) weightFile->Get("Pi");
  eta_spectra=(TH1D *) weightFile->Get("Eta");
  
   // wightFunction=(TF1 *) weightFile->Get("fun");
   // TError=(TGraphErrors *) weightFile->Get("gra");

}
 StPhotonicReEff::~StPhotonicReEff()
{

  //..
  //
}
void StPhotonicReEff::bookHistogram()
{ 
  cout<<"book histogram"<<endl;
  //track
  mSingle_Track=new TH1F("mSingle_Track","mSingle_Track",2000,0,20);
  //Pair
  mPair=new TH1F("mPair","mPir",2000,0,20);
  mInv_massVspt=new TH1F("mInv_massVspt","",2000,0,20);

  mh2Pair=new TH2F("mh2Pair","mPir",2000,0,20,2000,1e-12,2);
  mh2Single_Track=new TH2F("mh2Single_Track","mSingle_Track",2000,0,20,2000,1e-12,2);
  

}
void StPhotonicReEff::read(TString fileName)
{ 
  //  TH1F::SetDefaultSumw2();
  TFile * infile=new TFile(fileName.Data());
  // if( TNtuple* ntTracks=(TNtuple* ) infile->Get("nt_sngl"));
  // cout<<"   Get File"<<endl;
  TNtuple* ntTracks=(TNtuple* ) infile->Get("nt_sngl");
  TNtuple * ntPair=(TNtuple* )infile->Get("nt_pair");
    
  StPheRecoTrack * tracksMC=new StPheRecoTrack; 
  tracksMC->Init(ntTracks);
  StPheRecoPair *Pair=new StPheRecoPair;
  Pair->Init(ntPair);
  
  for(Int_t iTrk=0;iTrk<tracksMC->GetEntries();iTrk++)
    {
 
      tracksMC->GetEntry(iTrk);
      if(fabs(tracksMC->vz)<35 
	 && tracksMC->dca<1.5 
	 && tracksMC->nfit>20 
	 && tracksMC->nDedxPts>15 
	 && tracksMC->nfit/(Float_t)tracksMC->nmax > 0.52 
	 && sqrt(tracksMC->stpcx*tracksMC->stpcx+tracksMC->stpcy*tracksMC->stpcy)<73 
	 &&0.2<tracksMC->rpt
	  &&10<tracksMC->ncom
	 && abs(tracksMC->rpt-tracksMC->pt)/tracksMC->pt<0.15
	 && abs(tracksMC->eta)<0.7)
	{
	  Float_t weight= Weight(tracksMC->ppt,tracksMC->pgeantId);
	  mSingle_Track->Fill(tracksMC->rpt,weight);
	  mh2Single_Track->Fill(tracksMC->rpt,weight,weight);    
	  mInv_massVspt->Fill(tracksMC->rpt,weight);//Pair->p1rpt,Pair->massDCA,weight);
	}
    }
  for(Int_t iPair=0;iPair<Pair->GetEntries();iPair++)
    {
      Pair->GetEntry(iPair);
      if(fabs(Pair->vz)<35
	 && fabs(Pair->p1eta)<0.7 
	 &&0.3<Pair->p2rpt 
	 && Pair->p1nfit>20 
	 && Pair->p1nfit/(Float_t)Pair->p1nmax>0.52 
	 && Pair->p1ndedx>15 
	 && Pair->p1ncom>10
	 && Pair->p2ncom>10 

	 && 0.2<Pair->p1rpt
	 && abs(Pair->p1rpt-Pair->p1pt)/Pair->p1pt<0.15
	 && sqrt(Pair->p1stpcx*Pair->p1stpcx+Pair->p1stpcy*Pair->p1stpcy)<73
	 && Pair->p1dca<1.5 
	 && Pair->p2nfit>20 && (Float_t)Pair->p2nfit/(Float_t)Pair->p2nmax>0.52
	 &&Pair->massDCA<0.24
	 && Pair->pairDCA<1)	
	{ 
	  Float_t weight= Weight(Pair->p1ppt,Pair->p1pgeantId);
	  mPair->Fill(Pair->p1rpt,weight);
	  mh2Pair->Fill(Pair->p1rpt,weight,weight);
	

	}
    }
}
void StPhotonicReEff::WriteHistogram()
{ 
  mOutputFile->cd();
  mSingle_Track->Write();
  mPair->Write();
  mInv_massVspt->Write();
  mh2Pair->Write();
  mh2Single_Track->Write();

}
Double_t StPhotonicReEff::Weight(Float_t pt,Float_t parent_geantID)
{
  Double_t weight=0; 
  
  if(parent_geantID==1)
    weight=Gamma_spectra->GetBinContent(Gamma_spectra->FindBin(pt));

  if(parent_geantID==10007)
    weight=pi0_spectra->GetBinContent(pi0_spectra->FindBin(pt));
  
  if(parent_geantID==10017)
    weight=eta_spectra->GetBinContent(eta_spectra->FindBin(pt));   

  
 
  return weight;
}
