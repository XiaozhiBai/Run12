#include <iostream>
#include "StSingleElectron_EmbeddingQA.h"
#include "StSingle_Electron_tracks.h"
#include "StSingle_electron_EventCount.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TNtuple.h"
#include "cmath"
using namespace std;

ClassImp(StSingleElectron_EmbeddingQA)
  
StSingleElectron_EmbeddingQA::StSingleElectron_EmbeddingQA(const char* outName)
{
  
  mOutputFile=new TFile(outName,"RECREATE");

}
 StSingleElectron_EmbeddingQA::~StSingleElectron_EmbeddingQA()
{
  //..
  //
}
void StSingleElectron_EmbeddingQA::bookHistogram()
{ 
  //  Float_t Pt_xx[9]={2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,16.5};
  // Float_t Pt_xx[]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,2,2.5,3,3.5,4};
  mNTrack=new TH1F("mNTrack","mNTrack",1000,0,20);
  // mFitPos_pt=new TH2F("mFitPos_pt","mFitPos_pt",200,0,20,100,0,1);
  //  mNSMDPhi=new TH2F("mNSMDPhi","mNSMDPhi",200,0,20,50,0,5);
  //cut
  DSMAdc_over_adc0=new TH1F("DSMAdc_over_adc0","DSMAdc_over_adc0",100,0,0.3);

  low_DSMAdc=new TH1F("low_DSMAdc","low_DSMAdc",200,0,20);
  DSM_adc0=new TH2F("DSM_adc0","DSM_adc0",2000,0,2000,100,0,100);
  DSMadc_pt=new TH2F("DSMadc_pt","DSMadc",100,0,100,200,0,20);
  DSMadc_pt_cut=new TH2F("DSMadc_pt_cut","DSMadc",100,0,100,200,0,20);
  adc0_pt=new TH2F("adc0_pt","adc0_pt",2000,0,2000,200,0,20);
  adc0_pt_cut=new TH2F("adc0_pt_cut","adc0_pt_cut",2000,0,2000,200,0,20);


  TowerID_cut=new TH1F("TowerID_cut","",4800,0.5,4800.5);
  TowerID=new TH1F("TowerID","",4800,0.5,4800.5);
  TowerID_all=new TH1F("TowerID_all","",4800,0.5,4800.5);

  TowerID_cut_noHOt=new TH1F("TowerID_cut_noHOt","",4800,0.5,4800.5);
  TowerID_noHOt=new TH1F("TowerID_noHOt","",4800,0.5,4800.5);
  TowerID_all_noHOt=new TH1F("TowerID_all_noHOt","",4800,0.5,4800.5);


  mPhi_pt_cut=new TH2F ("mPhi_pt_cut","mcPhi_cut",200,0,20,400,-4,4);
  mEta_pt_cut=new TH2F("mEta_pt_cut","mcEta_cut",200,0,20,300,-1.5,1.5);
  mHitFit_pt_cut =new TH2F("mHitFit_pt_cut","mHitFit_pt_cut",200,0,20,60,0,60);
  mNSMDEta_pt_cut=new TH2F("mNSMDEta_pt_cut","mNSMDEta_pt_cut",200,0,20,10,0,10);
  mHitsDedx_pt_cut=new TH2F("mHitsDedx_pt_cut","mHitsDedx_pt_cut",200,0,20,60,0,60);
  mNTrack_cut=new TH1F("mNTrack_cut","mNTrack_cut",1000,0,20);
  mNTrack_cut_25=new TH1F("mNTrack_cut_25","mNTrack_cut_25",1000,0,20);
  mNTrack_cut_20=new TH1F("mNTrack_cut_20","mNTrack_cut_20",1000,0,20);
  mFitPos_pt_cut=new TH2F("mFitPos_pt_cut","mFitPos_pt_cut",200,0,20,100,0,1);
  mNSMDPhi_pt_cut=new TH2F("mNSMDPhi_pt_cut","mNSMDPhi_cut",200,0,20,10,0,10);
  mgDca_pt_cut=new TH2F("mgDca_pt_cut","mgDca",200,0,20,50,0,5);
 
  mNsigE_pt_cut=new TH2F("mNsigE_pt_cut","mNsigE_pt_cut",200,0,20,50,-5,5);
  mDedx_pt_cut=new TH2F("mDedx_pt_cut","mDedx_pt_cut",200,0,20,60,1,6); 
  mPoe_pt_cut= new TH2F ("mPoe_pt_cut","",200,0,20,100,0,2.5);
}
void StSingleElectron_EmbeddingQA::read(TString fileName)
{
  //  cout<<"!!!!!!!!!!!!!!!!!read"<<endl;
  TFile * infile=new TFile(fileName.Data());

  TNtuple* ntEventCount=(TNtuple*) infile->Get("eventCount");
  TNtuple* ntTracks=(TNtuple* ) infile->Get("tracksMC");
 

  StSingle_electron_EventCount* eventCount=new StSingle_electron_EventCount;
  StSingle_Electron_tracks * tracksMC=new StSingle_Electron_tracks; 
  tracksMC->Init(ntTracks);
  eventCount->Init(ntEventCount);
  
  cout<<"  Entries"<<eventCount->GetEntries() <<endl;
  cout<<"  Entries"<<tracksMC->GetEntries() <<endl; 
  
  for (Int_t iEvent=0;iEvent<eventCount->GetEntries();iEvent++)//
    {
      eventCount->GetEntry(iEvent);
      
      // mDelta_Vz->Fill(abs(eventCount->vz-eventCount->mcVz));
    }

  for(Int_t iTrk=0;iTrk<tracksMC->GetEntries();iTrk++)
    {
      // cout<<"sss"<<endl;
      tracksMC->GetEntry(iTrk);
      if(abs(tracksMC->eta)<0.7 )
	{
	  mNTrack->Fill(tracksMC->pt);

	}
      
      if(tracksMC->dca<1.5 && tracksMC->nfit>20 && tracksMC->nDedxPts>15 && tracksMC->nfit/(Float_t)tracksMC->nmax > 0.52 && sqrt(tracksMC->stpcx*tracksMC->stpcx+tracksMC->stpcy*tracksMC->stpcy)<73 && 0<tracksMC->rpt && abs(tracksMC->eta)<0.7 )
	{
	  mNTrack_cut->Fill(tracksMC->pt);
	  if(tracksMC->nfit>25)
	    {	
	      mNTrack_cut_25->Fill(tracksMC->pt);
	    }
	}

                  if(tracksMC->dca<1.5 && tracksMC->nfit>20 && tracksMC->nDedxPts>15 && tracksMC->nfit/(Float_t)tracksMC->nmax > 0.52 && sqrt(tracksMC->stpcx*tracksMC->stpcx+tracksMC->stpcy*tracksMC->stpcy)<71 && 0<tracksMC->rpt && abs(tracksMC->eta)<0.7 &&tracksMC->bsmdNEta>1 && tracksMC->bsmdNPhi>1 && 0.3<tracksMC->p/tracksMC->btowE0 && tracksMC->p/tracksMC->btowE0<1.5 && abs(tracksMC->bemcDistZ)<3&& abs(tracksMC->bemcDistPhi)<0.015)
    
	{
	  TowerID_all->Fill(tracksMC->btowId);
	  if(!isHotTower(tracksMC))
	    TowerID_all_noHOt->Fill(tracksMC->btowId);
	  
	  DSMadc_pt->Fill(tracksMC->btowDsmAdc0,tracksMC->pt);
	 
	  adc0_pt->Fill(tracksMC->btowAdc0,tracksMC->pt);
	  DSMAdc_over_adc0->Fill(tracksMC->btowDsmAdc0/(Float_t)tracksMC->btowAdc0);
	  if(tracksMC->btowDsmAdc0<=11)
	  low_DSMAdc->Fill(tracksMC->btowDsmAdc0);
	  if(tracksMC->btowDsmAdc0 > tracksMC->btowAdc0*0.1)
	    {
	    TowerID->Fill(tracksMC->btowId);
	    if(!isHotTower(tracksMC))
	      TowerID_noHOt->Fill(tracksMC->btowId);
	    }
	  if(tracksMC->btowDsmAdc0 < tracksMC->btowAdc0*0.1 && !isHotTower(tracksMC) )
	    {
	      TowerID_cut->Fill(tracksMC->btowId);
	      mPhi_pt_cut->Fill(tracksMC->pt,tracksMC->rphi);
	      //	      mEta_pt_cut->Fill(tracksMC->pt,tracksMC->reta);
	  
	      DSMadc_pt_cut->Fill(tracksMC->btowDsmAdc0,tracksMC->pt);
	      DSM_adc0->Fill(tracksMC->btowAdc0,tracksMC->btowDsmAdc0);
	      adc0_pt_cut->Fill(tracksMC->btowAdc0,tracksMC->pt);
	      mHitFit_pt_cut->Fill(tracksMC->pt,tracksMC->nfit);
	      mNSMDEta_pt_cut->Fill(tracksMC->pt,tracksMC->bsmdNEta);
	      mFitPos_pt_cut->Fill(tracksMC->pt,tracksMC->nfit/(Float_t)tracksMC->nmax);
	      mHitsDedx_pt_cut->Fill(tracksMC->pt,tracksMC->nDedxPts);
	      mNSMDPhi_pt_cut->Fill(tracksMC->pt,tracksMC->bsmdNPhi);
	      mPoe_pt_cut->Fill(tracksMC->pt,tracksMC->p/tracksMC->btowE0);
	      mDedx_pt_cut->Fill(tracksMC->pt,tracksMC->dedx);
	      mgDca_pt_cut->Fill(tracksMC->pt,tracksMC->dca);
	      mNsigE_pt_cut->Fill(tracksMC->pt,tracksMC->nSigE);
	      if(!isHotTower(tracksMC))
		TowerID_cut_noHOt->Fill(tracksMC->btowId);
	    }
	}

   }
  
}
void StSingleElectron_EmbeddingQA::WriteHistogram()
{
  cout<<"  write"<<endl;
  mOutputFile->cd();
  //  mDelta_Vz->Write();
    mPhi_pt_cut->Write();
    mEta_pt_cut->Write();
    //  mHitFit_pt_cut->Write();
  //  mNSMDEta_pt->Write();
  // mFitPos_pt->Write();
  // mHitsDedx->Write();
   mNTrack->Write();
   mPoe_pt_cut->Write();
  // mNSMDPhi->Write();
  // mgDca->Write();
  // mNsigE->Write();
  
    //  mHitFit_pt_cut->Write();
    // mNSMDEta_pt_cut->Write();
    // mFitPos_pt_cut->Write();
    // mHitsDedx_pt_cut->Write();
     mNTrack_cut->Write();
     mNTrack_cut_25->Write(); 
     //  mNTrack_cut_20->Write();
    // mNSMDPhi_pt_cut->Write();
    //  mgDca_pt_cut->Write();
    //  mNsigE_pt_cut->Write();
 //
 mHitFit_pt_cut->Write();
 // mHitMax_pt_cut->Write();
 mFitPos_pt_cut->Write();
 mDedx_pt_cut->Write();
  mHitsDedx_pt_cut->Write();
 mNSMDEta_pt_cut->Write();
 mNSMDPhi_pt_cut->Write();


 mgDca_pt_cut->Write();
 mNsigE_pt_cut->Write();
 DSMadc_pt->Write();
 DSMadc_pt_cut->Write();
 adc0_pt->Write();
 adc0_pt_cut->Write();
 DSM_adc0->Write();
 DSMAdc_over_adc0->Write();
 TowerID_cut->Write();
 TowerID->Write();
 TowerID_all->Write();

 TowerID_cut_noHOt->Write();
 TowerID_noHOt->Write();
 TowerID_all_noHOt->Write();
 low_DSMAdc->Write();
 // mDedx_pt_cut->Write();
}
bool StSingleElectron_EmbeddingQA::isHotTower(StSingle_Electron_tracks * trk)
{

  Int_t  Hot_towerlist[] ={32,52,115,246,268,276,294,386,510,562,682,750,773,894,987,994,1043,1064,1143,1233,1264,1285,1307,1487,1593,1710,1733,1823,1824,1851,1946,2022,2044,2064,2110,2146,2163,2165,2203,2291,2314,2522,2530,2634,2653,2835,2864,2866,2973,3006,3062,3533,3545,3727,3862,3949,4051,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767,32,52,115,268,276,294,295,510,562,682,750,987,994,1064,1143,1233,1264,1285,1307,1487,1593,1733,1824,1851,1946,2044,2163,2203,2634,2653,2835,2864,2866,2973,3006,3693,3727,3862,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};

  for(Int_t i=0;i<sizeof(Hot_towerlist)/sizeof(Hot_towerlist[0]);i++ )
    {
      if((Hot_towerlist[i]-1)==trk->btowId)

	return kTRUE;
    }
  return kFALSE;
}
