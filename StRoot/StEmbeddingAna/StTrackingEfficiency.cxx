/******************************************************************************************
 *Wed Mar 23 15:24:09 PDT 2016                                                            *
 * This class for the tracking efficiency                                                 *
 * By xiaozhi Bai xiaozhi@uic.edu                                                         *
 ******************************************************************************************/


#include <iostream>
#include "StTrackingEfficiency.h"
#include "StSingle_Electron_tracks.h"
#include "StSingle_electron_EventCount.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TNtuple.h"
#include "cmath"
using namespace std;

ClassImp(StTrackingEfficiency)
  
StTrackingEfficiency::StTrackingEfficiency(const char* outName)
{
  

  mOutputFile=new TFile(outName,"RECREATE");
  weightFile=new TFile("/global/homes/x/xiao00/work/run12/Embedding/fonll200gev.root","READ");
  gFONLLc=(TGraphErrors *) weightFile->Get("gFONLLc");  

}
 StTrackingEfficiency::~StTrackingEfficiency()
{
  //..
  //
}
void StTrackingEfficiency::bookHistogram()
{ 
 
  mNTrack_nocuts=new TH1F("mNTrack_nocuts","mNTrack_nocuts",1000,0,20);
  mNTrack_cut=new TH1F("mNTrack_cut","mNTrack_cut",1000,0,20);
  mNTrack_cut_25=new TH1F("mNTrack_cut_25","mNTrack_cut_25",1000,0,20);

  mNTrack_nocuts->Sumw2();
  mNTrack_cut->Sumw2();
  mNTrack_cut_25->Sumw2();
}
void StTrackingEfficiency::read(TString fileName)
{
  //  cout<<"!!!!!!!!!!!!!!!!!read"<<endl;
  TFile * infile=new TFile(fileName.Data());
  TNtuple* ntEventCount=(TNtuple*) infile->Get("eventCount");
  TNtuple* ntTracks=(TNtuple* ) infile->Get("tracksMC");

  StSingle_electron_EventCount* eventCount=new StSingle_electron_EventCount;
  StSingle_Electron_tracks * tracksMC=new StSingle_Electron_tracks; 
  tracksMC->Init(ntTracks);
  eventCount->Init(ntEventCount);

  cout<<" working on Tracking efficiency!!!! "<<endl;
  for (Int_t iEvent=0;iEvent<eventCount->GetEntries();iEvent++)//
    {
      eventCount->GetEntry(iEvent);
      // mDelta_Vz->Fill(abs(eventCount->vz-eventCount->mcVz));
    }
  for(Int_t iTrk=0;iTrk<tracksMC->GetEntries();iTrk++)
    {
      if(isHotTower(tracksMC)) continue;
      tracksMC->GetEntry(iTrk);
  
      Double_t weight =Weight(tracksMC->pt,tracksMC->geantId);
      
      if(abs(tracksMC->eta)<0.7 )
	{
	  mNTrack_nocuts->Fill(tracksMC->pt,weight);
	}
      if(tracksMC->dca<1.5 && tracksMC->nfit>20 && tracksMC->nDedxPts>15 && tracksMC->nfit/(Float_t)tracksMC->nmax > 0.52 && sqrt(tracksMC->stpcx*tracksMC->stpcx+tracksMC->stpcy*tracksMC->stpcy)<73 && 0<tracksMC->rpt && abs(tracksMC->eta)<0.7 )
	{
	  mNTrack_cut->Fill(tracksMC->pt,weight);
	  if(tracksMC->nfit>25)
	    {	
	      mNTrack_cut_25->Fill(tracksMC->pt,weight);
	    }
	}
    }
}
void StTrackingEfficiency::WriteHistogram()
{
  cout<<"  write"<<endl;

  mOutputFile->cd();
  mNTrack_nocuts->Write();
  mNTrack_cut->Write();
  mNTrack_cut_25->Write(); 
}
bool StTrackingEfficiency::isHotTower(StSingle_Electron_tracks * trk)
{
  Int_t  Hot_towerlist[] ={32,52,115,246,268,276,294,386,510,562,682,750,773,894,987,994,1043,1064,1143,1233,1264,1285,1307,1487,1593,1710,1733,1823,1824,1851,1946,2022,2044,2064,2110,2146,2163,2165,2203,2291,2314,2522,2530,2634,2653,2835,2864,2866,2973,3006,3062,3533,3545,3727,3862,3949,4051,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767,32,52,115,268,276,294,295,510,562,682,750,987,994,1064,1143,1233,1264,1285,1307,1487,1593,1733,1824,1851,1946,2044,2163,2203,2634,2653,2835,2864,2866,2973,3006,3693,3727,3862,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};
  
  for(Int_t i=0;i<sizeof(Hot_towerlist)/sizeof(Hot_towerlist[0]);i++ )
    {
      if((Hot_towerlist[i]-1)==trk->btowId)
	
	return kTRUE;
    }
  return kFALSE;
}
Double_t StTrackingEfficiency::Weight(Float_t pt,Float_t geantID)
{
  Double_t weight=0;
  if(abs(geantID)==2||abs(geantID)==3)
    weight=gFONLLc->Eval(pt);
  return weight;
}
