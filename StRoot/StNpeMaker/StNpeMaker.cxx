/******************************************************************************************
 *Mon Mar  7 19:33:04 EST 2016                                                            *
 * This class for read  DemesonTree Fill all the No-phoronic electron realted histogram   *
 * By xiaozhi Bai xiaozhi@uic.edu                                                         *
 ******************************************************************************************/
// StNpemaker.cxx

#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include "StNpeMaker.h"
#include "mBadRunList.h"
#include "StCuts.h"

#include "StTRIGGERS.h"

#include "prescales.h"
#include "StDmesonMaker/StDmesonEvent.h"
#include "StDmesonMaker/StDmesonTrack.h"
#include "StDmesonMaker/StElectronPair.h"

#include "StLorentzVectorF.hh"
#include "StMessMgr.h"

#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TClonesArray.h"

ClassImp(StNpeMaker)
using namespace std;
StNpeMaker::StNpeMaker(const char* outName){
  mOutputFile = new TFile(outName, "RECREATE");
  mPrescales=prescales::Instance();

  TH1F::SetDefaultSumw2();
  TH2F::SetDefaultSumw2();
  TH3F::SetDefaultSumw2();

  mNBadRuns_MB = sizeof(mBadRuns_MB)/sizeof(int);
  mNBadRuns_HT0 = sizeof(mBadRuns_HT0)/sizeof(int);
  mNBadRuns_HT2 = sizeof(mBadRuns_HT2)/sizeof(int);
  LOG_INFO <<"Numbber of Bad run MB : HT0 : HT2   "<< mNBadRuns_MB <<" :" <<mNBadRuns_HT0<< " :"<< mNBadRuns_HT2<<endm;
}
StNpeMaker::~StNpeMaker(){
  /*  */
}
void StNpeMaker::Clear(Option_t *opt) {
  /**/
}
void StNpeMaker::read(TString fileName){
  TFile* inFile = new TFile(fileName.Data(),"READ");
  TTree* tree = (TTree*)inFile->Get("T");
  tree->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  tree->SetBranchAddress("dEvent", &mEvent);

  if(!mEvent){
    LOG_WARN << " No DmesonEvent! Skip! " << endm;
    return ;  
  } 
  TClonesArray*   aTracks = 0;
  TClonesArray* aPairs=0;
  LOG_INFO <<"Read File "<<fileName<< " Total number of events"<< tree->GetEntriesFast()<<endm;
  for (UInt_t i = 0; i < tree->GetEntriesFast(); i++)
    {
      tree->GetEntry(i);
      int runId = mEvent->runId();
      if(!mEvent->isVPDMB()&&!mEvent->isHT0_BBCMB_TOF0()&&!mEvent->isHT2_BBCMB()) continue; // No intriesting trigger 
      if(mEvent->isVPDMB())   // reject bad run i MB
	{
	  for(int i=0;i<mNBadRuns_MB;i++){
	    if(runId==mBadRuns_MB[i])
	      {
		LOG_INFO<< " This is bad run !! rejected "<<endm;
		return ;
	      }
	  }
	}
      if(mEvent->isHT0_BBCMB_TOF0())   // reject bad run i HT0
	{
	  for(int i=0;i<mNBadRuns_HT0;i++){
	    if(runId==mBadRuns_HT0[i])
	      {
		LOG_INFO<< " This is bad run !! rejected "<<endm;
		return ;
	      }
	  }
	}
      if(mEvent->isHT2_BBCMB())   // reject bad run i HT2
	{
	  for(int i=0;i<mNBadRuns_HT2;i++){
	    if(runId==mBadRuns_HT2[i])
	      {
		LOG_INFO<< " This is bad run !! rejected "<<endm;
		return ;
	      }
	  }
	}


      // Run_QA(mEvent); 
      Events_Cuts_and_Counting(mEvent); // after Events cuts and trigger selection;
      aPairs=mEvent->electronPair();
      aTracks=mEvent->tracks();    
      if(mEvent->isHT0_BBCMB_TOF0())
	{ 
	  Int_t  bTrg=0;
	  Double_t  ps=mPrescales->GetPrescale(mEvent->runId(),HT0BBCMBTOF0);
	  if(ps<0) continue;
	  Fill_PhotonicE_hist ( bTrg, mEvent, ps);
	  Fill_Inclusive_hist ( bTrg, mEvent ,ps);
	}
      if(mEvent->isHT2_BBCMB())
        {
	  Int_t  bTrg=1;
	  Double_t  ps=mPrescales->GetPrescale(mEvent->runId(),HT2BBCMB);
          if(ps<0) continue;
          Fill_PhotonicE_hist ( bTrg, mEvent , ps );
          Fill_Inclusive_hist ( bTrg, mEvent , ps );
	} 
      
      if(mEvent->isVPDMB())
	{
	  continue; // MB trigger will add soon
	  Int_t       bTrg=2;
	  Double_t ps=mPrescales->GetPrescale(mEvent->runId(),VPDMB);
	  if(ps<0) continue;
	  // Fill_PhotonicE_hist ( bTrg, mNpeEvent , ps );
	  // Fill_Inclusive_hist ( bTrg, mNpeEvent , ps );
	  // Fill_Kaon_Kaon(mNpeEvent,bTrg);
	}       
      if(i%1000==0) cout<< " Working on event "<<i<<endl;
    }
}
Bool_t StNpeMaker:: Events_Cuts_and_Counting(StDmesonEvent* evt)
{
  if(!evt) return kFALSE;
  Float_t vz = evt->primaryVertex().z();
  if (!(evt->ranking() >cuts::mRanking )) return kFALSE;
  if(fabs(vz) > cuts::Tpc_vz) return kFALSE;
  if(evt->isVPDMB() && fabs(vz-evt->vzVpd())>cuts::diff_vzVpdVz) return kFALSE;

  // counting the #MB and #HT for Events counting
  if(evt->isBBCMB_TOF0() ||evt->isHT0_BBCMB_TOF0())
    {
      Int_t bTrg=0;
      Double_t ps_HT0_BBCMB_TOF0=mPrescales->GetPrescale(evt->runId(),HT0BBCMBTOF0);
      Double_t ps_BBCMB_TOF0=mPrescales->GetPrescale(evt->runId(),BBCMBTOF0);
      if(ps_BBCMB_TOF0>0 && ps_HT0_BBCMB_TOF0>0)
	{
	  if(evt->isBBCMB_TOF0())     
	    {  
	      mh1MB_Nevents[bTrg]->Fill(1.5);
	      mh1MB_Nevents_ps[bTrg]->Fill(1.5,ps_BBCMB_TOF0);
	    }  
	
	  if(evt->isHT0_BBCMB_TOF0())
	    {
	      mh1HT_Nevents[bTrg]->Fill(2.5);
	      mh1HT_Nevents_ps[bTrg]->Fill(2.5,ps_BBCMB_TOF0);
	      mh1hVz[bTrg]->Fill(evt->primaryVertex().z());
	      mh1hVz_ps[bTrg]->Fill(evt->primaryVertex().z(),ps_HT0_BBCMB_TOF0);
	    }

	}
    }
  // 
  if(evt->isBBCMB() ||evt->isHT2_BBCMB())
    { 
      Int_t bTrg=1;
      Double_t ps_BBCMB=mPrescales->GetPrescale(evt->runId(),BBCMB);
      Double_t ps_HT2BBCMB=mPrescales->GetPrescale(evt->runId(),HT2BBCMB);
      if( ps_BBCMB>0 && ps_HT2BBCMB>0)                
	{
	  if(evt->isBBCMB())
	    {
	      mh1MB_Nevents[bTrg]->Fill(1.5);
	      mh1MB_Nevents_ps[bTrg]->Fill(1.5,ps_BBCMB);
	    }   
	  if(evt->isHT2_BBCMB())
	    {
	      mh1HT_Nevents[bTrg]->Fill(2.5);
	      mh1HT_Nevents_ps[bTrg]->Fill(2.5,ps_HT2BBCMB);
	      mh1hVz[bTrg]->Fill(evt->primaryVertex().z());
	      mh1hVz_ps[bTrg]->Fill(evt->primaryVertex().z(),ps_HT2BBCMB);
	    }
	}
    }
  // VPDMB
  if(evt->isVPDMB())
    {
      Int_t bTrg=2;
      Double_t ps_VPDMB=mPrescales->GetPrescale(evt->runId(),VPDMB);
      if(ps_VPDMB>0)
	{               
	  mh1MB_Nevents[bTrg]->Fill(1.5);
	  mh1MB_Nevents_ps[bTrg]->Fill(1.5,ps_VPDMB);
	  mh1hVz[bTrg]->Fill(evt->primaryVertex().z());
	  mh1hVz[bTrg]->Fill(evt->primaryVertex().z(),ps_VPDMB);
	  mh2hVzVsDiff[bTrg]->Fill(evt->primaryVertex().z(),evt->primaryVertex().z()-evt->vzVpd());
	  mh2hVzVsVpdVz[bTrg]->Fill(evt->primaryVertex().z(),evt->vzVpd());
	}
    }
  if(!mEvent->isVPDMB()&&!mEvent->isHT0_BBCMB_TOF0()&&!mEvent->isHT2_BBCMB()) return kFALSE; // No intriesting trigger )
  else return kTRUE;
}
void StNpeMaker::Fill_Inclusive_hist (Int_t bTrg,StDmesonEvent * mEvent ,Double_t ps )
{
  TClonesArray*   aTracks = 0;
  aTracks=mEvent->tracks();
  for(Int_t it=0;it<mEvent->nTracks();it++)
    {
      StDmesonTrack* Trk = (StDmesonTrack*)aTracks->At(it);
      if(!pass_cut_acceptance(Trk)) continue; // apply pt eta cuts on aprimary electron
      if(!pass_cut_Track_qaulity(Trk)) continue; // pass track quality on primary electron 
      // below three for the hadron nsigmae calibratio
      //------------------------------------------------------------------------------------------------------
      mh2Pion_nSigmaElec[bTrg]->Fill(Trk->gMom().perp(),Trk->nSigmaElectron()-Trk->nSigmaPion());
      mh2Kaon_nSigmaElec[bTrg]->Fill(Trk->gMom().perp(),Trk->nSigmaElectron()-Trk->nSigmaKaon());
      mh2Proton_nSigmaElec[bTrg]->Fill(Trk->gMom().perp(),Trk->nSigmaElectron()-Trk->nSigmaProton());
      //------------------------------------------------------------------------------------------------------
      if(!pass_cut_Trigger_electron(Trk,bTrg )) continue; // pass triggere cuts     
      if((bTrg==0||bTrg==1)&&isHotTower(Trk,bTrg)) continue; // reject hot tower 
      if(Trk->trgTowDsmAdc()>Trk->adc0()*0.1) continue; // since embedding reject those band,keep same with embedding 
      if(pass_cut_poe(Trk)&&pass_cut_Match_BEMC(Trk)&&pass_cut_BSMD(Trk)) // without nsigma e cut ,fit purity
	{
	  mh2nSigmaElec[bTrg]->Fill(Trk->nSigmaElectron(),Trk->gMom().perp());
	  mh2nSigmaElec_ps[bTrg]->Fill(Trk->nSigmaElectron(),Trk->gMom().perp(),ps);
	  if(pass_cut_nsigmaE(Trk)) // now all the cuts applied 
	    {
	      mh1electronPt[bTrg]->Fill(Trk->gMom().perp());
	      mh1electronPt_ps[bTrg]->Fill(Trk->gMom().perp(),ps); // final spectra
	      //    cout<< " psss includeive "<<endl;
	    }
	}
    }
}
void StNpeMaker::Fill_PhotonicE_hist (Int_t bTrg,StDmesonEvent * mEvent ,Double_t ps )
{
  TClonesArray*   aTracks = 0;
  TClonesArray* aPairs=0;
  aPairs=mEvent->electronPair();
  aTracks=mEvent->tracks();
  for(Int_t ip=0;ip<mEvent->nElectronPair();ip++)
    {
      StElectronPair* pair = (StElectronPair*)aPairs->At(ip);
      StDmesonTrack* eTrk = (StDmesonTrack*)aTracks->At(pair->electronId()); // primary electron 
      StDmesonTrack* pTrk = (StDmesonTrack*)aTracks->At(pair->partnerId()); // partner electron 
      if(!pass_cut_eePair(pair)) continue; // pass pair cuts  
      if(!pass_cut_Trigger_electron(eTrk,bTrg )) continue; // pass triggere cuts     
      if((bTrg==0||bTrg==1)&&isHotTower(eTrk,bTrg)) continue; // reject hot tower 
      if(eTrk->trgTowDsmAdc()>eTrk->adc0()*0.1) continue; // since embedding reject those band,keep same with embedding 
      if(!pass_cut_acceptance(eTrk)) continue; // apply pt eta cuts on aprimary electron
      if(!pass_cut_Track_qaulity(eTrk)) continue; // pass track quality on primary electron 
      if(bTrg==0||bTrg==1) // HT0 and HT2
	{
	  Fill_pair_hist_HT(bTrg,pair,eTrk,pTrk,ps); // eID will applied
	}
    }
}
void StNpeMaker::Fill_pair_hist_HT(Int_t iTrg,StElectronPair * pair, StDmesonTrack* eTrk,StDmesonTrack * pTrk ,Double_t ps ) 
{
  // all the eid cuts applied 
  if(pass_cut_nsigmaE(eTrk)&&pass_cut_poe(eTrk)&&pass_cut_Match_BEMC(eTrk)&&pass_cut_BSMD(eTrk))
    {
      if(eTrk->charge()!=pTrk->charge())
	{ 
	  mHitFit_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit());
	  mNSMDEta_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nEta());
	  mHitsDedxUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsDedx());
	  mNTrack_cutUnlike[iTrg]->Fill(eTrk->gMom().perp());
	  mFitPos_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax());
	  mgDcaUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dca());
  	  mPhi_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().phi());
	  mEta_ptUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().pseudoRapidity());
	  mDedxUnlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dEdx());
	  mPoeUnlike[iTrg]->Fill(eTrk->gMom().mag()/eTrk->e0(),eTrk->gMom().perp());

	  mh2InvMassUnlike_ps[iTrg]->Fill(pair->m(),eTrk->gMom().perp(),ps);
	  mh2InvMassUnlike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  
	  mh3EMC_PartUnlike[iTrg]->Fill(pTrk->nSigmaElectron(),eTrk->gMom().perp(),pair->m());
	  mh3EtaPhiUnlike[iTrg]->Fill(eTrk->gMom().pseudoRapidity(),eTrk->gMom().phi(),eTrk->gMom().perp());
	  mh3EMC_ADCPartUnlike[iTrg]->Fill(pTrk->nSigmaElectron(),eTrk->trgTowDsmAdc(),pair->m());
	}
      if(eTrk->charge()==pTrk->charge()) 
	{
	  mHitFit_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit());
	  mHitsDedxlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsDedx());
	  mNTrack_cutlike[iTrg]->Fill(eTrk->gMom().perp());
	  mFitPos_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->nHitsFit()/(Float_t)eTrk->nHitsMax());
	  mgDcalike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dca()); 
	  mPhi_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().phi());
	  mEta_ptlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->gMom().pseudoRapidity());  
	  mDedxlike[iTrg]->Fill(eTrk->gMom().perp(),eTrk->dEdx()); 
	  mPoelike[iTrg]->Fill(eTrk->gMom().mag()/eTrk->e0(),eTrk->gMom().perp());
	  
	  mh2InvMasslike[iTrg]->Fill(pair->m(),eTrk->gMom().perp());
	  mh2InvMasslike_ps[iTrg]->Fill(pair->m(),eTrk->gMom().perp(),ps);
	  mh3EMC_Partlike[iTrg]->Fill(pTrk->nSigmaElectron(),eTrk->gMom().perp(),pair->m());
	  mh3EtaPhilike[iTrg]->Fill(eTrk->gMom().pseudoRapidity(),eTrk->gMom().phi(),eTrk->gMom().perp());
	  mh3EMC_ADCPartlike[iTrg]->Fill(pTrk->nSigmaElectron(),eTrk->trgTowDsmAdc(),pair->m());
	}
      //    cout<<" pass !!!"<<endl;
    }
}
Bool_t StNpeMaker::pass_cut_eePair(StElectronPair * Pair)
{
  Float_t pairMass=Pair->m();
  Float_t Pairdca=Pair->pairDca();
  // cout<< pairMass<< " "<< Pairdca<<endl;
  if(pairMass<cuts::InvMass&&Pairdca<cuts::pairDca)
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeMaker::pass_cut_Trigger_electron(StDmesonTrack *trk, int bTrg )
{
  if(bTrg==0&&trk->trgTowDsmAdc()>11 && trk->adc0()>180) return kTRUE;
  if(bTrg==1&&trk->trgTowDsmAdc()>18 && trk->adc0()>300) return kTRUE;
  if(bTrg==2) return kTRUE; // MB
  return kFALSE;
}
Bool_t StNpeMaker::pass_cut_acceptance(StDmesonTrack * trk)
{
  Float_t Pt=trk->gMom().perp();
  Float_t eTa=trk->gMom().pseudoRapidity();
 
  if(fabs(eTa)<cuts::eta_max && Pt>cuts::pt_min)
    return kTRUE;
  else return kFALSE;
}

Bool_t StNpeMaker::pass_cut_Track_qaulity(StDmesonTrack * trk)
{
  Int_t nhitDedx=trk->nHitsDedx();
  Int_t nhitsFit=trk->nHitsFit();
  Int_t nhitsMax=trk->nHitsMax();
  Float_t gDca=trk->dca();  
  Float_t first_point= trk->firstTpcPointR();
  if(gDca<cuts::gDca && cuts::nhitDedx<nhitDedx && nhitsFit>cuts::nhitsFit && (Float_t) nhitsFit/nhitsMax>cuts::nhit_Possmax  &&  first_point<cuts::first_point_cut)
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeMaker::pass_cut_nsigmaE(StDmesonTrack* trk)
{
  if(!trk) return kFALSE;
  Float_t nSigmaE=trk->nSigmaElectron();
  
  if(nSigmaE>cuts::nsigmae_low && nSigmaE<cuts::nsigmae_high)  
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeMaker::pass_cut_poe(StDmesonTrack * trk)
{
  Float_t P=trk->gMom().mag();
  Float_t E=trk->e0();
  
  if(P/E>cuts::poe_low && P/E<cuts::poe_high)
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeMaker::pass_cut_BSMD(StDmesonTrack * trk)
{
  Int_t nEta=trk->nEta();
  Int_t nPhi=trk->nPhi();
  if(nEta>cuts::NEta && nPhi>cuts::NPhi && cuts::Pt_cut<=trk->gMom().perp())
    return kTRUE;
  else return kFALSE;
}
Bool_t StNpeMaker::pass_cut_Match_BEMC(StDmesonTrack * trk )
{
  if(cuts::Dz>fabs(trk->zDist()) && cuts::Dphi>fabs(trk->phiDist())&&cuts::Pt_cut<=trk->gMom().perp())
     return kTRUE;
     else return kFALSE;
}
Bool_t StNpeMaker::isHotTower(StDmesonTrack *trk,Int_t bTrg)
{
  // Int_t  HT0_Hot_towerlist[] ={32,52,115,246,268,276,294,386,510,562,682,750,773,894,987,994,1043,1064,1143,1233,1264,1285,1307,1487,1593,1710,1733,1823,1824,1851,1946,2022,2044,2064,2110,2146,2163,2165,2203,2291,2314,2522,2530,2634,2653,2835,2864,2866,2973,3006,3062,3533,3545,3727,3862,3949,4051,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};

  // Int_t HT2_Hot_towerlist[]={32,52,115,268,276,294,295,510,562,682,750,987,994,1064,1143,1233,1264,1285,1307,1487,1593,1733,1824,1851,1946,2044,2163,2203,2634,2653,2835,2864,2866,2973,3006,3727,3862,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};

  // Int_t  Hot_towerlist[] ={32,52,115,246,268,276,294,295,386,510,562,682,750,773,894,987,994,1043,1064,1143,1233,1264,1285,1307,1487,1593,1710,1733,1823,1824,1851,1946,2022,2044,2064,2110,2146,2163,2165,2203,2291,2314,2522,2530,2634,2653,2835,2864,2866,2973,3006,3062,3533,3545,3691,3727,3862,3949,4051,4131,4170,4263,4431,4459,4684,4685,4686,4705,4767};

   Int_t  Hot_towerlist[]={22,30,31,114,251,275,308,509,552,555,681,682,691,749,772,1063,1263,1268,1284,1304,1306,1329,1486,1592,1709,1768,1770,1823,1882,1904,1909,1945,2022,2042,2043,2048,2051,2067,2145,2146,2162,2171,2190,2202,2272,2288,2290,2493,2504,2522,2529,2549,2706,2723,2712,2750,2863,2865,2868,2874,2952,3061,3007,3061,3092,3112,3154,3264,3166,3269,3271,3307,3326,3330,3331,3349,3365,3373,3532,3544,3626,3692,3821,3861,3932,3964,4130,4169,4226,4232,4249,4252,4262,4353,4430,4546,4749,4727,4766};

   for(Int_t i=0;i<sizeof(Hot_towerlist)/sizeof(Hot_towerlist[0]);i++ )
     {
       if(Hot_towerlist[i]==trk->btowId())
	 return kTRUE;
     }
   return kFALSE;
}

void StNpeMaker::bookObjects(){
  mEvent = new StDmesonEvent;
  LOG_INFO <<"\tInitialising Histograms"<<endm;
  
  for(Int_t iTrg=0;iTrg<3;iTrg++)  
    {
      mh1MB_Nevents[iTrg]=new TH1F(Form("mh1MB_NeventsTrg%i",iTrg),"mh1MB_Nevents",10,0,10);
      mh1HT_Nevents[iTrg]=new TH1F(Form("mh1HT_NeventsTrg%i",iTrg),"mh1HT_Nevents",10,0,10);
      mh1MB_Nevents_ps[iTrg]=new TH1F(Form("mh1MB_Nevents_psTrg%i",iTrg),"mh1MB_Nevents_ps",10,0,10);
      mh1HT_Nevents_ps[iTrg]=new TH1F(Form("mh1HT_Nevents_psTrg%i",iTrg),"mh1HT_Nevents_ps",10,0,10);
      
      mh1hVz[iTrg] = new TH1F(Form("mh1hVzTrg%i",iTrg), ";TPC Vz (cm)", 1000, -40+1e-6, 40+1e-6);
      mh1hVz_ps[iTrg] = new TH1F(Form("mh1hVz_psTrg%i",iTrg), ";TPC Vz (cm)", 400, -40+1e-6, 40+1e-6);

      mh2hVzVsDiff[iTrg] = new TH2F(Form("mh2hVzVsDiffTrg%i",iTrg), ";TPC Vz (cm); Tpc Vz-VPD Vz (cm)", 400, -40+1e-6, 40+1e-6, 200, -10+1e-6, 10+1e-6);
      mh2hVzVsVpdVz[iTrg] = new TH2F(Form("mh2hVzVsVpdVzTrg%i",iTrg), ";TPC Vz (cm);VPD Vz (cm)", 400, -40+1e-6, 40+1e-6, 400, -40+1e-6, 40+1e-6);

      // all the cuts appllied
      mPhi_ptUnlike[iTrg]=new TH2F(Form("mPhi_ptUnlikeTrg%i",iTrg)," ",200,0,20,400,-4,4);
      mEta_ptUnlike[iTrg]=new TH2F(Form("mEta_ptUnlikeTrg%i",iTrg)," ",200,0,20,300,-1.5,1.5);
      mPhi_ptlike[iTrg]=new TH2F(Form("mPhi_ptlikeTrg%i",iTrg)," ",200,0,20,400,-4,4);
      mEta_ptlike[iTrg]=new TH2F(Form("mEta_ptlikeTrg%i",iTrg)," ",200,0,20,300,-1.5,1.5);
      mHitFit_ptUnlike[iTrg] =new TH2F(Form("mHitFit_ptUnlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
      mHitFit_ptlike[iTrg] =new TH2F(Form("mHitFit_ptlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
      mNSMDEta_ptUnlike[iTrg]=new TH2F(Form("mNSMDEta_ptUnlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
      mNSMDEta_ptlike[iTrg]=new TH2F(Form("mNSMDEta_ptlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
      mNSMDPhi_ptUnlike[iTrg]=new TH2F(Form("mNSMDPhi_ptUnlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
      mNSMDPhi_ptlike[iTrg]=new TH2F(Form("mNSMDPhi_ptlikeTrg%i",iTrg)," ",200,0,20,10,0,10);
      mHitsDedxUnlike[iTrg]=new TH2F(Form("mHitsDedxUnlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
      mHitsDedxlike[iTrg]=new TH2F(Form("mHitsDedxlikeTrg%i",iTrg)," ",200,0,20,60,0,60);
      mNTracklike[iTrg]=new TH1F(Form("mNTracklikeTrg%i",iTrg)," ",200,0,20);
      mNTrackUnlike[iTrg]=new TH1F(Form("mNTrackUnlikeTrg%i",iTrg)," ",200,0,20);
      mNTrack_cutUnlike[iTrg]=new TH1F(Form("mNTrack_cutUnlikeTrg%i",iTrg)," ",200,0,20);
      mNTrack_cutlike[iTrg]=new TH1F(Form("mNTrack_cutlikeTrg%i",iTrg)," ",200,0,20);     
      mNTrack_cut25Unlike[iTrg]=new TH1F(Form("mNTrack_cut25UnlikeTrg%i",iTrg)," ",200,0,20);
      mNTrack_cut25like[iTrg]=new TH1F(Form("mNTrack_cut25likeTrg%i",iTrg)," ",200,0,20);
      mNTrack_cut20Unlike[iTrg]=new TH1F(Form("mNTrack_cut20UnlikeTrg%i",iTrg)," ",200,0,20);
      mNTrack_cut20like[iTrg]=new TH1F(Form("mNTrack_cut20likeTrg%i",iTrg)," ",200,0,20);
      mFitPos_ptlike[iTrg]=new TH2F(Form("mFitPos_ptlikeTrg%i",iTrg)," ",200,0,20,100,0,1);
      mFitPos_ptUnlike[iTrg]=new TH2F(Form("mFitPos_ptUnlikeTrg%i",iTrg)," ",200,0,20,100,0,1);
      mgDcalike[iTrg]=new TH2F(Form("mgDcalikeTrg%i",iTrg),"",200,0,20,50,0,5);
      mgDcaUnlike[iTrg]=new TH2F(Form("mgDcaUnlikeTrg%i",iTrg),"",200,0,20,50,0,5);
      mNsigElike[iTrg]=new TH2F(Form("mNsigElikeTrg%i",iTrg),"",200,0,20,100,-5+1e-6,5+1e-6);
      mNsigEUnlike[iTrg]=new TH2F(Form("mNsigEUnlikeTrg%i",iTrg),"",200,0,20,50,-5+1e-6,5+1e-6);
      mDedxlike[iTrg]=new TH2F(Form("mDedxlikeTrg%i",iTrg),"",200,0,20,60,1,6);
      mDedxUnlike[iTrg]=new TH2F(Form("mDedxUnlikeTrg%i",iTrg),"",200,0,20,60,1,6);
      mPoeUnlike[iTrg] = new TH2F(Form("mPoeUnlikeTrg%i",iTrg),"",100,0,2.5,200,0,20.);
      mPoelike[iTrg] = new TH2F(Form("mPoelikeTrg%i",iTrg),"",100,0,2.5,200,0,20.);


      mh2InvMassUnlike[iTrg]=new TH2F(Form("mh2InvMassUnlikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
      mh2InvMasslike[iTrg]=new TH2F(Form("mh2InvMasslikeTrg%i",iTrg),"",40,0,0.4,200 ,0,20);
      mh2InvMassUnlike_ps[iTrg]=new TH2F(Form("mh2InvMassUnlike_psTrg%i",iTrg),"",40,0,0.4,200 ,0,20.);
      mh2InvMasslike_ps[iTrg]=new TH2F(Form("mh2InvMasslike_psTrg%i",iTrg),"",40,0,0.4,200 ,0,20);

      mh3EMC_PartUnlike[iTrg]=new TH3F(Form("mh3EMC_PartUnlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20.,40,0,0.4);
      mh3EMC_Partlike[iTrg]=new TH3F(Form("mh3EMC_PartlikeTrg%i",iTrg),"",140,-7,7,200 ,0,20,40,0,0.4);
      mh3EMC_ADCPartUnlike[iTrg]=new TH3F(Form("mh3EMC_ADCPartUnlikeTrg%i",iTrg),"",140,-7,7,40 ,0,40.,40,0,0.4);
      mh3EMC_ADCPartlike[iTrg]=new TH3F(Form("mh3EMC_ADCPartlikeTrg%i",iTrg),"",140,-7,7,40 ,0,40,40,0,0.4);     
      mh3EtaPhiUnlike[iTrg]=new TH3F(Form("mh3EtaPhiUnlikeTrg%i",iTrg),"",20 ,-1,1,80,-4,4,200,0,20);
      mh3EtaPhilike[iTrg]=new TH3F(Form("mh3EtaPhilikeTrg%i",iTrg),"",20 ,-1,1,80,-4,4,200,0,20);

      //inclusive electron 

      mh2nSigmaElec[iTrg] = new TH2F(Form("mh2nSigmaElecTrg%i",iTrg),"",200,-10+1e-6,10+1e-6,200 ,0,20.);
      mh2nSigmaElec_ps[iTrg] = new TH2F(Form("mh2nSigmaElec_psTrg%i",iTrg),"",200,-10+1e-6,10+1e-6,200 ,0,20.);
      mh1electronPt[iTrg]=new TH1F(Form("mh1electronPtTrg%i",iTrg),"",200,0,20.);
      mh1electronPt_ps[iTrg]=new TH1F(Form("mh1electronPt_psTrg%i",iTrg),"",200,0,20.);

      mh2Pion_nSigmaElec[iTrg] =new TH2F(Form("mh2Pion_nSigmaElecTrg%i",iTrg),"",200,0,20, 700,-35+1e-6,35+1e-6);
      mh2Kaon_nSigmaElec[iTrg] =new TH2F(Form("mh2Kaon_nSigmaElecTrg%i",iTrg),"",200,0,20,700,-35+1e-6,35+1e-6);
      mh2Proton_nSigmaElec[iTrg] =new TH2F(Form("mh2Proton_nSigmaElecTrg%i",iTrg),"",200,0,20,700,-35+1e-6,35+1e-6);

      
      mh1MB_Nevents[iTrg]->Sumw2();
      mh1HT_Nevents[iTrg]->Sumw2();
      mh1MB_Nevents_ps[iTrg]->Sumw2();
      mh1HT_Nevents_ps[iTrg]->Sumw2();
      
      mh2hVzVsDiff[iTrg]->Sumw2();
      mh2hVzVsVpdVz[iTrg]->Sumw2();
      mh1hVz[iTrg]->Sumw2();
      mh1hVz_ps[iTrg]->Sumw2();
    }
}
void StNpeMaker::writeObjects(){
  
  mOutputFile->cd();
  for(Int_t iTrg=0;iTrg<nTrg;iTrg++)  
    {
      mh1MB_Nevents[iTrg]->Write();
      mh1HT_Nevents[iTrg]->Write();
      mh1MB_Nevents_ps[iTrg]->Write();
      mh1HT_Nevents_ps[iTrg]->Write();

      mh2hVzVsDiff[iTrg]->Write();
      mh2hVzVsVpdVz[iTrg]->Write();
      mh1hVz[iTrg]->Write();
      mh1hVz_ps[iTrg]->Write();

      // all cuts applied 
      mPhi_ptUnlike[iTrg]->Write();
      mEta_ptUnlike[iTrg]->Write();
      mPhi_ptlike[iTrg]->Write();
      mEta_ptlike[iTrg]->Write();
       
      mHitFit_ptUnlike[iTrg]->Write();
      mNSMDEta_ptUnlike[iTrg]->Write();
      mHitsDedxUnlike[iTrg]->Write();
      mHitFit_ptlike[iTrg]->Write();
      mNSMDEta_ptlike[iTrg]->Write();
      mHitsDedxlike[iTrg]->Write();

      mNTrackUnlike[iTrg]->Write();
      mNTrack_cutUnlike[iTrg]->Write();
      mNTracklike[iTrg]->Write();
      mNTrack_cutlike[iTrg]->Write();
      mFitPos_ptlike[iTrg]->Write();
      mgDcalike[iTrg]->Write();
      mFitPos_ptUnlike[iTrg]->Write();
      mgDcaUnlike[iTrg]->Write();
      mNSMDPhi_ptUnlike[iTrg]->Write();
      mNSMDPhi_ptlike[iTrg]->Write();
      mNTrack_cut25Unlike[iTrg]->Write();
      mNTrack_cut25like[iTrg]->Write();
      mNTrack_cut20Unlike[iTrg]->Write();
      mNTrack_cut20like[iTrg]->Write();

      mNsigElike[iTrg]->Write();
      mNsigEUnlike[iTrg]->Write();

      mDedxlike[iTrg]->Write();
      mDedxUnlike[iTrg]->Write();
      mPoeUnlike[iTrg] ->Write();
      mPoelike[iTrg] ->Write();

      mh2InvMassUnlike[iTrg] ->Write();
      mh2InvMasslike[iTrg] ->Write();
      mh2InvMassUnlike_ps[iTrg] ->Write();
      mh2InvMasslike_ps[iTrg] ->Write();

      mh3EMC_PartUnlike[iTrg] ->Write();
      mh3EMC_Partlike[iTrg] ->Write();
      mh3EMC_ADCPartUnlike[iTrg] ->Write();
      mh3EMC_ADCPartlike[iTrg] ->Write();
      mh3EtaPhiUnlike[iTrg] ->Write();
      mh3EtaPhilike[iTrg] ->Write();

      // inclusive e
      mh2nSigmaElec[iTrg]->Write();
      mh2nSigmaElec_ps[iTrg]->Write();
      mh1electronPt[iTrg]->Write(); 
      mh1electronPt_ps[iTrg]->Write();

      mh2Pion_nSigmaElec[iTrg]->Write();
      mh2Kaon_nSigmaElec[iTrg]->Write();
      mh2Proton_nSigmaElec[iTrg]->Write(); 
    }
}