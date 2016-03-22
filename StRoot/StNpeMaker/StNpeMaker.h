 /******************************************************************************************
 *Mon Mar  7 19:33:04 EST 2016                                                            *
 * This class for read  DemesonTree Fill all the No-phoronic electron realted histogram   *
 * By xiaozhi Bai xiaozhi@uic.edu                                                         *
 ******************************************************************************************/
 
// StNpeMaker.h 

#ifndef StNpeMaker_h
#define StNpeMaker_h

#include <exception>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include "TProfile.h"

class TString;
class TF1;
class TH1D;
class TH1F;
class TH2F;
class TH3F;
class TFile;
class TNtuple;
class StDmesonTrack;
class StDmesonEvent;
class prescales;
class StElectronPair;

class StNpeMaker 
{
 public:
  StNpeMaker(const char *outName);
  virtual ~StNpeMaker();
  virtual void  Clear(Option_t *opt="");
  void bookObjects();
  void read(TString fileName);
  void writeObjects();

  void   Run_QA(StDmesonEvent * ); 
  Bool_t isGoodEvent(StDmesonEvent*);
  Bool_t Events_Cuts_and_Counting(StDmesonEvent *);
  void   Fill_PhotonicE_hist (Int_t ,StDmesonEvent * ,Double_t );
  void   Fill_Inclusive_hist (Int_t ,StDmesonEvent * ,Double_t );
  void   Fill_pair_hist_HT(Int_t,StElectronPair*,StDmesonTrack*,StDmesonTrack*,Double_t);
  void   Fill_pair_hist_MB(Int_t,StElectronPair*,StDmesonTrack*,StDmesonTrack*,Double_t);


  Bool_t pass_cut_eePair(StElectronPair* );
  Bool_t isHotTower(StDmesonTrack *,Int_t ); 
  Bool_t pass_cut_Trigger_electron(StDmesonTrack*,Int_t);
  Bool_t pass_cut_acceptance(StDmesonTrack *);
  Bool_t pass_cut_Track_qaulity(StDmesonTrack *);
  Bool_t pass_cut_nsigmaE(StDmesonTrack* );
  Bool_t pass_cut_poe(StDmesonTrack *);
  Bool_t pass_cut_BSMD(StDmesonTrack *);      
  Bool_t pass_cut_Match_BEMC(StDmesonTrack * );

  //  Bool_t pass_cut_Match_EMC_Dz(StDmesonTrack * );

 private:
  TFile* mOutputFile;
  StDmesonEvent * mEvent;
  prescales *mPrescales;
  ifstream file_runNumber;
  ofstream  outfile;
  Int_t mNBadRuns_MB;
  Int_t mNBadRuns_HT0;
  Int_t mNBadRuns_HT2;

  const static  int nTrg=3;

  // Event 
  TH1F *mh1MB_Nevents[nTrg];
  TH1F *mh1HT_Nevents[nTrg];
  TH1F *mh1MB_Nevents_ps[nTrg];
  TH1F *mh1HT_Nevents_ps[nTrg];

  TH1F *mh1hVz[nTrg];
  TH1F *mh1hVz_ps[nTrg];

  TH2F *mh2hVzVsVpdVz[nTrg];
  TH2F *mh2hVzVsDiff[nTrg];

  //all the cuts applied 
  TH2F*   mPhi_ptUnlike[nTrg];
  TH2F*   mEta_ptUnlike[nTrg];
  TH2F*   mPhi_ptlike[nTrg];
  TH2F*   mEta_ptlike[nTrg];
        
  TH2F*   mHitFit_ptUnlike[nTrg];
  TH2F*   mNSMDEta_ptUnlike[nTrg];
  TH2F*   mHitsDedxUnlike[nTrg];
  TH2F*   mHitFit_ptlike[nTrg];
  TH2F*   mNSMDEta_ptlike[nTrg];
  TH2F*   mHitsDedxlike[nTrg];
        
  TH1F*   mNTrackUnlike[nTrg];
  TH1F*   mNTrack_cutUnlike[nTrg];
  TH1F*   mNTracklike[nTrg];
  TH1F*   mNTrack_cutlike[nTrg];
  TH2F*   mFitPos_ptlike[nTrg];
  TH2F*   mgDcalike[nTrg];
  TH2F*   mFitPos_ptUnlike[nTrg];
  TH2F*   mgDcaUnlike[nTrg];
        
  TH2F*   mNSMDPhi_ptUnlike[nTrg];
  TH2F*   mNSMDPhi_ptlike[nTrg];
  TH1F*   mNTrack_cut25Unlike[nTrg];
  TH1F*   mNTrack_cut25like[nTrg];

  TH1F*   mNTrack_cut20Unlike[nTrg];
  TH1F*   mNTrack_cut20like[nTrg];
        
  TH2F*   mNsigElike[nTrg];
  TH2F*   mNsigEUnlike[nTrg];
        
  TH2F*   mDedxlike[nTrg];
  TH2F*   mDedxUnlike[nTrg];
  TH2F*   mPoelike[nTrg];
  TH2F*  mPoeUnlike[nTrg];   


  TH3F*   mh3EMC_PartUnlike[nTrg];
  TH3F*   mh3EMC_Partlike[nTrg];
  TH3F*   mh3EMC_ADCPartUnlike[nTrg];
  TH3F*   mh3EMC_ADCPartlike[nTrg];
  
  TH3F * mh3EtaPhiUnlike[nTrg];
  TH3F * mh3EtaPhilike[nTrg]; 
  TH2F * mh2InvMasslike[nTrg];
  TH2F * mh2InvMassUnlike[nTrg];
  TH2F * mh2InvMasslike_ps[nTrg];
  TH2F * mh2InvMassUnlike_ps[nTrg];
  
  // inclusive electron 


  TH2F* mh2nSigmaElec[nTrg];
  TH2F* mh2nSigmaElec_ps[nTrg];
  TH1F * mh1electronPt[nTrg];
  TH1F * mh1electronPt_ps[nTrg];
  //
  TH2F * mh2Pion_nSigmaElec[nTrg];
  TH2F * mh2Kaon_nSigmaElec[nTrg];
  TH2F * mh2Proton_nSigmaElec[nTrg];
  
  
  ClassDef(StNpeMaker, 1)
};
#endif
