# ifndef StSingleElectron_EmbeddingQA_h
# define StSingleElectron_EmbeddingQA_h

#include "StSingle_Electron_tracks.h"
#include "StSingle_electron_EventCount.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

class TString;
class TFile;
class TF1;

class TH1F;
class TH2F;




class StSingle_Electron_tracks;
class StSingle_electron_EventCount;

class StSingleElectron_EmbeddingQA
{
 public:
  StSingleElectron_EmbeddingQA(const char *outName);
  virtual ~StSingleElectron_EmbeddingQA();

  void bookHistogram() ;
  void read(TString filename);
  void WriteHistogram();
  bool isHotTower(StSingle_Electron_tracks *);

 private:
  TFile * mOutputFile;
  StSingle_Electron_tracks *AcTracks;
  // StMcAccpEventCount *mAcEvent;   

  // Event

  TH1F * mDelta_Vz;
  //track
  TH2F * mPhi_pt_cut;
  TH2F * mEta_pt_cut;
  TH2F* mHitFit_pt;
  TH2F*  mHitMax_pt;
  TH2F* mFitPos_pt;
  TH2F* mPoe;
  TH2F *  mHitsDedx;
  TH2F * mNSMDEta_pt;
  TH2F * mNSMDPhi_pt;
  //
 TH1F *  mNTrack;
  TH1F *  mNTrack_cut;
  TH1F *mNTrack_cut_25;
  TH1F *mNTrack_cut_20;
  //cuts
  TH2F* mHitFit_pt_cut;
  TH2F*  mHitMax_pt_cut;
  TH2F* mFitPos_pt_cut;
 
  TH2F *  mHitsDedx_pt_cut;
  TH2F * mNSMDEta_pt_cut;
  TH2F * mNSMDPhi_pt_cut;
  TH2F * mPoe_pt_cut; 

  TH2F* mgDca_pt_cut;
  TH2F* mNsigE_pt_cut;
  TH2F * mDedx_pt_cut;
  TH2F*  DSMadc_pt;
  TH2F*  DSMadc_pt_cut;
  TH2F *  adc0_pt;
  TH2F *  adc0_pt_cut;
  TH2F * DSM_adc0;
  TH1F *DSMAdc_over_adc0;
  TH1F *low_DSMAdc;
  TH1F * TowerID_cut;
  TH1F * TowerID;
  TH1F * TowerID_all;
  TH1F * TowerID_cut_noHOt;
  TH1F * TowerID_noHOt;
  TH1F * TowerID_all_noHOt;

  //  Float_t Pt_xx[9]={2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,16.5};
 
  ClassDef(StSingleElectron_EmbeddingQA,1)

};

#endif
