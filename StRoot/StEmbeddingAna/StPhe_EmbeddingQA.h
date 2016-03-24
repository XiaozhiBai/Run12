# ifndef StPhe_EmbeddingQA_h
# define  StMakeDataQa_Gamma_h

/* #include "StMcAccpTracks.h" */
/* #include "StMcAccpEventCount.h" */
//#include "StMcAccpTracks.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

class TString;
class TF1;
class TH1D;
class TH2D;

class  TH3F;

class StPheTrack;
class StPhePair;

class StPhe_EmbeddingQA
{
 public:
  StPhe_EmbeddingQA(const char *outName);
  virtual ~StPhe_EmbeddingQA();

  void bookHistogram() ;
  void read(TString filename);
  void WriteHistogram();

 private:
  TFile * mOutputFile;
  StPheTrack *AcTracks;
  //  StMcAccpEventCount *mAcEvent;   

  // Event

  TH1F * mDelta_Vz;
  //track
  TH2F * mPhi_pt;
  TH2F * mEta_pt;
  TH2F* mHitFit_pt;
  TH2F*  mHitMax_pt;
  TH2F* mFitPos_pt;
  TH2F* mPoe;
  TH2F *  mHitsDedx;
  TH1F *  mNTrack;
  TH1F *  mNTrack_cut;

  //cuts
  TH2F* mHitFit_pt_cut;
   TH2F* mFitPos_pt_cut;
   TH2F *  mHitsDedx_cut;
  TH2F*  mNsigE_pt_cut;
 TH2F *  mPhi_pt_cut;
 TH2F* mEta_pt_cut;
 TH2F * mNSMDEta_pt_cut;
 TH2F * mNSMDPhi_pt_cut;
 TH2F * mgDca_pt_cut;
 TH2F * mDedx_pt_cut;
 TH2F* mInvMass_pt_cut;
  ClassDef(StPhe_EmbeddingQA,1)

};

#endif
