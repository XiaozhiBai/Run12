#ifndef StDmesonTrack_hh
#define StDmesonTrack_hh

class StPicoTrack;

#include "TObject.h"
#include "StThreeVectorF.hh"
#include <stdio.h>
#include <math.h>

class StDmesonTrack : public TObject {
 public:
  StDmesonTrack();
  ~StDmesonTrack();
  StDmesonTrack(StPicoTrack*,Short_t,UShort_t,UShort_t);

  void    Clear(const Option_t *opt="");
  //virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  StThreeVectorF pMom() const    { return mPMomentum; }
  StThreeVectorF gMom() const    { return mGMomentum; }
  Float_t dca() const            { return (Float_t)mGDca/1000.; }
  Short_t charge() const         { return (mNHitsFit>0) ? +1 : -1; }
  Int_t   nHitsFit() const       { return (mNHitsFit>0) ? (Int_t)mNHitsFit : (Int_t)(-1*mNHitsFit); }
  Int_t   nHitsMax() const       { return (Int_t)mNHitsMax; }
  Int_t   nHitsDedx() const      { return (Int_t)mNHitsDedx; }
  Float_t dEdx() const           { return (Float_t)mDedx/1000.; }
  Float_t nSigmaPion() const     { return (Float_t)mNSigmaPion/100.; }
  Float_t nSigmaKaon() const     { return (Float_t)mNSigmaKaon/100.; }
  Float_t nSigmaProton() const   { return (Float_t)mNSigmaProton/100.; }
  Float_t nSigmaElectron() const { return (Float_t)mNSigmaElectron/100.; }
  Float_t firstTpcPointR() const { return (Float_t)mFirstTpcPointR/100.;}

  Int_t   btofMatchFlag() const  { return (Int_t)mBTofMatchFlag; }
  Float_t btofBeta() const       { return (Float_t)mBTofBeta/20000.; }
  Float_t btofYLocal() const     { return (Float_t)mBTofYLocal/1000.; }

  Int_t   adc0() const           { return (Int_t)mBTOWADC0; }
  Float_t e0() const             { return (Float_t)mBTOWE0/1000.; }
  Float_t e() const              { return (Float_t)mBTOWE/1000.; }
  Float_t zDist() const          { return (Float_t)mBEMCDistZ/100.; }
  Float_t phiDist() const        { return (Float_t)mBEMCDistPhi/10000.; }
  Int_t   nEta() const           { return (Int_t)mBSMDNEta; }
  Int_t   nPhi() const           { return (Int_t)mBSMDNPhi; }  

  Int_t   btowId() const         { return (Int_t)mBTOWId; }
  Int_t   trgTowId() const       { return (Int_t)mTrgTowId;}
  Int_t   trgTowDsmAdc() const   { return (Int_t)mTrgTowDSMADC;}
  Int_t   trgTowAdc() const      { return (Int_t)mTrgTowADC;}

 protected:
  StThreeVectorF mPMomentum;  // primary momentum, (0.,0.,0.) if none
  StThreeVectorF mGMomentum;  // global  momentum
  UShort_t mGDca;             // global dca*1000
  UShort_t mDedx;             // dEdx*1000
  Char_t   mNHitsFit;         // q*nHitsFit
  UChar_t  mNHitsMax;         // nHitsMax
  UChar_t  mNHitsDedx;        // nHitsDedx
  Short_t  mNSigmaPion;       // nsigmaPi * 100
  Short_t  mNSigmaKaon;       // nsigmaK * 100
  Short_t  mNSigmaProton;     // nsigmaP * 100
  Short_t  mNSigmaElectron;   // nsigmaE * 100
  UShort_t mFirstTpcPointR;       // R=sqrt(x*x+y*y) of first TPC point * 100; mm precision

  // BTOF variables  
  UChar_t  mBTofMatchFlag;    // 0 - no match, 1 - one-to-one, 2 - one-to-multiple
  UShort_t mBTofBeta;         // beta * 20000
  Short_t  mBTofYLocal;       // ylocal * 1000

  // these variables are extracted from the standard BEMC cluster algorithm
  Short_t  mBTOWADC0;         // adc0 higest adc in the cluster
  Short_t  mBTOWE0;           // E0*1000 highest tower in the cluster
  Short_t  mBTOWE;            // EMC point E*1000 
  Short_t  mBEMCDistZ;        // z*100
  Short_t  mBEMCDistPhi;      // phi*10000
  UChar_t  mBSMDNEta;         // # of hits in eta
  UChar_t  mBSMDNPhi;         // # of hits in phi

  // these variables are purely from single tower or nearby towers  
  Short_t  mBTOWId;           // projected tower Id 1-4800

  // trigger tower associated with the track has |trkProjTwo φ  -trgTow φ| < 0.06 && |trkProjη - trgProjη| < 0.06 
  // see StDmesonMaker::getTrgBTOW 
  Short_t   mTrgTowId;       
  UShort_t  mTrgTowDSMADC;       // DSM adc of trigger tower
  UShort_t  mTrgTowADC;       // adc of trigger tower

  ClassDef(StDmesonTrack, 1)
};

#endif
