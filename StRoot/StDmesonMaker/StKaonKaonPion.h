#ifndef StKaonKaonPion_hh
#define StKaonKaonPion_hh
#ifdef __ROOT__

#include "TObject.h"
#include "TClonesArray.h"
#include "StLorentzVectorF.hh"

class StPicoTrack;
class StPicoEvent;
class StKaonKaon;

class StKaonKaonPion : public TObject {
 public:
  StKaonKaonPion();
  StKaonKaonPion(StKaonKaonPion*);
  StKaonKaonPion(StKaonKaon*, StPicoTrack*,Int_t,Int_t);
  virtual ~StKaonKaonPion(); 
  void clear();

  Int_t   kaonKaonId() const { return (Int_t)mKaonKaonId;}
  Int_t   pionId() const { return (Int_t)mPionId;}
  Float_t m()    const { return (Float_t)mM/10000.;}
  Float_t pt()   const { return (Float_t)mPt/1000.;}
  Float_t eta()   const { return (Float_t)mEta/1000.;}
  Float_t phi()   const { return (Float_t)mPhi/1000.;}
  Float_t cosThetaStarPion() const { return (Float_t)mCosThetaStarPion/100.;}
  StLorentzVectorF lorentzVector();

 protected:
  UShort_t mKaonKaonId;
  UShort_t mPionId;

  UShort_t mM; // mass * 10000. : KeV precision	
  UShort_t mPt; // pt * 1000.   : MeV precision
  Short_t mEta; // eta * 1000.
  Short_t mPhi; // phi * 1000.
  Char_t mCosThetaStarPion; // cosThetaStarPion * 100.

  ClassDef(StKaonKaonPion,1)
};


#endif
#endif

