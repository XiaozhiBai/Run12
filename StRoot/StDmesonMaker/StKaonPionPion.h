#ifndef StKaonPionPion_hh
#define StKaonPionPion_hh
#ifdef __ROOT__

#include "TObject.h"
#include "TClonesArray.h"
#include "StLorentzVectorF.hh"

class StPicoTrack;
class StPicoEvent;
class StKaonPion;

class StKaonPionPion : public TObject {
 public:
  StKaonPionPion();
  StKaonPionPion(StKaonPionPion*);
  StKaonPionPion(StKaonPion*, StPicoTrack*,Int_t,Int_t);
  virtual ~StKaonPionPion(); 
  void clear();

  Int_t   kaonPionId() const { return (Int_t)mKaonPionId;}
  Int_t   pionId() const { return (Int_t)mPionId;}
  Float_t m()    const { return (Float_t)mM/10000.;}
  Float_t pt()   const { return (Float_t)mPt/1000.;}
  Float_t eta()   const { return (Float_t)mEta/1000.;}
  Float_t phi()   const { return (Float_t)mPhi/1000.;}
  StLorentzVectorF lorentzVector();

 protected:
  UShort_t   mKaonPionId;
  UShort_t   mPionId;

  UShort_t mM; // mass * 10000. : KeV precision
  UShort_t mPt; // pt * 1000. : MeV precision
  Short_t mEta; // eta * 1000.
  Short_t mPhi; // phi * 1000.

    ClassDef(StKaonPionPion,1)
};


#endif
#endif

