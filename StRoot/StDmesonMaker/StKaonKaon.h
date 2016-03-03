#ifndef StKaonKaon_hh
#define StKaonKaon_hh
#ifdef __ROOT__

#include "TObject.h"
#include "TClonesArray.h"
#include "StLorentzVectorF.hh"

class StPicoTrack;
class StPicoEvent;

class StKaonKaon : public TObject {
 public:
  StKaonKaon();
  StKaonKaon(StKaonKaon*);
  StKaonKaon(StPicoTrack*, StPicoTrack*,Int_t k1Id,Int_t k2Id);
  virtual ~StKaonKaon(); 
  void clear();

  Int_t   kaon1Id() const { return (Int_t)mKaon1Id;}
  Int_t   kaon2Id() const { return (Int_t)mKaon2Id;}
  Float_t m()    const { return (Float_t)mM/10000.;}
  Float_t pt()   const { return (Float_t)mPt/1000.;}
  Float_t eta()   const { return (Float_t)mEta/1000.;}
  Float_t phi()   const { return (Float_t)mPhi/1000.;}
  Float_t openAngle() const { return (Float_t)mOpenAngle/1000.;}
  StLorentzVectorF lorentzVector();
  
 protected:
  UShort_t mKaon1Id;
  UShort_t mKaon2Id;
  
  UShort_t mM; // mass * 10000. : keV precision
  UShort_t mPt; // pt * 1000.   : MeV precision
  Short_t mEta; // eta * 1000.
  Short_t mPhi; // phi * 1000.
  Short_t mOpenAngle; // openAngle * 1000.

  ClassDef(StKaonKaon,1)
};


#endif
#endif

