#ifndef StElectronPair_hh
#define StElectronPair_hh
#ifdef __ROOT__

#include "TObject.h"
#include "TClonesArray.h"

class StPicoTrack;
class StPicoEvent;

class StElectronPair : public TObject {
 public:
  StElectronPair();
  StElectronPair(StElectronPair*);
  StElectronPair(StPicoTrack*, StPicoTrack*, Int_t,Int_t,Float_t);
  virtual ~StElectronPair(); 
  void clear();

  Int_t   electronId() const { return (Int_t)mElectronId;}
  Int_t   partnerId() const { return (Int_t)mPartnerId;}
  Float_t m() const { return (Float_t)mM/100000.;}
  Float_t pairDca() const { return mPairDca>0 ? 1.*mPairDca/100. : -1.*mPairDca/100.;}
  Short_t partnerCharge() const { return (mPairDca>0) ? +1 : -1; }

          
 protected:
  UShort_t  mElectronId; // StDmesonTrack id; not picoDst Id
  UShort_t  mPartnerId; // StDmesonTrack id; not picoDst Id

  UShort_t mM; // mass * 100000. : 100eV precision
  Char_t mPairDca; // partnerQ*pairDca*100. : mm precision

  ClassDef(StElectronPair,1)
};


#endif
#endif

