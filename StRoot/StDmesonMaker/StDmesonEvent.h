#ifndef StDmesonEvent__h
#define StDmesonEvent__h

class StPicoEvent;
class StPicoTrack;

#include "StDmesonTrack.h"
#include "StKaonPion.h"
#include "StKaonPionPion.h"
#include "StKaonKaon.h"
#include "StKaonKaonPion.h"
#include "StElectronPair.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "StThreeVectorF.hh"

class StDmesonEvent : public TObject {
  
 public:
  StDmesonEvent();
  virtual       ~StDmesonEvent() { clear("C"); }
  void          clear(const char *option ="");
  void addPicoEvent(const StPicoEvent& picoEvent);
  void          addTrack(StPicoTrack*,Short_t,UShort_t,UShort_t);
  void          addKaonPion(StKaonPion*);
  void          addKaonPionPion(StKaonPionPion*);
  void          addKaonKaon(StKaonKaon*);
  void          addKaonKaonPion(StKaonKaonPion*);
  void          addElectronPair(StElectronPair*);

  TClonesArray* tracks()      const { return mTracks; }
  TClonesArray* kaonPion()    const { return mKaonPion; }
  TClonesArray* kaonPionPion()    const { return mKaonPionPion; }
  TClonesArray* kaonKaon()    const { return mKaonKaon; }
  TClonesArray* kaonKaonPion()    const { return mKaonKaonPion; }
  TClonesArray* electronPair()    const { return mElectronPair;}

  Int_t        nTracks()    const { return mNTracks;}
  Int_t        nKaonPion()  const { return mNKaonPion;}
  Int_t        nKaonPionPion()  const { return mNKaonPionPion;}
  Int_t        nKaonKaon()  const { return mNKaonKaon;}
  Int_t        nKaonKaonPion()  const { return mNKaonKaonPion;}
  Int_t        nElectronPair()  const { return mNElectronPair;}

  // for variables from StPicoEvent
  Int_t    runId() const            { return mRunId; }
  Int_t    eventId() const          { return mEventId; }
  StThreeVectorF primaryVertex() const { return mPrimaryVertex; }
  Int_t    triggerWord() const      { return mTriggerWord; }
  Int_t    refMultPos() const       { return (Int_t)mRefMultPos; }
  Int_t    refMultNeg() const       { return (Int_t)mRefMultNeg; }
  Int_t    refMult() const          { return (Int_t)(mRefMultPos+mRefMultNeg); }
  		
  Float_t  vzVpd() const            { return (Float_t)mVzVpd/100.; }

  Float_t  ZDCx() const             { return mZDCx; }
  Float_t  BBCx() const             { return mBBCx; }
  
  UShort_t btofTrayMultiplicity() const { return mbTofTrayMultiplicity ; }
  UShort_t numberOfGlobalTracks() const { return mNumberOfGlobalTracks ; }

  Float_t ranking() const { return mRanking ; }
  UShort_t nBEMCMatch() const { return mNBEMCMatch ; }

  // triggers. See StPicoConstants.cxx
  Bool_t isMB() { return mTriggerWord & 0x2007; }
  Bool_t isVPDMB() { return mTriggerWord & 0x3; }
  Bool_t isBBCMB() { return mTriggerWord & 0x4; }
  Bool_t isBBCMB_TOF0() { return mTriggerWord & 0x2000; }

  Bool_t isHT()      { return mTriggerWord & 0x1f8;}
  Bool_t isHT0()     { return mTriggerWord & 0x48;}
  Bool_t isHT0_VPDMB()     { return mTriggerWord & 0x8;}
  Bool_t isHT0_BBCMB_TOF0()     { return mTriggerWord & 0x40;}
  Bool_t isHT1()     { return mTriggerWord & 0x90;}
  Bool_t isHT1_VPDMB()     { return mTriggerWord & 0x10;}
  Bool_t isHT1_BBCMB_TOF0()     { return mTriggerWord & 0x80;}
  Bool_t isHT2()     { return mTriggerWord & 0x120;}
  Bool_t isHT2_bare()     { return mTriggerWord & 0x20;}
  Bool_t isHT2_BBCMB()     { return mTriggerWord & 0x100;}

  private:

  // variables from StPicoEvent
  Int_t          mRunId;           // run number
  Int_t          mEventId;         // event number
  StThreeVectorF mPrimaryVertex;   // primary Vertex (1st)
  UShort_t       mTriggerWord;     // self-defined trigger word - see code for details
  UShort_t       mRefMultNeg;      // TPC refMult neg
  UShort_t       mRefMultPos;      // TPC refMult pos
  Short_t        mVzVpd;           // VzVpd*100.
  UShort_t        mZDCx;           // zdcX
  UShort_t        mBBCx;
  UShort_t mbTofTrayMultiplicity ; // BTOF tray multiplicity
  UShort_t mNumberOfGlobalTracks ; // # of global tracks
  Float_t mRanking ;
  UShort_t mNBEMCMatch ;

  //....
  Int_t         mNTracks;
  Int_t         mNKaonPion;
  Int_t         mNKaonPionPion;
  Int_t         mNKaonKaon;
  Int_t         mNKaonKaonPion;
  Int_t         mNElectronPair;

  TClonesArray*        mTracks;
  static TClonesArray* fgTracks;
  TClonesArray*        mKaonPion;
  static TClonesArray* fgKaonPion;
  TClonesArray*        mKaonPionPion;
  static TClonesArray* fgKaonPionPion;

  TClonesArray*        mKaonKaon;
  static TClonesArray* fgKaonKaon;
  TClonesArray*        mKaonKaonPion;
  static TClonesArray* fgKaonKaonPion;

  TClonesArray*        mElectronPair;
  static TClonesArray* fgElectronPair;

  ClassDef(StDmesonEvent,1)
};

#endif
