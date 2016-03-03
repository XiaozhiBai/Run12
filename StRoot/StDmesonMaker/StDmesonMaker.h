#ifndef StDmesonMaker_h
#define StDmesonMaker_h

#include "StMaker.h"
#include "StDmesonEvent.h"

class StEmcGeom;
class StPicoDst;
class StPicoDstMaker;
class StPicoTrack;
class TString;
class TH1F;
class TFile;

class StDmesonMaker : public StMaker {
  public:
    StDmesonMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StDmesonMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
  private:
    Bool_t  isGoodEvent();
    Bool_t  isGoodTrack(StPicoTrack*);
    Bool_t  isPion(StPicoTrack*);
    Bool_t  isKaon(StPicoTrack*);
    Bool_t  isElectron(StPicoTrack*);
    Bool_t  isPartnerElectron(StPicoTrack*);

    void getTrgBTOW(StPicoTrack* trk, Short_t& f_tow_id,UShort_t& f_DSM_adc, UShort_t& f_adc);
    
  private:
    StPicoDstMaker* mPicoDstMaker;
    StPicoDst*      mPicoDst;
    StPicoEvent*    mPicoEvent;
    StEmcGeom*      mEmcGeom;
    
    TFile* mOutputFile;
    TTree* mTree;
    StDmesonEvent* mDmesonEvent;

    ClassDef(StDmesonMaker, 1)
};

#endif
