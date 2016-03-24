#ifndef ST_MCACCPMAKER_H
#define ST_MCACCPMAKER_H

#include <vector>

class TFile;
class TNtuple;

class StMcTrack;
class StTrack;
class StGlobalTrack;
class StAssociationMaker;
class StEvent;
class StMcEvent;
class StMcPixelHit;
class StMcTpcHit;
class StEmcCollection;
class StEmcPosition;
class StEmcGeom;
class StTriggerSimuMaker;

class StEmcRawHit;
class StEmcCluster;
class StSPtrVecEmcCluster;
class StEmcCollection;
class StDedxPidTraits;

#include "StMaker.h"
#include "StThreeVectorD.hh"


class StMcAccpMaker : public StMaker
{
  private:
    char* outfile;
    int collisionEnergy;
    /* std::vector<unsigned int> interesting_triggers; */
    /* std::vector<float> firedTriggersIndices; */

    TFile* mFile;

    TNtuple* mTupleMC; //.. RC is a subset of MC
    TNtuple* mEventCount; //.. For counting purposes

    StAssociationMaker* mAssoc;
    StTriggerSimuMaker* mTrgSimMaker;
    const StTrack* findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits);
    const StMcTrack* findPartner(StGlobalTrack* rcTrack, int& maxCommonTpcHits);
    StDedxPidTraits* findDedxPidTraits(const StTrack* track);

    StEmcCollection* mEmcCollection;
    StEmcPosition* mEmcPosition;
    StEmcGeom*       mEmcGeom[4];
    StEmcRawHit*     mEmcIndex[4800];
    double mField; //.. magnetic field

    void initEmc();
    void buildEmcIndex();
    bool getBEMC(const StTrack *t, int *id, int *adc, float *ene, int *id0,float *d, int *nep, int *towid);
    bool getBtowDSMAdc(int btowId,int* adc,int* id1);
    bool getBtowDSMAdc(int btowId,int btowId0,int* adc,int* id1);

    bool pass_trigger(StEvent* evt);
    int fill_eventCount(StMcEvent* mcEvt, StEvent* evt,float nRTracks=-1,float nMcTracks=-1);
    int fill_mTuple_MC(StMcEvent* mcev,int& nRTracks,int& nMcTracks);

    int get_ndedx(const StTrack* t);

  public:
    StMcAccpMaker(StTriggerSimuMaker* trgSim,char* in,int energy);
    // ,       const std::vector<unsigned int> triggers);
    int Init();
    int Make();
    int Finish();

  ClassDef(StMcAccpMaker, 1)
};

#endif
