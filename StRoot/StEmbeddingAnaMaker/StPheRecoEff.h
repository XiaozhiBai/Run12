/////////////////////////////////////////////////////////////////////
// to calculation photonic electron reconstruction efficiency      //
/////////////////////////////////////////////////////////////////////

#ifndef StPheRecoEff_H
#define StPheRecoEff_H


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

class StEmcRawHit;
class StEmcCluster;
class StSPtrVecEmcCluster;
class StEmcCollection;
class StDedxPidTraits;

class ana_maker;


#include "StMaker.h"
#include "StThreeVectorD.hh"

class StPheRecoEff : public StMaker
{
  private:
      float t1rp;
      float t1rpt;
      float t1rpx;
      float t1rpy;
      float t1rpz;
      float t1reta;
      float t1rphi;
      float t1nfit;
      float t1ndedx;
      float t1nmax;
      float t1ncom;
      float t1key;
      float t1dca;
      float t1chi2;
      float t2rp;
      float t2rpt;
      float t2rpx;
      float t2rpy;
      float t2rpz;
      float t2reta;
      float t2rphi;
      float t2nfit;
      float t2nmax;
      float t2ncom;
      float t2key;
      float t2dca;
      float t2chi2;
      float p1bemcId;
      float p1btowAdc0;
      float p1btowE0;
      float p1btowE;
      float p1bemcDistZ;
      float p1bemcDistPhi;
      float p1bsmdNEta;
      float p1bsmdNPhi;
      float p1btowId;

      float p2bemcId;
      float p2btowAdc0;
      float p2btowE0;
      float p2btowE;
      float p2bemcDistZ;
      float p2bemcDistPhi;
      float p2bsmdNEta;
      float p2bsmdNPhi;
      float p2btowId;
      float t1pr_rp;
      float t1pr_rpT;
      float t2pr_rp;
      float t2pr_rpT;
      float t1_stpcx; //.. first TPC point x
      float t1_stpcy; //.. first TPC point y
      float t1_stpcz; //.. first TPC point z

      float pairDCA;
      float thetaDCA;
      float massPair;
      float mcMassPair;
      float pairPT;
      float dist;

      unsigned int mRefMultPos;
      unsigned int mRefMultNeg;

  private:
    int parent; 
    int collisionEnergy;
    /* std::vector<unsigned int> interesting_triggers; */
    /* std::vector<float> firedTriggersIndices; */

    char* outfile;
    TFile* mFile;

    TNtuple* nt_sngl; //.. single track
    TNtuple* nt_pair; //.. pair
    TNtuple* nt_mc; //.. this should contain mc gamma, pi0 and eta (not only embedded)
    TNtuple* mEventCount; //.. For counting purposes

    StAssociationMaker* mAssoc;
    const StTrack* findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits);

    bool is_correct_decay(StMcTrack* t);

    bool pass_trigger(StEvent* evt);
    int fill_eventCount(StMcEvent* mcEvt, StEvent* evt,float nMcTracks=-1);
    int fill_mTuple_MC(StMcEvent* mcev,int& nMcTracks);

    //.. emc ..
    StEmcCollection* mEmcCollection;
    StEmcPosition* mEmcPosition;
    StEmcGeom*       mEmcGeom[4];
    StEmcRawHit*     mEmcIndex[4800];
    double mField; //.. magnetic field

    void initEmc();
    void buildEmcIndex();
    bool getBEMC(const StTrack *t, int *id, int *adc, float *ene, float *d, int *nep, int *towid);

    StDedxPidTraits* findDedxPidTraits(const StTrack* track);

    void Reset_t1();
    void Reset_t2();

  public:
    StPheRecoEff(char* in,int parentGeantId,const int energy);
    // , const std::vector<unsigned int> triggers);
    int Init();
    int Make();
    int Finish();

  ClassDef(StPheRecoEff, 1)
};

#endif
