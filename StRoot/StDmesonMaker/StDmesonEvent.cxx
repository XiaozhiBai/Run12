#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "TMath.h"

#include "StDmesonEvent.h"
#include "StDmesonTrack.h"
#include "StKaonPion.h"
#include "StKaonPionPion.h"

ClassImp(StDmesonEvent)

TClonesArray *StDmesonEvent::fgTracks = 0;
TClonesArray *StDmesonEvent::fgKaonPion = 0;
TClonesArray *StDmesonEvent::fgKaonPionPion = 0;
TClonesArray *StDmesonEvent::fgKaonKaon = 0;
TClonesArray *StDmesonEvent::fgKaonKaonPion = 0;
TClonesArray *StDmesonEvent::fgElectronPair = 0;

//-----------------------------------------------------------------------
StDmesonEvent::StDmesonEvent()
{
  if (!fgTracks) fgTracks = new TClonesArray("StDmesonTrack");
  mTracks = fgTracks; mNTracks = 0;

  if (!fgKaonPion) fgKaonPion = new TClonesArray("StKaonPion");
  mKaonPion = fgKaonPion; mNKaonPion = 0;

  if (!fgKaonPionPion) fgKaonPionPion = new TClonesArray("StKaonPionPion");
  mKaonPionPion = fgKaonPionPion; mNKaonPionPion = 0;

  if (!fgKaonKaon) fgKaonKaon = new TClonesArray("StKaonKaon");
  mKaonKaon = fgKaonKaon; mNKaonKaon = 0;

  if (!fgKaonKaonPion) fgKaonKaonPion = new TClonesArray("StKaonKaonPion");
  mKaonKaonPion = fgKaonKaonPion; mNKaonKaonPion = 0;

  if (!fgElectronPair) fgElectronPair = new TClonesArray("StElectronPair");
  mElectronPair = fgElectronPair; mNElectronPair = 0;
}

//-----------------------------------------------------------------------
void StDmesonEvent::addPicoEvent(const StPicoEvent& picoEvent)
{
  // StPicoEvent variables
  mRunId = picoEvent.runId();
  mEventId = picoEvent.eventId();

  mPrimaryVertex = picoEvent.primaryVertex();
  mTriggerWord = picoEvent.triggerWord();
  mRefMultNeg = picoEvent.refMultNeg();
  mRefMultPos = picoEvent.refMultPos();
  mVzVpd = TMath::Nint(picoEvent.vzVpd() * 100.);
  mZDCx = picoEvent.ZDCx();
  mBBCx = picoEvent.BBCx();
  mbTofTrayMultiplicity = picoEvent.btofTrayMultiplicity();
  mNumberOfGlobalTracks = picoEvent.numberOfGlobalTracks();
  mRanking = picoEvent.ranking();
  mNBEMCMatch = picoEvent.nBEMCMatch();
}

//-----------------------------------------------------------------------
void StDmesonEvent::clear(const char *option)
{
  mTracks->Clear(option);
  mNTracks = 0;
  mKaonPion->Clear(option);
  mNKaonPion = 0;
  mKaonPionPion->Clear(option);
  mNKaonPionPion = 0;
  mKaonKaon->Clear(option);
  mNKaonKaon = 0;
  mKaonKaonPion->Clear(option);
  mNKaonKaonPion = 0;
  mElectronPair->Clear(option);
  mNElectronPair = 0;

  mPrimaryVertex.set(-999.,-999.,-999.);
}
//---------------------------------------------------------------------
void StDmesonEvent::addTrack(StPicoTrack* t,Short_t trg_tow,UShort_t dsm_adc,UShort_t adc) {
  TClonesArray &tracks = *mTracks;
  new(tracks[mNTracks++]) StDmesonTrack(t,trg_tow,dsm_adc,adc);
}
//---------------------------------------------------------------------
void StDmesonEvent::addKaonPion(StKaonPion* t) {
  TClonesArray &KaonPion = *mKaonPion;
  new(KaonPion[mNKaonPion++]) StKaonPion(t);
}
//---------------------------------------------------------------------
void StDmesonEvent::addKaonPionPion(StKaonPionPion* t) {
  TClonesArray &KaonPionPion = *mKaonPionPion;
  new(KaonPionPion[mNKaonPionPion++]) StKaonPionPion(t);
}
//---------------------------------------------------------------------
void StDmesonEvent::addKaonKaon(StKaonKaon* t) {
  TClonesArray &KaonKaon = *mKaonKaon;
  new(KaonKaon[mNKaonKaon++]) StKaonKaon(t);
}
//---------------------------------------------------------------------
void StDmesonEvent::addKaonKaonPion(StKaonKaonPion* t) {
  TClonesArray &KaonKaonPion = *mKaonKaonPion;
  new(KaonKaonPion[mNKaonKaonPion++]) StKaonKaonPion(t);
}
//---------------------------------------------------------------------
void StDmesonEvent::addElectronPair(StElectronPair* t) {
  TClonesArray &ElectronPair = *mElectronPair;
  new(ElectronPair[mNElectronPair++]) StElectronPair(t);
}
