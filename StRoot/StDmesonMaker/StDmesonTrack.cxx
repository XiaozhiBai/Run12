#include "StDmesonTrack.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoConstants.h"
#include "TMath.h"

ClassImp(StDmesonTrack)

//----------------------------------------------------------------------------------
StDmesonTrack::StDmesonTrack()
{
  Clear();
}

StDmesonTrack::StDmesonTrack(StPicoTrack* t,Short_t trgTowId, UShort_t trgTowDsmAdc,UShort_t trgTowAdc)
{
  mPMomentum = t->pMom();
  mGMomentum = t->gMom();
  mGDca      = (UShort_t)(t->dca()*1000.);
  mDedx      = (UShort_t)(t->dEdx()*1000.);
  int q      = t->charge();
  mNHitsFit  = (Char_t)(t->nHitsFit()*q);
  mNHitsMax  = (UChar_t)(t->nHitsMax());
  mNHitsDedx = (UChar_t)(t->nHitsDedx());
  mNSigmaPion     = (Short_t)(t->nSigmaPion()*100.);
  mNSigmaKaon     = (Short_t)(t->nSigmaKaon()*100.);
  mNSigmaProton   = (Short_t)(t->nSigmaProton()*100.);
  mNSigmaElectron = (Short_t)(t->nSigmaElectron()*100.);
  
  Float_t r = sqrt(t->origin().x()*t->origin().x()+t->origin().y()*t->origin().y());
  mFirstTpcPointR = (100.*r)>Pico::USHORTMAX ? 0 : (UShort_t)(100.*r);

  mBTofMatchFlag = (UChar_t)t->btofMatchFlag();
  mBTofBeta = (UShort_t)(t->btofBeta()*20000.);
  mBTofYLocal = (Short_t)(t->btofYLocal()*1000.);

  mBTOWADC0     = (Short_t)t->adc0();
  mBTOWE0       = (Short_t)(t->e0()*1000.);
  mBTOWE        = (Short_t)(t->e()*1000.);
  mBEMCDistZ    = (Short_t)(t->zDist()*100.0);
  mBEMCDistPhi  =  (Short_t)(t->phiDist()*10000.);
  mBSMDNEta = (Char_t)(t->nEta());
  mBSMDNPhi = (Char_t)(t->nPhi());

  mBTOWId       = (Short_t)t->btowId();

  mTrgTowId = trgTowId;
  mTrgTowDSMADC = trgTowDsmAdc;
  mTrgTowADC = trgTowAdc;
}

//----------------------------------------------------------------------------------
StDmesonTrack::~StDmesonTrack()
{ /* noop */ }

//----------------------------------------------------------------------------------
void StDmesonTrack::Clear(const Option_t* opt)
{
  mPMomentum.set(0.,0.,0.);
  mGDca = Pico::USHORTMAX;
  mDedx = 0;
  mNHitsFit  = 0;
  mNHitsMax  = 0;
  mNHitsDedx = 0;
  mNSigmaPion     = Pico::SHORTMAX;
  mNSigmaKaon     = Pico::SHORTMAX;
  mNSigmaProton   = Pico::SHORTMAX;
  mNSigmaElectron = Pico::SHORTMAX;
  mFirstTpcPointR = 0;

  mBTofMatchFlag = 0;
  mBTofBeta    = 0;
  mBTofYLocal  = Pico::SHORTMAX;

  mBTOWADC0     = 0;
  mBTOWE0       = -1*Pico::SHORTMAX;
  mBTOWE        = -1*Pico::SHORTMAX;

  mBTOWId       = -1;
  mTrgTowDSMADC = 0;
  mTrgTowADC = 0;
}
//----------------------------------------------------------------------------------
/*void StDmesonTrack::Print(const Char_t *option) const {
  if(strcmp(option,"tpc")==0 || strcmp(option,"")==0) {
    LOG_INFO << "id=" << id() << " flowflag=" << flowFlag() << " chi2=" << chi2() << " dca=" << dca() << endm;
    LOG_INFO << "gMom=" << gMom() << endm;
    LOG_INFO << "pMom=" << pMom() << endm;
    LOG_INFO << "Origin=" << origin() << endm;
//    LOG_INFO << "Q vector=" << mQXi << " " << mQYi << endm;
    LOG_INFO << " nHitsFit = " << nHitsFit() << " nHitsdEdx = " << nHitsDedx() << endm;
    LOG_INFO << " nSigma Pi/K/P/E = " << nSigmaPion() << "/" << nSigmaKaon() << "/" << nSigmaProton() << "/" << nSigmaElectron() << endm;
  }
  if(strcmp(option,"tof")==0 || strcmp(option,"")==0) {
    LOG_INFO << " BTOF cellId = " << btofCellId() << " tof = " << btof() << " beta = " << btofBeta() << endm;
    LOG_INFO << " BTOF match = " << btofMatchFlag() << " yLocal/zLocal " << btofYLocal() << " " << btofZLocal() << endm;
  }
  if(strcmp(option,"emc")==0 || strcmp(option,"")==0) {
    LOG_INFO << " BEMC Id = " << bemcId() << " BTOW ADC0 = " << adc0() << " energy0 = " << e0() << " e = " << e() << endm;
    LOG_INFO << " BEMC distz = " << zDist() << " distphi = " << phiDist() << endm;
    LOG_INFO << " BSMD nEta/nPhi = " << nEta() << "/" << nPhi() << endm;
    LOG_INFO << " BTOW Id = " << btowId() << " tower Id 2/3 = " << btowId2() << " " << btowId3() << endm;
    LOG_INFO << " BTOW energy = " << e1() << " " << e2() << " " << e3() << endm;
    LOG_INFO << " BTOW position to center = " << etaTowDist() << " " << phiTowDist() << endm;
  }
}*/
