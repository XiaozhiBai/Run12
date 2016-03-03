#ifdef __ROOT__
#include "StKaonKaonPion.h"
#include "StKaonKaon.h"

#include "TMath.h"

#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoConstants.h"
#include "StKaonPion.h"

ClassImp(StKaonKaonPion)
  


////////////////////////////////////////////////////////////////////
StKaonKaonPion::StKaonKaonPion() {
  clear();
}
////////////////////////////////////////////////////////////////////
StKaonKaonPion::StKaonKaonPion(StKaonKaonPion* t)
{
  clear();

  mKaonKaonId = t->mKaonKaonId;
  mPionId = t->mPionId;
  mM    = t->mM;
  mPt = t->mPt;
  mEta = t->mEta;
  mPhi = t->mPhi;
  mCosThetaStarPion = t->mCosThetaStarPion;
}

////////////////////////////////////////////////////////////////////
StKaonKaonPion::StKaonKaonPion(StKaonKaon* kaonKaon, StPicoTrack* pion, Int_t kkId,Int_t pId)
{
  clear();
  if (!kaonKaon || !pion) return;

  StLorentzVectorF kkFourMom = kaonKaon->lorentzVector();

  StLorentzVectorF pFourMom(pion->pMom(),pion->pMom().massHypothesis(.13957));
  StLorentzVectorF kkpFourMom = kkFourMom + pFourMom;

  StLorentzVectorF kkpFourMomReverse(-kkpFourMom.px(), -kkpFourMom.py(), -kkpFourMom.pz(), kkpFourMom.e());
  StLorentzVectorF pMomStar= pFourMom.boost(kkpFourMomReverse);

  mKaonKaonId = kkId;
  mPionId = pId;

  mM = (10000.*kkpFourMom.m())>Pico::USHORTMAX ? Pico::USHORTMAX : (UShort_t)(TMath::Nint(kkpFourMom.m()*10000.));
  mPt = (1000.*kkpFourMom.perp())>Pico::USHORTMAX ? Pico::USHORTMAX : (UShort_t)(TMath::Nint(kkpFourMom.perp()*1000.));
  mEta = (Short_t)(TMath::Nint(kkpFourMom.pseudoRapidity()*1000.));
  mPhi = (Short_t)(TMath::Nint(kkpFourMom.phi()*1000.));
  mCosThetaStarPion = (Char_t)(TMath::Nint(TMath::Cos(pMomStar.vect().angle(kkpFourMom.vect())) * 100. ));
}
//////////////////////////
StLorentzVectorF StKaonKaonPion::lorentzVector()
{
  Float_t px = pt()*TMath::Cos(phi());
  Float_t py = pt()*TMath::Sin(phi());
  Float_t pz = pt()*TMath::SinH(eta());
  Float_t e = TMath::Sqrt(px*px+py*py+pz*pz+m()*m());

  StLorentzVectorF kkpFourMom(px,py,pz,e);

  return kkpFourMom;
}
////////////////////////////////////////////////////////////////////
StKaonKaonPion::~StKaonKaonPion() {
}

////////////////////////////////////////////////////////////////////
void StKaonKaonPion::clear() 
{
  mKaonKaonId=Pico::USHORTMAX;
  mPionId=Pico::USHORTMAX;

  mM=Pico::USHORTMAX;
  mPt=0;
  mEta=0.;
  mPhi=0.;
  mCosThetaStarPion=0.;
};
#endif // __ROOT__
