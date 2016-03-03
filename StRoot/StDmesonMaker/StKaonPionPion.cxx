#ifdef __ROOT__
#include "StKaonPionPion.h"

#include "TMath.h"

#include "StLorentzVectorF.hh"
#include "StThreeVectorF.hh"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoConstants.h"
#include "StKaonPion.h"

ClassImp(StKaonPionPion)
  


////////////////////////////////////////////////////////////////////
StKaonPionPion::StKaonPionPion() {
  clear();
}
////////////////////////////////////////////////////////////////////
StKaonPionPion::StKaonPionPion(StKaonPionPion* t)
{
  clear();
  mKaonPionId = t->mKaonPionId;
  mPionId = t->mPionId;
  mM    = t->mM;
  mPt = t->mPt;
  mEta = t->mEta;
  mPhi = t->mPhi;
}

////////////////////////////////////////////////////////////////////
StKaonPionPion::StKaonPionPion(StKaonPion* kaonPion, StPicoTrack* pion, Int_t kpId, Int_t pId)
{
  clear();
  if (!kaonPion || !pion) return;

  StLorentzVectorF kpFourMom = kaonPion->lorentzVector();

  //cout<<kpFourMom.m()<<endl;
  //cout<<kpFourMom.perp()<<endl;
  //cout<<kpFourMom.pseudoRapidity()<<endl;
  //cout<<kpFourMom.phi()<<endl;

  StLorentzVectorF pFourMom(pion->pMom().x(),pion->pMom().y(),pion->pMom().z(),pion->pMom().massHypothesis(.13957));
  StLorentzVectorF kppFourMom = kpFourMom + pFourMom;

  mKaonPionId = (UShort_t)kpId;
  mPionId = (UShort_t)pId;

  mM = (10000.*kppFourMom.m())>Pico::USHORTMAX ? Pico::USHORTMAX : (UShort_t)(TMath::Nint(10000.*kppFourMom.m()));
  mPt = (1000.*kppFourMom.perp())>Pico::USHORTMAX ? Pico::USHORTMAX : (UShort_t)(TMath::Nint(1000.*kppFourMom.perp()));
  mEta = (Short_t)(TMath::Nint(kppFourMom.pseudoRapidity()*1000.));
  mPhi = (Short_t)(TMath::Nint(kppFourMom.phi()*1000.));
}
//////////////////////////
StLorentzVectorF StKaonPionPion::lorentzVector()
{
  Float_t px = pt()*TMath::Cos(phi());
  Float_t py = pt()*TMath::Sin(phi());
  Float_t pz = pt()*TMath::SinH(eta());
  Float_t e = TMath::Sqrt(px*px+py*py+pz*pz+m()*m());

  StLorentzVectorF kppFourMom(px,py,pz,e);

  return kppFourMom;
}
////////////////////////////////////////////////////////////////////
StKaonPionPion::~StKaonPionPion() {
}

////////////////////////////////////////////////////////////////////
void StKaonPionPion::clear() {

  mKaonPionId=Pico::USHORTMAX;
  mPionId=Pico::USHORTMAX;

  mM=Pico::USHORTMAX;
  mPt=0;
  mEta=0.;
  mPhi=0.;
};
#endif // __ROOT__
