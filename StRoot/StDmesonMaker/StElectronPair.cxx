#ifdef __ROOT__
#include "StElectronPair.h"

#include "TMath.h"

#include "StLorentzVectorD.hh"
#include "StThreeVectorD.hh"
#include "StPhysicalHelixD.hh"
#include "SystemOfUnits.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoConstants.h"

ClassImp(StElectronPair)

const Double_t ELECMASS=0.00051099906;

////////////////////////////////////////////////////////////////////
StElectronPair::StElectronPair() {
  clear();
}
////////////////////////////////////////////////////////////////////
StElectronPair::StElectronPair(StElectronPair* t) {
  clear();
  mElectronId = t->mElectronId;
  mPartnerId = t->mPartnerId;
  mM    = t->mM;
  mPairDca = t->mPairDca;
}

////////////////////////////////////////////////////////////////////
StElectronPair::StElectronPair(StPicoTrack* e1Trk, StPicoTrack* e2Trk,Int_t eId,Int_t pId,Float_t B)
{
  clear();
  if (!e1Trk || !e2Trk) return;

  StPhysicalHelixD e1Helix(e1Trk->gMom(), e1Trk->origin(), B*kilogauss, e1Trk->charge());
  StPhysicalHelixD e2Helix(e2Trk->gMom(), e2Trk->origin(), B*kilogauss, e2Trk->charge());

  pair<double,double> s = e1Helix.pathLengths(e2Helix);
  StThreeVectorD e1PosAtDca = e1Helix.at(s.first);
  StThreeVectorD e2PosAtDca = e2Helix.at(s.second);
  Double_t pairDCA = (e1PosAtDca - e2PosAtDca).mag();//pairDCA

  StThreeVectorD p1DCA=e1Helix.momentumAt(s.first,B*kilogauss);
  StThreeVectorD p2DCA=e2Helix.momentumAt(s.second,B*kilogauss);

  StLorentzVectorD LV1,LV2;
  LV1.setPx(p1DCA.x());
  LV1.setPy(p1DCA.y());
  LV1.setPz(p1DCA.z());
  LV1.setE(p1DCA.massHypothesis(ELECMASS));

  LV2.setPx(p2DCA.x());
  LV2.setPy(p2DCA.y());
  LV2.setPz(p2DCA.z());
  LV2.setE(p2DCA.massHypothesis(ELECMASS));

  Double_t massDCA=(LV1+LV2).m();

  mElectronId = (UShort_t)eId;
  mPartnerId = (UShort_t)pId;
  mM = (100000.*massDCA)>Pico::USHORTMAX ? Pico::USHORTMAX : (UShort_t)(TMath::Nint(100000.*massDCA));
  mPairDca = (100.*pairDCA)>Pico::CHARMAX ? Pico::CHARMAX : (Char_t)(e2Trk->charge()*TMath::Nint(100.*pairDCA));
}
////////////////////////////////////////////////////////////////////
StElectronPair::~StElectronPair() {
}

////////////////////////////////////////////////////////////////////
void StElectronPair::clear() {

  mElectronId=Pico::USHORTMAX;

  mM=Pico::USHORTMAX;
  mPairDca = Pico::CHARMAX;
};
#endif // __ROOT__
