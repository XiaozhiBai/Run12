#include <vector>
#include "StDmesonMaker.h"
#include "StCuts.h"

// EMC
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"

#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoDstMaker/StPicoTrigger.h"
#include "StRoot/StPicoDstMaker/StPicoBTOWHit.h"
#include "TTree.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"
#include "TH1F.h"
#include "TFile.h"

ClassImp(StDmesonMaker)

//-----------------------------------------------------------------------------
StDmesonMaker::StDmesonMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
   : StMaker(name)
{
   mPicoDstMaker = picoMaker;
   mPicoDst = 0;
   mPicoEvent = 0;
   mDmesonEvent = new StDmesonEvent();
   mEmcGeom = StEmcGeom::getEmcGeom(detname[0].Data());

   mOutputFile = new TFile(outName, "RECREATE");
   mOutputFile->SetCompressionLevel(1);
   int BufSize = (int)pow(2., 16.);
   int Split = 1;
   mTree = new TTree("T", "T", BufSize);
   mTree->SetAutoSave(1000000); // autosave every 10 Mbytes
   mTree->Branch("dEvent", "StDmesonEvent", &mDmesonEvent, BufSize, Split);
}

//-----------------------------------------------------------------------------
StDmesonMaker::~StDmesonMaker()
{
   /*  */
}

//-----------------------------------------------------------------------------
Int_t StDmesonMaker::Init()
{
   return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StDmesonMaker::Finish()
{
   mOutputFile->cd();
   return kStOK;
}
//-----------------------------------------------------------------------------
void StDmesonMaker::Clear(Option_t *opt)
{
   mDmesonEvent->clear();
}

//-----------------------------------------------------------------------------
Int_t StDmesonMaker::Make()
{
   if (!mPicoDstMaker)
   {
      LOG_WARN << " No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   mPicoDst = mPicoDstMaker->picoDst();
   if (!mPicoDst)
   {
      LOG_WARN << " No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if (isGoodEvent())
   {

      UInt_t nTracks = mPicoDst->numberOfTracks();

      std::vector<Int_t> idPicoDstKaons;
      std::vector<Int_t> idPicoDstPions;
      std::vector<Int_t> idDmesonKaons;
      std::vector<Int_t> idDmesonPions;

      std::vector<Int_t> idPrimaryElectrons;
      std::vector<Int_t> idGlobalElectrons;
      std::vector<Int_t> idPicoPrimaryElectrons;
      std::vector<Int_t> idPicoGlobalElectrons;

      Bool_t bKaon = 0;
      Bool_t bPion = 0;
      Bool_t bElectron = 0;
      Bool_t bGlobalElectron = 0;

      for (UInt_t iTrack = 0; iTrack < nTracks; iTrack++)
      {
         StPicoTrack* trk = mPicoDst->track(iTrack);

         if (!trk || !isGoodTrack(trk)) continue;

         bKaon = 0;
         bPion = 0;
	 bElectron = 0;
	 bGlobalElectron = 0;

	 Short_t  trg_tow_id = -1;
	 UShort_t dsm_adc = 0;
	 UShort_t adc = 0;
	 //if(mDmesonEvent->isHT()) getTrgBTOW(trk, trg_tow_id,dsm_adc, adc);
	 if(trk->gMom().perp()>1)  getTrgBTOW(trk, trg_tow_id,dsm_adc, adc);

	 if(isPartnerElectron(trk))
	 {
		 idPicoGlobalElectrons.push_back(iTrack);
		 idGlobalElectrons.push_back(mDmesonEvent->nTracks());
		 bGlobalElectron = kTRUE;
	 }

	 //if(isElectron(trk) && (mDmesonEvent->isVPDMB() || dsm_adc>0))
	 if(isElectron(trk))
	 {
		 idPicoPrimaryElectrons.push_back(iTrack);
		 idPrimaryElectrons.push_back(mDmesonEvent->nTracks());
		 bElectron = kTRUE;
	 }

         if (isPion(trk))
         {
            idPicoDstPions.push_back(iTrack);
            idDmesonPions.push_back(mDmesonEvent->nTracks());
	    bPion = kTRUE;
         }

         if (isKaon(trk))
         {
            idPicoDstKaons.push_back(iTrack);
            idDmesonKaons.push_back(mDmesonEvent->nTracks());
	    bKaon = kTRUE;
         }

	 // Shusu needs all hadrons for e-h correlation study
         if ((trk->gMom().perp()>cuts::hadronsGlobalPt && trk->dca()<cuts::hadronsDca)
			 || bKaon || bPion || bElectron || bGlobalElectron)
         {
            mDmesonEvent->addTrack(trk, trg_tow_id,dsm_adc,adc);
         }
      } // .. end tracks loop

      // make photonic electron
      for(UInt_t ie1 = 0; ie1 < idPicoPrimaryElectrons.size(); ie1++)
      {
	      StPicoTrack* e1Trk = (StPicoTrack*)mPicoDst->track(idPicoPrimaryElectrons[ie1]);

	      for(UInt_t ie2 = 0; ie2 < idPicoGlobalElectrons.size(); ie2++)
	      {
		      if(idPicoPrimaryElectrons[ie1] == idPicoGlobalElectrons[ie2]) continue;

		      StPicoTrack* e2Trk = (StPicoTrack*)mPicoDst->track(idPicoGlobalElectrons[ie2]);

		      StElectronPair* pair = new StElectronPair(e1Trk,e2Trk,idPrimaryElectrons[ie1],idGlobalElectrons[ie2],mPicoEvent->bField());

		      if(pair->m()<0.3 && pair->pairDca()<1.0)
		      {
			      mDmesonEvent->addElectronPair(pair);
		      }

		      delete pair;
	      }
      } // ... end electron pair loop

      for (UInt_t ik = 0; ik < idPicoDstKaons.size(); ik++)
      {
         StPicoTrack* kaon = mPicoDst->track(idPicoDstKaons[ik]);
         // make Kπ pairs
         for (UInt_t ip = 0; ip < idPicoDstPions.size(); ip++)
         {
            if (idPicoDstKaons[ik] == idPicoDstPions[ip]) continue;

            StPicoTrack* pion = mPicoDst->track(idPicoDstPions[ip]);
	    if(pion->pMom().perp()<=0.2 || fabs(pion->pMom().pseudoRapidity())>=cuts::pionEta) continue; // for the Kπ pair we need both to have pT>0.2 

            StKaonPion* kaonPion = new StKaonPion(kaon, pion, idDmesonKaons[ik], idDmesonPions[ip]);

            if (kaonPion->m() <= 0.48 || kaonPion->m() > 2.5)
            {
               delete kaonPion;
               continue;
            }

            if ((kaonPion->m() >= 1.7 && kaonPion->m() <= 2.1) && kaon->charge()*pion->charge() < 0)
            {
               // make Kππ
               for (UInt_t ip2 = ip + 1; ip2 < idPicoDstPions.size(); ip2++)
               {
                  StPicoTrack* pion2 = mPicoDst->track(idPicoDstPions[ip2]);
                  StKaonPionPion* kaonPionPion = new StKaonPionPion(kaonPion, pion2, mDmesonEvent->nKaonPion(), idDmesonPions[ip2]);

                  if (kaonPionPion->m() - kaonPion->m() < 0.16)
                  {
                     mDmesonEvent->addKaonPionPion(kaonPionPion);
                  }

                  delete kaonPionPion;
               }
            }

            mDmesonEvent->addKaonPion(kaonPion);
            delete kaonPion;
         } // .. end make Kπ pairs


         // make KK pairs
         for (UInt_t ik2 = ik + 1; ik2 < idPicoDstKaons.size(); ik2++)
         {
            StPicoTrack* kaon2 = mPicoDst->track(idPicoDstKaons[ik2]);
            StKaonKaon* kaonKaon = new StKaonKaon(kaon, kaon2, idDmesonKaons[ik], idDmesonKaons[ik2]);

            // store φ mass only
            if (kaonKaon->m() < 0.98 || kaonKaon->m() > 1.08)
            {
               delete kaonKaon;
               continue;
            }

	    // fill only for vpdMB events
	    if(!mDmesonEvent->isVPDMB()) continue;

            // make KKπ
            for (UInt_t ip = 0; ip < idPicoDstPions.size(); ip++)
            {
               if (idPicoDstPions[ip] == idPicoDstKaons[ik] || idPicoDstPions[ip] == idPicoDstKaons[ik2]) continue;
               StPicoTrack* pion = mPicoDst->track(idPicoDstPions[ip]);

               StKaonKaonPion* kaonKaonPion = new StKaonKaonPion(kaonKaon, pion, mDmesonEvent->nKaonKaon(), idDmesonPions[ip]);
               if (kaonKaonPion->m() > 1.6 && kaonKaonPion->m() < 2.5)
               {
                  mDmesonEvent->addKaonKaonPion(kaonKaonPion);
               }
               delete kaonKaonPion;
            }

            mDmesonEvent->addKaonKaon(kaonKaon);
            delete kaonKaon;
         }
      } // .. end of kaons loop


      idPicoDstKaons.clear();
      idPicoDstPions.clear();
      idDmesonKaons.clear();
      idDmesonPions.clear();
      idPrimaryElectrons.clear();
      idPicoPrimaryElectrons.clear();
      idPicoGlobalElectrons.clear();

   } //.. end of good event fill

   mTree->Fill();
   mDmesonEvent->clear("C");


   return kStOK;
}

//-----------------------------------------------------------------------------
Bool_t StDmesonMaker::isGoodEvent()
{
   mPicoEvent = mPicoDst->event();
   mDmesonEvent->addPicoEvent(*mPicoEvent);

   // these are the triggers we will use
   // reconstruct only (VPDMB || HT0 || HT1 || HT2)
   // see StPicoConstants.cxx for triggerWord definition
   if(!(mPicoEvent->triggerWord() & 0x1FB)) return kFALSE;

   // cuts
   Float_t vz = mPicoEvent->primaryVertex().z();
   if (!(mPicoEvent->ranking() > 0)) return kFALSE;
   if (fabs(vz) > cuts::vz) return kFALSE;

   // this cut is only for minBias events (i.e not HT events)
   if (!(mPicoEvent->triggerWord() & 0x1F8) && fabs(vz - mPicoEvent->vzVpd()) > cuts::vzVpdVz) return kFALSE;

   return kTRUE;
}
//-----------------------------------------------------------------------------
Bool_t StDmesonMaker::isGoodTrack(StPicoTrack* trk)
{
   if (trk->nHitsFit() >= cuts::nHitsFit &&
         (Float_t)trk->nHitsFit() / (Float_t)trk->nHitsMax() > cuts::nHitsFitnHitsMax) return kTRUE;
   
   return kFALSE;
}
//-----------------------------------------------------------------------------
Bool_t StDmesonMaker::isPion(StPicoTrack* trk)
{
	// no cut on Eta because the soft pion can have any eta
   if (trk->pMom().perp() >= cuts::pionPt && fabs(trk->nSigmaPion()) < cuts::nSigmaPion && trk->dca()<cuts::dca) return kTRUE;
   
   return kFALSE;
}
//-----------------------------------------------------------------------------
Bool_t StDmesonMaker::isKaon(StPicoTrack* trk)
{
   if (trk->pMom().perp() >= cuts::kaonPt  && fabs(trk->nSigmaKaon()) < cuts::nSigmaKaon 
		   && fabs(trk->pMom().pseudoRapidity()) < cuts::kaonEta && trk->dca()<cuts::dca) return kTRUE;

   return kFALSE;
}
//------------------------------------------------------------------
Bool_t StDmesonMaker::isElectron(StPicoTrack* trk)
{
   if (trk->gMom().perp() >= cuts::electronPt  && trk->dca()<cuts::dca && 
		   trk->nSigmaElectron()>cuts::nSigmaElectronLow && trk->nSigmaElectron()<cuts::nSigmaElectronHigh
		   && fabs(trk->gMom().pseudoRapidity())< cuts::electronEta
		   && (trk->gMom().perp()>1.0 || (trk->btofBeta()>0 && fabs(1/trk->btofBeta()-1)<cuts::electronsBeta))) return kTRUE;
   
   return kFALSE;
}
//------------------------------------------------------------------
Bool_t StDmesonMaker::isPartnerElectron(StPicoTrack* trk)
{
	return trk->gMom().perp()>cuts::partnerElecGlobalPt && fabs(trk->nSigmaElectron())<cuts::partnerElecNSigmaElectron;
}
//-----------------------------------------------------
void StDmesonMaker::getTrgBTOW(StPicoTrack* trk, Short_t& f_tow_id,UShort_t& f_DSM_adc, UShort_t& f_adc)
{
   f_tow_id = -1;
   f_DSM_adc = 0;
   f_adc = 0;

   Float_t f_dPhi = -999.;
   Float_t f_dEta = -999.;

   if (!(mPicoEvent->triggerWord() & 0x1f8)) return; // HT events only

   Int_t triggerTowerId = 0;
   for (UInt_t i = 0; i < mPicoDst->numberOfTriggers(); i++)
   {
      // triggers for pp200GeV Run12
      if ((mPicoEvent->triggerWord() & 0x48 && mPicoDst->trigger(i)->flag() & 0x1) ||
            (mPicoEvent->triggerWord() & 0x90 && mPicoDst->trigger(i)->flag() & 0x2) ||
            (mPicoEvent->triggerWord() & 0x120 && mPicoDst->trigger(i)->flag() & 0x4))
      {
         triggerTowerId = mPicoDst->trigger(i)->id();

         if (triggerTowerId >= 1 && triggerTowerId <= 4800 && trk->btowId() >= 1 && trk->btowId() <= 4800)
         {
            Float_t triggerTowerPhi, triggerTowerEta;
            Float_t trackProjectionPhi, trackProjectionEta;
            mEmcGeom->getEtaPhi(triggerTowerId, triggerTowerEta, triggerTowerPhi);
            mEmcGeom->getEtaPhi(trk->btowId(), trackProjectionEta, trackProjectionPhi);

            Float_t dPhi = triggerTowerPhi - trackProjectionPhi;
            while (dPhi > TMath::Pi()) dPhi -= 2 * TMath::Pi();
            while (dPhi < -TMath::Pi()) dPhi += 2 * TMath::Pi();
            Float_t dEta = triggerTowerEta - trackProjectionEta;

            if (fabs(dPhi) < 0.06 && fabs(dEta) < 0.06 && fabs(dPhi) < fabs(f_dPhi) && fabs(dEta) < fabs(f_dEta))
            {
               f_dPhi = dPhi;
               f_dEta = dEta;
               f_DSM_adc = mPicoDst->trigger(i)->adc();
	       f_tow_id = mPicoDst->trigger(i)->id();
            }
         }
      }
   }

   if(f_tow_id<=0) return; 

   for(unsigned int i=0; i<mPicoDst->numberOfBTOWHits(); i++)
   {
	   StPicoBTOWHit* hit = (StPicoBTOWHit*)mPicoDst->btowHit(i);

	   if(f_tow_id == hit->id())
	   {
		   f_adc = hit->adc();
		   break;
	   }
   }
}
