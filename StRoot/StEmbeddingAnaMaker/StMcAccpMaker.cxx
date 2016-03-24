#include "TFile.h"
#include "TNtuple.h"
#include "TVector2.h"

#include "StParticleDefinition.hh"
#include "StTrack.h"
#include "StPrimaryTrack.h"
#include "StTrackGeometry.h"
#include "StMcEvent/StMcEventTypes.hh"
#include "StMcEvent/StMcPixelHit.hh"
#include "StMcEvent/StMcContainers.hh"
#include "SystemOfUnits.h"
#include "StEvent.h"
#include "StTrackNode.h"
#include "StGlobalTrack.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcPoint.h"
#include "StEmcCluster.h"
#include "StEmcCollection.h"

#include "StTriggerUtilities/StTriggerSimuMaker.h"
#include "StTriggerUtilities/Bemc/StBemcTriggerSimu.h"

#include "StMcEvent.hh"
#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"
#include "TVector3.h"
#include "StMcEvent/StMcTpcHit.hh"
#include "StMcAccpMaker.h"
#include "StMcTpcHitCollection.hh"

#include "StEvent/StEventTypes.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StEvent/StEventSummary.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTpcDedxPidAlgorithm.h"
#include "StarClassLibrary/StParticleTypes.hh"

#include "StDetectorDbMaker/StDetectorDbTriggerID.h"

const int nadc_max = 4;
std::vector<unsigned int> Interesting_triggers;
std::vector<float> FiredTriggersIndices;
ClassImp(StMcAccpMaker);

StMcAccpMaker::StMcAccpMaker(StTriggerSimuMaker* trgSim,char* in,int energy)
//  ,      const std::vector<unsigned int> triggers)
{
    mTrgSimMaker = trgSim;
    outfile = in;
    collisionEnergy = energy;
    // Interesting_triggers = triggers;

    Interesting_triggers.push_back(370542); //
    Interesting_triggers.push_back(370522); //
    for(int i=0;i<6;i++) {FiredTriggersIndices.push_back(-1);};
    
    cout<<"StMcAccpMaker::StMcAccpMaker - DONE"<<endl;
}
//__________________________________
int StMcAccpMaker::Init()
{
    initEmc();

    mFile = new TFile(outfile, "recreate");
    assert(mFile && !mFile->IsZombie());
        
    const char* varlist = "geantId:p:pt:svx:svy:svz:pgeantId:pp:ppt:phi:y:eta:peta:label:"
    "rp:rpt:reta:rphi:nfit:ncom:nmax:"
    "nDedxPts:dedx:dedx_2:nSigPi:nSigK:nSigP:nSigE:"
    "vx:vy:vz:"
    "bemcId:btowId0:btowAdc0:btowE0:btowE:bemcDistZ:bemcDistPhi:bsmdNEta:bsmdNPhi:btowId:"
    "btowDsmAdcProj:btowDsmAdcNext:btowIdNext:btowDsmAdc0:" // dsm adc of projection tower, dsm adc of "next" highest tower in the nearby 8 towers, 
                                                            // id of "next" tower, dsm ADC of btowId0 (highest tower in cluster)
    "dca:pr_rp:pr_rpT:stpcx:stpcy:stpcz";
    
    const char* evtVarList = "runId:eventId:mcVx:mcVy:mcVz:vx:vy:vz:vzVpd:posRefMult:negRefMult:"
    "zdc:bbc:nMcTracks:nRTracks:magField:"
    "t0:t1:t2:t3:t4:t5";

    mTupleMC= new TNtuple("tracksMC", "tracks", varlist);
    mEventCount = new TNtuple("eventCount", "eventCount", evtVarList);

    mAssoc = (StAssociationMaker*)GetMaker("StAssociationMaker");

    cout<<"StMcAccpMaker::Init - DONE"<<endl;
    return StMaker::Init();
}

//__________________________________
int StMcAccpMaker::Make()
{
    cout<<"StMcAccpMaker::Make() - call"<<endl;
    StMcEvent* mcEvent = (StMcEvent*)GetDataSet("StMcEvent");
    
    if (!mcEvent) 
    {
        LOG_WARN << "No StMcEvent" << endm;
        return kStWarn;
    }

    StEvent* event = (StEvent*)GetDataSet("StEvent");
    if (!event) {
        LOG_WARN << "No StEvent" << endm;
        return kStWarn;
    }
    
    cout<<"StMcAccpMaker::Make() : event: "<<event->id()<<endl;

    mEmcCollection = event->emcCollection();
    if(mEmcCollection) buildEmcIndex();

    mField = (float)event->summary()->magneticField();

    
    // Check if event fired any interesting triggers
    bool event_fired_Interesting_triggers = pass_trigger(event);
    
    // Fill
    int nRTracks = -1;
    int nMcTracks = -1;
    if(event_fired_Interesting_triggers)
    {
        cout<<"Good event."<<endl;
        fill_mTuple_MC(mcEvent,nRTracks,nMcTracks);
    }
    else
    {
        cout<<"No interesting triggers. Counting event then skipping."<<endl;
    }
    
    fill_eventCount(mcEvent,event,(float)nRTracks,(float)nMcTracks);

    return kStOk;
}
//____________________________________
int StMcAccpMaker::fill_eventCount(StMcEvent* mcEvt, StEvent* evt,float nRTracks,float nMcTracks)
{   
    cout<<"Taking event counts..."<<endl;
    float vars[30];
    
    float vpdVz = -999;
    StBTofHeader* tofheader = 0;
    if(evt->btofCollection())  tofheader = evt->btofCollection()->tofHeader();
    if(tofheader) vpdVz = tofheader->vpdVz();
    
    int i = 0;
    vars[i++] = (float)evt->runId();
    vars[i++] = (float)evt->id();
    vars[i++] = (float)mcEvt->primaryVertex()->position().x();
    vars[i++] = (float)mcEvt->primaryVertex()->position().y();
    vars[i++] = (float)mcEvt->primaryVertex()->position().z();
    vars[i++] = (float)evt->primaryVertex()->position().x();
    vars[i++] = (float)evt->primaryVertex()->position().y();
    vars[i++] = (float)evt->primaryVertex()->position().z();
    vars[i++] = vpdVz;
    vars[i++] = (float)uncorrectedNumberOfPositivePrimaries(*evt);
    vars[i++] = (float)uncorrectedNumberOfNegativePrimaries(*evt);
    vars[i++] = (float)evt->runInfo()->zdcCoincidenceRate();
    vars[i++] = (float)evt->runInfo()->bbcCoincidenceRate();
    vars[i++] = nMcTracks; 
    vars[i++] = nRTracks;
    vars[i++] = (float)evt->summary()->magneticField()/10;
    vars[i++] = FiredTriggersIndices.at(0);
    vars[i++] = FiredTriggersIndices.at(1);
    vars[i++] = FiredTriggersIndices.at(2);
    vars[i++] = FiredTriggersIndices.at(3);
    vars[i++] = FiredTriggersIndices.at(4);
    vars[i++] = FiredTriggersIndices.at(5);
    
    mEventCount->Fill(vars);
    
    return kStOk;
}
//____________________________________
bool StMcAccpMaker::pass_trigger(StEvent* evt)
{
    cout<<"Checking triggers..."<<endl;
    bool interesting_event = false;
    
    const StTriggerId* st_trgid = evt->triggerIdCollection()->nominal();
    
    for(unsigned int i=0; i<FiredTriggersIndices.size(); i++){FiredTriggersIndices[i]=-1;}
    
    // Fill interesting triggers
    cout<<"Interesting fired triggers: ";
    for(unsigned int i = 0; i<st_trgid->maxTriggerIds(); i++) 
    {
        unsigned int id = st_trgid->triggerId(i);
        
        int trgIndex = -1;
        
        for (unsigned int j = 0; j < Interesting_triggers.size(); j++) 
        {
            if (Interesting_triggers[j] == id) 
            {
                trgIndex = j;
                interesting_event = true;
                cout<<id<<" ";
                break;
            }
        }
        
        if(trgIndex!=-1) FiredTriggersIndices[trgIndex] = 1.0;
    }
    
    return interesting_event;
}
//____________________________________
int StMcAccpMaker::fill_mTuple_MC(StMcEvent* mcEvent,int& nRTracks, int& nMcTracks)
{
    StSPtrVecMcTrack& trks = mcEvent->tracks();
    cout<<"Filling "<<trks.size()<<" mcTracks..."<<endl;   
    
    nRTracks = 0;
    nMcTracks = 0;

    for (unsigned int i = 0;  i<trks.size(); i++) 
    {
        StMcTrack* t1 = trks[i];

        if(t1->geantId()!=2 && t1->geantId()!=3) continue;

        float t1svx, t1svy, t1svz;
        
        if(t1->startVertex()) 
        {
            t1svx = t1->startVertex()->position().x();
            t1svy = t1->startVertex()->position().y();
            t1svz = t1->startVertex()->position().z();
        }

        if(sqrt(t1svx*t1svx + t1svy*t1svy)>3) continue;

	++nMcTracks;

        float vtx = mcEvent->primaryVertex()->position().x();
        float vty = mcEvent->primaryVertex()->position().y();
        float vtz = mcEvent->primaryVertex()->position().z();

        int t1geantId = t1->geantId();

        float t1phi = t1->momentum().phi();
        float t1y = t1->rapidity();
        float t1eta = t1->pseudoRapidity();

        float t1p = t1->momentum().mag();

        float pgeantId = 0, pmom = 0, pt = 0,peta = 0, label = 0;
        
        if(t1->parent()) 
        {
            pgeantId = t1->parent()->geantId();
            pmom = t1->parent()->momentum().mag();
            pt= t1->parent()->pt();
            peta = t1->parent()->pseudoRapidity();
            label = t1->parent()->eventGenLabel();
        }

        float t1px = t1->momentum().x();
        float t1py = t1->momentum().y();
        //float t1pz = t1->momentum().z();

        float t1pt = sqrt(t1px*t1px+t1py*t1py);

        float rp = -9999;
        float rpt = -9999;
        float rpx = -9999;
        float rpy = -9999;
        float rpz = -9999;
        float reta = -9999;
        float rphi = -9999;
        float ncom = -9999;
        float key = -9999;

	int bemcId = -1; 
	int adc0 = -9999; 
	int btowId0 = -9999;
	float e[5]; 
	float dist[4]; 
	int nhit[2]; 
	int ntow[3];
	for(int i=0;i<5;i++) { e[i] = 0.; }
	for(int i=0;i<4;i++) { dist[i] = 1.e9; }
	for(int i=0;i<2;i++) { nhit[i] = 0; }
	for(int i=0;i<3;i++) { ntow[i] = -1; }

        int nfit = -9999;
        int nmax = -9999;
        int ncommonhits = -9999;
        
        float nDedxPts = -9999;
        float dedx = -9999;
	float dedx_2 = -9999;
        float nSigPi = -9999;
        float nSigK = -9999;
        float nSigP = -9999;
        float nSigE = -9999;
        
        float dca = -9999;
        float t1pr_rp = -9999;
        float t1pr_rpT = -9999;
        float t1_stpcx = -9999;
        float t1_stpcy = -9999;
        float t1_stpcz = -9999;

        const StTrack* rcTrack = findPartner(t1, ncommonhits);

        if(rcTrack) 
        {
            ++nRTracks;
            rp = rcTrack->geometry()->momentum().mag();
            rpt = rcTrack->geometry()->momentum().perp();
            rpx = rcTrack->geometry()->momentum().x();
            rpy = rcTrack->geometry()->momentum().y();
            rpz = rcTrack->geometry()->momentum().z();
            reta = rcTrack->geometry()->momentum().pseudoRapidity();
            rphi = rcTrack->geometry()->momentum().phi();
            nfit = rcTrack->fitTraits().numberOfFitPoints(kTpcId);
            nmax = rcTrack->numberOfPossiblePoints(kTpcId);
            ncom =ncommonhits; 
            key = rcTrack->key();
            
            // dedx info
            static StTpcDedxPidAlgorithm aplus;
            static StPionPlus* Pion = StPionPlus::instance();
            static StKaonPlus* Kaon = StKaonPlus::instance();
            static StProton* Proton = StProton::instance();
            static StElectron* Electron = StElectron::instance();
            const StParticleDefinition* prtcl = rcTrack->pidTraits(aplus);
            if (prtcl) 
            {
                nDedxPts = aplus.traits()->numberOfPoints();
                dedx = aplus.traits()->mean();
                nSigPi = aplus.numberOfSigma(Pion);
                nSigK = aplus.numberOfSigma(Kaon);
                nSigP = aplus.numberOfSigma(Proton);
                nSigE = aplus.numberOfSigma(Electron);
            }

	    StDedxPidTraits* tpc_pid = findDedxPidTraits(rcTrack);
	    if(tpc_pid) dedx_2 = tpc_pid->mean();

	    getBEMC(rcTrack, &bemcId, &adc0, e, &btowId0, dist, nhit, ntow);

            //...find the DCA of the track w.r.t. the zvtx 
            StThreeVectorF origin(vtx,vty,vtz);
            dca = rcTrack->geometry()->helix().distance(origin);

	    //.. check if this track is a primary track
	    StPrimaryTrack* t1_pr_trk = (StPrimaryTrack*) rcTrack->node()->track(primary);
	    if(t1_pr_trk) 
            {
		    t1pr_rp = t1_pr_trk->geometry()->momentum().mag();
		    t1pr_rpT = t1_pr_trk->geometry()->momentum().perp();
            }

	    t1_stpcx = rcTrack->detectorInfo()->firstPoint().x();
	    t1_stpcy = rcTrack->detectorInfo()->firstPoint().y();
	    t1_stpcz = rcTrack->detectorInfo()->firstPoint().z();

        }
        
	float btowId = (ntow[0]<=0 || ntow[0]>4800) ? -1 : (float)ntow[0];
	// dsm adc info
	int dsmAdc[3] = {-9999,-9999,-9999};
	int id1 = -9999;
	getBtowDSMAdc((int)btowId,btowId0,dsmAdc,&id1);

        float tuple[] = {t1geantId,t1p,t1pt,t1svx,t1svy,t1svz,pgeantId,pmom,pt,t1phi,t1y,t1eta,peta,label,
        rp, rpt, reta, rphi,nfit, ncom, nmax, 
        nDedxPts,dedx,dedx_2,nSigPi,nSigK,nSigP,nSigE,
        vtx, vty, vtz, 
	bemcId,btowId0,adc0,e[0],e[1],dist[0],dist[1],nhit[0],nhit[1],btowId,
	dsmAdc[0],dsmAdc[2],(float)id1,dsmAdc[1],
        dca, t1pr_rp, t1pr_rpT, t1_stpcx, t1_stpcy, t1_stpcz};

        mTupleMC->Fill(tuple);
    }

    return kStOk;
}


//________________________________________________
const StTrack* StMcAccpMaker::findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits) 
{//..StMcTrack find partner from the StTracks
    if(!mAssoc) {
        cout<<" empty StAssociateMaker, stop!!"<<endl;
        exit(1);
    }
    pair<mcTrackMapIter, mcTrackMapIter> p = mAssoc->mcTrackMap()->equal_range(mcTrack);

    const StTrack* maxTrack = 0;
    maxCommonTpcHits = 0;
    for (mcTrackMapIter k = p.first; k != p.second; ++k) {
        int commonTpcHits = k->second->commonTpcHits();
        const StTrack* track = k->second->partnerTrack()->node()->track(global);//should be global
        if (track && commonTpcHits > maxCommonTpcHits) {
            maxTrack = track;
            maxCommonTpcHits = commonTpcHits;
        }
    }
    return maxTrack;
}


//________________________________________________
const StMcTrack* StMcAccpMaker::findPartner(StGlobalTrack* rcTrack, int& maxCommonTpcHits) 
{//.. StGlobalTracks find partner from StMcTracks. 
    //.. See example from StRoot/StMcAnalysisMaker

    if(!mAssoc) {
        cout<<" empty StAssociateMaker, stop!!"<<endl;
        exit(1);
    }

    pair<rcTrackMapIter,rcTrackMapIter> p = mAssoc->rcTrackMap()->equal_range(rcTrack);

    const StMcTrack* maxTrack = 0;
    maxCommonTpcHits = 0;
    for (rcTrackMapIter k = p.first; k != p.second; ++k) {
        int commonTpcHits = k->second->commonTpcHits();

        const StMcTrack* track = k->second->partnerMcTrack();

        if (track && commonTpcHits > maxCommonTpcHits) {
            maxTrack = track;
            maxCommonTpcHits = commonTpcHits;
        }
    }
    return maxTrack;
}
//______________________________________________________________________________
StDedxPidTraits* StMcAccpMaker::findDedxPidTraits(const StTrack* track)
{
  StDedxPidTraits* pid=0;
  StPtrVecTrackPidTraits traits = track->pidTraits(kTpcId);
  
  for (UInt_t i = 0; i < traits.size(); i++) {
    pid = dynamic_cast<StDedxPidTraits*>(traits[i]);
    if (pid && pid->method() == kTruncatedMeanId) break;
  }
  return pid;
} 
//______________________________________________________
int StMcAccpMaker::Finish()
{
    mFile->Write();
    mFile->Close();
    return kStOk;
}
//___________________________________________
int StMcAccpMaker::get_ndedx(const StTrack* t)
{
    int ndedx = -1;
    StPtrVecTrackPidTraits pidTraits = t->pidTraits(kTpcId);   // t is defined as StPrimaryTrack* t
    if(pidTraits.size()){
        StDedxPidTraits* pid;
        for(unsigned int i = 0; i < pidTraits.size(); i++){
            pid = dynamic_cast<StDedxPidTraits*>(pidTraits[i]);
            if(pid && (pid->method() == kTruncatedMeanId)){
                //this->dEdx  = pid->mean();                        //dedx
                ndedx = pid->numberOfPoints();            //number of dedx hits
                break;
            }
        }
    }

    return ndedx;
}
//-----------------------------------------------------------------------
bool StMcAccpMaker::getBEMC(const StTrack *t, int *id, int *adc, float *ene, int* id0, float *d, int *nep, int *towid) {
  *id = -1; *adc = 0;
  for(int i=0;i<5;i++) { ene[i] = 0.; }
  for(int i=0;i<4;i++) { d[i] = 1.e9; }
  for(int i=0;i<2;i++) { nep[i] = 0; }
  for(int i=0;i<3;i++) { towid[i] = -1; }

  if(!mEmcCollection) {
    LOG_WARN << " No Emc Collection for this event " << endm;
    return kFALSE;
  }

  StThreeVectorD position, momentum;
  StThreeVectorD positionBSMDE, momentumBSMDE;
  StThreeVectorD positionBSMDP, momentumBSMDP;
  Double_t bFld = mField*kilogauss/tesla; // bFld in Tesla
  bool ok       = false;
  bool okBSMDE  = false;
  bool okBSMDP  = false;
  if(mEmcPosition) {
    ok      = mEmcPosition->projTrack(&position, &momentum, t, bFld, mEmcGeom[0]->Radius());
    okBSMDE = mEmcPosition->projTrack(&positionBSMDE, &momentumBSMDE, t, bFld, mEmcGeom[2]->Radius());
    okBSMDP = mEmcPosition->projTrack(&positionBSMDP, &momentumBSMDP, t, bFld, mEmcGeom[3]->Radius());
  }
//  if(!ok || !okBSMDE || !okBSMDP) {
  if(!ok) {
    LOG_WARN << " Projection failed for this track ... " << endm;
    return kFALSE;
  }

  if(ok && okBSMDE && okBSMDP) {

  Int_t mod, eta, sub;
  StSPtrVecEmcPoint& bEmcPoints = mEmcCollection->barrelPoints();
  int index=0;
  float mindist=1.e9;
  mEmcGeom[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub); //project on SMD plan
  for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it != bEmcPoints.end(); it++, index++) {
    bool associated=false;
    StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
    if(bEmcClusters.size()==0 ) continue;
    if(bEmcClusters[0]==NULL) continue;
    for(StPtrVecEmcClusterIterator cIter = bEmcClusters.begin(); cIter != bEmcClusters.end(); cIter++) {
      StPtrVecEmcRawHit& bEmcHits = (*cIter)->hit();
      for(StPtrVecEmcRawHitIterator hIter = bEmcHits.begin(); hIter != bEmcHits.end(); hIter++) {
        if(mod == (Int_t)(*hIter)->module() && eta == (Int_t)(*hIter)->eta() && sub == (Int_t)(*hIter)->sub()) {
          associated=true;
          break;
        }
      }
      if(associated) {
        for(StPtrVecEmcRawHitIterator hitit=bEmcHits.begin(); hitit!=bEmcHits.end();hitit++) {
          if((*hitit)->energy()>ene[0]) ene[0]=(*hitit)->energy();
          if((int)(*hitit)->adc()>(*adc)) 
	  {
		  *adc=(*hitit)->adc();
		  *id0=(*hitit)->softId(1);
	  }
        }
      }
    }

    StPtrVecEmcCluster& smdeClusters = (*it)->cluster(kBarrelSmdEtaStripId);
    StPtrVecEmcCluster& smdpClusters = (*it)->cluster(kBarrelSmdPhiStripId);

    if(associated) {
      *id = index;
      ene[1] = ene[1] + (*it)->energy(); //use point's energy, not tower cluster's energy

      float deltaphi=(*it)->position().phi()-positionBSMDP.phi();
      if(deltaphi>=TMath::Pi()) deltaphi=deltaphi-TMath::TwoPi();
      if(deltaphi<-TMath::Pi()) deltaphi=deltaphi+TMath::TwoPi();

      float rsmdp=mEmcGeom[3]->Radius();
      float pointz=(*it)->position().z();
      float deltaz=pointz-positionBSMDE.z();
      if(sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz)<mindist) {
        d[1]=deltaphi;
        d[0]=deltaz;
        if(smdeClusters.size()>=1) nep[0]=smdeClusters[0]->nHits();
        if(smdpClusters.size()>=1) nep[1]=smdpClusters[0]->nHits();
        mindist=sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz);
      }
    }//associated
  }

  } // end if (ok && okBSMDE && okBSMDP)

  //Get BEMC tower energy from matched tower + 2 nearest towers

  int towerId;
  int localTowerId = -1;
  int localId1 = -1;
  int localId2 = -1;
  double energy1 = 0, energy2 = 0;
  double energyTemp;
  double dist1 = 1000, dist2 = 1000;
  double distTemp;
  Float_t etaTemp, phiTemp;

  mEmcGeom[0]->getId(position.phi(),position.pseudoRapidity(),towerId);
  for(int ieta=-1;ieta<2;ieta++){
    for(int iphi=-1;iphi<2;iphi++){
      localTowerId++;//loops from 0 to 8
      int nextTowerId = mEmcPosition->getNextTowerId(towerId, ieta, iphi);
      if(nextTowerId < 1 || nextTowerId > 4800) continue;
      StEmcRawHit* emcHit = mEmcIndex[nextTowerId-1];
      if (emcHit==0) continue;
      if (emcHit->energy()<0.2) continue; // don't include any noise tower
      if(ieta==0&&iphi==0) {
        mEmcGeom[0]->getEta(nextTowerId, etaTemp);
        mEmcGeom[0]->getPhi(nextTowerId, phiTemp);
        ene[2] = emcHit->energy();
        d[2] = position.pseudoRapidity() - etaTemp; 
        d[3] = position.phi() - phiTemp; 
      }
      else {
        energyTemp = emcHit->energy();
        mEmcGeom[0]->getEta(nextTowerId, etaTemp);
        mEmcGeom[0]->getPhi(nextTowerId, phiTemp);
        distTemp = sqrt((etaTemp-position.pseudoRapidity())*(etaTemp-position.pseudoRapidity()) + (phiTemp-position.phi())*(phiTemp-position.phi()));
        if(distTemp < dist1) {
          dist2 = dist1;
          dist1 = distTemp;
	  energy2 = energy1;
          energy1 = energyTemp;
          localId1 = localTowerId;
        }
        else if(distTemp < dist2){
          dist2 = distTemp;
          energy2 = energyTemp;
          localId2 = localTowerId;
        }
      }
    }
  }
  towid[0] = towerId;
  ene[3] = energy1;//closest tower
  towid[1] = localId1;
  ene[4] = energy2;//2nd closest tower
  towid[2] = localId2;

  if(Debug()) {
//  if(1) {
    cout << " ====== BEMC results ====== " << endl;
    cout << " Energy = " << ene[0] << " " << ene[1] << " " << ene[2] << " " << ene[3] << " " << ene[4] << endl;
    cout << " BSMD = " << nep[0] << " " << nep[1] << endl;
    cout << " TowerId = " << towid[0] << " " << towid[1] << " " << towid[2] << endl;
  }

  return kTRUE;
}
//-----------------------------------------------------------------------
bool StMcAccpMaker::getBtowDSMAdc(int btowId,int btowId0,int* adc,int* id1)
{
	if(!mTrgSimMaker || !mEmcPosition) return kFALSE;

	adc[0] = -9999;
	adc[1] = -9999;
	adc[2] = -9999;
	*id1 = -9999;

	if(btowId>0 && btowId<=4800) adc[0] = mTrgSimMaker->bemc->barrelHighTowerAdc(btowId); // dsm adc of track projection tower
	if(btowId0>0 && btowId0<=4800) adc[1] = mTrgSimMaker->bemc->barrelHighTowerAdc(btowId0);

	if(btowId<0 || btowId>4800) return false;

	for(int ieta=-1;ieta<2;ieta++)
	{
		for(int iphi=-1;iphi<2;iphi++)
		{
			if(ieta==0 && iphi==0) continue;

			int nextTowerId = mEmcPosition->getNextTowerId(btowId, ieta, iphi);
			if(nextTowerId < 1 || nextTowerId > 4800) continue;
			int tmpAdc = mTrgSimMaker->bemc->barrelHighTowerAdc(nextTowerId);

			if(tmpAdc >adc[2])
			{
				// adc and id of the highest tower in the nearby eight towers
				adc[2] = tmpAdc; 
				*id1 = nextTowerId;
			}
		}
	}

	return kTRUE;
}
//-----------------------------------------------------------------------
void StMcAccpMaker::initEmc() {
  mEmcPosition = new StEmcPosition();
  for(int i=0;i<4;i++) {
    mEmcGeom[i] = StEmcGeom::getEmcGeom(detname[i].Data());
  }
}
//-----------------------------------------------------------------------
void StMcAccpMaker::buildEmcIndex() {
  StEmcDetector *mEmcDet = mEmcCollection->detector(kBarrelEmcTowerId);
  memset(mEmcIndex, 0, sizeof(mEmcIndex));
  if(!mEmcDet) return;
  for (size_t iMod=1; iMod<=mEmcDet->numberOfModules(); ++iMod) {
    StSPtrVecEmcRawHit& modHits = mEmcDet->module(iMod)->hits();
    for (size_t iHit=0; iHit<modHits.size(); ++iHit) {
      StEmcRawHit* rawHit = modHits[iHit];
      if(!rawHit) continue;
      unsigned int softId = rawHit->softId(1); 
      if (mEmcGeom[0]->checkId(softId)==0) { // OK
        mEmcIndex[softId-1] = rawHit;
      }
    }
  }
}
