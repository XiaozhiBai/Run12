/******************************************************************************************
 *Mon Mar  7 19:33:04 EST 2016                                                            *
 * This class for read  DemesonTree Fill all the No-phoronic electron realted histogram   *
 * By xiaozhi Bai xiaozhi@uic.edu                                                         *
 ******************************************************************************************/
// StNpemaker.cxx

#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include "StNpeMaker.h"
//#include "StCuts.h"

//#include "StTRIGGERS.h"

//#include "prescales.h"
#include "StDmesonMaker/StDmesonEvent.h"
#include "StDmesonMaker/StDmesonTrack.h"
#include "StDmesonMaker/StElectronPair.h"

#include "StLorentzVectorF.hh"
#include "StMessMgr.h"

#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TClonesArray.h"

ClassImp(StNpeMaker)
using namespace std;
StNpeMaker::StNpeMaker(const char* outName){
  mOutputFile = new TFile(outName, "RECREATE");
  TH1F::SetDefaultSumw2();
  TH2F::SetDefaultSumw2();
  
}
StNpeMaker::~StNpeMaker(){
  /*  */
}
void StNpeMaker::Clear(Option_t *opt) {
  /**/
}
void StNpeMaker::read(TString fileName){
  TFile* inFile = new TFile(fileName.Data(),"READ");
  TTree* tree = (TTree*)inFile->Get("T");
  tree->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  tree->SetBranchAddress("dEvent", &mEvent);

  if(!mEvent){
    LOG_WARN << " No DmesonEvent! Skip! " << endm;
    return kStOK;  
  } 
    
  TClonesArray*   aTracks = 0;
  TClonesArray* aPairs=0;
  LOG_INFO <<"Read File "<<fileName<< " Total number of events"<< tree->GetEntriesFast()<<endm;
  for (UInt_t i = 0; i < tree->GetEntriesFast(); i++)
    {
      tree->GetEntry(i);
      Run_QA(mEvent); 
      
      //     LOG_INFO <<mEvent->primaryVertex().z()<<endm;
      // aPairs=mEvent->electronPair();
      // aTracks=mEvent->tracks();    
    }
}

void StNpeMaker::bookObjects(){
  mEvent = new StDmesonEvent;
  HT0_HT2=new TH2F("HT0_HT2","",2,0,2,2,0,2);
  cout<< " Book"<<endl;
}
void StNpeMaker::writeObjects(){
  HT0_HT2->Write();
}
