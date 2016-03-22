#include "StMcAccpEventCount.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

StMcAccpEventCount::StMcAccpEventCount()
{
}

StMcAccpEventCount::~StMcAccpEventCount()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t StMcAccpEventCount::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t StMcAccpEventCount::GetEntries()
{
    if(fChain)return fChain->GetEntriesFast();
    return 0;
}

void StMcAccpEventCount::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runId", &runId, &b_runId);
   fChain->SetBranchAddress("eventId", &eventId, &b_eventId);
   fChain->SetBranchAddress("mcVx", &mcVx, &b_mcVx);
   fChain->SetBranchAddress("mcVy", &mcVy, &b_mcVy);
   fChain->SetBranchAddress("mcVz", &mcVz, &b_mcVz);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("vzVpd", &vzVpd, &b_vzVpd);
   fChain->SetBranchAddress("posRefMult", &posRefMult, &b_posRefMult);
   fChain->SetBranchAddress("negRefMult", &negRefMult, &b_negRefMult);
   fChain->SetBranchAddress("zdc", &zdc, &b_zdc);
   fChain->SetBranchAddress("bbc", &bbc, &b_bbc);
   fChain->SetBranchAddress("nMcTracks", &nMcTracks, &b_nMcTracks);
   fChain->SetBranchAddress("nRTracks", &nRTracks, &b_nRTracks);
   fChain->SetBranchAddress("magField", &magField, &b_magField);
   fChain->SetBranchAddress("t0", &t0, &b_t0);
   fChain->SetBranchAddress("t1", &t1, &b_t1);
   fChain->SetBranchAddress("t2", &t2, &b_t2);
   fChain->SetBranchAddress("t3", &t3, &b_t3);
   fChain->SetBranchAddress("t4", &t4, &b_t4);
   fChain->SetBranchAddress("t5", &t5, &b_t5);
}
