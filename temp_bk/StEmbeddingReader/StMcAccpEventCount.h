//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 22 17:13:19 2012 by ROOT version 5.32/00
// from TTree eventCount/eventCount
// found on file: all_2.root
//////////////////////////////////////////////////////////

#ifndef StMcAccpEventCount_h
#define StMcAccpEventCount_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class StMcAccpEventCount {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         runId;
   Float_t         eventId;
   Float_t         mcVx;
   Float_t         mcVy;
   Float_t         mcVz;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         vzVpd;
   Float_t         posRefMult;
   Float_t         negRefMult;
   Float_t         zdc;
   Float_t         bbc;
   Float_t         nMcTracks;
   Float_t         nRTracks;
   Float_t         magField;
   Float_t         t0;
   Float_t         t1;
   Float_t         t2;
   Float_t         t3;
   Float_t         t4;
   Float_t         t5;

   // List of branches
   TBranch        *b_runId;   //!
   TBranch        *b_eventId;   //!
   TBranch        *b_mcVx;   //!
   TBranch        *b_mcVy;   //!
   TBranch        *b_mcVz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vzVpd;   //!
   TBranch        *b_posRefMult;   //!
   TBranch        *b_negRefMult;   //!
   TBranch        *b_zdc;   //!
   TBranch        *b_bbc;   //!
   TBranch        *b_nMcTracks;   //!
   TBranch        *b_nRTracks;   //!
   TBranch        *b_magField;   //!
   TBranch        *b_t0;   //!
   TBranch        *b_t1;   //!
   TBranch        *b_t2;   //!
   TBranch        *b_t3;   //!
   TBranch        *b_t4;   //!
   TBranch        *b_t5;   //!

   StMcAccpEventCount();
   virtual ~StMcAccpEventCount();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Int_t GetEntries();
};

#endif
