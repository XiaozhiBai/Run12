//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 15 11:11:24 2012 by ROOT version 5.32/00
// from TTree nt_pair/pairs
// found on file: st_physics_gamma.root
//////////////////////////////////////////////////////////

#ifndef StPheRecoPair_h
#define StPheRecoPair_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class StPheRecoPair {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         p1geantId;
   Float_t         p1p;
   Float_t         p1pt;
   Float_t         p1svx;
   Float_t         p1svy;
   Float_t         p1svz;
   Float_t         p1pgeantId;
   Float_t         p1pp;
   Float_t         p1ppt;
   Float_t         p1eta;
   Float_t         p1peta;
   Float_t         p1label;
   Float_t         p1rp;
   Float_t         p1rpt;
   Float_t         p1reta;
   Float_t         p1rphi;
   Float_t         p1ndedx;
   Float_t         p1nfit;
   Float_t         p1ncom;
   Float_t         p1nmax;
   Float_t         p1dca;
   Float_t         p1chi2;
   Float_t         p2geantId;
   Float_t         p2p;
   Float_t         p2pt;
   Float_t         p2svx;
   Float_t         p2svy;
   Float_t         p2svz;
   Float_t         p2pgeantId;
   Float_t         p2pp;
   Float_t         p2ppt;
   Float_t         p2eta;
   Float_t         p2peta;
   Float_t         p2label;
   Float_t         p2rp;
   Float_t         p2rpt;
   Float_t         p2reta;
   Float_t         p2rphi;
   Float_t         p2nfit;
   Float_t         p2ncom;
   Float_t         p2nmax;
   Float_t         p2dca;
   Float_t         p2chi2;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         massDCA;
   Float_t         mcMassPair;
   Float_t         pairDCA;
   Float_t         pairPT;
   Float_t         theta;
   Float_t         xpair;
   Float_t         ypair;
   Float_t         zpair;
   Float_t         mult;
   Float_t         dist;
   Float_t         p1bemcId;
   Float_t         p1btowAdc0;
   Float_t         p1btowE0;
   Float_t         p1btowE;
   Float_t         p1bemcDistZ;
   Float_t         p1bemcDistPhi;
   Float_t         p1bsmdNEta;
   Float_t         p1bsmdNPhi;
   Float_t         p1btowId;
   Float_t         p2bemcId;
   Float_t         p2btowAdc0;
   Float_t         p2btowE0;
   Float_t         p2btowE;
   Float_t         p2bemcDistZ;
   Float_t         p2bemcDistPhi;
   Float_t         p2bsmdNEta;
   Float_t         p2bsmdNPhi;
   Float_t         p2btowId;
   Float_t         p1pr_rp;
   Float_t         p1pr_rpT;
   Float_t         p2pr_rp;
   Float_t         p2pr_rpT;
   Float_t         p1stpcx;
   Float_t         p1stpcy;
   Float_t         p1stpcz;
   Float_t         p1gpgeantId;
   Float_t         p1gppt;
   Float_t         p2gpgeantId;
   Float_t         p2gppt;

   // List of branches
   TBranch        *b_p1geantId;   //!
   TBranch        *b_p1p;   //!
   TBranch        *b_p1pt;   //!
   TBranch        *b_p1svx;   //!
   TBranch        *b_p1svy;   //!
   TBranch        *b_p1svz;   //!
   TBranch        *b_p1pgeantId;   //!
   TBranch        *b_p1pp;   //!
   TBranch        *b_p1ppt;   //!
   TBranch        *b_p1eta;   //!
   TBranch        *b_p1peta;   //!
   TBranch        *b_p1label;   //!
   TBranch        *b_p1rp;   //!
   TBranch        *b_p1rpt;   //!
   TBranch        *b_p1reta;   //!
   TBranch        *b_p1rphi;   //!
   TBranch        *b_p1ndedx;   //!
   TBranch        *b_p1nfit;   //!
   TBranch        *b_p1ncom;   //!
   TBranch        *b_p1nmax;   //!
   TBranch        *b_p1dca;   //!
   TBranch        *b_p1chi2;   //!
   TBranch        *b_p2geantId;   //!
   TBranch        *b_p2p;   //!
   TBranch        *b_p2pt;   //!
   TBranch        *b_p2svx;   //!
   TBranch        *b_p2svy;   //!
   TBranch        *b_p2svz;   //!
   TBranch        *b_p2pgeantId;   //!
   TBranch        *b_p2pp;   //!
   TBranch        *b_p2ppt;   //!
   TBranch        *b_p2eta;   //!
   TBranch        *b_p2peta;   //!
   TBranch        *b_p2label;   //!
   TBranch        *b_p2rp;   //!
   TBranch        *b_p2rpt;   //!
   TBranch        *b_p2reta;   //!
   TBranch        *b_p2rphi;   //!
   TBranch        *b_p2nfit;   //!
   TBranch        *b_p2ncom;   //!
   TBranch        *b_p2nmax;   //!
   TBranch        *b_p2dca;   //!
   TBranch        *b_p2chi2;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_massDCA;   //!
   TBranch        *b_mcMassPair;
   TBranch        *b_pairDCA;   //!
   TBranch        *b_pairPT;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_xpair;   //!
   TBranch        *b_ypair;   //!
   TBranch        *b_zpair;   //!
   TBranch        *b_mult;   //!
   TBranch        *b_dist;   //!
   TBranch        *b_p1bemcId;   //!
   TBranch        *b_p1btowAdc0;   //!
   TBranch        *b_p1btowE0;   //!
   TBranch        *b_p1btowE;   //!
   TBranch        *b_p1bemcDistZ;   //!
   TBranch        *b_p1bemcDistPhi;   //!
   TBranch        *b_p1bsmdNEta;   //!
   TBranch        *b_p1bsmdNPhi;   //!
   TBranch        *b_p1btowId;   //!
   TBranch        *b_p2bemcId;   //!
   TBranch        *b_p2btowAdc0;   //!
   TBranch        *b_p2btowE0;   //!
   TBranch        *b_p2btowE;   //!
   TBranch        *b_p2bemcDistZ;   //!
   TBranch        *b_p2bemcDistPhi;   //!
   TBranch        *b_p2bsmdNEta;   //!
   TBranch        *b_p2bsmdNPhi;   //!
   TBranch        *b_p2btowId;   //!
   TBranch        *b_p1pr_rp;   //!
   TBranch        *b_p1pr_rpT;   //!
   TBranch        *b_p2pr_rp;   //!
   TBranch        *b_p2pr_rpT;   //!
   TBranch        *b_p1stpcx;   //!
   TBranch        *b_p1stpcy;   //!
   TBranch        *b_p1stpcz;   //!
   TBranch        *b_p1gpgeantId;   //!
   TBranch        *b_p1gppt;   //!
   TBranch        *b_p2gpgeantId;   //!
   TBranch        *b_p2gppt;   //!

   StPheRecoPair();
   virtual ~StPheRecoPair();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual void     Init(TTree *tree);
};

#endif
