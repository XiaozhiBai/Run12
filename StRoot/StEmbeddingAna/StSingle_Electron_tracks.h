//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 22 17:13:16 2012 by ROOT version 5.32/00
// from TTree tracksMC/tracks
// found on file: all_2.root
//////////////////////////////////////////////////////////

#ifndef StSingle_Electron_tracks_h
#define StSingle_Electron_tracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class StSingle_Electron_tracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         geantId;
   Float_t         p;
   Float_t         pt;
   Float_t         svx;
   Float_t         svy;
   Float_t         svz;
   Float_t         pgeantId;
   Float_t         pp;
   Float_t         ppt;
   Float_t         phi;
   Float_t         y;
   Float_t         eta;
   Float_t         peta;
   Float_t         label;
   Float_t         rp;
   Float_t         rpt;
   Float_t         reta;
   Float_t         rphi;
   Float_t         nfit;
   Float_t         ncom;
   Float_t         nmax;
   Float_t         nDedxPts;
   Float_t         dedx;
   Float_t         dedx_2;
   Float_t         nSigPi;
   Float_t         nSigK;
   Float_t         nSigP;
   Float_t         nSigE;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         bemcId;
   Float_t         btowId0;
   Float_t         btowAdc0;
   Float_t         btowE0;
   Float_t         btowE;
   Float_t         bemcDistZ;
   Float_t         bemcDistPhi;
   Float_t         bsmdNEta;
   Float_t         bsmdNPhi;
   Float_t         btowId;
   Float_t         btowDsmAdcProj;
   Float_t         btowDsmAdcNext;
   Float_t         btowIdNext;
   Float_t         btowDsmAdc0;
   Float_t         dca;
   Float_t         pr_rp;
   Float_t         pr_rpT;
   Float_t         stpcx;
   Float_t         stpcy;
   Float_t         stpcz;
   /*   
//
   Float_t         geantId;
   Float_t         p;  //momentum
   Float_t         pt; // pt
   Float_t         svx;//?
   Float_t         svy;
   Float_t         svz;
   Float_t         pgeantId; //primary track
   Float_t         pp; 
   Float_t         ppt;
   Float_t         phi;  //phi
   Float_t         y; //rapidity
   Float_t         eta;   // pseudoRapidity
   Float_t         peta;  //primary pseudora
   Float_t         label;  //?
   Float_t         rp;   //    " r" means reconstruction track 
   Float_t         rpt;   //    
   Float_t         reta;  //
   Float_t         rphi;  //
   Float_t         nfit;  // number of TPC hit for Fit track 
   Float_t         ncom;  //  the number of common point of Reconstruction track and MCtrack     
   Float_t         nmax;  //maxmum hit  
   Float_t         nDedxPts;  //?
   Float_t         dedx;  // energy lose
   Float_t         nSigPi; // nsigma Pi or nsigma Pi of electron 
   Float_t         nSigK;
   Float_t         nSigP;
   Float_t         nSigE;
   Float_t         vx;   // primary vertex { startVertex} 
   Float_t         vy;
   Float_t         vz;
   Float_t         bemcId;   //    // difference between tower id
   Float_t         btowAdc0;  //  offline adc
   Float_t         btowE0;  //the highest energy of the cluster 
   Float_t         btowE;  //energy of the cluter
   Float_t         bemcDistZ;  //tower Dz
   Float_t         bemcDistPhi;  //tower Dphi
   Float_t         bsmdNEta;  // the strip number of etaPlane 
   Float_t         bsmdNPhi;  //phi
   Float_t         btowId; //tower id
   Float_t         dca; //global dca  
   Float_t         pr_rp;  // ? primary momentum  
   Float_t         pr_rpT;//?   primari pt
   Float_t         stpcx; //  first  tpc point x
   Float_t         stpcy;//               y
   Float_t         stpcz;//               z
   */
   // List of branches
   TBranch        *b_geantId;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_svx;   //!
   TBranch        *b_svy;   //!
   TBranch        *b_svz;   //!
   TBranch        *b_pgeantId;   //!
   TBranch        *b_pp;   //!
   TBranch        *b_ppt;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_y;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_peta;   //!
   TBranch        *b_label;   //!
   TBranch        *b_rp;   //!
   TBranch        *b_rpt;   //!
   TBranch        *b_reta;   //!
   TBranch        *b_rphi;   //!
   TBranch        *b_nfit;   //!
   TBranch        *b_ncom;   //!
   TBranch        *b_nmax;   //!
   TBranch        *b_nDedxPts;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_dedx_2;   //!
   TBranch        *b_nSigPi;   //!
   TBranch        *b_nSigK;   //!
   TBranch        *b_nSigP;   //!
   TBranch        *b_nSigE;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_bemcId;   //!
   TBranch        *b_btowId0;   //!
   TBranch        *b_btowAdc0;   //!
   TBranch        *b_btowE0;   //!
   TBranch        *b_btowE;   //!
   TBranch        *b_bemcDistZ;   //!
   TBranch        *b_bemcDistPhi;   //!
   TBranch        *b_bsmdNEta;   //!
   TBranch        *b_bsmdNPhi;   //!
   TBranch        *b_btowId;   //!
   TBranch        *b_btowDsmAdcProj;   //!
   TBranch        *b_btowDsmAdcNext;   //!
   TBranch        *b_btowIdNext;   //!
   TBranch        *b_btowDsmAdc0;   //!
   TBranch        *b_dca;   //!
   TBranch        *b_pr_rp;   //!
   TBranch        *b_pr_rpT;   //!
   TBranch        *b_stpcx;   //!
   TBranch        *b_stpcy;   //!
   TBranch        *b_stpcz;   //!
   /*
   //
   TBranch        *b_geantId;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_svx;   //!
   TBranch        *b_svy;   //!
   TBranch        *b_svz;   //!
   TBranch        *b_pgeantId;   //!
   TBranch        *b_pp;   //!
   TBranch        *b_ppt;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_y;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_peta;   //!
   TBranch        *b_label;   //!
   TBranch        *b_rp;   //!
   TBranch        *b_rpt;   //!
   TBranch        *b_reta;   //!
   TBranch        *b_rphi;   //!
   TBranch        *b_nfit;   //!
   TBranch        *b_ncom;   //!
   TBranch        *b_nmax;   //!
   TBranch        *b_nDedxPts;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_nSigPi;   //!
   TBranch        *b_nSigK;   //!
   TBranch        *b_nSigP;   //!
   TBranch        *b_nSigE;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_bemcId;   //!
   TBranch        *b_btowAdc0;   //!
   TBranch        *b_btowE0;   //!
   TBranch        *b_btowE;   //!
   TBranch        *b_bemcDistZ;   //!
   TBranch        *b_bemcDistPhi;   //!
   TBranch        *b_bsmdNEta;   //!
   TBranch        *b_bsmdNPhi;   //!
   TBranch        *b_btowId;   //!
   TBranch        *b_dca;   //!
   TBranch        *b_pr_rp;   //!
   TBranch        *b_pr_rpT;   //!
   TBranch        *b_stpcx;   //!
   TBranch        *b_stpcy;   //!
   TBranch        *b_stpcz;   //!
   */
   StSingle_Electron_tracks();
   virtual ~StSingle_Electron_tracks();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual void     Init(TTree *tree);
};

#endif
