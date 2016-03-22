#define StPheRecoTrack_cxx
#include "StPheRecoTrack.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

StPheRecoTrack::StPheRecoTrack() 
{
}

StPheRecoTrack::~StPheRecoTrack()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t StPheRecoTrack::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t StPheRecoTrack::GetEntries()
{
    if(fChain) return fChain->GetEntriesFast();
    return 0;
}

void StPheRecoTrack::Init(TTree *tree)
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

   fChain->SetBranchAddress("geantId", &geantId, &b_geantId);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("svx", &svx, &b_svx);
   fChain->SetBranchAddress("svy", &svy, &b_svy);
   fChain->SetBranchAddress("svz", &svz, &b_svz);
   fChain->SetBranchAddress("pgeantId", &pgeantId, &b_pgeantId);
   fChain->SetBranchAddress("pp", &pp, &b_pp);
   fChain->SetBranchAddress("ppt", &ppt, &b_ppt);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("pphi", &pphi, &b_pphi);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("peta", &peta, &b_peta);
   fChain->SetBranchAddress("label", &label, &b_label);
   fChain->SetBranchAddress("plabel", &plabel, &b_plabel);
   fChain->SetBranchAddress("gplabel", &gplabel, &b_gplabel);
   fChain->SetBranchAddress("rp", &rp, &b_rp);
   fChain->SetBranchAddress("rpt", &rpt, &b_rpt);
   fChain->SetBranchAddress("reta", &reta, &b_reta);
   fChain->SetBranchAddress("rphi", &rphi, &b_rphi);
   fChain->SetBranchAddress("nfit", &nfit, &b_nfit);
   fChain->SetBranchAddress("ncom", &ncom, &b_ncom);
   fChain->SetBranchAddress("nmax", &nmax, &b_nmax);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("dca", &dca, &b_dca);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("nDedxPts", &nDedxPts, &b_nDedxPts);
   fChain->SetBranchAddress("dedx", &dedx, &b_dedx);
   fChain->SetBranchAddress("dedx_2", &dedx_2, &b_dedx_2);
   fChain->SetBranchAddress("nSigPi", &nSigPi, &b_nSigPi);
   fChain->SetBranchAddress("nSigK", &nSigK, &b_nSigK);
   fChain->SetBranchAddress("nSigP", &nSigP, &b_nSigP);
   fChain->SetBranchAddress("nSigE", &nSigE, &b_nSigE);
   fChain->SetBranchAddress("bemcId", &bemcId, &b_bemcId);
   fChain->SetBranchAddress("btowAdc0", &btowAdc0, &b_btowAdc0);
   fChain->SetBranchAddress("btowE0", &btowE0, &b_btowE0);
   fChain->SetBranchAddress("btowE", &btowE, &b_btowE);
   fChain->SetBranchAddress("bemcDistZ", &bemcDistZ, &b_bemcDistZ);
   fChain->SetBranchAddress("bemcDistPhi", &bemcDistPhi, &b_bemcDistPhi);
   fChain->SetBranchAddress("bsmdNEta", &bsmdNEta, &b_bsmdNEta);
   fChain->SetBranchAddress("bsmdNPhi", &bsmdNPhi, &b_bsmdNPhi);
   fChain->SetBranchAddress("btowId", &btowId, &b_btowId);
   fChain->SetBranchAddress("pr_rp", &pr_rp, &b_pr_rp);
   fChain->SetBranchAddress("pr_rpT", &pr_rpT, &b_pr_rpT);
   fChain->SetBranchAddress("stpcx", &stpcx, &b_stpcx);
   fChain->SetBranchAddress("stpcy", &stpcy, &b_stpcy);
   fChain->SetBranchAddress("stpcz", &stpcz, &b_stpcz);
   fChain->SetBranchAddress("gpgeantId", &gpgeantId, &b_gpgeantId);
   fChain->SetBranchAddress("gppt", &gppt, &b_gppt);
}
