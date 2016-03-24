#define StPhePair_cxx
#include "StPhePair.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


StPhePair::StPhePair()
{

}

StPhePair::~StPhePair()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t StPhePair::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t StPhePair::GetEntries()
{
  if(fChain)return fChain->GetEntriesFast();
  return 0;
}

void StPhePair::Init(TTree *tree)
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

   fChain->SetBranchAddress("p1geantId", &p1geantId, &b_p1geantId);
   fChain->SetBranchAddress("p1p", &p1p, &b_p1p);
   fChain->SetBranchAddress("p1pt", &p1pt, &b_p1pt);
   fChain->SetBranchAddress("p1svx", &p1svx, &b_p1svx);
   fChain->SetBranchAddress("p1svy", &p1svy, &b_p1svy);
   fChain->SetBranchAddress("p1svz", &p1svz, &b_p1svz);
   fChain->SetBranchAddress("p1pgeantId", &p1pgeantId, &b_p1pgeantId);
   fChain->SetBranchAddress("p1pp", &p1pp, &b_p1pp);
   fChain->SetBranchAddress("p1ppt", &p1ppt, &b_p1ppt);
   fChain->SetBranchAddress("p1eta", &p1eta, &b_p1eta);
   fChain->SetBranchAddress("p1peta", &p1peta, &b_p1peta);
   fChain->SetBranchAddress("p1label", &p1label, &b_p1label);
   fChain->SetBranchAddress("p1rp", &p1rp, &b_p1rp);
   fChain->SetBranchAddress("p1rpt", &p1rpt, &b_p1rpt);
   fChain->SetBranchAddress("p1reta", &p1reta, &b_p1reta);
   fChain->SetBranchAddress("p1rphi", &p1rphi, &b_p1rphi);
   fChain->SetBranchAddress("p1ndedx", &p1ndedx, &b_p1ndedx);
   fChain->SetBranchAddress("p1nfit", &p1nfit, &b_p1nfit);
   fChain->SetBranchAddress("p1ncom", &p1ncom, &b_p1ncom);
   fChain->SetBranchAddress("p1nmax", &p1nmax, &b_p1nmax);
   fChain->SetBranchAddress("p1dca", &p1dca, &b_p1dca);
   fChain->SetBranchAddress("p1chi2", &p1chi2, &b_p1chi2);
   fChain->SetBranchAddress("p2geantId", &p2geantId, &b_p2geantId);
   fChain->SetBranchAddress("p2p", &p2p, &b_p2p);
   fChain->SetBranchAddress("p2pt", &p2pt, &b_p2pt);
   fChain->SetBranchAddress("p2svx", &p2svx, &b_p2svx);
   fChain->SetBranchAddress("p2svy", &p2svy, &b_p2svy);
   fChain->SetBranchAddress("p2svz", &p2svz, &b_p2svz);
   fChain->SetBranchAddress("p2pgeantId", &p2pgeantId, &b_p2pgeantId);
   fChain->SetBranchAddress("p2pp", &p2pp, &b_p2pp);
   fChain->SetBranchAddress("p2ppt", &p2ppt, &b_p2ppt);
   fChain->SetBranchAddress("p2eta", &p2eta, &b_p2eta);
   fChain->SetBranchAddress("p2peta", &p2peta, &b_p2peta);
   fChain->SetBranchAddress("p2label", &p2label, &b_p2label);
   fChain->SetBranchAddress("p2rp", &p2rp, &b_p2rp);
   fChain->SetBranchAddress("p2rpt", &p2rpt, &b_p2rpt);
   fChain->SetBranchAddress("p2reta", &p2reta, &b_p2reta);
   fChain->SetBranchAddress("p2rphi", &p2rphi, &b_p2rphi);
   fChain->SetBranchAddress("p2nfit", &p2nfit, &b_p2nfit);
   fChain->SetBranchAddress("p2ncom", &p2ncom, &b_p2ncom);
   fChain->SetBranchAddress("p2nmax", &p2nmax, &b_p2nmax);
   fChain->SetBranchAddress("p2dca", &p2dca, &b_p2dca);
   fChain->SetBranchAddress("p2chi2", &p2chi2, &b_p2chi2);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("massDCA", &massDCA, &b_massDCA);
   fChain->SetBranchAddress("mcMassPair", &mcMassPair, &b_mcMassPair);
   fChain->SetBranchAddress("pairDCA", &pairDCA, &b_pairDCA);
   fChain->SetBranchAddress("pairPT", &pairPT, &b_pairPT);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("xpair", &xpair, &b_xpair);
   fChain->SetBranchAddress("ypair", &ypair, &b_ypair);
   fChain->SetBranchAddress("zpair", &zpair, &b_zpair);
   fChain->SetBranchAddress("mult", &mult, &b_mult);
   fChain->SetBranchAddress("dist", &dist, &b_dist);
   fChain->SetBranchAddress("p1bemcId", &p1bemcId, &b_p1bemcId);
   fChain->SetBranchAddress("p1btowAdc0", &p1btowAdc0, &b_p1btowAdc0);
   fChain->SetBranchAddress("p1btowE0", &p1btowE0, &b_p1btowE0);
   fChain->SetBranchAddress("p1btowE", &p1btowE, &b_p1btowE);
   fChain->SetBranchAddress("p1bemcDistZ", &p1bemcDistZ, &b_p1bemcDistZ);
   fChain->SetBranchAddress("p1bemcDistPhi", &p1bemcDistPhi, &b_p1bemcDistPhi);
   fChain->SetBranchAddress("p1bsmdNEta", &p1bsmdNEta, &b_p1bsmdNEta);
   fChain->SetBranchAddress("p1bsmdNPhi", &p1bsmdNPhi, &b_p1bsmdNPhi);
   fChain->SetBranchAddress("p1btowId", &p1btowId, &b_p1btowId);
   fChain->SetBranchAddress("p2bemcId", &p2bemcId, &b_p2bemcId);
   fChain->SetBranchAddress("p2btowAdc0", &p2btowAdc0, &b_p2btowAdc0);
   fChain->SetBranchAddress("p2btowE0", &p2btowE0, &b_p2btowE0);
   fChain->SetBranchAddress("p2btowE", &p2btowE, &b_p2btowE);
   fChain->SetBranchAddress("p2bemcDistZ", &p2bemcDistZ, &b_p2bemcDistZ);
   fChain->SetBranchAddress("p2bemcDistPhi", &p2bemcDistPhi, &b_p2bemcDistPhi);
   fChain->SetBranchAddress("p2bsmdNEta", &p2bsmdNEta, &b_p2bsmdNEta);
   fChain->SetBranchAddress("p2bsmdNPhi", &p2bsmdNPhi, &b_p2bsmdNPhi);
   fChain->SetBranchAddress("p2btowId", &p2btowId, &b_p2btowId);
   fChain->SetBranchAddress("p1pr_rp", &p1pr_rp, &b_p1pr_rp);
   fChain->SetBranchAddress("p1pr_rpT", &p1pr_rpT, &b_p1pr_rpT);
   fChain->SetBranchAddress("p2pr_rp", &p2pr_rp, &b_p2pr_rp);
   fChain->SetBranchAddress("p2pr_rpT", &p2pr_rpT, &b_p2pr_rpT);
   fChain->SetBranchAddress("p1stpcx", &p1stpcx, &b_p1stpcx);
   fChain->SetBranchAddress("p1stpcy", &p1stpcy, &b_p1stpcy);
   fChain->SetBranchAddress("p1stpcz", &p1stpcz, &b_p1stpcz);
   fChain->SetBranchAddress("p1gpgeantId", &p1gpgeantId, &b_p1gpgeantId);
   fChain->SetBranchAddress("p1gppt", &p1gppt, &b_p1gppt);
   fChain->SetBranchAddress("p2gpgeantId", &p2gpgeantId, &b_p2gpgeantId);
   fChain->SetBranchAddress("p2gppt", &p2gppt, &b_p2gppt);
}
