/******************************************************************************************
 *Wed Mar 23 14:03:09 PDT 2016                                                            *
 * This class to read embedding analysis tree, extract trigger efficiency                 *
 * By xiaozhi Bai xiaozhi@uic.edu                                                         *
 ******************************************************************************************/
// StTrigEfficiency.cxx

#include <iostream>
#include "StTrigEfficiency.h"
#include "StSingle_Electron_tracks.h"
#include "TFile.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "cmath"

using namespace std;

ClassImp(StTrigEfficiency)
StTrigEfficiency::StTrigEfficiency(const char* outName)
{
  TH1F::SetDefaultSumw2();  
  mOutputFile=new TFile(outName,"RECREATE");
  weightFile=new TFile("NPE_spectra.root","READ");
}
StTrigEfficiency::~StTrigEfficiency()
{
  //..
  //
}
void StTrigEfficiency::bookHistogram()
{ 
  cout<<"gbook histogram"<<endl;
  Float_t Pt_xx[9]={2.5,3.5,4.5,5.5,6.5,7.5,8.5,10.5,16.5};
  
  //track
  m_DSMadc11_cut=new TH1F("m_DSMadc11_cut","m_DSMadc11_cut",200,0,20);
  m_DSMadc18_cut=new TH1F("m_DSMadc18_cut","m_DSMadc18_cut",200,0,20); 
  m_DSMadc=new TH1F("m_DSMadc","m_DSMadc",200,0,20);
}
void StTrigEfficiency::read(TString fileName)
{ 
  TFile * infile=new TFile(fileName.Data());
  TNtuple* ntTracks=(TNtuple* ) infile->Get("tracksMC");
  StSingle_Electron_tracks * tracksMC=new StSingle_Electron_tracks; 
  tracksMC->Init(ntTracks);
  weightFile->cd();
  
  TH1F *N_npe=(TH1F *) weightFile->Get("N_npe");
  TF1* f1 = new TF1("f1","[0]/(pow(x,[1])+[2])",0,20);
  f1->SetParameter(0,10);
  f1->SetParameter(1,8);
  f1->SetParameter(2,0.7);
  N_npe->Fit("f1","R0Q","" ,2.5,14);
  
  for(Int_t iTrk=0;iTrk<tracksMC->GetEntries();iTrk++)
    {
      tracksMC->GetEntry(iTrk);
      if(tracksMC->dca<1.5 && tracksMC->nfit>20 && tracksMC->nDedxPts>15 && tracksMC->nfit/(Float_t)tracksMC->nmax > 0.52 && sqrt(tracksMC->stpcx*tracksMC->stpcx+tracksMC->stpcy*tracksMC->stpcy)<73 && abs(tracksMC->eta)<0.7 && tracksMC->btowDsmAdc0 < tracksMC->btowAdc0*0.1 && 0<tracksMC->rpt  && tracksMC->bsmdNEta>1 && tracksMC->bsmdNPhi>1 && 0.3<tracksMC->p/tracksMC->btowE0 && tracksMC->p/tracksMC->btowE0<1.5 && abs(tracksMC->bemcDistZ)<3 && abs(tracksMC->bemcDistPhi)<0.015)
	{
	  Float_t weight =f1->Eval(tracksMC->pt);
	  cout<<"*********"<<weight<<endl;
	  if(11<tracksMC->btowDsmAdc0 )
	    m_DSMadc11_cut->Fill(tracksMC->pt,weight);
	  if(18<tracksMC->btowDsmAdc0 )
	    m_DSMadc18_cut->Fill(tracksMC->pt,weight);	
	  m_DSMadc->Fill(tracksMC->pt,weight);
	}
    }
}
void StTrigEfficiency::WriteHistogram()
{ 
  cout<<"  write"<<endl;
  mOutputFile->cd();
  m_DSMadc11_cut->Write();
  m_DSMadc18_cut->Write();
  m_DSMadc->Write();
}
/*
Double_t StTrigEfficiency::Weight(Float_t pt)
{
    Double_t weight; 
    //    TFile * weightFile=new TFile("pi0_gamma_fit.root","READ");
    weightFile->cd();

    TH1F *N_npe=(TH1F *) weightFile->Get("N_npe");

    TF1* f1 = new TF1("f1","[0]/(pow(x,[1])+[2])",0,20);
    f1->SetParameter(0,10);
    f1->SetParameter(1,8);
    f1->SetParameter(2,0.7); 
    N_npe->Fit("f1","R0Q","" ,2.5,14);

    weight=f1->Eval(pt);

    cout<<weight<<"!!!!!!!!!!!!!!!!!"<<" pt  "<<pt<<endl;

    return weight;


   
}
*/
