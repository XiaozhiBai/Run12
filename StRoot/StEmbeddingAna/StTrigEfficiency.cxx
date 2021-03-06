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
  weightFile=new TFile("/global/homes/x/xiao00/work/run12/Embedding/fonll200gev.root","READ");
  gFONLLc=(TGraphErrors *) weightFile->Get("gFONLLc");
}
StTrigEfficiency::~StTrigEfficiency()
{
  //..
  //
}
void StTrigEfficiency::bookHistogram()
{ 
  
  //track
  m_DSMadc11_cut=new TH1F("m_DSMadc11_cut","m_DSMadc11_cut",2000,0,20);
  m_DSMadc18_cut=new TH1F("m_DSMadc18_cut","m_DSMadc18_cut",2000,0,20); 
  m_DSMadc=new TH1F("m_DSMadc","m_DSMadc",2000,0,20);

  m_DSMadc11_cut->Sumw2();
  m_DSMadc18_cut->Sumw2();
  m_DSMadc->Sumw2();

}
void StTrigEfficiency::read(TString fileName)
{ 
  TFile * infile=new TFile(fileName.Data());
  TNtuple* ntTracks=(TNtuple* ) infile->Get("tracksMC");
  StSingle_Electron_tracks * tracksMC=new StSingle_Electron_tracks; 
  tracksMC->Init(ntTracks);
  weightFile->cd();
  

  cout<< "Working on Trigger efficiency!!!!"<<endl; 

  for(Int_t iTrk=0;iTrk<tracksMC->GetEntries();iTrk++)
    {
      tracksMC->GetEntry(iTrk);
     
      if(tracksMC->dca<1.5 && tracksMC->nfit>20 && tracksMC->nDedxPts>15 && tracksMC->nfit/(Float_t)tracksMC->nmax > 0.52 && sqrt(tracksMC->stpcx*tracksMC->stpcx+tracksMC->stpcy*tracksMC->stpcy)<73 && abs(tracksMC->eta)<0.7 && tracksMC->btowDsmAdc0 < tracksMC->btowAdc0*0.1 && 0<tracksMC->rpt  && tracksMC->bsmdNEta>1 && tracksMC->bsmdNPhi>1 && 0.3<tracksMC->p/tracksMC->btowE0 && tracksMC->p/tracksMC->btowE0<1.5 && abs(tracksMC->bemcDistZ)<3 && abs(tracksMC->bemcDistPhi)<0.015&&10<tracksMC->ncom)
	{
	  Float_t weight =Weight(tracksMC->pt,tracksMC->geantId);

	  m_DSMadc->Fill(tracksMC->pt,weight);
	  if(11<tracksMC->btowDsmAdc0&&tracksMC->btowAdc0>180)
	    m_DSMadc11_cut->Fill(tracksMC->pt,weight);
	  if(18<tracksMC->btowDsmAdc0 &&tracksMC->btowAdc0>300)
	    m_DSMadc18_cut->Fill(tracksMC->pt,weight);	
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
Double_t StTrigEfficiency::Weight(Float_t pt,Float_t geantID)
{
  Double_t weight=0;

  if(abs(geantID)==2||abs(geantID)==3)
    weight=gFONLLc->Eval(pt);
  return weight;
}

