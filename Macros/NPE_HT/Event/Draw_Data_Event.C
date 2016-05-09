/*--------------------Xiaozhi------------------------------------ 


  by xiaozhi
  Wed Dec 16 14:56:01 CST 2015
  
  ---------------------------------------------------------------
*/
#include<fstream>
#include <iostream>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TNtuple.h"

using namespace std;
const int nTrg=3;


TH1F *mh1Vz[nTrg];
TH1F *mh1Vz_ps[nTrg];

TH1F *mh1Nevents[nTrg];
TH1F *mh1Nevents_ps[nTrg];

TH1F *mh1Nevents_MB[nTrg];
TH1F *mh1Nevents_MBps[nTrg];

TString mh1Name_Vz[nTrg]={"mh1hVzTrg0","mh1hVzTrg1","mh1hVzTrg2"};
TString mh1Name_Vzps[nTrg]={"mh1hVz_psTrg0","mh1hVz_psTrg1","mh1hVz_psTrg2"};

TString mh1Name_Nevents[nTrg]={"mh1HT_NeventsTrg0","mh1HT_NeventsTrg1","mh1HT_NeventsTrg2"};
TString mh1Name_Nevents_ps[nTrg]={"mh1HT_Nevents_psTrg0","mh1HT_Nevents_psTrg1","mh1HT_Nevents_psTrg2"};

TString mh1Name_Nevents_MB[nTrg]={"mh1MB_NeventsTrg0","mh1MB_NeventsTrg1","mh1MB_NeventsTrg2"};
TString mh1Name_Nevents_MBps[nTrg]={"mh1MB_Nevents_psTrg0","mh1MB_Nevents_psTrg1","mh1MB_Nevents_psTrg2"};

TH1F *nEvents=new TH1F("nEvents","",8,1,9);
TH1F *nEvents_ps=new TH1F("nEvents_ps","",8,1,9);

TH1F *nEvents_MB=new TH1F("nEventsMB","",8,1,9);
TH1F *nEvents_MBps=new TH1F("nEvents_MBps","",8,1,9);

TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
void Draw_NEvents(TH1F *, TH1F *,TH1F *, TH1F *);
char buf[1024];

int Draw_Data_Event(){
  TH1F:: SetDefaultSumw2();
  TFile * inFile= new TFile("../RootFile/Root_File_5_2/hist_5_2.root","READ");
  for(int iTrg=0;iTrg<nTrg;iTrg++)
    {
      mh1Vz[iTrg]=(TH1F *) inFile->Get(mh1Name_Vz[iTrg]);	 
      mh1Vz_ps[iTrg]=(TH1F *) inFile->Get(mh1Name_Vzps[iTrg]);	 
      
      mh1Vz[iTrg]->GetXaxis()->SetTitle("TPc Vz (cm)");
      mh1Vz[iTrg]->GetYaxis()->SetTitle("Counts");
      mh1Vz_ps[iTrg]->GetXaxis()->SetTitle("TPc Vz (cm)");
      mh1Vz_ps[iTrg]->GetYaxis()->SetTitle("Counts");

      mh1Vz[iTrg]->Rebin(10);
      mh1Vz_ps[iTrg]->Rebin(10);

      mh1Nevents[iTrg]=(TH1F *) inFile->Get(mh1Name_Nevents[iTrg]);	 
      mh1Nevents_ps[iTrg]=(TH1F *) inFile->Get(mh1Name_Nevents_ps[iTrg]);	 
      
      mh1Nevents_MB[iTrg]=(TH1F *) inFile->Get(mh1Name_Nevents_MB[iTrg]);	 
      mh1Nevents_MBps[iTrg]=(TH1F *) inFile->Get(mh1Name_Nevents_MBps[iTrg]);	 
    
    }
  nEvents->SetBinContent(2,mh1Nevents[0]->GetBinContent(3));
  nEvents->SetBinContent(4,mh1Nevents[1]->GetBinContent(3));
  nEvents->SetBinContent(6,mh1Nevents[2]->GetBinContent(3));

  nEvents_ps->SetBinContent(2,mh1Nevents_ps[0]->GetBinContent(3));
  nEvents_ps->SetBinContent(4,mh1Nevents_ps[1]->GetBinContent(3));
  nEvents_ps->SetBinContent(6,mh1Nevents_ps[2]->GetBinContent(3));

  
  nEvents_MB->SetBinContent(2,mh1Nevents_MB[0]->GetBinContent(2));
  nEvents_MB->SetBinContent(4,mh1Nevents_MB[1]->GetBinContent(2));
  nEvents_MB->SetBinContent(6,mh1Nevents_MB[2]->GetBinContent(2));
  
  nEvents_MBps->SetBinContent(2,mh1Nevents_MBps[0]->GetBinContent(2));
  nEvents_MBps->SetBinContent(4,mh1Nevents_MBps[1]->GetBinContent(2));
  nEvents_MBps->SetBinContent(6,mh1Nevents_MBps[2]->GetBinContent(2));

  for(int iTrg=0;iTrg<nTrg;iTrg++)
    {
      TCanvas *c2=new TCanvas("c2","",1200,1000);
      c2->Divide(2,1);

      c2->cd(1);

      mh1Vz[iTrg]->Draw();
      c2->cd(2);
      mh1Vz_ps[iTrg]->Draw();
      sprintf(buf,"mh1Vz%i.pdf",iTrg);
      c2->SaveAs(buf);
    }
  Draw_NEvents(nEvents,nEvents_ps,nEvents_MB,nEvents_MBps);
  return 0;
}
void Draw_NEvents(TH1F *nEvent,TH1F *nEvent_ps,TH1F *nEvents_MB,TH1F *nEvents_MBps)
{
  gStyle->SetOptStat(00000);

  TCanvas *c2=new TCanvas("c2","",800,600);
  gPad->cd();
  gPad->SetLogy(1);
  gPad->SetBottomMargin(0.13);
  
  nEvent->SetFillColor(4);
  nEvent->SetFillStyle(3001);

  nEvent_ps->SetFillColor(4);
  nEvent_ps->SetFillStyle(3001);

  nEvents_MB->SetFillColor(4);
  nEvents_MB->SetFillStyle(3001);

  nEvents_MBps->SetFillColor(4);
  nEvents_MBps->SetFillStyle(3001);
  
  TH2F *h1=new TH2F("h1","",8,1,9,100,1,2e9);
  TH2F *h2=new TH2F("h2","",8,1,9,100,1,1e13);
  TH2F *h3=new TH2F("h3","",8,1,9,100,1,2e9);
  TH2F *h4=new TH2F("h4","",8,1,9,100,1,2e15);
  
  h1->GetXaxis()->SetBinLabel(2,"HT0BBCMBTOF0");
  h2->GetXaxis()->SetBinLabel(2,"HT0BBCMBTOF0");
  h3->GetXaxis()->SetBinLabel(2,"BBCMBTOF0");
  h4->GetXaxis()->SetBinLabel(2,"BBCMBTOF0");
  
  h1->GetXaxis()->SetBinLabel(4,"HT2BBCMB");
  h2->GetXaxis()->SetBinLabel(4,"HT2BBCMB");
  h3->GetXaxis()->SetBinLabel(4,"BBCMB");
  h4->GetXaxis()->SetBinLabel(4,"BBCMB");
  
  h1->GetXaxis()->SetBinLabel(6,"VPDMB");
  h2->GetXaxis()->SetBinLabel(6,"VPDMB");
  h3->GetXaxis()->SetBinLabel(6,"VPDMD");
  h4->GetXaxis()->SetBinLabel(6,"VPDMB");
  
  h1->GetYaxis()->SetTitle("Counts");
  h2->GetYaxis()->SetTitle("Counts");
  h3->GetYaxis()->SetTitle("Counts");
  h4->GetYaxis()->SetTitle("Counts");
  
  h1->GetXaxis()->SetLabelSize(0.06);
  h2->GetXaxis()->SetLabelSize(0.06);
  h3->GetXaxis()->SetLabelSize(0.06);
  h4->GetXaxis()->SetLabelSize(0.06);



  h1->GetYaxis()->SetTitleOffset(1.3);
  h2->GetYaxis()->SetTitleOffset(1.3);
  h3->GetYaxis()->SetTitleOffset(1.3);
  h4->GetYaxis()->SetTitleOffset(1.3);

  // h1->GetXaxis()->SetLabelOffset(0.002);
  // h2->GetXaxis()->SetLabelOffset(0.002);
  // h3->GetXaxis()->SetLabelOffset(0.002);
  // h4->GetXaxis()->SetLabelOffset(0.002);

  
  h1->SetTitle("EMC triggered #Events after prescle");
  h2->SetTitle("EMC triggered #Events before prescle");
  h3->SetTitle("MinBias triggered #Events after prescle");
  h4->SetTitle("MinBias triggered #Events before prescle");

  h1->Draw();
  nEvent->Draw("SAMEBAR");
  sprintf(buf,"HT0:%E",nEvent->GetBinContent(2));
  drawLatex(0.2,0.85,buf,62,0.035,2);

  sprintf(buf,"HT2:%E",nEvent->GetBinContent(4));
  drawLatex(0.2,0.8,buf,62,0.035,2);

  // sprintf(buf,"MB:%E",nEvents_ps->GetBinContent(6));
  // drawLatex(0.2,0.75,buf,62,0.035,2);
  c2->SaveAs("nEvents.pdf");

  h2->Draw();
  nEvent_ps->Draw("SAMEBAR");
  sprintf(buf,"HT0:%E",nEvent_ps->GetBinContent(2));
  drawLatex(0.2,0.85,buf,62,0.035,2);

  sprintf(buf,"HT2:%E",nEvent_ps->GetBinContent(4));
  drawLatex(0.2,0.8,buf,62,0.035,2);

  // sprintf(buf,"MB:%E",nEvents_ps->GetBinContent(6));
  // drawLatex(0.2,0.75,buf,62,0.035,2);
  c2->SaveAs("nEvent_ps.pdf");

  h3->Draw();
  nEvents_MB->Draw("SAMEBAR");
  sprintf(buf,"HT0:%E",nEvents_MB->GetBinContent(2));
  drawLatex(0.2,0.85,buf,62,0.035,2);

  sprintf(buf,"HT2:%E",nEvents_MB->GetBinContent(4));
  drawLatex(0.2,0.8,buf,62,0.035,2);

  sprintf(buf,"MB:%E",nEvents_MB->GetBinContent(6));
  drawLatex(0.2,0.75,buf,62,0.035,2);
  c2->SaveAs("nEvent_MB.pdf");

  h4->Draw();
  nEvents_MBps->Draw("SAMEBAR");
  sprintf(buf,"HT0:%E",nEvents_MBps->GetBinContent(2));
  drawLatex(0.2,0.85,buf,62,0.035,2);

  sprintf(buf,"HT2:%E",nEvents_MBps->GetBinContent(4));
  drawLatex(0.2,0.8,buf,62,0.035,2);

  sprintf(buf,"MB:%E",nEvents_MBps->GetBinContent(6));
  drawLatex(0.2,0.75,buf,62,0.035,2);
  c2->SaveAs("nEvent_MBps.pdf");
}
TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
  {
    TLatex *latex = new TLatex(x,y,text); 
    latex->SetNDC();
     latex->SetTextFont(textFont);
    latex->SetTextSize(textSize);   
    latex->SetTextColor(colorIndex); 
    latex->Draw("same");
    return latex;
  }


