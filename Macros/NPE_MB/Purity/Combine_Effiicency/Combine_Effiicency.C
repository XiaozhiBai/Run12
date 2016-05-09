#include<fstream>
#include <iostream>
#include<iomanip>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "../../mBinning_MB.h"

using namespace std;
void Combine_Effiicency()
{

  TFile *file_puirty_low=new TFile("Inclusive_purity_MB_EMC.root","READ");
  TFile *file_puirty_high=new TFile("purity_MB.root","READ");

  TH1F *purity_low=(TH1F *) file_puirty_low->Get("mh1purity_44")->Clone("purity_low");
  TH1F *purity_high=(TH1F *) file_puirty_high->Get("purity_MB")->Clone("purity_high");
  
  purity_low->Draw();
  purity_high->Draw("sameP");
  
  TH1F *purity_MB_full=new TH1F("purity_MB_full","",NpT_bins_run12_MB,pt_run12_MB);

  for(int ipt=0;ipt<NpT_bins_run12_MB;ipt++)
    {
      if(ipt<20){
      purity_MB_full->SetBinContent(ipt+1,purity_low->GetBinContent(ipt+1));
      purity_MB_full->SetBinError(ipt+1,purity_low->GetBinError(ipt+1));
      cout<<purity_low->GetBinLowEdge(ipt+1)<<endl;;
      }
      else{
	purity_MB_full->SetBinContent(ipt+1,purity_high->GetBinContent(ipt+1-20));
	purity_MB_full->SetBinError(ipt+1,purity_high->GetBinError(ipt+1-20));
	cout<<purity_high->GetBinLowEdge(ipt+1-20)<<endl;;
      }
    }
  purity_MB_full->Draw();
  TFile *file_combine_purity=new TFile("purity_MB_com.root","RECREATE");
  purity_MB_full->Write();
  file_combine_purity->Close();
  
}
