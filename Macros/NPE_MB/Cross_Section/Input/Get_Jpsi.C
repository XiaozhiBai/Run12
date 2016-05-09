#include <iostream>
#include <fstream>
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
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "binning_util.h"
#include "../../mBinning_HT.h"
//define Nbins_HT 13
using namespace std;

//Double_t Pt_bin_HT[14]={2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.5,10,14};

TH1F * JPsi=new TH1F("JPsi","",Nbins_HT,Pt_bin_HT);
TH1F *Unit_fector=new TH1F("Unit_fector","",Nbins_HT,Pt_bin_HT);
void Get_Jpsi()
{

  Float_t a,b,c,d,e;
  TH1F *Up=new TH1F("Up","",15,2.5,10);
  TH1F *Central=new TH1F("Central","",15,2.5,10);
  TH1F *Down=new TH1F("Down","",15,2.5,10);
 ifstream file;
  file.open("Jpisi.dat",ios::in);

while(!file.eof())
    { 
      file>>a>>b>>c>>d>>e;
 
  cout<<a<<"  "<<b<<"  "<<c<<endl;
  Up->SetBinContent(Up->FindBin((a+b)/2.0),c);
  Down->SetBinContent(Down->FindBin((a+b)/2.0),e);
  Central->SetBinContent(Central->FindBin((a+b)/2.0),d);
  Up->SetBinError(Up->FindBin((a+b)/2.0),abs(c-e)/2.0);
  Down->SetBinError(Down->FindBin((a+b)/2.0),abs(c-e)/2.0);
  Central->SetBinError(Central->FindBin((a+b)/2.0),abs(c-e)/2.0);

}

 TF1* f1 = new TF1("f1","[0]/(pow(x,[1])+[2])",2,16.5);
 f1->SetParameter(0,10);
 f1->SetParameter(1,8);
 f1->SetParameter(2,0.7); 

TF1* f2 = new TF1("f2","[0]/(pow(x,[1])+[2])",2,16.5);
 f2->SetParameter(0,10);
 f2->SetParameter(1,8);
 f2->SetParameter(2,0.7); 

TF1* f3 = new TF1("f3","[0]/(pow(x,[1])+[2])",2,16.5);
 f3->SetParameter(0,10);
 f3->SetParameter(1,8);
 f3->SetParameter(2,0.7); 

 TCanvas *c1=new TCanvas("c1","",800,600); 

 Up->Fit("f1","R","" ,2,16.5);
Central->Fit("f2","R","" ,2,16.5);
Down->Fit("f3","R","" ,2,16.5);
  
 for(Int_t i=0;i<Nbins_HT;i++)
   {
     if(i<11)
       {
	 Double_t err=fabs(Up->GetBinContent(i+1)-Down->GetBinContent(i+1))/2.0;

       JPsi->SetBinContent(i+1,Central->GetBinContent(i+1));
       JPsi->SetBinError(i+1,err);
       }
     if(11<=i)
       {
     Double_t Pt_bin=(Pt_bin_HT[i]+Pt_bin_HT[i+1])/2.0;
     Double_t jpsi=(f1->Eval(Pt_bin)+f3->Eval(Pt_bin))/2.0;
     Double_t jpsi_err=fabs(f1->Eval(Pt_bin)-f3->Eval(Pt_bin))/2.0;
     JPsi->SetBinContent(i+1,jpsi);
     JPsi->SetBinError(i+1,jpsi_err);
       }
     Unit_fector->SetBinContent(i+1,1E6);
     Unit_fector->SetBinError(i+1,0);
   }


 JPsi->Sumw2();
 JPsi->Divide(Unit_fector);
 JPsi->Draw("PE1");

 TFile *File=new TFile("Jpsi.root","recreate");
 JPsi->Write();
 File->Close();

 // Unit_fector->Draw();
// Up->Draw();
// Central->Draw("");
// Down->Draw("same");

}

