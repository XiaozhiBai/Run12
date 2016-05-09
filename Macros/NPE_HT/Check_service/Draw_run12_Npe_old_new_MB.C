#include <iostream>
#include<iomanip>
#include <fstream>
#include "TLatex.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

// #include "../mBinning_HT.h"
// #include "../../NPE_MB_wide_bin/mBinning_MB.h"

using namespace std;

Int_t gColor08=1;
Int_t gColor08_fill=1;
Int_t gMStyle08=20;


Int_t gColor12HT=2;
Int_t gColor12HT_fill=42;
Int_t gMStyle12HT=20;

Int_t gColor12MB=2;
Int_t gColor12MB_fill=42;
Int_t gMStyle12MB=24;


Int_t lWidth=3;
Int_t lStyle=2;
Int_t lStyle_c=1;
Int_t lColor=4;

TLatex* drawLatex(Double_t, Double_t, char* , Int_t , Double_t , Int_t);
char buf[1024];
void Draw_run12_Npe_old_new_MB()

{

  char buf[1024];
  TH1::SetDefaultSumw2();
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleOffset(1,"X");
  gStyle->SetTitleOffset(1,"Y");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13); 
  //  gStyle->SetEndErrorSize(6);
  gStyle->SetOptStat(000);


  //
  // FONLL
  TFile *file_data_fonll=new TFile("Data_Model/fonll200gev.root","READ");
  TGraphErrors *gFONLLu=( TGraphErrors *) file_data_fonll->Get("gFONLLu");
  TGraphErrors *gFONLLc=( TGraphErrors *) file_data_fonll->Get("gFONLLc");
  TGraphErrors *gFONLLl=( TGraphErrors *) file_data_fonll->Get("gFONLLl");

  TGraphErrors *gR_FONLLu= (TGraphErrors *) file_data_fonll->Get("gR_FONLLu");
  TGraphErrors *gR_FONLLc= (TGraphErrors *) file_data_fonll->Get("gR_FONLLc");
  TGraphErrors *gR_FONLLl= (TGraphErrors *) file_data_fonll->Get("gR_FONLLl");


  gFONLLu->SetLineStyle(lStyle);
  gFONLLc->SetLineStyle(lStyle_c);
  gFONLLl->SetLineStyle(lStyle);

  gFONLLu->SetLineWidth(lWidth);
  gFONLLc->SetLineWidth(lWidth);
  gFONLLl->SetLineWidth(lWidth);

  gFONLLu->SetLineColor(lColor);
  gFONLLc->SetLineColor(lColor);
  gFONLLl->SetLineColor(lColor);

  // run12 HT
  TFile *file_run12_HT=new TFile("Data_Model/run12_Npe_MB.root","READ");

  TH1F *run12_HT_NPE_sts=(TH1F *) file_run12_HT->Get("run12_MB_NPE_sts");
  TH1F *run12_HT_NPE_sys=(TH1F *) file_run12_HT->Get("run12_MB_NPE_sys");
  
  TH1F *run12_HT_NPE_FONLL_sts=(TH1F *)  file_run12_HT->Get("run12_MB_NPE_FONLL_sts"); 
  TH1F *run12_HT_NPE_FONLL_sys=(TH1F *)  file_run12_HT->Get("run12_MB_NPE_FONLL_sys"); 
  
  // run12
  TGraphErrors* gr_run12_HT_sts = new TGraphErrors(run12_HT_NPE_sts);
  TGraphErrors* gr_run12_HT_sys = new TGraphErrors(run12_HT_NPE_sys);

  TGraphErrors* gr_run12_HT_FONLL_sts = new TGraphErrors(run12_HT_NPE_FONLL_sts);
  TGraphErrors* gr_run12_HT_FONLL_sys = new TGraphErrors(run12_HT_NPE_FONLL_sys);



  gr_run12_HT_sys->SetLineColor(gColor12HT);
  gr_run12_HT_sts->SetMarkerStyle(gMStyle12HT);
  gr_run12_HT_sts->SetMarkerColor(gColor12HT);
  gr_run12_HT_sts->SetLineColor(gColor12HT);
  gr_run12_HT_sts->SetMarkerSize(1);
  gr_run12_HT_sys->SetMarkerSize(1);




  gr_run12_HT_FONLL_sys->SetLineColor(gColor12HT);
  gr_run12_HT_FONLL_sts->SetMarkerStyle(gMStyle12HT);
  gr_run12_HT_FONLL_sts->SetMarkerColor(gColor12HT);
  gr_run12_HT_FONLL_sts->SetLineColor(gColor12HT);
  gr_run12_HT_FONLL_sts->SetMarkerSize(1);
  gr_run12_HT_FONLL_sys->SetMarkerSize(1);



  for(Int_t ipt=0;ipt<run12_HT_NPE_sts->GetNbinsX();ipt++)
    {
      gr_run12_HT_sys->SetPointError(ipt,0.04,gr_run12_HT_sys->GetErrorY(ipt));
      gr_run12_HT_sts->SetPointError(ipt,0,gr_run12_HT_sts->GetErrorY(ipt));
      
      gr_run12_HT_FONLL_sys->SetPointError(ipt,0.04,gr_run12_HT_FONLL_sys->GetErrorY(ipt));
      gr_run12_HT_FONLL_sts->SetPointError(ipt,0,gr_run12_HT_FONLL_sts->GetErrorY(ipt));
    }

  gr_run12_HT_sys->SetLineColor(gColor12HT_fill);
  gr_run12_HT_sys->SetFillStyle(1001);
  gr_run12_HT_sys->SetFillColor(gColor12HT_fill);

  gr_run12_HT_FONLL_sys->SetLineColor(gColor12HT_fill);
  gr_run12_HT_FONLL_sys->SetFillStyle(1001);
  gr_run12_HT_FONLL_sys->SetFillColor(gColor12HT_fill);


  TGraphErrors* gr_run08_HT_sts = new TGraphErrors("Data_Model/NPE_Run08_subtract_Jpsi.dat","%lg %lg %lg");
  TGraphErrors* gr_run08_HT_sys = new TGraphErrors("Data_Model/NPE_Run08_subtract_Jpsi.dat","%lg %lg %*lg %lg");

  TFile *Run08_data_fonll=new TFile("Data_Model/Run08_data_fonll.root","READ");
  TGraphErrors *gr_run08_HT_FONLL_sts=(TGraphErrors *) Run08_data_fonll->Get("Graph_from_run08_sts_ratio");
  TGraphErrors *gr_run08_HT_FONLL_sys=(TGraphErrors *) Run08_data_fonll->Get("Graph_from_run08_sys_ratio");
  
  gr_run08_HT_sts->SetMarkerStyle(gMStyle08);
  gr_run08_HT_sts->SetMarkerColor(gColor08);
  gr_run08_HT_sts->SetLineColor(gColor08);
  gr_run08_HT_sys->SetMarkerColor(gColor08);
  gr_run08_HT_sys->SetLineColor(gColor08);
  gr_run08_HT_sts->SetMarkerSize(1);
  gr_run08_HT_sys->SetMarkerSize(1);


   gr_run08_HT_FONLL_sts->SetMarkerSize(1);
   gr_run08_HT_FONLL_sys->SetMarkerSize(1);
   gr_run08_HT_FONLL_sts->SetMarkerStyle(gMStyle08);
   gr_run08_HT_FONLL_sys->SetMarkerStyle(gMStyle08);
   gr_run08_HT_FONLL_sts->SetMarkerColor(gColor08);
   gr_run08_HT_FONLL_sys->SetMarkerColor(gColor08);
   gr_run08_HT_FONLL_sts->SetLineColor(gColor08);
   gr_run08_HT_FONLL_sys->SetLineColor(gColor08);
   gr_run08_HT_FONLL_sys->SetLineWidth(1);
   gr_run08_HT_FONLL_sys->SetFillStyle(3001);
   gr_run08_HT_FONLL_sys->SetFillColor(gColor08_fill);


  for(Int_t ipt=0;ipt<9;ipt++)
    {
      gr_run08_HT_sts->SetPointError(ipt,0,gr_run08_HT_sts->GetErrorY(ipt));
      gr_run08_HT_sys->SetPointError(ipt,0.1,gr_run08_HT_sys->GetErrorY(ipt));
    }


   for(Int_t ipt=0;ipt<9;ipt++)
   {
     gr_run08_HT_FONLL_sts->SetPointError(ipt,0,gr_run08_HT_FONLL_sts->GetErrorY(ipt));
     gr_run08_HT_FONLL_sys->SetPointError(ipt,0.1,gr_run08_HT_FONLL_sys->GetErrorY(ipt));
   }

   // for(Int_t ipt=0;ipt<9;ipt++)
   //   {
   //     Double_t x=0,y=0;
   //     gr_run08_HT_sts->GetPoint(ipt,x,y);
   //     gr_run08_HT_sts->SetPoint(ipt,x-0.12,y);
   //     gr_run08_HT_sys->SetPoint(ipt,x-0.12,y);
   //   }

   for(Int_t ipt=0;ipt<9;ipt++)
     {
       Double_t x=0,y=0;
       gr_run08_HT_FONLL_sts->GetPoint(ipt,x,y);
       gr_run08_HT_FONLL_sts->SetPoint(ipt,x-0.12,y);
       gr_run08_HT_FONLL_sys->SetPoint(ipt,x-0.12,y);
     }


   gr_run08_HT_sys->SetFillStyle(3001);
   gr_run08_HT_sys->SetFillColor(gColor08_fill);

   //phex
  //  TGraphErrors* PhNEX_sts = new TGraphErrors("Data_Model/Npe_Phenix_cross_section.dat","%lg %lg %lg");
  // TGraphErrors* PhNEX_sys = new TGraphErrors("Data_Model/Npe_Phenix_cross_section.dat","%lg %lg %*lg %lg");


  TGraphErrors* PhNEX_sts = new TGraphErrors("Data_Model/NPE_cross_section_MB.dat","%lg %lg %lg");
  TGraphErrors* PhNEX_sys = new TGraphErrors("Data_Model/NPE_cross_section_MB.dat","%lg %lg %*lg %lg");

  
  PhNEX_sts->SetMarkerStyle(20);
  PhNEX_sys->SetMarkerStyle(20);
  
  PhNEX_sts->SetName("PhNEX_sts");
  PhNEX_sys->SetName("PhNEX_sys");

  PhNEX_sts->SetMarkerColor(3);
  PhNEX_sys->SetMarkerColor(3);

  PhNEX_sts->SetLineColor(1);
  PhNEX_sys->SetLineColor(3);

 


  //  TGraphErrors* PhNEX_sts_fonll = new TGraphErrors("Data_Model/Npe_Phenix_cross_section_fonll_ratio.dat","%lg %lg %lg");
  // TGraphErrors* PhNEX_sys_fonll = new TGraphErrors("Data_Model/Npe_Phenix_cross_section_fonll_ratio.dat","%lg %lg %*lg %lg");

  //    TGraphErrors* PhNEX_sts_fonll = new TGraphErrors("Data_Model/NPE_Ratio12.dat","%lg %lg %lg");
  // TGraphErrors* PhNEX_sys_fonll = new TGraphErrors("Data_Model/NPE_Ratio12.dat","%lg %lg %*lg %lg");

  TFile *File=new TFile("Data_Model/data_fonll.root","READ");

  TGraphErrors *PhNEX_sts_fonll=(TGraphErrors *) File->Get("Ratio_sts");
  TGraphErrors *PhNEX_sys_fonll=(TGraphErrors *) File->Get("Ratio_sys");

  
  PhNEX_sts_fonll->SetMarkerStyle(20);
  PhNEX_sys_fonll->SetMarkerStyle(20);
  
  PhNEX_sts_fonll->SetName("PhNEX_sts_fonll");
  PhNEX_sys_fonll->SetName("PhNEX_sys_fonll");

  PhNEX_sts_fonll->SetMarkerColor(3);
  PhNEX_sys_fonll->SetMarkerColor(3);

  PhNEX_sts_fonll->SetLineColor(1);
  PhNEX_sys_fonll->SetLineColor(3);

  

  TH2D *h_1 = new TH2D("h_1","",150,0.2,4,100,0.5e-8,1);
  h_1->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  h_1->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})"); 
  // h_1->GetXaxis()->CenterTitle();  
  // h_1->GetYaxis()->CenterTitle();  

  h_1->GetYaxis()->SetTitleSize(0.06);
  h_1->GetYaxis()->SetTitleOffset(0.8);

  h_1->GetYaxis()->SetLabelSize(0.05);

  for(Int_t i=1;i<gFONLLu->GetN();i++) {
    gFONLLu->SetPointError(i,0,0);
    gFONLLc->SetPointError(i,0,0);
    gFONLLl->SetPointError(i,0,0);
  }

  for(Int_t i=1;i<gR_FONLLu->GetN();i++) {
    gR_FONLLu->SetPointError(i,0,0);
    gR_FONLLc->SetPointError(i,0,0);
    gR_FONLLl->SetPointError(i,0,0);
  }


  TCanvas *cc1 = new TCanvas("cc1","cc1",900,1000);
  cc1->Divide(1,2);
  cc1->cd(1)->SetPad(0.01,0.4,0.99,0.99);
  cc1->cd(1)->SetBottomMargin(0.0);
  cc1->cd(1)->SetTopMargin(0.02);
  cc1->cd(1)->SetLeftMargin(0.1);
  gPad->SetLogy();
  h_1->DrawCopy();


   // gr_run08_HT_sys->Draw("sameE2");
   // gr_run08_HT_sts->Draw("sameP");


   // gr_run12_MB_sys->Draw("sameE2");
   // gr_run12_MB_sts->Draw("sameP");
   //   run12_MB_NPE_sts_all_new->Draw("samePE");
     
   gr_run12_HT_sys->Draw("sameE2");
   gr_run12_HT_sts->Draw("sameP");


   PhNEX_sts->Draw("sameP");
   PhNEX_sys->Draw("p[]:same");


  gFONLLu->Draw("sameL");
  gFONLLc->Draw("sameL");
  gFONLLl->Draw("sameL");

  TLegend *legend_2  = new TLegend(0.7,0.5,0.95,0.95);
  //    legend_2 ->AddEntry(gr_run12_MB_sts,"Run12 MinBias","lp");
  legend_2 ->AddEntry(gr_run12_HT_sts,"Run12 MB new","lp");
  // legend_2 ->AddEntry(data09_stat_err,"Run09","lp");
  //  legend_2 ->AddEntry(gr_run08_HT_sts,"Run05+Run08","lp");
   legend_2 ->AddEntry(PhNEX_sts,"Run12 MB old","lp");
  legend_2 ->AddEntry(gFONLLc,"FONLL","lp");
  //  legend_2 ->AddEntry(gFONLLl,"FONLL Uncertainty","lp");
  legend_2 ->SetBorderSize(0);
  legend_2 ->SetTextSize(0.05);
  legend_2 ->SetFillStyle(0);
  legend_2->SetTextFont(62);
  legend_2->Draw();



  sprintf(buf,"STAR Preliminary");
  drawLatex(0.25,0.85,buf,62,0.09,2);

  sprintf(buf,"p+p@#sqrt{s}=200 GeV");


  sprintf(buf,"|#eta|<0.7");
  drawLatex(0.76,0.43,buf,62,0.077,1);

 
  cc1->cd(2)->SetTopMargin(0.0);
  cc1->cd(2)->SetBottomMargin(0.14);
  cc1->cd(2)->SetLeftMargin(0.1);
  gPad->SetLogy();
  TLine *line1 = new TLine(0.2,1,13,1);
  line1->SetLineStyle(lStyle_c);
  line1->SetLineColor(lColor);
  line1->SetLineWidth(lWidth);
  

  TH2D *h_2 = new TH2D("h_2","",150,0.2,4,100,0.25,5.999);
  // h_frame->GetYaxis()->SetRangeUser(0,6);
  // h_frame->SetTitle("ratio");
  h_2->GetYaxis()->SetTitle("Data/FONLL");
  h_2->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");

  h_2->GetYaxis()->SetTitleSize(0.08);
  h_2->GetXaxis()->SetTitleSize(0.08);
  h_2->GetYaxis()->SetTitleOffset(0.6);
  h_2->GetXaxis()->SetTitleOffset(0.8);

  h_2->GetXaxis()->SetLabelSize(0.06);
  h_2->GetYaxis()->SetLabelSize(0.06);





  // h_2->GetYaxis()->CenterTitle();
  // h_2->GetXaxis()->CenterTitle();


  h_2->DrawCopy();

  gR_FONLLu->SetLineWidth(lWidth);
  gR_FONLLl->SetLineWidth(lWidth);

  gR_FONLLu->SetLineStyle(lStyle);
  gR_FONLLl->SetLineStyle(lStyle);

  gR_FONLLu->SetLineColor(lColor);
  gR_FONLLl->SetLineColor(lColor);

  // gr_run08_HT_FONLL_sts->Draw("sameP");
  // gr_run08_HT_FONLL_sys->Draw("sameE2");


   // gr_run12_MB_FONLL_sys->Draw("sameE2");
   // gr_run12_MB_FONLL_sts->Draw("sameP");



   //   run12_MB_NPE_FONLL_sts_new->Draw("samePE");


   gr_run12_HT_FONLL_sys->Draw("sameE2");
   gr_run12_HT_FONLL_sts->Draw("sameP");


  gR_FONLLu->Draw("sameL");
  gR_FONLLl->Draw("sameL");
  line1->Draw("same");

  PhNEX_sts_fonll->Draw("samePE");
  PhNEX_sys_fonll->Draw("p[]:same");

  cc1->SaveAs("run12_cross_section_MB.pdf");
  // cc1->SaveAs("run12_cross_section.png");
  // cc1->SaveAs("run12_cross_section.gif");
}

  
  

//   TCanvas *c2=new TCanvas("c2","",800,600);
//   gPad->SetLogy();
//   TH2F *h2=new TH2F("h2","",100,0.2,4,100,1e-8,0.2);
//   h2->Draw();

//   h2->SetTitle("NPE Cross Section");
//   h2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//   h2->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})");
//   h2->GetXaxis()->CenterTitle(true);
//   h2->GetYaxis()->CenterTitle(true);
//   h2->Draw();



   
//    TFile *file_MB2012=new TFile("Data_Model/file_MB2012.root","READ");

//    TH1F *run12MB_sts_comb=(TH1F *) file_MB2012->Get("run12MB_sts_comb");
//    TH1F *run12MB_sys_comb=(TH1F *) file_MB2012->Get("run12MB_sys_comb");



//    TGraphErrors *gr_run12_Comb_sts=new TGraphErrors(run12MB_sts_comb);
//    TGraphErrors *gr_run12_Comb_sys=new TGraphErrors(run12MB_sys_comb);
   

//    gr_run12_Comb_sts->SetMarkerStyle(gMStyle12MB);
//    gr_run12_Comb_sys->SetMarkerStyle(gMStyle12MB);
//    gr_run12_Comb_sts->SetMarkerColor(gColor12MB);
//    gr_run12_Comb_sys->SetMarkerColor(gColor12MB);
//    gr_run12_Comb_sts->SetLineColor(gColor12MB);
//    gr_run12_Comb_sys->SetLineColor(gColor12MB);
   
//    for(Int_t ipt=0;ipt<Nbins_AA;ipt++)
//     {
//       gr_run12_Comb_sts->SetPointError(ipt,0, gr_run12_Comb_sts->GetErrorY(ipt));
//       gr_run12_Comb_sys->SetPointError(ipt,0.09, gr_run12_Comb_sys->GetErrorY(ipt));

//     }
//    gr_run12_Comb_sys->SetLineColor(gColor12MB_fill);
//    gr_run12_Comb_sys->SetFillStyle(1001);
//    gr_run12_Comb_sys->SetFillColor(gColor12MB_fill);





//   gr_run12_Comb_sts->Draw("sameP");
//   gr_run12_Comb_sys->Draw("p:[]same");




// fstream  outfile("NPE_cross_section_MB_combine.dat",ios::out|ios::trunc);
 
//  for(Int_t ipt=0;ipt<Nbins_AA;ipt++){
//    outfile<<run12MB_sts_comb ->GetBinCenter(ipt+1)<<"  " <<run12MB_sts_comb ->GetBinContent(ipt+1)<<"  "<<run12MB_sts_comb->GetBinError(ipt+1)<< "  "<<run12MB_sys_comb->GetBinError(ipt+1)<<endl;
//  }


//   //..............PHEX

//   TGraphErrors* PhNEX_sts = new TGraphErrors("phenixresult.dat","%lg %lg %lg");
//   TGraphErrors* PhNEX_sys = new TGraphErrors("phenixresult.dat","%lg %lg %*lg %lg");
  
//   PhNEX_sts->SetMarkerStyle(20);
//   PhNEX_sys->SetMarkerStyle(20);
  
//   PhNEX_sts->SetName("PhNEX_sts");
//   PhNEX_sys->SetName("PhNEX_sys");

//   PhNEX_sts->SetMarkerColor(3);
//   PhNEX_sys->SetMarkerColor(3);

//   PhNEX_sts->SetLineColor(3);
//   PhNEX_sys->SetLineColor(4);

//   PhNEX_sts->Draw("sameP");
//   PhNEX_sys->Draw("p[]:same");


//   TLegend *legend  = new TLegend(0.5,0.5,0.8,0.85);
//   legend ->AddEntry(gr_run12MB_sts,"Run12 MB xiaozhi","lp");
//   legend ->AddEntry(gr_run12MB_1sts,"Run12 MB shenghui","lp");
//   legend ->AddEntry(gr_run12_Comb_sts,"Combine","lp");

//   legend ->AddEntry(PhNEX_sts,"PHENIX","lp");
//   legend ->AddEntry(gFONLLc,"FONLL","lp");
//   legend ->AddEntry(gFONLLl,"FONLL Uncertainty","lp");
  
//   legend ->SetBorderSize(0);
//   legend ->SetTextSize(0.04);
//   legend ->SetFillStyle(0);
//   legend->SetTextFont(62);
//   legend ->Draw("same");





//   c2->SaveAs("run12MBcrossSection.pdf");




//   PhNEX_sts->GetPoint(2,a,b);
//   cout<< a<<" "<<b<<endl;

//   TCanvas *c_Ratio=new TCanvas("c_Ratio","",800,600);

//   TH2F *h3=new TH2F("h3","",100,0,4,100,0,6);
//   h3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//   h3->GetYaxis()->SetTitle("Data/FONLL");
  
//   h3->GetXaxis()->CenterTitle();
//   h3->GetYaxis()->CenterTitle();  
//   h3->Draw();

//   run12MB_sts_comb_datafonll->Add(run12MB_sts_comb);
//   run12MB_sys_comb_datafonll->Add(run12MB_sys_comb);


//   run12MB_sts_comb_datafonll->Divide(FONLL);
//   run12MB_sys_comb_datafonll->Divide(FONLL);

//   cout<< "Ratio"<< run12MB_sts_comb_datafonll->GetNbinsX()<<endl;
// //   for(Int_t ipt=0;ipt<Nbins_AA;ipt++)
// //     {
// //       cout << run12MB_sts_comb_datafonll->GetBinCenter(ipt+1)<< " " << run12MB_sts_comb_datafonll->GetBinContent(ipt+1)<< "  "<< run12MB_sts_comb_datafonll->GetBinError(ipt+1)<<"  "<<
// // 	run12MB_sys_comb_datafonll->GetBinError(ipt+1) <<"  "<<
// // 	run12MB_sys_comb->GetBinError(ipt+1)<<" "<<
// // 	FONLL->GetBinContent(ipt+1)<<
// // endl;;
// //     }
// //   return;
  

//   TGraphErrors *gR_FONLLu= (TGraphErrors *) file_data_fonll->Get("gR_FONLLu");
//   TGraphErrors *gR_FONLLc= (TGraphErrors *) file_data_fonll->Get("gR_FONLLc");
//   TGraphErrors *gR_FONLLl= (TGraphErrors *) file_data_fonll->Get("gR_FONLLl");

//   for(Int_t i=1;i<gR_FONLLu->GetN();i++) {
//     gR_FONLLu->SetPointError(i,0,0);
//     gR_FONLLc->SetPointError(i,0,0);
//     gR_FONLLl->SetPointError(i,0,0);
//   }
  

//   gR_FONLLu->SetLineStyle(2);
//   gR_FONLLc->SetLineStyle(2);
//   gR_FONLLl->SetLineStyle(2);
  
//   gR_FONLLu->SetLineWidth(1);
//   gR_FONLLc->SetLineWidth(1.5);
//   gR_FONLLl->SetLineWidth(1);
  
//   gR_FONLLu->SetLineColor(1);
//   gR_FONLLc->SetLineColor(1);
//   gR_FONLLl->SetLineColor(1);
  

//    TLine *line=new TLine(0,1,10,1);
//    line->SetLineStyle(7);
//    gR_FONLLu->Draw("same");
//    gR_FONLLl->Draw("same");
//    line->Draw("same");


//    TFile *file_DataFONLL=new TFile("Data/DataFONLL.root","READ");
//    TGraph *phenix_ratio=(TGraph *)file_DataFONLL->Get("phenix");
//    TGraph *phenix_ratio_sys=(TGraph *)file_DataFONLL->Get("phenixSys");

//    phenix_ratio->SetLineColor(3);
//    phenix_ratio->SetMarkerColor(3);
//    phenix_ratio_sys->SetLineColor(4);

//    phenix_ratio->Draw("sameP");
//    phenix_ratio_sys->Draw("p[]:same");
   

//    //--------------------------run08--------------------------------------------

//    //--------------------------run08--------------------------------------------


   
//    //------------------run09------------------------------
//    TGraphErrors* Ratio09_stat_err = new TGraphErrors("Ratio09.dat","%lg %lg %lg");
//    TGraphErrors* Ratio09_sys_err = new TGraphErrors("Ratio09.dat","%lg %lg %*lg %lg");
//    Ratio09_stat_err->SetMarkerStyle(20);
//    Ratio09_sys_err->SetLineColor(4);
//    Ratio09_stat_err->SetLineColor(7);
//    Ratio09_stat_err->SetMarkerColor(7);
//    Ratio09_sys_err->SetMarkerColor(4);
//    Ratio09_stat_err->SetMarkerSize(1);
//    Ratio09_sys_err->SetMarkerSize(1);
//    Ratio09_stat_err->Draw("samePE1");
//    Ratio09_sys_err->Draw("P[]:same");
//    TFile *file_09=new TFile("NPE_pp09_04_10_GeV.root","READ");
//    TH1F * ratio_low_pt= ( TH1F *)file_09->Get("FONLL_ratio_control");
//    ratio_low_pt->SetBinContent(1,-999);
//    ratio_low_pt->SetBinContent(2,-999);
//    TGraphErrors *gr_Run09MB_data_fonll=new TGraphErrors(ratio_low_pt);
//    gr_Run09MB_data_fonll->SetMarkerStyle(20);
//    gr_Run09MB_data_fonll->SetMarkerColor(7);
//    gr_Run09MB_data_fonll->SetLineColor(7);
//    gr_Run09MB_data_fonll->SetMarkerSize(1);
//    // gr_Run09MB_data_fonll->Draw("samePE1");

//    //------------------run09------------------------------


//   //----------------------------run12-----------------------------------------

//    // for(Int_t ipt=0;ipt<Nbins_AA;ipt++)
//    //   {


//    //   }

//   TGraphErrors *  gr_run12MB_sts_comb_datafonll =new TGraphErrors(run12MB_sts_comb_datafonll); 
//   TGraphErrors *  gr_run12MB_sys_comb_datafonll =new TGraphErrors(run12MB_sys_comb_datafonll); 
//   gr_run12MB_sts_comb_datafonll->SetMarkerStyle(gMStyle12MB);
//   gr_run12MB_sys_comb_datafonll->SetMarkerStyle(gMStyle12MB);

//   gr_run12MB_sts_comb_datafonll->SetMarkerColor(gColor12MB);
//   gr_run12MB_sys_comb_datafonll->SetMarkerColor(gColor12MB);

//   gr_run12MB_sts_comb_datafonll->SetLineColor(gColor12MB);
//   gr_run12MB_sys_comb_datafonll->SetLineColor(gColor12MB);

//   gr_run12MB_sts_comb_datafonll->SetMarkerSize(1);
//   gr_run12MB_sys_comb_datafonll->SetMarkerSize(1);


//   for(Int_t ipt=0;ipt<Nbins_AA;ipt++)
//     {
//       gr_run12MB_sts_comb_datafonll->SetPointError(ipt,0, gr_run12MB_sts_comb_datafonll->GetErrorY(ipt));
//       gr_run12MB_sys_comb_datafonll->SetPointError(ipt,0.08, gr_run12MB_sys_comb_datafonll->GetErrorY(ipt));

//     }

//   gr_run12MB_sys_comb_datafonll->SetLineColor(gColor12MB_fill);
//   gr_run12MB_sys_comb_datafonll->SetFillStyle(1001);
//   gr_run12MB_sys_comb_datafonll->SetFillColor(gColor12MB_fill);
  




//   gr_run12MB_sts_comb_datafonll->Draw("sameP");
//   gr_run12MB_sys_comb_datafonll->Draw("p[]:same");
   
//   TFile *Run12_data_fonll=new TFile("Run12HT_data_fonll.root","READ");
//   // TGraphErrors *Graph_from_run12_sts_ratio=(TGraphErrors *) Run12_data_fonll->Get("Graph_from_run12_sts_ratio");
//   // TGraphErrors *Graph_from_run12_sys_ratio=(TGraphErrors *) Run12_data_fonll->Get("Graph_from_run12_sys_ratio");

//   TGraphErrors* Graph_from_run12_sts_ratio = new TGraphErrors("Data/NPE_Ratio12.dat","%lg %lg %lg");
//   TGraphErrors* Graph_from_run12_sys_ratio  = new TGraphErrors("Data/NPE_Ratio12.dat","%lg %lg %*lg %lg");
   

//    Graph_from_run12_sts_ratio->SetMarkerStyle(gMStyle12HT);
//    Graph_from_run12_sys_ratio->SetMarkerStyle(gMStyle12HT);


//    Graph_from_run12_sts_ratio->SetMarkerSize(1);
//    Graph_from_run12_sys_ratio->SetMarkerSize(1);
//    Graph_from_run12_sts_ratio->SetMarkerColor(gColor12HT);
//    Graph_from_run12_sys_ratio->SetMarkerColor(gColor12HT);
//    Graph_from_run12_sts_ratio->SetLineColor(gColor12HT);
//    Graph_from_run12_sys_ratio->SetLineColor(gColor12HT);
//    // Graph_from_run12_sts_ratio->Draw("sameP");
//    // Graph_from_run12_sys_ratio->Draw("p[]:same");
//   //----------------------------run12-----------------------------------------





//   TLegend *legend_1  = new TLegend(0.5,0.5,0.8,0.85);
//   legend_1 ->AddEntry(gr_run12MB_sts_comb_datafonll,"Run12 MB","lp");
//   legend_1 ->AddEntry(Graph_from_run12_sys_ratio,"Run12 HT","lp");
//   // legend_1 ->AddEntry(Ratio09_stat_err,"Run09 HT","lp");
//   legend_1 ->AddEntry(gr_Run09MB_data_fonll,"Run09","lp");
//   legend_1 ->AddEntry(Graph_from_run08_sys_ratio,"Run08 HT","lp");
//   legend_1 ->AddEntry(phenix_ratio,"PHNEX","lp");

//   legend_1 ->AddEntry(line,"FONLL","l");
//   legend_1 ->AddEntry(gFONLLl,"FONLL Uncertainty","l");
  
//   legend_1 ->SetBorderSize(0);
//   legend_1 ->SetTextSize(0.04);
//   legend_1 ->SetFillStyle(0);
//   legend_1->SetTextFont(62);
//   legend ->Draw("same");
//   c_Ratio->SaveAs("Data_fonll_full.pdf");



//   // run12
//   TGraphErrors* data12_stat_err_HT = new TGraphErrors("Data/NPE_cross_section.dat","%lg %lg %lg");
//   TGraphErrors* data12_sys_err_HT = new TGraphErrors("Data/NPE_cross_section.dat","%lg %lg %*lg %lg");
  
//   data12_stat_err_HT->SetLineColor(1);
//   data12_sys_err_HT->SetLineColor(gColor12HT);
  
//   data12_stat_err_HT->SetMarkerStyle(gMStyle12HT);
//   data12_stat_err_HT->SetMarkerColor(gColor12HT);
  
//   data12_stat_err_HT->SetMarkerSize(1);
//   data12_sys_err_HT->SetMarkerSize(1);
  

//   for(Int_t ipt=0;ipt<13;ipt++)
//     {

//       data12_sys_err_HT->SetPointError(ipt,0.12,data12_sys_err_HT->GetErrorY(ipt));
   
//         Graph_from_run12_sys_ratio->SetPointError(ipt,0.12,Graph_from_run12_sys_ratio->GetErrorY(ipt));
//       //   Graph_from_run12_sts_ratio->SetPointError(ipt,0,Graph_from_run12_sts_ratio->GetErrorY(ipt));
      
//     }

//   data12_sys_err_HT->SetLineColor(gColor12HT_fill);
//   data12_sys_err_HT->SetFillStyle(1001);
//   data12_sys_err_HT->SetFillColor(gColor12HT_fill);

//   Graph_from_run12_sys_ratio->SetLineColor(gColor12HT_fill);
//   Graph_from_run12_sys_ratio->SetFillStyle(1001);
//   Graph_from_run12_sys_ratio->SetFillColor(gColor12HT_fill);








//   //09

//   TGraphErrors* data09_stat_err = new TGraphErrors("Data/NPERun09_corss_section.dat","%lg %lg %lg");
//   TGraphErrors* data09_sys_err = new TGraphErrors("Data/NPERun09_corss_section.dat","%lg %lg %*lg %lg");
  
//   data09_stat_err->SetMarkerStyle(20);
//   data09_stat_err->SetMarkerColor(7);
//   data09_stat_err->SetLineColor(7);
//   data09_sys_err->SetLineColor(4);

//   data09_stat_err->SetMarkerSize(1);
//   data09_sys_err->SetMarkerSize(1);

//   TFile *file_run09MB=new TFile("Data/Run09_cross_section.root");

//   TH1F *Run09MB=(TH1F *) file_09->Get("yield_kor");

//   TGraphErrors * gr_Run09MB=new TGraphErrors(Run09MB);
  
//   gr_Run09MB->SetMarkerStyle(20);
//   gr_Run09MB->SetMarkerColor(2);
//   gr_Run09MB->SetLineColor(2);
//   gr_Run09MB->SetMarkerSize(1);

//   TGraphErrors * gr_run09mb=(TGraphErrors *) file_run09MB->Get("Graph_from_yield_kor");

//   gr_run09mb->SetMarkerStyle(20);
//   gr_run09mb->SetMarkerColor(7);
//   gr_run09mb->SetLineColor(7);
//   gr_run09mb->SetMarkerSize(1);



  
//   //run08 

//   TGraphErrors* gr_stat_err = new TGraphErrors("Data/NPE_Run08_subtract_Jpsi.dat","%lg %lg %lg");
//   TGraphErrors* gr_sys_err = new TGraphErrors("Data/NPE_Run08_subtract_Jpsi.dat","%lg %lg %*lg %lg");
//   gr_stat_err->SetMarkerStyle(gMStyle08);
//   gr_stat_err->SetMarkerColor(gColor08);
//   gr_stat_err->SetLineColor(gColor08);

//   gr_sys_err->SetMarkerColor(gColor08);
//   gr_sys_err->SetLineColor(gColor08);
 

//   gr_stat_err->SetMarkerSize(1);
//   gr_sys_err->SetMarkerSize(1);



//    for(Int_t ipt=0;ipt<9;ipt++)
//    {
//      //dddd
//      gr_stat_err->SetPointError(ipt,0,gr_stat_err->GetErrorY(ipt));
//      gr_sys_err->SetPointError(ipt,0.12,gr_sys_err->GetErrorY(ipt));



//    }
//    gr_sys_err->SetFillStyle(3001);
//    gr_sys_err->SetFillColor(gColor08_fill);


//   //--------------------------------------Draw--------------------------


//   TLegend *legend_2  = new TLegend(0.7,0.5,0.95,0.95);
//   legend_2 ->AddEntry(gr_run12_Comb_sts,"Run12 MinBias","lp");
//   legend_2 ->AddEntry(data12_stat_err_HT,"Run12 BEMC ","lp");
//   // legend_2 ->AddEntry(data09_stat_err,"Run09","lp");
//   legend_2 ->AddEntry(gr_stat_err,"Run05+Run08","lp");
//   // legend_2 ->AddEntry(PhNEX_sts,"PHENIX ^{[2]}","lp");
//   legend_2 ->AddEntry(gFONLLc,"FONLL","lp");
//   //  legend_2 ->AddEntry(gFONLLl,"FONLL Uncertainty","lp");
//   legend_2 ->SetBorderSize(0);
//   legend_2 ->SetTextSize(0.07);
//   legend_2 ->SetFillStyle(0);
//   legend_2->SetTextFont(62);


  
//   TH2D *h_1 = new TH2D("h_1","",150,0.2,13,100,0.5e-13,0.1);
//   h_1->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
//   h_1->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mb Gev^{-2}c^{3})"); 
//   // h_1->GetXaxis()->CenterTitle();  
//   // h_1->GetYaxis()->CenterTitle();  

//   h_1->GetYaxis()->SetTitleSize(0.06);
//   h_1->GetYaxis()->SetTitleOffset(0.8);

//   h_1->GetYaxis()->SetLabelSize(0.05);

//   // h_1->GetXaxis()->SetLabelSize(0.08);




//   for(Int_t i=1;i<gFONLLu->GetN();i++) {
//     gFONLLu->SetPointError(i,0,0);
//     gFONLLc->SetPointError(i,0,0);
//     gFONLLl->SetPointError(i,0,0);
//   }


//   TCanvas *cc1 = new TCanvas("cc1","cc1",900,1000);
//   cc1->Divide(1,2);
//   cc1->cd(1)->SetPad(0.01,0.4,0.99,0.99);
//   cc1->cd(1)->SetBottomMargin(0.0);
//   cc1->cd(1)->SetTopMargin(0.02);
//   cc1->cd(1)->SetLeftMargin(0.1);

  
//   gPad->SetLogy();
//   h_1->DrawCopy();

//   gFONLLugFONLLu->SetLineStyle(lStyle);
//   gFONLLc->SetLineStyle(lStyle_c);
//   gFONLLl->SetLineStyle(lStyle);

//   gFONLLu->SetLineWidth(lWidth);
//   gFONLLc->SetLineWidth(lWidth);
//   gFONLLl->SetLineWidth(lWidth);

//   gFONLLu->SetLineColor(lColor);
//   gFONLLc->SetLineColor(lColor);
//   gFONLLl->SetLineColor(lColor);




//   gFONLLu->Draw("sameL");
//   gFONLLc->Draw("sameL");
//   gFONLLl->Draw("sameL");

//   //PHEIX

//   // PhNEX_sts->Draw("sameP");
//   // PhNEX_sys->Draw("p[]:same");
  

  
//   // 08
//   gr_sys_err->Draw("sameE2");
//   gr_stat_err->Draw("SAMEPE1");

  
//   // //09 
//   // data09_stat_err->Draw("SAMEP");
//   // data09_sys_err->Draw("P[]:same");  

//   //  gr_run09mb->Draw("samePE1");
//   //cc1->cd(2)->SetPad(0.01,0.01,0.99,0.4);
//   // run12
//   gr_run12_Comb_sys->Draw("sameE2");
//   gr_run12_Comb_sts->Draw("sameP");


//   data12_sys_err_HT->Draw("sameE2");
//   data12_stat_err_HT->Draw("sameP");


//   //  data12_sys_err_HT->Draw("P[]:same");


//   legend_2 ->Draw("same");



//   sprintf(buf,"STAR Preliminary");
//   drawLatex(0.25,0.85,buf,62,0.09,2);

//   sprintf(buf,"p+p@#sqrt{s}=200 GeV");
//   drawLatex(0.25,0.74,buf,62,0.08,1);

//   sprintf(buf,"|#eta|<0.7");
//   drawLatex(0.76,0.43,buf,62,0.077,1);

 
//   cc1->cd(2)->SetTopMargin(0.0);
//   cc1->cd(2)->SetBottomMargin(0.14);
//   cc1->cd(2)->SetLeftMargin(0.1);
//   gPad->SetLogy();
//   TLine *line1 = new TLine(0.2,1,13,1);
//   line1->SetLineStyle(lStyle_c);
//   line1->SetLineColor(lColor);
//   line1->SetLineWidth(lWidth);
  

//   TH2D *h_2 = new TH2D("h_2","",150,0.2,13,100,0.25,5.999);
//   // h_frame->GetYaxis()->SetRangeUser(0,6);
//   // h_frame->SetTitle("ratio");
//   h_2->GetYaxis()->SetTitle("Data/FONLL");
//   h_2->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");

//   h_2->GetYaxis()->SetTitleSize(0.08);
//   h_2->GetXaxis()->SetTitleSize(0.08);
//   h_2->GetYaxis()->SetTitleOffset(0.6);
//   h_2->GetXaxis()->SetTitleOffset(0.8);

//   h_2->GetXaxis()->SetLabelSize(0.06);
//   h_2->GetYaxis()->SetLabelSize(0.06);





//   // h_2->GetYaxis()->CenterTitle();
//   // h_2->GetXaxis()->CenterTitle();


//   h_2->DrawCopy();

//   gR_FONLLu->SetLineWidth(lWidth);
//   gR_FONLLl->SetLineWidth(lWidth);

//   gR_FONLLu->SetLineStyle(lStyle);
//   gR_FONLLl->SetLineStyle(lStyle);

//   gR_FONLLu->SetLineColor(lColor);
//   gR_FONLLl->SetLineColor(lColor);


//   gR_FONLLu->Draw("sameL");
//   gR_FONLLl->Draw("sameL");
  

//   // phenix_ratio->Draw("sameP");
//   // phenix_ratio_sys->Draw("p[]:same");

//   //08

//   for(Int_t ipt=0;ipt<9;ipt++)
//      {
//        Double_t x=0,y=0;
//        Graph_from_run08_sys_ratio->GetPoint(ipt,x,y);
//        Graph_from_run08_sys_ratio->SetPoint(ipt,x-0.1,y);
//        Graph_from_run08_sts_ratio->SetPoint(ipt,x-0.1,y);
//      }

//   Graph_from_run08_sys_ratio->Draw("sameE2");
//   Graph_from_run08_sts_ratio->Draw("sameP");

 
//  //09
//    // Ratio09_stat_err->Draw("samePE1");
//    // Ratio09_sys_err->Draw("P[]:same");
//    // gr_Run09MB_data_fonll->Draw("samePE1");
//    //run12
//   //gr_run12MB_sts_comb_datafonll->SetPoint();
//    gr_run12MB_sys_comb_datafonll->Draw("sameE2");    
//    gr_run12MB_sts_comb_datafonll->Draw("sameP");

//    for(Int_t ipt=0;ipt<13;ipt++)
//      {
//        Double_t x=0,y=0;
//        Graph_from_run12_sys_ratio->GetPoint(ipt,x,y);
//        Graph_from_run12_sys_ratio->SetPoint(ipt,x+0.1,y);
//        Graph_from_run12_sts_ratio->SetPoint(ipt,x+0.1,y);
//      }
   

//    Graph_from_run12_sys_ratio->Draw("sameE2");  
//    Graph_from_run12_sts_ratio->Draw("samePE");

//    //dev/  Graph_from_run12_sys_ratio->Draw("p[]:same");
 

//    line1->Draw("same");
//    cc1->SaveAs("pprun12.pdf");
//    cc1->SaveAs("pprun12.gif");
//    cc1->SaveAs("pprun12.png");
// }

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
