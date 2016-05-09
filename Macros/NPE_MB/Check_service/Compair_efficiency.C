#include <iostream>
#include<iomanip>
#include <fstream>
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
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

using namespace std;

void Compair_efficiency()
{

  const int nhist=6;
  
  TH1F *efficiency_new[nhist];
  TH1F *efficiency_old[nhist];

  TFile *file_new[nhist];
  TFile *file_new_old[nhist];

  

  TString histname_new[nhist]={"efficiency_tof_match","nsigmaE_MB","Tof_beta_cut_MB","Tracking_efficiency_MB","PHE_re_efficiency","purity_MB"};
  TString histname_old[nhist]={"Tof_Match_E","nsigmaE_cut_efficiency","cut_eff","Tracking_Efficiency","Combine","purity"};

  TString file_name_new[nhist]={"TOF_Match_efficiency.root","nSigma_Cut_efficiency.root","Tof_Beta_cuts_efficiency.root","Tracking_efficiency_MB.root","Photonic_re_Efficiency.root","purity_MB.root"};
  TString file_name_old[nhist]={"ToF_Match_efficiency.root","Nsigma_Cut_efficiency.root","TOF_Cut_efficiency.root","Tracking_efficiency.root","Photonic_Electron_reconstruction_efficiency_combine_MB.root","Purity.root"};



  
  for(Int_t ihist=0;ihist<nhist;ihist++)
    {
      file_new[ihist]=new TFile(TString("../Cross_Section/Input/")+file_name_new[ihist],"READ");
      efficiency_new[ihist]=(TH1F *) file_new[ihist]->Get(histname_new[ihist]);
      
      file_new[ihist]=new TFile(TString("/Users/xiaozhi/NPE_work/NPE_low_pt/NPE_Spectra_applied_EMC_cut_9_18/")+file_name_old[ihist],"READ");
      //     cout<< file_new[ihist]<<endl;
      //efficiency_new[ihist]=(TH1F *) 
      
      
    }

  efficiency_new[2]->Draw();
  return ;
  
  // TFile *infile_data =new TFile("../Cross_Section/Input/hist_4_10.root","read");
  TFile *infile_data =new TFile("../Cross_Section/Input/hist_4_13.root","read");
  TFile *infile_Tof_match =new TFile("../Cross_Section/Input/TOF_Match_efficiency.root","read");
  TFile *infile_dEdx_cut =new TFile("../Cross_Section/Input/nSigma_Cut_efficiency.root","read");
  TFile *infile_Tof_cut =new TFile("../Cross_Section/Input/Tof_Beta_cuts_efficiency.root","read");
  TFile *infile_Tracking =new TFile("../Cross_Section/Input/Tracking_efficiency_MB.root","read");
  TFile *infile_PHE_re =new TFile("../Cross_Section/Input/Photonic_re_Efficiency.root","read");
  TFile *infile_purity =new TFile("../Cross_Section/Input/purity_MB.root","read");
  
  TH1F *eff_Tof_match=(TH1F *) infile_Tof_match->Get("efficiency_tof_match");
  TH1F *eff_dedx=(TH1F *) infile_dEdx_cut->Get("nsigmaE_MB");
  TH1F *eff_Tof_cut=(TH1F *) infile_Tof_cut->Get("Tof_beta_cut_MB");
  TH1F *eff_Tracking=(TH1F *) infile_Tracking->Get("Tracking_efficiency_MB");
  TH1F *eff_PHE_re=(TH1F *) infile_PHE_re->Get("PHE_re_efficiency");
  TH1F *eff_purity=(TH1F *) infile_purity->Get("purity_MB");

  eff_Tof_match->Draw();
  eff_dedx->Draw();
  eff_Tof_cut->Draw();
  eff_Tracking->Draw();
  eff_PHE_re->Draw();
  eff_purity->Draw();




  // old 
  
  TFile *infile_data_old =new TFile("/Users/xiaozhi/NPE_work/NPE_high_pt/Get_Inclsive_Photonic/hist_4_27.root","read");

  TFile *infile_Tof_match_old =new TFile("../../../NPE_high_pt_7_30/NPE_Cross_section/EMC_cut_efficiency.root","read");
  TFile *infile_dEdx_cut_old =new TFile("../../../NPE_high_pt_7_30/NPE_Cross_section/nSigma_Cut_efficiency.root","read");
  TFile *infile_Tof_cut_old =new TFile("../../../NPE_high_pt_7_30/NPE_Cross_section/TrigEfficiency.root","read");
  TFile *infile_Tracking_old =new TFile("../../../NPE_high_pt_7_30/NPE_Cross_section/Tracking_efficiency.root","read");
  TFile *infile_PHE_re_old =new TFile("../../../NPE_high_pt_7_30/NPE_Cross_section/Photonic_ReconStructionEfficiency.root","read");
  TFile *infile_purity_old =new TFile("../../../NPE_high_pt_7_30/NPE_Cross_section/Purity.root","read");
  
  TH1F *eff_Tof_match_old=(TH1F *) infile_Tof_match_old->Get("EMC_cut_effciency");
  TH1F *eff_dedx_old=(TH1F *) infile_dEdx_cut_old->Get("nsigma_cut_Efficiency");
  TH1F *eff_Tof_cut_HT0_old=(TH1F *) infile_Tof_cut_old->Get("hh_11_0");
  TH1F *eff_Tof_cut_HT2_old=(TH1F *) infile_Tof_cut_old->Get("hh_18_0");
  TH1F *eff_Tracking_old=(TH1F *) infile_Tracking_old->Get("Tracking_Efficiency");
  TH1F *eff_PHE_re_old=(TH1F *) infile_PHE_re_old->Get("ReConEffienc");
  TH1F *eff_purity_old=(TH1F *) infile_purity_old->Get("purity");

}
/*
  //old
 eff_Tof_match_old->SetMarkerColor(4);
 eff_Tof_match_old->SetMarkerStyle(20);
 eff_Tof_match_old->SetTitle("Tof_match");

 eff_dedx_old->SetMarkerColor(4);
 eff_dedx_old->SetMarkerStyle(20);
 eff_dedx_old->SetTitle("nsigma E"); 
 
 eff_Tracking_old->SetMarkerStyle(20);
 eff_Tracking_old->SetMarkerColor(4);
 eff_Tracking_old->SetTitle("Tracking");

 eff_PHE_re_old->SetMarkerStyle(20);
 eff_PHE_re_old->SetMarkerColor(4);
 eff_PHE_re_old->SetTitle("PHE reconstruction");

 eff_purity_old->SetMarkerStyle(20);
 eff_purity_old->SetMarkerColor(4);
 eff_purity_old->SetTitle("purity");

eff_Tof_cut_HT0_old->SetMarkerStyle(20);
eff_Tof_cut_HT0_old->SetMarkerColor(4);
//eff_Tof_cut_HT0_old->SetMarkerColor(4);

eff_Tof_cut_HT2_old->SetMarkerStyle(20);
eff_Tof_cut_HT2_old->SetMarkerColor(4);


//new
 eff_Tracking->SetMarkerColor(2);
 eff_Tracking->SetMarkerStyle(20);


 eff_PHE_re->SetMarkerStyle(20);
 eff_PHE_re->SetMarkerColor(2);

 eff_purity->SetMarkerStyle(20);
 eff_purity->SetMarkerColor(2);

eff_Tof_cut_HT0->SetMarkerStyle(20);
eff_Tof_cut_HT0->SetMarkerColor(2);

eff_Tof_cut_HT2->SetMarkerStyle(20);
eff_Tof_cut_HT2->SetMarkerColor(2);

 // eff_Tof_match_old->SetMarkerColor(4);
 // eff_Tof_match_old->SetMarkerStyle(20);


 TCanvas *c2=new TCanvas("c2","",1200,1000);
 c2->Divide(2,2);
 
 TCanvas *c3=new TCanvas("c3","",1200,1000);
 c3->Divide(2,2);
 
 
 TH2F *h2=new TH2F("h2","",100,2,10,100,0,1.2);
 h2->GetXaxis()->SetTitle("Pt");
 h2->GetYaxis()->SetTitle("Tof_match efficiency");
 //h2->SetTitle("Tof_match");
 // sprintf(buf,"Tof_match efficiency");
 //drawLatex(0.45,0.75,buf,70,0.035,2);
 c2->cd(1);
 h2->Draw(); 
 
 eff_Tof_match_old->Draw("sameP");
 eff_Tof_match->Draw("sameP");
 
 c2->cd(2);
 h2->SetTitle("nisgma E");
 h2->Draw(); 
 eff_dedx_old->Draw("same");
 eff_dedx->Draw("sameP");
 
 c2->cd(3);
 h2->Draw(); 
 h2->SetTitle("Tracking");
 eff_Tracking_old->Draw("same");
 eff_Tracking->Draw("sameP");
 
 c2->cd(4);
 h2->SetTitle("PHE recontruction");
 h2->Draw(); 

 eff_PHE_re_old->Draw("same");
 eff_PHE_re->Draw("sameP");

 c3->cd(1);
 h2->Draw();
 h2->SetTitle("purity");
 eff_purity_old->Draw("sameP");
 eff_purity->Draw("sameP");

 c3->cd(2);
 // gPad->SetLogy();
 eff_Tof_cut_HT0_old->Draw();
 eff_Tof_cut_HT0->Draw("samePE");

 c3->cd(3);
 // gPad->SetLogy();
 eff_Tof_cut_HT2_old->Draw();
 eff_Tof_cut_HT2->Draw("samePE");

 c2->SaveAs("efficiency_check_c2.pdf");
 c3->SaveAs("efficiency_check_c3.pdf");

 //TH1F *Inclusive[nTrg];
 char buf[1024];
 const int nTrg=2;
 TH1F *Inclusive_ps[nTrg];
 TH1F *Inclusive_ps_old[nTrg];
 
  TH2F *Photonic_ps_unlike[nTrg];
  TH2F *Photonic_ps_like[nTrg];
  TH2F *Photonic_ps_unlike_like[nTrg];


  TH1F *Photonic_unlike_pt_ps[nTrg];
  TH1F *Photonic_like_pt_ps[nTrg];
  TH1F *Photonic_unlike_like_pt_ps[nTrg];

 
  TH2F *Photonic_ps_unlike_old[nTrg];
  TH2F *Photonic_ps_like_old[nTrg];
  TH2F *Photonic_ps_unlike_like_old[nTrg];
  TH1F *Photonic_unlike_pt_ps_old[nTrg];
  TH1F *Photonic_like_pt_ps_old[nTrg];
  TH1F *Photonic_unlike_like_pt_ps_old[nTrg];


 TString Inclusive_name_ps[nTrg]={"mh1electronPt_psTrg0","mh1electronPt_psTrg1"};
  TString Photonic_unlikename_ps[nTrg]={"mh2InvMassUnlike_psTrg0","mh2InvMassUnlike_psTrg1"};
  TString Photonic_likename_ps[nTrg]={"mh2InvMasslike_psTrg0","mh2InvMasslike_psTrg1"};

  TString Photonic_unlikename_ps_old[nTrg]={"mh2InvMassUnlikeTrg0","mh2InvMassUnlikeTrg2"};
  TString Photonic_likename_ps_old[nTrg]={"mh2InvMasslikeTrg0","mh2InvMasslikeTrg2"};

for(Int_t iTrg=0;iTrg<nTrg;iTrg++)
    {
      Inclusive_ps[iTrg]=(TH1F *) infile_data->Get(Inclusive_name_ps[iTrg]);
      Inclusive_ps[iTrg]->SetLineColor(2);

      Photonic_ps_unlike[iTrg]=(TH2F *) infile_data->Get(Photonic_unlikename_ps[iTrg]);
      Photonic_ps_like[iTrg]=(TH2F *) infile_data->Get(Photonic_likename_ps[iTrg]);

      sprintf(buf,"unlike_ps%i",iTrg);
      Photonic_unlike_pt_ps[iTrg]=(TH1F *) Photonic_ps_unlike[iTrg]->ProjectionY(buf);
      sprintf(buf,"like_ps%i",iTrg);
      Photonic_like_pt_ps[iTrg]=(TH1F *) Photonic_ps_like[iTrg]->ProjectionY(buf);

      sprintf(buf,"unlike_like_ps%i",iTrg);
      Photonic_unlike_like_pt_ps[iTrg]=(TH1F *) Photonic_unlike_pt_ps[iTrg]->Clone(buf); 
      Photonic_unlike_like_pt_ps[iTrg]->Sumw2();
      Photonic_unlike_like_pt_ps[iTrg]->Add(Photonic_like_pt_ps[iTrg],-1);

      Photonic_unlike_like_pt_ps[iTrg]->SetLineColor(2);


      // old



      Photonic_ps_unlike_old[iTrg]=(TH2F *) infile_data_old->Get(Photonic_unlikename_ps_old[iTrg]);
      Photonic_ps_like_old[iTrg]=(TH2F *) infile_data_old->Get(Photonic_likename_ps_old[iTrg]);

      sprintf(buf,"unlike_ps_old%i",iTrg);
      Photonic_unlike_pt_ps_old[iTrg]=(TH1F *) Photonic_ps_unlike_old[iTrg]->ProjectionY(buf);
      sprintf(buf,"like_ps_old%i",iTrg);
      Photonic_like_pt_ps_old[iTrg]=(TH1F *) Photonic_ps_like_old[iTrg]->ProjectionY(buf);

      sprintf(buf,"unlike_like_ps_old%i",iTrg);
      Photonic_unlike_like_pt_ps_old[iTrg]=(TH1F *) Photonic_unlike_pt_ps_old[iTrg]->Clone(buf); 
      Photonic_unlike_like_pt_ps_old[iTrg]->Sumw2();
      Photonic_unlike_like_pt_ps_old[iTrg]->Add(Photonic_like_pt_ps_old[iTrg],-1);

           Photonic_unlike_like_pt_ps_old[iTrg]->SetLineColor(4);

    }

 Inclusive_ps_old[0]=(TH1F *) infile_data_old->Get("mh1electronPtTrg0");
 Inclusive_ps_old[1]=(TH1F *) infile_data_old->Get("mh1electronPtTrg2");

 Inclusive_ps_old[0]->SetLineColor(4);
 Inclusive_ps_old[1]->SetLineColor(4);

 cout<< Inclusive_ps_old[1]->Integral()<<endl;
 cout<< Inclusive_ps[1]->Integral()<<endl;


 TH1F *ratio_HT0_old=(TH1F *)Inclusive_ps_old[0]->Clone("ratio_HT0_old");
 TH1F *ratio_HT2_old=(TH1F *) Inclusive_ps_old[1]->Clone("ratio_HT2_old");

 TH1F *ratio_HT0=(TH1F *) Inclusive_ps[0]->Clone("ratio_HT0");
 TH1F *ratio_HT2=(TH1F *) Inclusive_ps[1]->Clone("ratio_HT2");

 ratio_HT0->Divide(Photonic_unlike_like_pt_ps[0]);
 ratio_HT2->Divide(Photonic_unlike_like_pt_ps[1]);

 ratio_HT0_old->Divide(Photonic_unlike_like_pt_ps_old[0]);
 ratio_HT2_old->Divide(Photonic_unlike_like_pt_ps_old[1]);




 TCanvas *c4=new TCanvas("c4","",1200,1000);
 c4->Divide(2,3);

 c4->cd(1);
 gPad->SetLogy();
 Inclusive_ps[0]->Draw();
 Inclusive_ps_old[0]->Draw("same");

 c4->cd(2);
 gPad->SetLogy();
 Inclusive_ps[1]->Draw();
 Inclusive_ps_old[1]->Draw("same");

 c4->cd(3);
 gPad->SetLogy();
 Photonic_unlike_like_pt_ps[0]->Draw();
 Photonic_unlike_like_pt_ps_old[0]->Draw("same");

 // Inclusive_ps_old[0]->Draw("same");
 c4->cd(4);
 gPad->SetLogy();
 Photonic_unlike_like_pt_ps[1]->Draw();
 Photonic_unlike_like_pt_ps_old[1]->Draw("same");


 c4->cd(5);

 ratio_HT0->GetYaxis()->SetRangeUser(0,5);
 ratio_HT0->Draw();
  ratio_HT0_old->Draw("same");

 c4->cd(6);
 ratio_HT2->GetYaxis()->SetRangeUser(0,5);
 ratio_HT2->Draw();
  ratio_HT2_old->Draw("same");

 c4->SaveAs("pt_spectra.pdf");

 return;


}

*/
