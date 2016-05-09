#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <fstream>
#include "TFile.h"
#include "TLatex.h"
#include "../mBinning_MB.h"
TH1D* histo(char *name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, char* xTitle, char* yTitle){
  TH1D *dd = new TH1D(name,"",100,xlow,xup);
  dd->SetMinimum(ylow);
  dd->SetMaximum(yup);
  dd->GetXaxis()->SetTitle(xTitle);
  dd->GetYaxis()->SetTitle(yTitle);

  dd->GetXaxis()->SetTitleSize(0.055);
  dd->GetXaxis()->SetTitleOffset(0.9);
  dd->GetXaxis()->SetLabelSize(0.045);
  dd->GetYaxis()->SetTitleSize(0.055);
  dd->GetYaxis()->SetTitleOffset(1);
  dd->GetYaxis()->SetLabelSize(0.045);
  //dd->GetXaxis()->CenterTitle(kTRUE);
  //dd->GetYaxis()->CenterTitle(kTRUE);
  dd->GetXaxis()->SetNdivisions(512);
  return dd;
}

TLatex* drawLatex(Double_t x,Double_t y,char* text,Int_t textFont,Double_t textSize,Int_t colorIndex){
        TLatex *latex =new TLatex(x,y,text);
        latex->SetNDC();
        latex->SetTextFont(textFont);
        latex->SetTextSize(textSize);
        latex->SetTextColor(colorIndex);
        latex->Draw("same");
        return latex;
}
void Tof_Beta_Cut_efficiency()
{

  //gStyle->SetOptStat(00000);
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleOffset(1.3,"X");
  gStyle->SetTitleOffset(1.3,"Y");
  gStyle->SetPadTopMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14); 
  
  char buf[1024];
  
  double pt[NpT_bins_run12_MB];
  double ptrange[NpT_bins_run12_MB];
  double mean[NpT_bins_run12_MB],meant;
  double meanerr[NpT_bins_run12_MB],meanerrt;
  double sigma[NpT_bins_run12_MB],sigmat;
  double sigmaerr[NpT_bins_run12_MB],sigmaerrt;
  double cov[NpT_bins_run12_MB];
  double cor[NpT_bins_run12_MB],cort;
  double eff[NpT_bins_run12_MB];
  double efferr[NpT_bins_run12_MB];
  
  TGraphErrors *gr1;
  
  TLatex tx;
  tx.SetTextSize(0.05);
  tx.SetNDC();
  
  sprintf(buf,"Tof_beta_mean_And_sigma.dat");
  ifstream indata(buf);
  //cout<<buf<<endl;
  if(!indata)    {
    cout<<buf<<"doesn't exist!"<<endl;
    return;
  }
  
  for(int ipt=0;ipt<NpT_bins_run12_MB;ipt++) {
    indata>>pt[ipt]>>ptrange[ipt]>>mean[ipt]>>meanerr[ipt]>>sigma[ipt]>>sigmaerr[ipt]>>cov[ipt];
    meanerr[ipt]=sqrt(meanerr[ipt]);
    sigmaerr[ipt]=sqrt(sigmaerr[ipt]);
    cor[ipt]=cov[ipt]/meanerr[ipt]/sigmaerr[ipt];
  }
    indata.close();
    
    TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
   
    ofstream outdata("nsigecuteffvspt.dat");
    TH1F *heff;
    TH1F *Tof_beta_cut_MB=new TH1F("Tof_beta_cut_MB","",NpT_bins_run12_MB,pt_run12_MB);
 
    for(int ipt=0;ipt<NpT_bins_run12_MB;ipt++) {
      
	char name[100],name1[10000],name2[100];
	sprintf(name,"twogaus%d",ipt);
	sprintf(name2,"eff%d",ipt);
	heff = new TH1F(name2,"",1000,0.8,1.1); 
	
	meant=mean[ipt];
	meanerrt=meanerr[ipt];
	sigmat=sigma[ipt];
	sigmaerrt=sigmaerr[ipt];
	cort=cor[ipt];
	
	TF2 *twogaus = new TF2(name,"1./2./3.14/[0]/[1]/sqrt(1-[2]*[2])*exp(-1./2./(1-[2]*[2])*(pow((x[0]-[3])/[0],2)-2*[2]*(x[0]-[3])*(x[1]-[4])/[0]/[1]+pow((x[1]-[4])/[1],2)))",meant-3*meanerrt,meant+3*meanerrt,sigmat-3*sigmaerrt,sigmat+3*sigmaerrt);
	twogaus->SetParameter(0,meanerrt);
	twogaus->SetParameter(1,sigmaerrt);
	twogaus->SetParameter(2,cort);
	twogaus->SetParameter(3,meant);
	twogaus->SetParameter(4,sigmat);

	twogaus->GetXaxis()->SetTitle("mean");
	twogaus->GetYaxis()->SetTitle("sigma");

	TF1 *gaus=new TF1("gaus","exp(-0.5*pow((x-[0])/[1],2))/sqrt(2.*TMath::Pi())/[1]",-4,4);
	for(int j=0;j<40000;j++){
	  if(j%5000==0) cout << "begin " << j << "th entry...." << endl;
	  double mean1,sigma1;
	  twogaus->GetRandom2(mean1,sigma1);
	  gaus->SetParameter(0,mean1);
	  gaus->SetParameter(1,sigma1);
	  heff->Fill(gaus->Integral(-0.03,0.03)/gaus->Integral(-1,1));
	}
	
	TCanvas *c2 = new TCanvas("c2", "",0,0,1200,600);
	gStyle->SetOptFit(1111);
	c2->Divide(2,1);
	c2->cd(1);
	heff->Sumw2();
	heff->SetMarkerStyle(8);
	heff->SetMarkerSize(1);
	TF1 *Gaus = new TF1("Gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))/sqrt(2.*TMath::Pi())/[2]",0.9,1.1);
	Gaus->SetParNames("counts","mean","sigma");
    	// // Gaus->SetParLimits(1,heff->GetMean()-1*heff->GetRMS(),heff->GetMean()+1*heff->GetRMS());
	// // Gaus->SetParLimits(2,heff->GetRMS()-100*heff->GetRMSError(),heff->GetRMS()+100*heff->GetRMSError());

    	 Gaus->SetParameter(1,heff->GetMean());
	 Gaus->SetParameter(2,heff->GetRMS());
	 
	 heff->GetXaxis()->SetTitle("#beta cut efficiency");
	 heff->GetYaxis()->SetTitle("Counts");
	 heff->SetTitle(mh1_pT_Title[ipt]);

	 heff->Fit(Gaus,"R","",heff->GetMean()-5*heff->GetRMS(),heff->GetMean()+5*heff->GetRMS());
	
	 Tof_beta_cut_MB->SetBinContent(ipt+1,Gaus->GetParameter(1));
	 Tof_beta_cut_MB->SetBinError(ipt+1,Gaus->GetParameter(2));

	heff->Draw("PE");
	
	c2->cd(2);
	
	twogaus->Draw("SURF2");
	
	char name5[100];
	sprintf(name5,"htmp%d.pdf",ipt);    
	
	c2->SaveAs(name5);
	
	
    }

    Tof_beta_cut_MB->SetMarkerStyle(20);
    Tof_beta_cut_MB->SetMarkerSize(1);
    Tof_beta_cut_MB->SetMarkerColor(2);
    
    TH2F *h2=new TH2F("h2","",100,1,4,100,0,1.1);
    h2->GetXaxis()->SetTitle("p_{T} GeV/c");
    h2->GetYaxis()->SetTitle("beta cut efficiency"); 
    
    gStyle->SetOptStat(0);     
    TCanvas *c3= new TCanvas("c3","",800,600);
    Tof_beta_cut_MB->Draw("P");
    h2->Draw();
    Tof_beta_cut_MB->Draw("same");
    c3->SaveAs("Tof_beta_cut.pdf");
    
    gStyle->SetOptStat(1111);     
    TFile *file=new TFile("Tof_Beta_cuts_efficiency.root","RECREATE");
    Tof_beta_cut_MB->Write();
    file->Close();
    

    // float x1,x2;
    // x1=heff->GetMean()-4*heff->GetRMS();
    // x2=heff->GetMean()+4*heff->GetRMS();

    // TF1 *Gaus = new TF1("Gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))/sqrt(2.*TMath::Pi())/[2]",x1,x2);
    // Gaus->SetLineWidth(3);
    // Gaus->SetParNames("counts","mean","Width");
    // 	Gaus->SetParLimits(1,heff->GetMean()-1*heff->GetRMS(),heff->GetMean()+1*heff->GetRMS());
    // Gaus->SetParLimits(2,heff->GetRMS()-100*heff->GetRMSError(),heff->GetRMS()+100*heff->GetRMSError());
    // Gaus->SetParLimits(0,0.00001, 3e2);
    // 	heff->Fit(Gaus,"INOR","",x1,x2);

    // double *par= Gaus->GetParameters();
    // Gaus->SetParameter(0,par[0]);
    // Gaus->SetParameter(1,par[1]);
    // Gaus->SetParameter(2,par[2]);

    // heff->Fit(Gaus,"INOR","",x1,x2);
    // Gaus->Draw("same");
    // Gaus->SetLineColor(kGreen);
    // Gaus->SetLineWidth(3);
    // twogaus->Draw("SURF1");

    // char txlb[100];
    // sprintf(txlb,"%4.1f < p_{T} < %4.1f GeV/c",pt_run12_MB[ipt],pt_run12_MB[i+1]);
    // tx.DrawLatex(0.6,0.87,txlb);
    
    // char chh1[50],chh5[50],chh6[50],chh7[50];
    // sprintf(chh1,"#chi^{2} / ndf = %6.2f / %d",Gaus->GetChisquare(),Gaus->GetNDF());
    // sprintf(chh5,"N = %5.3e #pm %4.2e",Gaus->GetParameter(0),Gaus->GetParError(0));
    // sprintf(chh6,"#mu = %5.4f #pm %5.4f",Gaus->GetParameter(1),Gaus->GetParError(1));
    // sprintf(chh7,"#sigma = %5.4f #pm %5.4f",Gaus->GetParameter(2),Gaus->GetParError(2));

    // TPaveStats *ptstats = new TPaveStats(0.15,0.62,0.45,0.92,"brNDC");
    // ptstats->SetName("stats");
    // ptstats->SetBorderSize(2);
    // ptstats->SetFillColor(10);
    // ptstats->SetTextAlign(12);
    // ptstats->SetTextSize(0.03);
    // TText *text = ptstats->AddText(chh1);
    // text = ptstats->AddText(chh5);
    // text = ptstats->AddText(chh6);
    // text = ptstats->AddText(chh7);
    // //ptstats->Draw();
    
    // char gifname[100];
    // //sprintf(gifname,"picHT/nsigeeff_pt_%3.2f_%3.2f.gif",pt_run12_MB[ipt],pt_run12_MB[i+1]);
    // sprintf(gifname,"picHT/nsigmatwogaus_pt_%3.2f_%3.2f.gif",pt_run12_MB[ipt],pt_run12_MB[i+1]);

    // mypdf->On();
    // cc->Update();
    // mypdf->NewPage();
    // cc->Modified();
    // mypdf->Off();
    // cc->SaveAs(gifname);
    
    // eff[ipt]=Gaus->GetParameter(1);
    // efferr[ipt]=Gaus->GetParameter(2);
    // outdata << pt[ipt] << "  " << ptrange[ipt] << "  "<< eff[ipt]<<" "<<efferr[ipt] << endl;
    
    // 	}
    // mypdf->On();
    // mypdf->Close();     
}
