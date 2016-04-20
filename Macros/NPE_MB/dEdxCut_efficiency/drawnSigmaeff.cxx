#include "TF2.h"
#include "TH1F.h"

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
void drawnSigmaeff()
{
    float xmin[22];//= 0.65;
    float xmax[22];//= 0.75;
    for(int i=0; i<2; i++) {xmin[i]=0.5; xmax[i]=0.7;}
    for(int i=2; i<3; i++) {xmin[i]=0.55; xmax[i]=0.7;}
    for(int i=3; i<6; i++) {xmin[i]=0.48; xmax[i]=0.7;}
    for(int i=6; i<11; i++) {xmin[i]=0.45; xmax[i]=0.7;}
    for(int i=11; i<20; i++) {xmin[i]=0.45; xmax[i]=0.7;}
    xmin[20]=0.0;xmax[20]=0.7;
	xmin[21]=0.0;xmax[21]=0.7;

    char buf[1024];
    const Int_t nPtBins = 21;
    float  ptbin[nPtBins+1] ;
    /*
	for(int i=0; i<14; i++)  ptbin[i] = 1.5 + 0.5*i;     //0.2 - 0.9
    ptbin[14] = 10;
    ptbin[15] = 14.0;
	ptbin[16] = 20.0;
	*/
	for(int i=0; i<11; i++)  ptbin[i] = 1.5 + 0.1*i;     //0.2 - 0.9
    for(int i=11; i<14; i++) ptbin[i] = 2.6 + 0.2*(i-11); // 1.0 - 1.5
	ptbin[14] = 3.5;
	ptbin[15] = 4.0;
	ptbin[16] = 5.0;
	ptbin[17] = 6.0;
	ptbin[18] = 7.0;
	ptbin[19] = 9.0;
	ptbin[20] = 13.0;
	ptbin[21] = 20.0;
    
    double pt[nPtBins];
    double ptrange[nPtBins];
    double mean[nPtBins],meant;
    double meanerr[nPtBins],meanerrt;
    double sigma[nPtBins],sigmat;
    double sigmaerr[nPtBins],sigmaerrt;
    double cov[nPtBins];
	double cor[nPtBins],cort;
    double eff[nPtBins];
    double efferr[nPtBins];
    TGraphErrors *gr1;

	TLatex tx;
    tx.SetTextSize(0.05);
    tx.SetNDC();

    sprintf(buf,"./DataHT/nsigemeanAndsigmaeff.dat");
    ifstream indata(buf);
    //cout<<buf<<endl;
    if(!indata)    {
      cout<<buf<<"doesn't exist!"<<endl;
      return;
    }

    for(int i=0;i<21;i++) {
    indata>>pt[i]>>ptrange[i]>>mean[i]>>meanerr[i]>>sigma[i]>>sigmaerr[i]>>cov[i];
    meanerr[i]=sqrt(meanerr[i]);
    sigmaerr[i]=sqrt(sigmaerr[i]);
    cor[i]=cov[i]/meanerr[i]/sigmaerr[i];
    //cout<<meanerr[i]<<" "<<sigmaerr[i]<<" "<<cov[i]<<endl;
    }
    indata.close();

    TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
    //TPDF *mypdf = new TPDF("./result/nsigeHTeff.pdf",111);
    TPDF *mypdf = new TPDF("./result/nsigeHTefftwogaus.pdf",111);
	mypdf->Off();
    ofstream outdata("./DataHT/nsigecuteffvspt.dat");

	TH1F *heff;

	for(int i=0;i<21;i++) {
    char name[100],name1[10000],name2[100];
    sprintf(name,"twogaus%d",i);
    sprintf(name2,"eff%d",i);
    
    if(i<3) heff = new TH1F(name2,"",1000,0.,1.); 
	else if(i<19) {delete heff; heff = new TH1F(name2,"",500,0.,1.);}
	else if(i<20) {delete heff; heff = new TH1F(name2,"",200,0.,1.);}
	else {delete heff; heff = new TH1F(name2,"",100,0.,1.);}
    
    meant=mean[i];
    meanerrt=meanerr[i];
    sigmat=sigma[i];
    sigmaerrt=sigmaerr[i];
    cort=cor[i];

    TF2 *twogaus = new TF2(name,"1./2./3.14/[0]/[1]/sqrt(1-[2]*[2])*exp(-1./2./(1-[2]*[2])*(pow((x[0]-[3])/[0],2)-2*[2]*(x[0]-[3])*(x[1]-[4])/[0]/[1]+pow((x[1]-[4])/[1],2)))",meant-3*meanerrt,meant+3*meanerrt,sigmat-3*sigmaerrt,sigmat+3*sigmaerrt);
    twogaus->SetParameter(0,meanerrt);
    twogaus->SetParameter(1,sigmaerrt);
    twogaus->SetParameter(2,cort);
    twogaus->SetParameter(3,meant);
    twogaus->SetParameter(4,sigmat);

	TF1 *gaus=new TF1("gaus","exp(-0.5*pow((x-[0])/[1],2))/sqrt(2.*TMath::Pi())/[1]",-4,4);
    for(int j=0;j<10000;j++){
    if(j%5000==0) cout << "begin " << j << "th entry...." << endl;
    double mean1,sigma1;
    twogaus->GetRandom2(mean1,sigma1);
    gaus->SetParameter(0,mean1);
    gaus->SetParameter(1,sigma1);
    heff->Fill(gaus->Integral(-1,3)/gaus->Integral(-3.5,3.5));
    }
   
    TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0.01);
    gStyle->SetPalette(1,0);
    cc->SetTickx();
    cc->SetTicky();
    cc->SetFillColor(10);
    cc->SetBorderMode(0);
    cc->SetBorderSize(2);
    cc->SetFrameFillColor(0);
    cc->SetFrameBorderMode(0);
    cc->SetLeftMargin(0.12);
    cc->SetBottomMargin(0.15);
    cc->SetTopMargin(0.06);
    cc->SetRightMargin(0.04);

    char name5[100];
    sprintf(name5,"htmp%d",i);    
    float xxmin=xmin[i];
    float xxmax=xmax[i];
    TH1D *htmp = new TH1D(name5,"",1,xxmin,xxmax);
    float ymin = 0;//hbeta->GetMinimum()*1.2;
    float ymax = heff->GetMaximum()*1.5;
    htmp->SetAxisRange(xxmin, xxmax, "X");
    htmp->GetYaxis()->SetTitle("Counts");
    htmp->GetYaxis()->SetTitleOffset(1.);
	htmp->GetYaxis()->SetTitleSize(0.06);
    htmp->SetMinimum(ymin);
    htmp->SetMaximum(ymax);
    htmp->GetXaxis()->SetNdivisions(208);
    htmp->GetXaxis()->SetTitle("n#sigma_{e} cut efficiency");
    htmp->GetXaxis()->SetTitleOffset(1);
    htmp->GetXaxis()->SetTitleSize(0.06);
    htmp->GetXaxis()->SetLabelOffset(0.01);
    htmp->GetXaxis()->SetLabelSize(0.045);
    htmp->GetYaxis()->SetLabelSize(0.045);
    htmp->GetXaxis()->SetLabelFont(42);
    htmp->GetYaxis()->SetNdivisions(505);
    htmp->Draw();

    heff->Sumw2();
    heff->SetMarkerStyle(8);
    heff->SetMarkerSize(1);
    heff->Draw("psame");

    float x1,x2;
    x1=heff->GetMean()-4*heff->GetRMS();
    x2=heff->GetMean()+4*heff->GetRMS();

    TF1 *Gaus = new TF1("Gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))/sqrt(2.*TMath::Pi())/[2]",x1,x2);
    Gaus->SetLineWidth(3);
    Gaus->SetParNames("counts","mean","Width");
	Gaus->SetParLimits(1,heff->GetMean()-1*heff->GetRMS(),heff->GetMean()+1*heff->GetRMS());
    Gaus->SetParLimits(2,heff->GetRMS()-100*heff->GetRMSError(),heff->GetRMS()+100*heff->GetRMSError());
    Gaus->SetParLimits(0,0.00001, 3e2);
	heff->Fit(Gaus,"INOR","",x1,x2);

    double *par= Gaus->GetParameters();
    Gaus->SetParameter(0,par[0]);
    Gaus->SetParameter(1,par[1]);
    Gaus->SetParameter(2,par[2]);

    heff->Fit(Gaus,"INOR","",x1,x2);
    Gaus->Draw("same");
    Gaus->SetLineColor(kGreen);
    Gaus->SetLineWidth(3);
    twogaus->Draw("SURF1");

    char txlb[100];
    sprintf(txlb,"%4.1f < p_{T} < %4.1f GeV/c",ptbin[i],ptbin[i+1]);
    tx.DrawLatex(0.6,0.87,txlb);
    
    char chh1[50],chh5[50],chh6[50],chh7[50];
    sprintf(chh1,"#chi^{2} / ndf = %6.2f / %d",Gaus->GetChisquare(),Gaus->GetNDF());
    sprintf(chh5,"N = %5.3e #pm %4.2e",Gaus->GetParameter(0),Gaus->GetParError(0));
    sprintf(chh6,"#mu = %5.4f #pm %5.4f",Gaus->GetParameter(1),Gaus->GetParError(1));
    sprintf(chh7,"#sigma = %5.4f #pm %5.4f",Gaus->GetParameter(2),Gaus->GetParError(2));

    TPaveStats *ptstats = new TPaveStats(0.15,0.62,0.45,0.92,"brNDC");
    ptstats->SetName("stats");
    ptstats->SetBorderSize(2);
    ptstats->SetFillColor(10);
    ptstats->SetTextAlign(12);
    ptstats->SetTextSize(0.03);
    TText *text = ptstats->AddText(chh1);
    text = ptstats->AddText(chh5);
    text = ptstats->AddText(chh6);
    text = ptstats->AddText(chh7);
    //ptstats->Draw();
    
    char gifname[100];
    //sprintf(gifname,"picHT/nsigeeff_pt_%3.2f_%3.2f.gif",ptbin[i],ptbin[i+1]);
    sprintf(gifname,"picHT/nsigmatwogaus_pt_%3.2f_%3.2f.gif",ptbin[i],ptbin[i+1]);
	mypdf->On();
    cc->Update();
    mypdf->NewPage();
    cc->Modified();
    mypdf->Off();
    cc->SaveAs(gifname);

    eff[i]=Gaus->GetParameter(1);
    efferr[i]=Gaus->GetParameter(2);
    outdata << pt[i] << "  " << ptrange[i] << "  "<< eff[i]<<" "<<efferr[i] << endl;

}
    mypdf->On();
    mypdf->Close();     
}
