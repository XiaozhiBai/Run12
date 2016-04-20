TLatex* drawLatex(Double_t x,Double_t y,char* text,Int_t textFont,Double_t textSize,Int_t colorIndex){
        TLatex *latex =new TLatex(x,y,text);
        latex->SetNDC();
        latex->SetTextFont(textFont);
        latex->SetTextSize(textSize);
        latex->SetTextColor(colorIndex);
        latex->Draw("same");
        return latex;
}
TH1* histo(char *name, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, char* xTitle, char* yTitle){
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
  dd->GetXaxis()->SetNdivisions(512);
  return dd;
}
void drawnSigma()
{
    TF1 *Gaus = new TF1("Gaus","[0]*exp(-0.5*pow((x-[1])/[2],2))/sqrt(2.*TMath::Pi())/[2]",-3.5,3.5);
    Gaus->SetLineWidth(3);
    Gaus->SetParNames("counts","mean","Width");
    Gaus->SetParLimits(2,0.5, 1.5);
    Gaus->SetParLimits(1,-1., 1.5);
    Gaus->SetParLimits(0,0.00001, 1e6);

    const Int_t nPtBins = 21;
    const float xmin = -6;
    const float xmax = 6;
    float  ptbin[nPtBins+1] ;
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
    /*for(int i=0; i<14; i++)  ptbin[i] = 1.5 + 0.5*i;     //0.2 - 0.9
    ptbin[14] = 10;
    ptbin[15] = 14.0;
    ptbin[16] = 20.0;
    */
    double pt[nPtBins];
    double yield[nPtBins];
    double eyield[nPtBins];
    double ept[nPtBins];
    double ycount[nPtBins];
    double eycount[nPtBins];
    double mean[nPtBins];
    double meanerr[nPtBins];
    double sigma[nPtBins];
    double sigmaerr[nPtBins];

    TLatex tx;
    tx.SetTextSize(0.05);
    tx.SetNDC();

    char listname[100],listname1[100];
    //sprintf(listname,"../DataMaker/pheHT.root");
    sprintf(listname,"../../NPE_HT/RootFile/Root_File_2_26/hist_2_26.root");
    cout << listname << endl;
    TFile *f=new TFile(listname);
    TH2F *hnSigEvspt=(TH2F *)f->Get("");
    TH2F *hnSigEvsptlike=(TH2F *)f->Get("nSigEPtPartLike");

    TH1D *htmp = new TH1D("htmp","",1,-6.,6.);
    TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
    TPDF *mypdf = new TPDF("./result/nsigepartHT.pdf",111);
    mypdf->Off();
    ofstream outdata("./DataHT/nsigemeanAndsigmaeff.dat");
	//ofstream outdata("./DataHT/meanAndsigma.dat");
    //ofstream outdata("./Data/nsigecuteffCountsSysUnc.dat"); 

    for(Int_t i=0;i<21;i++){
    char name[100],name1[100],name2[100];
    sprintf(name1,"nSigEUnlike_%d",i);
    sprintf(name2,"nSigELike_%d",i);
    sprintf(name,"nSigE_%d",i);

    Int_t ibin1,ibin2;
    int ptlw=(int)hnSigEvspt->GetXaxis()->FindBin(ptbin[i]+1e-6);
    int ptup=(int)hnSigEvspt->GetXaxis()->FindBin(ptbin[i+1]-1e-6);
    
    TH1F *hnSigEUnlike=(TH1F *)hnSigEvspt->ProjectionY(name1,ptlw,ptup);
    TH1F *hnSigELike=(TH1F *)hnSigEvsptlike->ProjectionY(name2,ptlw,ptup);
    TH1F *hnSigE = (TH1F *)hnSigEUnlike->Clone();
    hnSigE->Sumw2();
    hnSigE->Add(hnSigELike,-1);
    hnSigE->SetName(name);

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

    float ymin = hnSigE->GetMinimum()*2.2;
    float ymax = hnSigEUnlike->GetMaximum()*4;
    htmp->GetYaxis()->SetTitle("Counts");
    htmp->GetYaxis()->SetTitleOffset(0.95);
	htmp->GetYaxis()->SetTitleSize(0.06);
    htmp->SetMinimum(ymin);
    htmp->SetMaximum(ymax);
    htmp->GetXaxis()->SetNdivisions(208);
    htmp->GetXaxis()->SetTitle("n#sigma_{e}");
    htmp->GetXaxis()->SetTitleOffset(0.95);
    htmp->GetXaxis()->SetTitleSize(0.07);
    htmp->GetXaxis()->SetLabelOffset(0.01);
    htmp->GetXaxis()->SetLabelSize(0.045);
    htmp->GetYaxis()->SetLabelSize(0.045);
    htmp->GetXaxis()->SetLabelFont(42);
    htmp->GetYaxis()->SetNdivisions(505);
    htmp->Draw();

    ibin1=(int)hnSigE->FindBin(-3.5+1e-6);
    ibin2=(int)hnSigE->FindBin(3.0-1e-6);
    double nMc1=hnSigE->Integral(ibin1,ibin2);

    ibin1=(int)hnSigE->FindBin(-1.0+1e-6);
    ibin2=(int)hnSigE->FindBin(3.0-1e-6);
    double nData1=hnSigE->Integral(ibin1,ibin2);

    hnSigEUnlike->Rebin(2);
    hnSigELike->Rebin(2);
    hnSigE->Rebin(2);

    hnSigELike->SetAxisRange(-3.5,3.0,"X");
    hnSigELike->SetMarkerStyle(8);
    hnSigELike->SetMarkerSize(0.9);
    hnSigELike->SetMarkerColor(4);
    hnSigELike->Draw("psame");
  
    hnSigEUnlike->SetAxisRange(-3.5,3.0,"X");
    hnSigEUnlike->SetMarkerStyle(8);
    hnSigEUnlike->SetMarkerSize(0.9);
    hnSigEUnlike->SetMarkerColor(kRed);
    hnSigEUnlike->Draw("psame");

    hnSigE->SetMarkerStyle(8);
    hnSigE->SetMarkerSize(1.3);
    hnSigE->SetMarkerColor(kGreen);
    hnSigE->Draw("psame");
    hnSigE->Fit(Gaus,"INOR","",-3.5,3.0);

    double *par= Gaus->GetParameters();
    Gaus->SetParameter(0,par[0]);
    Gaus->SetParameter(1,par[1]);
    Gaus->SetParameter(2,par[2]);

    hnSigE->Fit(Gaus,"INOR","",-3.5,3.0);
    Gaus->Draw("same");

    char txlb[100];
    sprintf(txlb,"%4.1f < p_{T} < %4.1f GeV/c",ptbin[i],ptbin[i+1]);
    tx.DrawLatex(0.6,0.87,txlb);

    TLegend *lg = new TLegend(0.565,0.67,0.9,0.82);
    lg->SetFillStyle(0);
    lg->SetFillColor(10);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.05);
    lg->AddEntry(hnSigEUnlike,"Unlike Sign","p");
    lg->AddEntry(hnSigELike,"Like Sign","p");
    lg->AddEntry(hnSigE,"Unlike - Like Sign","p");
    lg->Draw();

     float dsige = (20.+20.)/hnSigE->GetNbinsX();
     pt[i]=(ptbin[i]+ptbin[i+1])/2.0;
     ept[i]=0.5*(ptbin[i+1]-ptbin[i]);

     double nMc=Gaus->Integral(-3.5,3.0)/dsige;
     double nData=Gaus->Integral(-1.0,3.0)/dsige;

     yield[i]=(nData+1)/(nMc+2);
     ycount[i]=(nData1+1)/(nMc1+2);
     eyield[i]=sqrt((nData+1)*(nData+2)/(nMc+2)/(nMc+3)-pow(nData+1,2)/pow(nMc+2,2));
     eycount[i]=sqrt((nData1+1)*(nData1+2)/(nMc1+2)/(nMc1+3)-pow(nData1+1,2)/pow(nMc1+2,2));
     cout<<yield[i]<<" "<<ycount[i]<<" "<<eyield[i]<<" "<<eycount[i]<<endl;

     mean[i]=Gaus->GetParameter(1);
     meanerr[i]=Gaus->GetParError(1);
     sigma[i]=Gaus->GetParameter(2);
     sigmaerr[i]=Gaus->GetParError(2);

     TVirtualFitter * fitter = TVirtualFitter::GetFitter();
     assert(fitter != 0);
     double * cov =fitter->GetCovarianceMatrix();

     cout<<cov[0]<<" "<<sqrt(cov[4])<<" "<<sqrt(cov[8])<<" "<<cov[7]<<" meanerr="<<Gaus->GetParError(1)<<"sigmaerr "<<Gaus->GetParError(2)<<endl;
     outdata << pt[i]       << "  " << ept[i]    << "  "
             << mean[i]     << "  " << cov[4]    << "  " 
             << sigma[i]    << "  " << cov[8]    << "  "<<cov[7]<< endl; 
     //outdata << pt[i]    << "  " << 0.5*(ptbin[i+1]-ptbin[i]) << "  "
     //        << mean[i]  << "  " << meanerr[i]                << "  "     
     //        << sigma[i] << "  " << sigmaerr[i]               << endl;
     //outdata << pt[i]    << "  " << 0.5*(ptbin[i+1]-ptbin[i]) << "  "
     //                    <<ycount[i]<<" "<<eycount[i] << endl;
    
    char chh1[50],chh2[50],chh3[50],chh4[50],chh5[50],chh6[50],chh7[50];

	sprintf(chh2,"Entries = %5.0f", hnSigE->Integral());
	sprintf(chh3,"mean = %5.3f #pm %5.3f",hnSigE->GetMean(),hnSigE->GetMeanError());
	sprintf(chh4,"RMS = %5.3f #pm %5.3f",hnSigE->GetRMS(),hnSigE->GetRMSError());
	sprintf(chh1,"#chi^{2} / ndf = %6.2f / %d",Gaus->GetChisquare(),Gaus->GetNDF());
    sprintf(chh5,"N = %5.3e #pm %4.2e",Gaus->GetParameter(0),Gaus->GetParError(0));
    sprintf(chh6,"#mu = %5.3f #pm %5.3f",Gaus->GetParameter(1),Gaus->GetParError(1));
	sprintf(chh7,"#sigma = %5.3f #pm %5.3f",Gaus->GetParameter(2),Gaus->GetParError(2));

     TPaveStats *ptstats = new TPaveStats(0.15,0.5,0.43,0.92,"brNDC");
     ptstats->SetName("stats");
     ptstats->SetBorderSize(2);
     ptstats->SetFillColor(10);
     ptstats->SetTextAlign(12);
     ptstats->SetTextSize(0.03);
     TText *text = ptstats->AddText(chh2);
	 text = ptstats->AddText(chh3);
	 text = ptstats->AddText(chh4);
	 text = ptstats->AddText(chh1);
     text = ptstats->AddText(chh5);
     text = ptstats->AddText(chh6);
     text = ptstats->AddText(chh7);
     ptstats->Draw();

     char gifname[100];
     sprintf(gifname,"picHT/nSigE_pt_%3.2f_%3.2f.gif",ptbin[i],ptbin[i+1]);
     mypdf->On();
     cc->Update();
     mypdf->NewPage();
     cc->Modified();
     mypdf->Off();
     cc->SaveAs(gifname);
    }
     mypdf->On();
     mypdf->Close();
}
