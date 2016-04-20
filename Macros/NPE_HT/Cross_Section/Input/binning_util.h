#include "TH1D.h"
#include "TF1.h"

TH1D* get_bin_shift(unsigned long int nSample,TF1* f,unsigned int nbins,double* xbins)
{
	// nSample is the number of times to sample the function 
	// f is the assumption function
	// nbins is the size of your array
	// xbins is the array of bins you desire
	
	TH1D* hist = new TH1D(Form("hCorr_%s",f->GetName()),"",nbins-1,xbins);
	hist->Sumw2();

	// sampe f and fill the histgram
	while(hist->GetEntries()<nSample) hist->Fill(f->GetRandom());


	// divide by bin width
	for(int iBin=1;iBin<=hist->GetNbinsX();iBin++)
	{
		float w = hist->GetBinWidth(iBin);
		hist->SetBinContent(iBin,hist->GetBinContent(iBin)/w);
		hist->SetBinError(iBin,hist->GetBinError(iBin)/w);
	}

	// scale by total number of counts
	double fIntg = f->Integral(xbins[0],xbins[nbins-1]);
	hist->Scale(fIntg/(double)hist->GetEntries());

	// make correction
	for(int iBin=1;iBin<=hist->GetNbinsX();iBin++)
	{
		float pt = hist->GetBinCenter(iBin);
		float fc = f->Eval(pt);
		hist->SetBinContent(iBin,hist->GetBinContent(iBin)/fc);
		hist->SetBinError(iBin,hist->GetBinError(iBin)/fc);
	}

	return hist;
}
