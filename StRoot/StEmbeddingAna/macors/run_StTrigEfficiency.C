void run_StTrigEfficiency(TString filelist,TString outFileName)
{
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StEmbeddingAna");

	StTrackingEfficiency* MC = new StTrackingEfficiency(outFileName.Data());

	MC->bookHistogram();
	char filename[1000];
	ifstream fstream(filelist.Data());
	int ifile = 0;
	while (fstream >> filename)
	{
		++ifile;
		MC->read(filename);
	}
	MC->WriteHistogram();
}
