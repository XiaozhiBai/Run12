void run_StPhe_EmbeddingQA(TString filelist,TString outFileName)
{   
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StEmbeddingAna");

	StPhe_EmbeddingQA* MC = new StPhe_EmbeddingQA(outFileName.Data());
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
