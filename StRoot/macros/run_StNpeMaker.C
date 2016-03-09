void run_StNpeMaker(TString filelist,TString outFileName)
{
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StDmesonMaker");
  gSystem->Load("StNpeMaker");

  StNpeMaker* NpeMaker = new StNpeMaker(outFileName.Data());
  NpeMaker->bookObjects();

  char filename[1000];
  ifstream fstream(filelist.Data());
  int ifile = 0;
  while (fstream >> filename)
    {
      ++ifile;
      //      cout <<"sngl_file: "<<ifile<<" : "<<filename<<endl;
      NpeMaker->read(filename);
    }
  NpeMaker->writeObjects();
}
