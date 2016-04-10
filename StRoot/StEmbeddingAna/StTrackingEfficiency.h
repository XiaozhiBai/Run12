# ifndef StTrackingEfficiency_h
# define StTrackingEfficiency_h

#include "StSingle_Electron_tracks.h"
#include "StSingle_electron_EventCount.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

class TString;
class TFile;
class TF1;

class TH1F;
class TH2F;

class StSingle_Electron_tracks;
class StSingle_electron_EventCount;
class TGraphErrors;
class StTrackingEfficiency
{
 public:
  StTrackingEfficiency(const char *outName);
  virtual ~StTrackingEfficiency();

  void bookHistogram() ;
  void read(TString filename);
  void WriteHistogram();
  bool isHotTower(StSingle_Electron_tracks *);
  Double_t Weight(Float_t ,Float_t);
  
 private:
  TFile * mOutputFile;
  StSingle_Electron_tracks *AcTracks;
  StSingle_electron_EventCount *mAcEvent;   
  TFile * weightFile;
  TGraphErrors * gFONLLc; 
  
 
  TH1F *  mNTrack_nocuts;
  TH1F *  mNTrack_cut;
  TH1F *  mNTrack_cut_25;
  TH1F *  mNTrack_cut_20;

 
  ClassDef(StTrackingEfficiency,1)

};

#endif
