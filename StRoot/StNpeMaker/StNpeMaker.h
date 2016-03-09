/******************************************************************************************
 *Mon Mar  7 19:33:04 EST 2016                                                            *
 * This class for read  DemesonTree Fill all the No-phoronic electron realted histogram   *
 * By xiaozhi Bai xiaozhi@uic.edu                                                         *
 ******************************************************************************************/
 
// StNpeMaker.h 

#ifndef StNpeMaker_h
#define StNpeMaker_h

#include <exception>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include "TProfile.h"

class TString;
class TF1;
class TH1D;
class TH1F;
class TH2F;
class TH3F;
class TFile;
class TNtuple;
class StDmesonTrack;
class StDmesonEvent;
class prescales;
class StElectronPair;

class StNpeMaker 
{
 public:
  StNpeMaker(const char *outName);
  virtual ~StNpeMaker();
  virtual void  Clear(Option_t *opt="");
  void bookObjects();
  void read(TString fileName);
  void writeObjects();

  void Run_QA(StDmesonEvent * ); 

 private:
  TFile* mOutputFile;
  StDmesonEvent * mEvent;
  prescales *mPrescales;
  ifstream file_runNumber;
  ofstream  outfile;
  // Event 
  TH2F *HT0_HT2;
  ClassDef(StNpeMaker, 1)
};
#endif
