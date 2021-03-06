#include "StPicoTrigger.h"
#include "StPicoConstants.h"
#include "StPicoDstMaker.h"

ClassImp(StPicoTrigger)

//----------------------------------------------------------------------------------
StPicoTrigger::StPicoTrigger()
{
  Clear();
}

//----------------------------------------------------------------------------------
StPicoTrigger::StPicoTrigger(int flag, int id, int adc)
{
  Clear();

  if(flag<0) mFlag = 0;
  if(id  <0) mId   = 0;
  if(adc <0) mAdc  = 0;

  mFlag = (flag>Pico::UCHARMAX)  ? Pico::UCHARMAX  : (UChar_t)flag;
  mId   = (id  >Pico::USHORTMAX) ? Pico::USHORTMAX : (UShort_t)id;
  mAdc  = (adc >Pico::USHORTMAX) ? Pico::USHORTMAX : (UShort_t)adc;
}

//----------------------------------------------------------------------------------
StPicoTrigger::~StPicoTrigger()
{ /* noop */ }

//----------------------------------------------------------------------------------
void StPicoTrigger::Clear(const Option_t* opt)
{
  mFlag = 0;
  mId = 0;
  mAdc = 0;
}
//----------------------------------------------------------------------------------
void StPicoTrigger::Print(const Char_t *option) const {
  LOG_INFO << " Flag = " << mFlag << " Id = " << mId << " Adc = " << mAdc << endm;
}
