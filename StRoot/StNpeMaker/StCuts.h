/*
 *
Wed Dec  2 09:31:44 PST 2015
by xiaozhi Bai xiaozhi@uic.edu
*
*
*/


#ifndef CUTS_H
#define CUTS_H
namespace cuts
{
  const Int_t adc0_HT1=256; //  
  const Int_t adc0_HT2=304; //  
  const Int_t adc0_HT3=416; //  

  const Int_t DSMadc_HT1=15; //  
  const Int_t DSMadc_HT2=18; //  
  const Int_t DSMadc_HT3=25; //  

 //event
  const Float_t Tpc_vz = 35.;
  const Float_t diff_vzVpdVz = 6.0;
  const Float_t mRanking=0; // >0 
 
  // track cuts


  const Float_t eta_max =0.7;
  const Float_t pt_min=0.2;
  const Float_t Pt_cut=2.0;
  
  const Float_t gDca = 1.5;
  const Float_t first_point_cut =73;
  const Int_t nhitDedx=15;
  const Int_t nhitsFit=20;
  const Int_t nhit_Possmax=0.52;    

  // ePID cuts


  const Int_t nsigmae_low=-1;
  const Int_t nsigmae_high=3;

  const Float_t Dphi =0.015;
  const Float_t Dz=3;
  const Float_t poe_low=0.3;
  const Float_t poe_high=1.5;
  const Int_t NEta =1;
  const Int_t NPhi=1;

  /* const Float_t eta_max=0.7; */
  /* const Float_t pt_min=0.2; */
  /* const Int_t nhitsFit = 20; */
  /* const Float_t nHitsRatio = 0.52; */
  /* const Int_t nhitDedx=15; */
  /* const Float_t gDca = 1.5; */
  /* const Float_t tofbeta=0.03;  */
  /* const Float_t tof_Y=1.8; */

  /* const Float_t Pt_cut_EMC=1.5; */

  // ee pair cuts
  const Float_t InvMass=0.05;
  const Float_t pairDca=0.5;

  /* // photonic electron partner cuts; */
  /* const Float_t nsigmaE_partner_low=-1.5; */
  /* const Float_t nsigmaE_partner_high=3; */
  /* const Float_t  pt_min_partner=0.2; */
  /* const Float_t eta_partner=1; */
  /* const Float_t  nhits_partner=15; // anatree */
  /* const Float_t Ratio_partner_low=0.52; // anasTree */
  /* const Float_t Ratio_Partner_high=1.02; // anaTree */
  /* const Float_t gDca_parner=3; // anaTree */

}
#endif
