#ifndef CUTS_H
#define CUTS_H

// pay attention that all the pairs invariant mass cuts are applied in the code
namespace cuts
{
	//event
	const Float_t vz = 100.0;
	const Float_t vzVpdVz = 6.0;

	//tracking
	const Int_t nHitsFit = 15;
	const Float_t nHitsFitnHitsMax = 0.52;
	const Float_t dca = 2.0; // this is not applied for global tracks, of course :-)

	// hadrons
	const Float_t hadronsGlobalPt = 0.2;
	const Float_t hadronsDca = 3.0;

	//pions 
	const Float_t nSigmaPion = 3.0;
	const Float_t pionPt = 0.15; // this is used for soft pions, for other pions 0.2 is used, see code
	const Float_t pionEta = 1.0; // this is used for all but soft pions

	//kaons 
	const Float_t kaonPt = 0.2;
	const Float_t kaonEta = 1.0;
	const Float_t nSigmaKaon = 3.0;

	// electrons
	const Float_t electronPt = 0.2;
	const Float_t electronEta = 1.0;
	const Float_t nSigmaElectronHigh = 5.0;
	const Float_t nSigmaElectronLow = -5.0;
	const Float_t electronsBeta = 0.035; // for tracks with pT< 1.0

	// partner electron cuts
	const Float_t partnerElecNSigmaElectron = 5.0;
	const Float_t partnerElecGlobalPt = 0.2;
}
#endif
