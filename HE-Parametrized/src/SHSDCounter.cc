#include "SHSDCounter.hh"
#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"

#include <iostream>

using namespace std;
using namespace CLHEP;

const Double_t SCINTYIELD_SCSN81 = 10000;

//	Constructor
//
SHSDCounter::SHSDCounter(G4String name, RunParams inData, TTree *tree)
	: G4VSensitiveDetector(name),
	runParams(inData)
{
	shTree = tree;
	shTree->Branch("numPhotons", &numPhotons, "numPhotons/D");
	rand.SetSeed(runParams.seed);
	
	
}

//	Destructor
//
SHSDCounter::~SHSDCounter()
{

}

//	Initialize
//
void SHSDCounter::Initialize(G4HCofThisEvent*)
{
	//	Clear
	//
	numPhotons = 0;
	eneDep = 0;
//	waveSpectrum.clear();
}

//	Process Hits
//
G4bool SHSDCounter::ProcessHits(G4Step *aStep, G4TouchableHistory*)
{
	//	Process
	//
	
	//cout << "Inside the Counter..." << endl;
	//cout << "PName: " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << endl;
	if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == 
			"opticalphoton")
	{
/*		if (runParams.verbosity>1)
		{
			G4cout << "### SD: Optical Photon..." << G4endl;
			G4cout << "### Creator Process Name: " 
				<<	aStep->GetTrack()->GetCreatorProcess()->GetProcessName() 
				<< G4endl;
		}

		G4double totE = aStep->GetPreStepPoint()->GetTotalEnergy();
		G4double totE1 = aStep->GetTrack()->GetTotalEnergy();
//		cout << totE/eV << "  " << totE1/eV << endl;
		Double_t lyambda = 1239.8/(totE/eV);

		numPhotons++;
//		waveSpectrum.push_back(lyambda);
*/
//		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		double d =0;
	}

	if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() !=
			"opticalphoton")
	{
		double thisDep = aStep->GetTotalEnergyDeposit()/MeV;
		eneDep += thisDep;
		int gen = ComputeOP(thisDep, SCINTYIELD_SCSN81);
		numPhotons += gen;
	}

	return true;
}

//	Get the PMT ID by name
//
int SHSDCounter::GetPMT(G4String name)
{
	int pos_ = name.find("_");
	string strID = name.substr(pos_+1, name.length());
	return atoi(strID.c_str());
}

//	Get the Number of Photoelectrons
//
double SHSDCounter::GetPENum(G4double lyambda)
{
/*	//	For wavelength ourside of our range...
	//
	if ((lyambda > QPE[numBins-1][0]*nm) or (lyambda < QPE[0][0]*nm))
		return 0;

	//	If lyambda is equal to lower edge of our spectrum
	//
	if (lyambda == QPE[0][0]*nm)
		return QPE[0][1];

	//	If lyamda is equal to the upper edge of our spectrum
	//
	if (lyambda == QPE[numBins-1][0]*nm)
		return QPE[numBins - 1][1];

	//	Proceed... Find the right bin for interpolation
	//
	int l2pos = 0;
	while (lyambda > QPE[l2pos][0]*nm)
		l2pos++;
	double pes = ComputePE(l2pos-1, l2pos, lyambda/nm);

	return pes;
	*/
	return 0;
}

//	Interpolate/Compute QPE
//
double SHSDCounter::ComputePE(int b1, int b2, double x)
{
/*	double x1 = QPE[b1][0];
	double x2 = QPE[b2][0];
	double y1 = QPE[b1][1];
	double y2 = QPE[b2][1];
	double y = y1 + (y2 - y1)*(x - x1)/(x2 - x1);
	return y;
	*/
}

//	Finalize the event
//
void SHSDCounter::EndOfEvent(G4HCofThisEvent*)
{
	cout << "### photons=" << numPhotons << endl;
	shTree->Fill();

	//	Print the pe array for the event
	//
/*	if (runParams.verbosity)
	{
		for (int i=0; i < numPMTs; i++)
			cout << "### PMT#" << i+1 << " PEs: " << pe[i] << endl;
	}
	*/
}

int SHSDCounter::ComputeOP(double ene, double opPerMeV)
{
	int genPhotons = 0;
	double mean = ene*opPerMeV;
	double sigma = sqrt(mean);
	genPhotons = floor(rand.Gaus(mean, sigma));

	return genPhotons;
}



















