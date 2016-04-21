#ifndef CCOUNTERSD_H
#define CCOUNTERSD_H 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "globals.hh"

#include <vector>
#include <string>

#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"

#include "SHDefs.hh"

using namespace std;

//	Define the SD Class
//
class SHSDCounter : public G4VSensitiveDetector
{
	public:
		SHSDCounter(G4String, RunParams, TTree*);
		virtual ~SHSDCounter();

		void Initialize(G4HCofThisEvent*);
		G4bool ProcessHits(G4Step*, G4TouchableHistory*);
		void EndOfEvent(G4HCofThisEvent*);
		
		//	Returns QE based on the lyambda 
		//	Interpolates the QE Array provided by Burak
		//
		double GetPENum(double lyambda);
		double ComputePE(int b1, int b2, double x);
		int ComputeOP(double ene, double opPerMeV);
		int GetPMT(G4String name);

		//	Members
		//
//		double pe[numPMTs];
		Double_t numPhotons;
		Double_t eneDep;
		TRandom rand;
		//vector<Double_t> waveSpectrum;

		TTree *shTree;
		RunParams runParams;
};

#endif
