#include "SHPrimaryGeneratorAction.hh"

#include <iostream>

#include "G4HEPEvtInterface.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "CLHEP/Random/RandFlat.h"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

using namespace CLHEP;
using namespace std;

//	0-7 are old energies(for EM)
//	8-12 are new energies(for HAD)
//
const int numEnergies = 13;
const int numPrims = 3;
G4double primEnergies[numEnergies] = {
	1.*GeV, 2.*GeV, 4.*GeV, 8.*GeV, 16.*GeV, 32.*GeV, 50.*GeV, 60.*GeV,
	20*GeV, 50*GeV, 100*GeV, 200*GeV, 300*GeV
};
G4String primNames[numPrims] = {"e-", "pi-", "mu-"};

//	Constructor
//
SHPrimaryGeneratorAction::SHPrimaryGeneratorAction(RunParams params, TTree *tree)
	: runParams(params),
	shTree(tree)
{
	//	Define a generator
	//
	primEnergy = primEnergies[runParams.iEnergy];
	primName = primNames[runParams.iPrim];
//	primName = "pi-";
//	primEnergy = 10*GeV;
//	primName = "proton";
	particleGun = new G4ParticleGun(1);
	particleGun->SetParticleEnergy(primEnergy);
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *particle = particleTable->FindParticle(primName);
	particleGun->SetParticleDefinition(particle);

	primPos = G4ThreeVector(0, 0, -2.*m);
	primDir = G4ThreeVector(0.*m, 0.*m, 1.0*m);
	particleGun->SetParticlePosition(primPos);
	particleGun->SetParticleMomentumDirection(primDir);
	particleGun->SetParticleEnergy(primEnergy);
}

//	Destructor
//
SHPrimaryGeneratorAction::~SHPrimaryGeneratorAction()
{
	delete particleGun;
}

//	Generate Primary Event
//
void SHPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
	//	All gun's settings have been set up in Constructor
	//
	//particleGun->SetParticleEnergy(10.*keV);
	G4cout << particleGun->GetParticleEnergy()/GeV << "  " 
		<< particleGun->GetParticleDefinition()->GetParticleName()
		<< G4endl; 
	particleGun->GeneratePrimaryVertex(anEvent);
}
