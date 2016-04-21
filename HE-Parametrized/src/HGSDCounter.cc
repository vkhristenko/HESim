#include "HGSDCounter.hh"
#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"

#include "HGSDCounter.hh"

#include <iostream>

using namespace std;
using namespace CLHEP;

//
//	Constructor
//
HGSDCounter::HGSDCounter(G4String  name, RunParams inData, TTree *tree)
	: G4VSensitiveDetector(name), 
	_runParams(inData)
{
	_tree = tree;
	_subDName = name;
//	_tree->Branch("dr", &_dr);
	_tree->Branch("response", &_response, "response/D");
}

//
//	Desctructor
//
HGSDCounter::~HGSDCounter()
{}

//
//	Initialize
//
void HGSDCounter::Initialize(G4HCofThisEvent*)
{
//	_dx.clear();
//	_dy.clear();
//	_dz.clear();
//	_dr.clear();
	_response=0;
}

//
//	Process Hits
//
G4bool HGSDCounter::ProcessHits(G4Step *aStep, G4TouchableHistory*)
{
	if (aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge()!=0)
	{
		G4ThreeVector dPos = aStep->GetDeltaPosition();
//		_dx.push_back(dPos.x()/um);
//		_dy.push_back(dPos.y()/um);
//		_dz.push_back(dPos.z()/um);
		_dr.push_back(dPos.mag()/um);

		_response += (dPos.mag()/um)*80;

		if (_runParams.verbosity>1)
		{
			G4String preName = 
				aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
			G4String postName =
				aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
	
			//	Get local/global position for pre Volume
			//
			G4ThreeVector preGlobalPos = 
				aStep->GetPreStepPoint()->GetPosition();
			G4TouchableHandle theTouchable = 
				aStep->GetPreStepPoint()->GetTouchableHandle();
			G4ThreeVector preLocalPos =
				theTouchable->GetHistory()->GetTopTransform().TransformPoint(
				preGlobalPos);


			G4ThreeVector postGlobalPos = 
				aStep->GetPostStepPoint()->GetPosition();
			theTouchable = 
				aStep->GetPostStepPoint()->GetTouchableHandle();
			G4ThreeVector postLocalPos =
				theTouchable->GetHistory()->GetTopTransform().TransformPoint(
				postGlobalPos);

			G4cout << "### A hit in Part: " << _subDName
				<< G4endl
				<< preName << "  " << postName << " " << dPos.mag()/um
				<< G4endl
				<< preGlobalPos/um << "  " << preLocalPos/um
				<< G4endl
				<< postGlobalPos/um << "  " << postLocalPos/um
				<< G4endl;
		}
	}

	return true;
}

//
//	End of Event
//
void HGSDCounter::EndOfEvent(G4HCofThisEvent*)
{
	_tree->Fill();
}
