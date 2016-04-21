#include "G4NistManager.hh"
#include "G4UnitsTable.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"

#include "SHDetectorConstruction.hh"
#include "HGSDCounter.hh"


/*
 *	Volume's dimensions
 */
//	World
//
G4double fullWorldX = 5.*m;
G4double fullWorldY = 5.*m;
G4double fullWorldZ = 5.*m;

//	Shashlik
//
G4double fullShashlikX = 14.*mm;
G4double fullShashlikY = 14.*mm;
G4double fullShashlikZ = 114.*mm;

//	Shashlik's container
//
G4double fullShashlikContainerX;
G4double fullShashlikContainerY;
G4double fullShashlikContainerZ = fullShashlikZ;

//	Absorber
//
G4double fullAbsX = 14.*mm;
G4double fullAbsY = 14.*mm;
G4double fullAbsZ = 2.5*mm;

//	Active Material
//
G4double fullActX = 14.*mm;
G4double fullActY = 14.*mm;
G4double fullActZ = 1.5*mm;

//	1 Layer
//
G4double fullLayerX = 14.*mm;
G4double fullLayerY = 14.*mm;
G4double fullLayerZ = fullAbsZ + fullActZ;

//	Fiber
//
G4double inRFiber = 0.*mm;
G4double outRFiber = 0.5*mm;
G4double fullAbsFiberZ = fullAbsZ;
G4double fullActFiberZ = fullActZ;

//	Gap
//
G4double fullGapX = 0.75*mm;
G4double fullGapY = 0.75*mm;

/*
 *	Member Functions
 */
SHDetectorConstruction::~SHDetectorConstruction()
{

}

//	
//	New Constructor: As of V5. To separate EM/FH/BH parts output.
//
SHDetectorConstruction::SHDetectorConstruction(RunParams params,
		TTree *emTree, TTree *fhTree, TTree *bhTree) :
	_emTree(emTree),
	_fhTree(fhTree),
	_bhTree(bhTree)
{
	MyConstructor(params, emTree);
}

//
//	Constructor
//
SHDetectorConstruction::SHDetectorConstruction(RunParams params,
		TTree *tree)
{
	MyConstructor(params, tree);
}

//	MyConstructor
//
void SHDetectorConstruction::MyConstructor(RunParams params, 
		TTree *tree)
{
	runParams = params;
	shTree = tree;
	_heTree = shTree;

	//	Set the Shashlik's Container Dimensions
	//
	fullShashlikContainerX = runParams.numModules*fullShashlikX + 
		(runParams.numModules-1)*fullGapX;
	fullShashlikContainerY = runParams.numModules*fullShashlikY + 
		(runParams.numModules-1)*fullGapY;

	//	Let's build materials
	//
	G4double a; G4double z; G4double density;
	G4double weightRatio; G4String name; G4String symbol;
	G4int nElem; G4int natoms;

	//	Elements go first
	//
	G4Element *eC = new G4Element(name="Carbon", symbol="C", z=6., 
			a=12.01*g/mole);
	G4Element *eN = new G4Element(name="Nitrogen", symbol="N", z=7.,
			a=14.01*g/mole);
	eO = new G4Element(name="Oxygen", symbol="O", z=8,
			a=16.00*g/mole);
	G4Element *eFe = new G4Element(name="Iron", symbol="Fe", z=26.,
			a=55.845*g/mole);
//	eH = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
	eH = new G4Element("Hydrogen", symbol="H", z=1, a=1.00794*g/mole);

	eLu = new G4Element(name="Lutetium", symbol="Lu", z=71.,
			a=174.97*g/mole);
	eSi = new G4Element(name="Silicium", symbol="Si", z=14.,
			a=28.09*g/mole);
	eY = new G4Element(name="Yttrium", symbol="Y", z=39.,
			a=88.91*g/mole);

	density = 2.700*g/cm3;
	a = 26.98*g/mole;
	G4Element *eAl = new G4Element(name="Aluminum", symbol="Al", z=13.,
			a);
	G4Element *eBe = new G4Element(name="Beryllium", symbol="Be", z=4, 
			a=9.012*g/mole);
	G4Element *eMg = new G4Element(name="Magnesium", symbol="Mg", z=12., 
			a=24.305*g/mole);
	G4Element *eTi = new G4Element(name="Titanium", symbol="Ti",  z=22., 
			a=47.867*g/mole);
	G4Element *eCs = new G4Element(name="Cesium", symbol="Cs", z=55., 
			a=132.90546*g/mole);
	G4Element *eSb = new G4Element(name="Antimony", symbol="Sb", z=51., 
			a=121.76*g/mole);
	G4Element *eGa = new G4Element(name="Gallium", symbol="Ga", z=31., 
			a=69.723*g/mole);
	G4Element *eP = new G4Element(name="Phosphorus", symbol="P", z=15., 
			a=30.97376*g/mole);
	G4Element *eAs = new G4Element(name="Arsenic", symbol="As", z=33., 
			a=74.9216*g/mole);
	G4Element *eIn = new G4Element(name="Indium", symbol="In", z=49., 
			a=114.818*g/mole);

	// Build the final compositions -> TODO
	//
	mAl2O3 = new G4Material(name="AluminumOxide", density=3.95*g/cm3, 
			nElem=2);
	mAl2O3->AddElement(eAl, natoms=2);
	mAl2O3->AddElement(eO, natoms=3);

	//mBeO = new G4Material(name="BerylliumOxide", density=3.02*g/cm3,
	//		nElem=2);
	mBeO = G4NistManager::Instance()->FindOrBuildMaterial("G4_BERYLLIUM_OXIDE");
//	mBeO->AddElement(eBe, natoms=1);
//	mBeO->AddElement(eO, natoms=1);
	
	mMgO = new G4Material(name="MagnesiumOxide", density=3.58*g/cm3,
			nElem=2);
	mMgO->AddElement(eMg, natoms=1);
	mMgO->AddElement(eO, natoms=1);

	mTiO = new G4Material(name="TitaniumOxide", density=4.23*g/cm3,
			nElem=2);
	mTiO->AddElement(eTi, natoms=1);
	mTiO->AddElement(eO, natoms=1);

	mCs3Sb = new G4Material(name="AntimonyTriCesium", density=3.076*g/cm3,
			nElem=2);
	mCs3Sb->AddElement(eCs, natoms=3);
	mCs3Sb->AddElement(eSb, natoms=1);
	
	mGaP = new G4Material(name="GalliumPhosphide", density=4.138*g/cm3,
			nElem=2);
	mGaP->AddElement(eGa, natoms=1);
	mGaP->AddElement(eP, natoms=1);

	mGaAsP = new G4Material(name="GalliumArsenicPhosphide", 
			density=4.482*g/cm3, nElem=3);
	mGaAsP->AddElement(eGa, natoms=1);
	mGaAsP->AddElement(eAs, natoms=1);
	mGaAsP->AddElement(eP, natoms=1);

	mGaPCs = new G4Material(name="GalliumCesiumPhosphide", 
			density=3.2*g/cm3, nElem=3);
	mGaPCs->AddElement(eGa, natoms=1);
	mGaPCs->AddElement(eP, natoms=1);
	mGaPCs->AddElement(eCs, natoms=1);

	mGaInP = new G4Material(name="GalliumIndiumPhosphide", 
			density=5.012*g/cm3, nElem=3);
	mGaInP->AddElement(eGa, natoms=1);
	mGaInP->AddElement(eIn, natoms=1);
	mGaInP->AddElement(eP, natoms=1);

	//mVacuum = new G4Material("Vacuum", z=1., a=1.008*g/mole, density=1.3-25*g/cm3, kStateGas, 2.73*kelvin, 3.e-18*pascal);	
//	mVacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
	
	density = universe_mean_density;
	G4double pressure = 1.e-19*pascal;
	G4double temperature = 0.1*kelvin;
	mVacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
		kStateGas, temperature, pressure);

	//	Tungsten
	//
	mW = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

	//	Lead
	//
	mPb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

	//	Copper
	//
	mCu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");

	//	Silicon Material
	//
	mSi = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

	//	Zink
	//
	mZn = G4NistManager::Instance()->FindOrBuildMaterial("G4_Zn");

	//	Brass
	//
	mBrass = new G4Material("Brass", density=8.5*g/cm3, 2);
	mBrass->AddMaterial(mCu, 70*perCent);
	mBrass->AddMaterial(mZn, 30*perCent);
	
	//	SiO2 or Quartz
	//
	density = 2.648*g/cm3;
	mSiO2 = new G4Material(name="SiO2", density, nElem=2);
	mSiO2->AddElement(eSi, natoms=1);
	mSiO2->AddElement(eO, natoms=2);

	//	Clm from DHCAL
	//
	mClm = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cl");

	//	Glass from DHCAL
	//
	mGlass = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");

	//	Epoxy from DHCAL
	//
	mEpoxy = new G4Material("Epoxy", density=1.3*g/cm3, 3);
	mEpoxy->AddElement(eH, 44);
	mEpoxy->AddElement(eC, 15);
	mEpoxy->AddElement(eO, 7);

	//	G10 from DHCAL code
	//
	mG10 = new G4Material("G10", density=1.9*g/cm3, 3);
	mG10->AddMaterial(mClm, 0.08);
	mG10->AddMaterial(mGlass, 0.773);
	mG10->AddMaterial(mEpoxy, 0.147);

	//	LYSO
	//
	mLYSO = new G4Material(name="LYSO", density=7.1*g/cm3, nElem=4);
	mLYSO->AddElement(eLu, 0.31101534);
	mLYSO->AddElement(eY, 0.368765605);
	mLYSO->AddElement(eSi, 0.083209699);
	mLYSO->AddElement(eO, 0.237009356);

	//
	//	Now, we have to build the Material's Properties Table
	//
	
	//	LYSO MPT
	//
	G4MaterialPropertiesTable *mptLYSO = new G4MaterialPropertiesTable();

	//	More realistic Emission Spectrum
	//
	const G4int NENTRIES = 27;
	G4double PHOTONENERGY_LYSO[NENTRIES] = {2.066333*eV, 2.101356*eV,
		2.137586*eV, 2.175088*eV, 2.213929*eV, 2.254182*eV, 2.295926*eV,
		2.339245*eV, 2.384231*eV, 2.43098*eV, 2.4796*eV, 2.530204*eV,
		2.582917*eV, 2.637872*eV, 2.695217*eV, 2.755111*eV, 2.817727*eV,
		2.883256*eV, 2.917176*eV, 2.951905*eV, 3.023902*eV, 3.0995*eV,
		3.178974*eV, 3.262632*eV, 3.350811*eV, 3.443889*eV, 3.542286*eV};
	G4double RINDEX_LYSO[NENTRIES] = {
		1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82,
		1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82, 1.82,
		1.82, 1.82, 1.82, 1.82, 1.82};
	G4double SLOW_LYSO[NENTRIES] = {
		0.000313, 0.000938, 0.003125, 0.00625, 0.009375, 0.01875, 0.025, 0.03125,
		0.046875, 0.0625, 0.0875, 0.125, 0.1875, 0.21875, 0.3125, 0.515625,
		0.6875, 0.84375, 0.94375, 0.9375, 0.9375, 1, 0.75, 0.5625, 0.0625,
		0.00625, 0.000313};
	//	For Debugging...
	//
	for (int i=0; i<NENTRIES; i++)
		cout << i << "  " << PHOTONENERGY_LYSO[i]/eV 
			<< "  " << SLOW_LYSO[i] << endl;

	mptLYSO->AddProperty("RINDEX", PHOTONENERGY_LYSO, RINDEX_LYSO, NENTRIES);
	mptLYSO->AddProperty("SLOWCOMPONENT", PHOTONENERGY_LYSO, SLOW_LYSO, NENTRIES);
	mptLYSO->AddConstProperty("SCINTILLATIONYIELD", 32000/MeV);
	mptLYSO->AddConstProperty("RESOLUTIONSCALE", 1.0);
	mptLYSO->AddConstProperty("SLOWTIMECONSTANT", 41.*ns);

	//	POLYETHYLENE + its MPT
	//	Properties for MPT are for SCSN81
	//
	mPOLYETHYLENE = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4MaterialPropertiesTable *mptSCSN = new G4MaterialPropertiesTable();
	const G4int NENTRIES_SCSN = 15;
	G4double PHOTONENERGY_SCSN[NENTRIES_SCSN] = {
		2.25418*eV, 2.36152*eV, 2.4796*eV, 2.53020*eV, 2.69522*eV,
		2.75511*eV, 2.81773*eV, 2.88326*eV, 2.95190*eV, 3.02390*eV,
		3.06123*eV, 3.09950*eV, 3.17897*eV, 3.26263*eV, 3.54229*eV};
	G4double RINDEX_SCSN[NENTRIES_SCSN] = {
		1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57,
		1.57, 1.57, 1.57, 1.57, 1.57};
	G4double SLOW_SCSN[NENTRIES_SCSN] = {
		0.05, 0.1, 0.16, 0.2, 0.6, 0.67, 0.8, 1., 0.8, 0.72, 0.75, 
		0.6, 0.1, 0.05, 0.001};
	mptSCSN->AddProperty("RINDEX", PHOTONENERGY_SCSN, RINDEX_SCSN, 
			NENTRIES_SCSN);
	mptSCSN->AddProperty("SLOWCOMPONENT", PHOTONENERGY_SCSN, SLOW_SCSN,
			NENTRIES_SCSN);
	mptSCSN->AddConstProperty("SCINTILLATIONYIELD", 10000/MeV);
	mptSCSN->AddConstProperty("RESOLUTIONSCALE", 1.0);
	mptSCSN->AddConstProperty("SLOWTIMECONSTANT", 2.5*ns);

	mPOLYETHYLENE->SetMaterialPropertiesTable(mptSCSN);
	mHEScintMain = mPOLYETHYLENE;

/*	
	const G4int nEntries = 2;
	G4double photonEnergy[nEntries] = {1.5*eV, 6.2*eV};
	G4double refractiveIndex[nEntries] = {1.82, 1.82};
	G4double fast[nEntries] = {1, 1};
	G4double slow[nEntries] = {1, 1};
	mptLYSO->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);
	mptLYSO->AddProperty("FASTCOMPONENT", photonEnergy, fast, nEntries);
	mptLYSO->AddProperty("SLOWCOMPONENT", photonEnergy, slow, nEntries);

	mptLYSO->AddConstProperty("SCINTILLATIONYIELD", 32000/MeV);
	mptLYSO->AddConstProperty("RESOLUTIONSCALE", 1.0);
	mptLYSO->AddConstProperty("FASTTIMECONSTANT", 41.*ns);
	mptLYSO->AddConstProperty("SLOWTIMECONSTANT", 41.*ns);
	mptLYSO->AddConstProperty("YIELDRATIO", 0.5);
*/

	mLYSO->SetMaterialPropertiesTable(mptLYSO);

	//	Choose Materials for EM Absorber
	//	1 - Lead
	//	2 - Copper
	//
	if (hgcParameters.em.iAbsMaterial == 1)
		mEMAbsMat = mPb;
	else
		mEMAbsMat = mCu;

	//	Choose Material for FH Absorber
	//	1 - Brass
	//	2 - Silicon??????????????????????? Doens't make much sense here...
	//
	if (hgcParameters.fh.iAbsMaterial == 1)
		mFHAbsMat = mBrass;
	else 
		mFHAbsMat = mSi;

	//	Choose a material for BH Absorber
	//	1 - Brass
	//	2 - Silicon??????????????????????????? Again, see above
	//
	if (hgcParameters.bh.iAbsMaterial == 1)
		mBHAbsMat = mBrass;
	else
		mBHAbsMat = mSi;

	//	Set Material for Electronics
	//
	mElectronicsMat = mG10;

	//	Set the HE Abs
	//
	mHEAbsMat = mBrass;

}

//	G4 Method
//
G4VPhysicalVolume* SHDetectorConstruction::Construct()
{
	return this->BuildGeometry();
}


//	Build all the Geometry
//	Include The EM Field, 
//
G4VPhysicalVolume* SHDetectorConstruction::BuildGeometry()
{
	//	Define World Dimensions
	//
//	G4double fullWorldZ = 1*m;
//	G4double fullWorldX = 1*m;
//	G4double fullWorldY = 1*m;

	//	Create the World iteself first
	//
	solidWorld = new G4Box("solidWorld", fullWorldX/2.0, fullWorldY/2.0, 
			fullWorldZ/2.0);
	logicWorld = new G4LogicalVolume(solidWorld, mVacuum, "logicWorld",
			0,0,0);
	physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "physWorld",
			0, 0, 0, false);

	//	Double checking...
	//
	G4cout << "### World: " << fullWorldX/cm << "  " << fullWorldY/cm
		<< "  " << fullWorldZ/cm
		<< G4endl;

	//	Build the Shashlik Set of Modules
	//
//	BuildShashlik();

	//	Build HGCal
	//
//	BuildHGCal();

	//	Build HE
	//
	BuildHE();

	//	Return pointer to the WOrld
	//
	return physWorld;
}

/*
 *	Building the Shashlik itself
 */
void SHDetectorConstruction::BuildShashlik()
{
	//	SD
	//
	G4SDManager *SDManager = G4SDManager::GetSDMpointer();
	SHSDCounter *SD = new SHSDCounter("data", runParams, shTree);
	SDManager->AddNewDetector(SD);

	//	Place the Container first
	//
	solidShashlikContainer = new G4Box("solidShashlikContainer",
			fullShashlikContainerX/2., fullShashlikContainerY/2.,
			fullShashlikContainerZ/2.);
	logicShashlikContainer = new G4LogicalVolume(solidShashlikContainer,
			mVacuum, "logicShashlikContainer", 0, 0, 0); 
	physShashlikContainer = new G4PVPlacement(0, G4ThreeVector(),
			logicShashlikContainer, "physShashlikContainer", logicWorld,
			0, 0, true);

	//	Place the Shashlik Module
	//
	solidShashlik = new G4Box("solidShashlik", fullShashlikX/2.0, 
			fullShashlikY/2., fullShashlikZ/2.);
	logicShashlik = new G4LogicalVolume(solidShashlik, mVacuum, "logicShashlik",
			0, 0, 0);
	for (int iModuleX=0; iModuleX<runParams.numModules; iModuleX++)
	{
		G4double xpos = -fullShashlikContainerX/2. + 
			fullShashlikX*(iModuleX + 0.5) + iModuleX*fullGapX;
		for (int iModuleY=0; iModuleY<runParams.numModules; iModuleY++)
		{
			G4double ypos = -fullShashlikContainerY/2. + 
				fullShashlikY*(iModuleY + 0.5) + iModuleY*fullGapY;
			physShashlik = new G4PVPlacement(0, G4ThreeVector(xpos, ypos, 0),
					logicShashlik, "physShashlik", logicShashlikContainer,
					0, 2*iModuleX + iModuleY, true);
		}
	}

	//	Build the Layer Up without placing it...
	//
	solidLayer = new G4Box("solidLayer", fullLayerX/2., fullLayerY/2.,
			fullLayerZ/2);
	logicLayer = new G4LogicalVolume(solidLayer, mVacuum,
			"logicLayer", 0, 0, 0);

	//	Build the Abs Part and place it inside the Layer
	//
	solidAbs = new G4Box("solidAbs", fullAbsX/2., fullAbsY/2.,
			fullAbsZ/2.);
	logicAbs = new G4LogicalVolume(solidAbs, mW, "logicAbs", 0, 0, 0);
	physAbs = new G4PVPlacement(0, G4ThreeVector(0, 0, 
				-fullLayerZ/2. + fullAbsZ/2), logicAbs, "physAbs", logicLayer, 
			0, 0, true);

	//	Build and place AbsFibers inside the Abs
	//	4 Placements
	//
	solidAbsFiber = new G4Tubs("solidAbsFiber", inRFiber, outRFiber, 
			fullAbsFiberZ/2., 0, 360*deg);
	logicAbsFiber = new G4LogicalVolume(solidAbsFiber, mSiO2,
			"logicAbsFiber", 0, SD, 0);
	for (int i=0; i<=1; i++)
	{
		G4double xpos = -0.5*fullAbsX + 0.25*fullAbsX + i*0.5*fullAbsX;
		for (int j=0; j<=1; j++)
		{
			G4double ypos = -0.5*fullAbsY + 0.25*fullAbsY + j*0.5*fullAbsY;
			physAbsFiber = new G4PVPlacement(0, G4ThreeVector(xpos, ypos, 0),
					logicAbsFiber, "physAbsFiber", logicAbs, 0, 2*i+j, true);
		}
	}

	//	Build the Act Part and place it inside the Layer
	//
	solidAct = new G4Box("solidAct", fullActX/2., fullActY/2., fullActZ/2.);
	logicAct = new G4LogicalVolume(solidAct, mLYSO,
			"logicAct", 0, SD, 0);
	physAct = new G4PVPlacement(0, G4ThreeVector(0, 0,
				-fullLayerZ/2. + fullAbsZ + fullActZ/2.),
			logicAct, "physAct", logicLayer, 0, 0, true);

	//	Build and place ActFibers inside the Act
	//	4 Placements
	//
	solidActFiber = new G4Tubs("solidActFiber", inRFiber, outRFiber,
			fullActFiberZ/2., 0, 360*deg);
	logicActFiber = new G4LogicalVolume(solidActFiber, mSiO2,
			"logicActFiber", 0, SD, 0);
	for (int i=0; i<=1; i++)
	{
		G4double xpos = -0.5*fullActX + 0.25*fullActX + i*0.5*fullActX;
		for (int j=0; j<=1; j++)
		{
			G4double ypos = -0.5*fullActY + 0.25*fullActY + j*0.5*fullActY;
			physActFiber = new G4PVPlacement(0, G4ThreeVector(xpos, ypos,0),
					logicActFiber, "physActFiber", logicAct, 0, 2*i+j, true);
		}
	}

	//	A single Layer has been Built, but NOT Placed!!!
	//	Place all of them
	//
	for (int iLayer=0; iLayer<runParams.numLayers; iLayer++)
	{
		G4double zpos = -fullShashlikZ/2.+ fullActZ + fullLayerZ*(iLayer + 0.5);
		physLayer = new G4PVPlacement(0, G4ThreeVector(0,0,zpos),
				logicLayer, "physLayer", logicShashlik, 0, iLayer, true);
	}

	//	Note: We put Abs first into the module, but the number of LYSO plates 
	//	should be more by 1 ===>>> place a LYSO just inside the Shashlik.
	//	NOTE: We arrange a spot for that just above
	//
	physAct = new G4PVPlacement(0, G4ThreeVector(0,0,
				-fullShashlikZ/2. + fullActZ/2.), logicAct, "physAct", 
			logicShashlik, 0, runParams.numLayers, true);	

	//	Optical Surgace for the Shashlik's LYSO
	//
/*	G4OpticalSurface *opShashlikSurface = new G4OpticalSurface("ShashlikSurface");
	opShashlikSurface->SetType(dielectric_metal);
	opShashlikSurface->SetFinish(polishedtyvekair);
	opShashlikSurface->SetModel(glisur);
	G4MaterialPropertiesTable *mptSurface = new G4MaterialPropertiesTable();
	G4double reflectivity[2] = {1.0, 1.0};
	G4double photonEnergy[2] = {1.5*eV, 6.2*eV};
	mptSurface->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, 2);
	opShashlikSurface->SetMaterialPropertiesTable(mptSurface);

	G4LogicalSkinSurface *skinSurface = new G4LogicalSkinSurface("Skin",
			logicAct, opShashlikSurface);
*/
	return;
}//	end of BuildShashlik


/*
 *	Build HGCAL
 */
void SHDetectorConstruction::BuildHGCal()
{
	//	SD
	//
	G4SDManager *SDManager = G4SDManager::GetSDMpointer();
	HGSDCounter *EMSD = new HGSDCounter("EMSD", runParams, _emTree);
	HGSDCounter *FHSD = new HGSDCounter("FHSD", runParams, _fhTree);
	HGSDCounter *BHSD = new HGSDCounter("BHSD", runParams, _bhTree);
	SDManager->AddNewDetector(EMSD);
	SDManager->AddNewDetector(FHSD);
	SDManager->AddNewDetector(BHSD);

	//	Limit the step size inside the SDs
	//
//	G4UserLimits *sdLimits = new G4UserLimits(
//			hgcParameters.em.fullEMPadXYZ[2]);

	//	Read in the Input Parameters for HGCAL
	//
	ReadHGConfigFile();

	//	Additional Dimensions for EM Part
	//
	G4double fullEMLayerX_1 = hgcParameters.em.fullEMAbsXYZ[0]; 
	G4double fullEMLayerY_1 = hgcParameters.em.fullEMAbsXYZ[1];
	G4double fullEMLayerZ_1 = hgcParameters.em.fullEMAbsXYZ[2] + 
		hgcParameters.em.fullEMPadXYZ[2] + 
		hgcParameters.em.fullEMPadReadoutXYZ[2];

	G4double fullEMLayerX_2 = hgcParameters.em.fullEMAbsXYZ[0]; 
	G4double fullEMLayerY_2 = hgcParameters.em.fullEMAbsXYZ[1];
	G4double fullEMLayerZ_2 = hgcParameters.em.fullEMAbsXYZ[3] + 
		hgcParameters.em.fullEMPadXYZ[2] + 
		hgcParameters.em.fullEMPadReadoutXYZ[2];

	G4double fullEMLayerX_3 = hgcParameters.em.fullEMAbsXYZ[0]; 
	G4double fullEMLayerY_3 = hgcParameters.em.fullEMAbsXYZ[1];
	G4double fullEMLayerZ_3 = hgcParameters.em.fullEMAbsXYZ[4] + 
		hgcParameters.em.fullEMPadXYZ[2] + 
		hgcParameters.em.fullEMPadReadoutXYZ[2];

	int nPadsX_1 = fullEMLayerX_1/hgcParameters.em.fullEMPadXYZ[0];
	int nPadsY_1 = fullEMLayerY_1/hgcParameters.em.fullEMPadXYZ[1];
	int nPadsX_2 = fullEMLayerX_1/hgcParameters.em.fullEMPadXYZ[3];
	int nPadsY_2 = fullEMLayerY_1/hgcParameters.em.fullEMPadXYZ[4];

	G4double fullEMPadLayerX = fullEMLayerX_1;
	G4double fullEMPadLayerY = fullEMLayerY_1;
	G4double fullEMPadLayerZ = hgcParameters.em.fullEMPadXYZ[2];

	G4double fullEMModuleX = fullEMLayerX_1;
	G4double fullEMModuleY = fullEMLayerY_1;
	G4double fullEMModuleZ = hgcParameters.em.nLayers_1*fullEMLayerZ_1 + 
		hgcParameters.em.nLayers_2*fullEMLayerZ_2 +
		hgcParameters.em.nLayers_3*fullEMLayerZ_3;

	//	Additional Dimensions for Front HCal Part
	//
	G4double fullFHLayerX = hgcParameters.fh.fullFHAbsXYZ[0];
	G4double fullFHLayerY = hgcParameters.fh.fullFHAbsXYZ[1];
	G4double fullFHLayerZ = hgcParameters.fh.fullFHAbsXYZ[2] +
		hgcParameters.fh.fullFHPadXYZ[2] + 
		hgcParameters.fh.fullFHPadReadoutXYZ[2];

	int nFHPadsX = fullFHLayerX/hgcParameters.fh.fullFHPadXYZ[0]; 
	int nFHPadsY = fullFHLayerY/hgcParameters.fh.fullFHPadXYZ[1];

	G4double fullFHPadLayerX = fullFHLayerX;
	G4double fullFHPadLayerY = fullFHLayerY;
	G4double fullFHPadLayerZ = hgcParameters.fh.fullFHPadXYZ[2];

	G4double fullFHModuleX = fullFHLayerX; 
	G4double fullFHModuleY = fullFHLayerY;
	G4double fullFHModuleZ = hgcParameters.fh.nLayers_Total*fullFHLayerZ; 

	//	Additional Dimensions for Back HCal part
	//
	G4double fullBHLayerX = hgcParameters.bh.fullBHAbsXYZ[0];
	G4double fullBHLayerY = hgcParameters.bh.fullBHAbsXYZ[1];
	G4double fullBHLayerZ = hgcParameters.bh.fullBHAbsXYZ[2] +
		hgcParameters.bh.fullBHPadXYZ[2] + 
		hgcParameters.bh.fullBHPadReadoutXYZ[2];

	int nBHPadsX = fullBHLayerX/hgcParameters.bh.fullBHPadXYZ[0]; 
	int nBHPadsY = fullBHLayerY/hgcParameters.bh.fullBHPadXYZ[1];

	G4double fullBHPadLayerX = fullBHLayerX;
	G4double fullBHPadLayerY = fullBHLayerY;
	G4double fullBHPadLayerZ = hgcParameters.bh.fullBHPadXYZ[2];

	G4double fullBHModuleX = fullBHLayerX; 
	G4double fullBHModuleY = fullBHLayerY;
	G4double fullBHModuleZ = hgcParameters.bh.nLayers_Total*fullBHLayerZ; 
	
	G4double fullHGCalX = fullBHModuleX;
	G4double fullHGCalY = fullBHModuleY;
//	G4double fullHGCalZ = 2.*m;
	G4double fullHGCalZ = fullEMModuleZ + fullFHModuleZ + fullBHModuleZ;

	//
	//	Section for Debug Purposes: Print all the Dimensions
	//
	G4cout << "### HGCAl: " << fullHGCalX/cm << "  " << fullHGCalY/cm
		<< "  " << fullHGCalZ/cm
		<< G4endl
		<< "### Parameters: " << hgcParameters.em.nLayers_Total
		<< "  " << hgcParameters.em.nLayers_1 << "  "
		<< hgcParameters.em.nLayers_2 << "  "
		<< hgcParameters.em.nLayers_3 << "  "
		<< hgcParameters.em.iAbsMaterial
		<< G4endl
		<< "### EM: " << fullEMModuleX/cm << "  " << fullEMModuleY/cm
		<< "  " << fullEMModuleZ/cm
		<< G4endl
		<< "### EM: Layer_1: " << fullEMLayerX_1/cm << "  "
		<< fullEMLayerY_1/cm << "  " << fullEMLayerZ_1/cm
		<< G4endl
		<< "### EM: Layer_1: Abs: " << hgcParameters.em.fullEMAbsXYZ[0]/cm
		<< "  " << hgcParameters.em.fullEMAbsXYZ[1]/cm
		<< "  " << hgcParameters.em.fullEMAbsXYZ[2]/cm
		<< G4endl
		<< "### EM: Layer_1: PadLayer: " << fullEMPadLayerX/cm
		<< "  " << fullEMPadLayerY/cm << "  " << fullEMPadLayerZ/cm
		<< G4endl
		<< "### EM: Layer_1: PadReadout: "
		<< hgcParameters.em.fullEMPadReadoutXYZ[0]/cm
		<< "  " << hgcParameters.em.fullEMPadReadoutXYZ[1]/cm << "  "
		<< hgcParameters.em.fullEMPadReadoutXYZ[2]/cm
		<< G4endl
		<< "### EM: Layer_2: " << fullEMLayerX_2/cm << "  "
		<< fullEMLayerY_2/cm << "  " << fullEMLayerZ_2/cm
		<< G4endl
		<< "### EM: Layer_3: " << fullEMLayerX_3/cm << "  "
		<< fullEMLayerY_3/cm << "  " << fullEMLayerZ_3/cm
		<< G4endl
		<< "### FH: " << fullFHModuleX/cm << "  " << fullFHModuleY/cm
		<< "  " << fullFHModuleZ/cm
		<< G4endl
		<< "### Parameters: " << hgcParameters.fh.nLayers_Total
		<< "  " << hgcParameters.fh.iAbsMaterial
		<< G4endl
		<< "### FH: Layer: " << fullFHLayerX/cm << "  " << fullFHLayerY/cm
		<< "  " << fullFHLayerZ/cm
		<< G4endl
		<< "### BH: " << fullBHModuleX/cm << "  " << fullBHModuleY/cm
		<< "  " << fullBHModuleZ/cm
		<< G4endl
		<< "### Parameters: " << hgcParameters.bh.nLayers_Total
		<< "  " << hgcParameters.bh.iAbsMaterial
		<< G4endl
		<< "### BH: Layer: " << fullBHLayerX/cm << "  " << fullBHLayerY/cm
		<< "  " << fullBHLayerZ/cm
		<< G4endl;

	cout << "### HERE: " << fullHGCalX/cm << "  " << fullHGCalY/cm << "  "
		<< fullHGCalZ/cm << endl;
	cout << "### HERE: " << fullWorldX/cm << "  " << fullWorldY/cm << "  "
		<< fullWorldZ/cm << endl;

	//	Place the Whole HGCAL first
	//
	solidHGCal = new G4Box("solidHGCal", fullHGCalX/2., fullHGCalY/2., 
			fullHGCalZ/2.);
	logicHGCal = new G4LogicalVolume(solidHGCal, mVacuum, "logicHGCal");
	physHGCal = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicHGCal, 
			"physHGCal", logicWorld, 0, 0, true);


	//	Build the EM Part
	//
	solidHG_EMModule = new G4Box("solidHG_EMModule", fullEMModuleX/2., 
			fullEMModuleY/2., fullEMModuleZ/2.);
	logicHG_EMModule = new G4LogicalVolume(solidHG_EMModule, mVacuum, 
			"logicHG_EMModule");
	G4double zpos = -fullHGCalZ/2. + fullEMModuleZ/2;
	if (hgcParameters.emOnOff)
		physHG_EMModule = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMModule, "physHG_EMModule", logicHGCal, 0, 0, true);

	//
	//	There are 3 sections for EM part
	//	See pdf for specifics
	//	
	
	//	First Part of EM;
	//
	solidHG_EMLayer[0] = new G4Box("solidHG_EMLayer_1",
			fullEMLayerX_1/2., fullEMLayerY_1/2., fullEMLayerZ_1/2.);
	logicHG_EMLayer[0] = new G4LogicalVolume(solidHG_EMLayer[0], mVacuum,
			"logicHG_EMLayer_1");

	solidHG_EMAbs[0] = new G4Box("solidHG_EMAbs_1",
			hgcParameters.em.fullEMAbsXYZ[0]/2., 
			hgcParameters.em.fullEMAbsXYZ[1]/2.,
			hgcParameters.em.fullEMAbsXYZ[2]/2.);
	logicHG_EMAbs[0] = new G4LogicalVolume(solidHG_EMAbs[0],
			mEMAbsMat, "logicHG_EMAbs_1");
	zpos = -fullEMLayerZ_1/2. + hgcParameters.em.fullEMAbsXYZ[2]/2.;
	physHG_EMAbs = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMAbs[0], "phys_HF_EMAbs", logicHG_EMLayer[0], 0, 0, true);

	solidHG_EMPadLayer = new G4Box("solidHG_EMPadLayer",
			fullEMPadLayerX/2., fullEMPadLayerY/2., fullEMPadLayerZ/2.);
	logicHG_EMPadLayer = new G4LogicalVolume(solidHG_EMPadLayer,
			mSi, "logicHG_EMPadLayer", 0, EMSD, 0);
//	logicHG_EMPadLayer->SetUserLimits(sdLimits);
	zpos = -fullEMLayerZ_1/2. + hgcParameters.em.fullEMAbsXYZ[2] +
		fullEMPadLayerZ/2.;
	physHG_EMPadLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMPadLayer, "physHG_EMPadLayer", logicHG_EMLayer[0], 0, 0,
			true);
	
	solidHG_EMPadReadout = new G4Box("solidHG_EMPadReadout",
			hgcParameters.em.fullEMPadReadoutXYZ[0]/2.,
			hgcParameters.em.fullEMPadReadoutXYZ[1]/2.,
			hgcParameters.em.fullEMPadReadoutXYZ[2]/2.);
	logicHG_EMPadReadout = new G4LogicalVolume(solidHG_EMPadReadout,
			mElectronicsMat, "logicHG_EMPadReadout");
	zpos = -fullEMLayerZ_1/2. + hgcParameters.em.fullEMAbsXYZ[2] +
		fullEMPadLayerZ + hgcParameters.em.fullEMPadReadoutXYZ[2]/2.;
	physHG_EMPadReadout = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMPadReadout, "physHG_EMPadReadout", 
			logicHG_EMLayer[0], 0, 0, true);

	for (int iLayer=0; iLayer<hgcParameters.em.nLayers_1; iLayer++)
	{
		zpos = -fullEMModuleZ/2. + (iLayer + 0.5)*fullEMLayerZ_1;
		physHG_EMLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
				logicHG_EMLayer[0], "physHG_EMLayer", logicHG_EMModule,
				0, iLayer, true);
	}

	//	Second Part of EM
	//
	solidHG_EMLayer[1] = new G4Box("solidHG_EMLayer_2",
			fullEMLayerX_2/2., fullEMLayerY_2/2., fullEMLayerZ_2/2.);
	logicHG_EMLayer[1] = new G4LogicalVolume(solidHG_EMLayer[1], mVacuum,
			"logicHG_EMLayer_2");

	solidHG_EMAbs[1] = new G4Box("solidHG_EMAbs_2",
			hgcParameters.em.fullEMAbsXYZ[0]/2., 
			hgcParameters.em.fullEMAbsXYZ[1]/2.,
			hgcParameters.em.fullEMAbsXYZ[3]/2.);
	logicHG_EMAbs[1] = new G4LogicalVolume(solidHG_EMAbs[1],
			mEMAbsMat, "logicHG_EMAbs_2");
	zpos = -fullEMLayerZ_2/2. + hgcParameters.em.fullEMAbsXYZ[3]/2.;
	physHG_EMAbs = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMAbs[1], "phys_HF_EMAbs", logicHG_EMLayer[1], 0, 0, true);

	solidHG_EMPadLayer = new G4Box("solidHG_EMPadLayer",
			fullEMPadLayerX/2., fullEMPadLayerY/2., fullEMPadLayerZ/2.);
	logicHG_EMPadLayer = new G4LogicalVolume(solidHG_EMPadLayer,
			mSi, "logicHG_EMPadLayer", 0, EMSD, 0);
//	logicHG_EMPadLayer->SetUserLimits(sdLimits);
	zpos = -fullEMLayerZ_2/2. + hgcParameters.em.fullEMAbsXYZ[3] +
		fullEMPadLayerZ/2.;
	physHG_EMPadLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMPadLayer, "physHG_EMPadLayer", logicHG_EMLayer[1], 0, 0,
			true);
	
	solidHG_EMPadReadout = new G4Box("solidHG_EMPadReadout",
			hgcParameters.em.fullEMPadReadoutXYZ[0]/2.,
			hgcParameters.em.fullEMPadReadoutXYZ[1]/2.,
			hgcParameters.em.fullEMPadReadoutXYZ[2]/2.);
	logicHG_EMPadReadout = new G4LogicalVolume(solidHG_EMPadReadout,
			mElectronicsMat, "logicHG_EMPadReadout");
	zpos = -fullEMLayerZ_2/2. + hgcParameters.em.fullEMAbsXYZ[3] +
		fullEMPadLayerZ + hgcParameters.em.fullEMPadReadoutXYZ[2]/2.;
	physHG_EMPadReadout = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMPadReadout, "physHG_EMPadReadout", 
			logicHG_EMLayer[1], 0, 0, true);

	for (int iLayer=0; iLayer<hgcParameters.em.nLayers_2; iLayer++)
	{
		zpos = -fullEMModuleZ/2. + hgcParameters.em.nLayers_1*fullEMLayerZ_1 + 
			(iLayer + 0.5)*fullEMLayerZ_2;
		physHG_EMLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
				logicHG_EMLayer[1], "physHG_EMLayer", logicHG_EMModule,
				0, iLayer+hgcParameters.em.nLayers_1, true);
	}

	//	Third part of EM
	//
	solidHG_EMLayer[2] = new G4Box("solidHG_EMLayer_3",
			fullEMLayerX_3/2., fullEMLayerY_3/2., fullEMLayerZ_3/2.);
	logicHG_EMLayer[2] = new G4LogicalVolume(solidHG_EMLayer[2], mVacuum,
			"logicHG_EMLayer_3");
	
	solidHG_EMAbs[2] = new G4Box("solidHG_EMAbs_3",
			hgcParameters.em.fullEMAbsXYZ[0]/2., 
			hgcParameters.em.fullEMAbsXYZ[1]/2.,
			hgcParameters.em.fullEMAbsXYZ[4]/2.);
	logicHG_EMAbs[2] = new G4LogicalVolume(solidHG_EMAbs[2],
			mEMAbsMat, "logicHG_EMAbs_3");
	zpos = -fullEMLayerZ_3/2. + hgcParameters.em.fullEMAbsXYZ[4]/2.;
	physHG_EMAbs = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMAbs[2], "phys_HF_EMAbs", logicHG_EMLayer[2], 0, 0, true);

	solidHG_EMPadLayer = new G4Box("solidHG_EMPadLayer",
			fullEMPadLayerX/2., fullEMPadLayerY/2., fullEMPadLayerZ/2.);
	logicHG_EMPadLayer = new G4LogicalVolume(solidHG_EMPadLayer,
			mSi, "logicHG_EMPadLayer", 0, EMSD, 0);
//	logicHG_EMPadLayer->SetUserLimits(sdLimits);
	zpos = -fullEMLayerZ_3/2. + hgcParameters.em.fullEMAbsXYZ[4] +
		fullEMPadLayerZ/2.;
	physHG_EMPadLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMPadLayer, "physHG_EMPadLayer", logicHG_EMLayer[2], 0, 0,
			true);
	
	solidHG_EMPadReadout = new G4Box("solidHG_EMPadReadout",
			hgcParameters.em.fullEMPadReadoutXYZ[0]/2.,
			hgcParameters.em.fullEMPadReadoutXYZ[1]/2.,
			hgcParameters.em.fullEMPadReadoutXYZ[2]/2.);
	logicHG_EMPadReadout = new G4LogicalVolume(solidHG_EMPadReadout,
			mElectronicsMat, "logicHG_EMPadReadout");
	zpos = -fullEMLayerZ_3/2. + hgcParameters.em.fullEMAbsXYZ[4] +
		fullEMPadLayerZ + hgcParameters.em.fullEMPadReadoutXYZ[2]/2.;
	physHG_EMPadReadout = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_EMPadReadout, "physHG_EMPadReadout", 
			logicHG_EMLayer[2], 0, 0, true);

	for (int iLayer=0; iLayer<hgcParameters.em.nLayers_3; iLayer++)
	{
		zpos = -fullEMModuleZ/2. + hgcParameters.em.nLayers_1*fullEMLayerZ_1 + 
			hgcParameters.em.nLayers_2*fullEMLayerZ_2 + 
			(iLayer + 0.5)*fullEMLayerZ_3;
		physHG_EMLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
				logicHG_EMLayer[2], "physHG_EMLayer", logicHG_EMModule,
				0, iLayer + hgcParameters.em.nLayers_1 + 
				hgcParameters.em.nLayers_2, true);
	}

	//
	//	Place the Front HCAL Part
	//
	solidHG_FHModule = new G4Box("solidHG_FHModule",
			fullFHModuleX/2., fullFHModuleY/2., fullFHModuleZ/2.);
	logicHG_FHModule = new G4LogicalVolume(solidHG_FHModule,
			mFHAbsMat, "logicHG_FHModule");
	zpos = -fullHGCalZ/2. + fullEMModuleZ + fullFHModuleZ/2.;

	if (hgcParameters.fhOnOff)
		physHG_FHModule = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_FHModule, "physHG_FHModule", logicHGCal, 0, 0, true);

	//	Just 1 section
	//
	solidHG_FHLayer = new G4Box("solidHG_FHLayer",
			fullFHLayerX/2., fullFHLayerY/2., fullFHLayerZ/2.);
	logicHG_FHLayer = new G4LogicalVolume(solidHG_FHLayer, mVacuum,
			"logicHG_FHLayer");

	solidHG_FHAbs = new G4Box("solidHG_FHAbs",
			hgcParameters.fh.fullFHAbsXYZ[0]/2.,
			hgcParameters.fh.fullFHAbsXYZ[1]/2.,
			hgcParameters.fh.fullFHAbsXYZ[2]/2.);
	logicHG_FHAbs = new G4LogicalVolume(solidHG_FHAbs, mFHAbsMat, 
			"logicHG_FHAbs");
	zpos = -fullFHLayerZ/2. + hgcParameters.fh.fullFHAbsXYZ[2]/2.;
	physHG_FHAbs = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_FHAbs, "physHG_FHAbs", logicHG_FHLayer, 0, 0, true);

	solidHG_FHPadLayer = new G4Box("solidHG_FHPadLayer",
			fullFHPadLayerX/2., fullFHPadLayerY/2., fullFHPadLayerZ/2.);
	logicHG_FHPadLayer = new G4LogicalVolume(solidHG_FHPadLayer,
			mSi, "logicHG_FHPadLayer", 0, FHSD, 0);
//	logicHG_FHPadLayer->SetUserLimits(sdLimits);
	zpos = -fullFHLayerZ/2. + hgcParameters.fh.fullFHAbsXYZ[2] + 
		fullFHPadLayerZ/2.;
	physHG_FHPadLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_FHPadLayer, "physHG_FHPadLayer", logicHG_FHLayer, 0,0,
			true);

	solidHG_FHPadReadout = new G4Box("solidHG_FHPadReadout",
			hgcParameters.fh.fullFHPadReadoutXYZ[0]/2.,
			hgcParameters.fh.fullFHPadReadoutXYZ[1]/2.,
			hgcParameters.fh.fullFHPadReadoutXYZ[2]/2.);
	logicHG_FHPadReadout = new G4LogicalVolume(solidHG_FHPadReadout,
			mElectronicsMat, "logicHG_FHPadReadout");
	zpos = -fullFHLayerZ/2. + hgcParameters.fh.fullFHAbsXYZ[2] + 
			fullFHPadLayerZ + hgcParameters.fh.fullFHPadReadoutXYZ[2]/2.;
	physHG_FHPadReadout = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_FHPadReadout, "physHG_FHPadRedout",
			logicHG_FHLayer, 0, 0, true);

	for (int iLayer=0; iLayer<hgcParameters.fh.nLayers_Total; iLayer++)
	{
		zpos = -fullFHModuleZ/2. + (iLayer + 0.5)*fullFHLayerZ;
		physHG_FHLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
				logicHG_FHLayer, "physHG_FHLayer", logicHG_FHModule,
				0, iLayer + hgcParameters.em.nLayers_Total, true);
	}

	//
	//	Place the back HCAL Part
	//
	solidHG_BHModule = new G4Box("solidHG_BHModule",
			fullBHModuleX/2., fullBHModuleY/2., fullBHModuleZ/2.);
	logicHG_BHModule = new G4LogicalVolume(solidHG_BHModule,
			mBHAbsMat, "logicHG_BHModule");
	zpos = -fullHGCalZ/2. + fullEMModuleZ + fullFHModuleZ +
		fullBHModuleZ/2.;

	if (hgcParameters.bhOnOff)
		physHG_BHModule = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_BHModule, "physHG_BHModule", logicHGCal, 0, 0, true);

	//	Just 1 section
	//
	solidHG_BHLayer = new G4Box("solidHG_BHLayer",
			fullBHLayerX/2., fullBHLayerY/2., fullBHLayerZ/2.);
	logicHG_BHLayer = new G4LogicalVolume(solidHG_BHLayer, mVacuum,
			"logicHG_BHLayer");

	solidHG_BHAbs = new G4Box("solidHG_BHAbs",
			hgcParameters.bh.fullBHAbsXYZ[0]/2.,
			hgcParameters.bh.fullBHAbsXYZ[1]/2.,
			hgcParameters.bh.fullBHAbsXYZ[2]/2.);
	logicHG_BHAbs = new G4LogicalVolume(solidHG_BHAbs, mBHAbsMat, 
			"logicHG_BHAbs");
	zpos = -fullBHLayerZ/2. + hgcParameters.bh.fullBHAbsXYZ[2]/2.;
	physHG_BHAbs = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_BHAbs, "physHG_BHAbs", logicHG_BHLayer, 0, 0, true);

	solidHG_BHPadLayer = new G4Box("solidHG_BHPadLayer",
			fullBHPadLayerX/2., fullBHPadLayerY/2., fullBHPadLayerZ/2.);
	logicHG_BHPadLayer = new G4LogicalVolume(solidHG_BHPadLayer,
			mSi, "logicHG_BHPadLayer", 0, BHSD, 0);
//	logicHG_BHPadLayer->SetUserLimits(sdLimits);
	zpos = -fullBHLayerZ/2. + hgcParameters.bh.fullBHAbsXYZ[2] + 
		fullBHPadLayerZ/2.;
	physHG_BHPadLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_BHPadLayer, "physHG_BHPadLayer", logicHG_BHLayer, 0,0,
			true);

	solidHG_BHPadReadout = new G4Box("solidHG_BHPadReadout",
			hgcParameters.bh.fullBHPadReadoutXYZ[0]/2.,
			hgcParameters.bh.fullBHPadReadoutXYZ[1]/2.,
			hgcParameters.bh.fullBHPadReadoutXYZ[2]/2.);
	logicHG_BHPadReadout = new G4LogicalVolume(solidHG_BHPadReadout,
			mElectronicsMat, "logicHG_BHPadReadout");
	zpos = -fullBHLayerZ/2. + hgcParameters.bh.fullBHAbsXYZ[2] + 
			fullBHPadLayerZ + hgcParameters.bh.fullBHPadReadoutXYZ[2]/2.;
	physHG_BHPadReadout = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHG_BHPadReadout, "physHG_BHPadRedout",
			logicHG_BHLayer, 0, 0, true);

	for (int iLayer=0; iLayer<hgcParameters.bh.nLayers_Total; iLayer++)
	{
		zpos = -fullBHModuleZ/2. + (iLayer + 0.5)*fullBHLayerZ;
		physHG_FHLayer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
				logicHG_BHLayer, "physHG_BHLayer", logicHG_BHModule,
				0, iLayer + hgcParameters.em.nLayers_Total +
				hgcParameters.fh.nLayers_Total, true);
	}	


	return;
}//	end of BuildHGCal

/*
 *	Read in HGCAL configuration data
 *	NOTE: 
 *	--	All the input sizes are in mm
 *	--	Also, the order of dimensios is exactly as in SHDefs.hh
 */
int SHDetectorConstruction::ReadHGConfigFile()
{
	cout << "### Reading in HGCAL Configuration File..." << endl;

	//	Init/Open/Check file
	//
	ifstream hgConfigFile(runParams.hgCalInputFileName);
	if (!hgConfigFile)
	{
		cout << "### ERROR: File " << runParams.hgCalInputFileName
			<< "  hasn't been found!!!" << endl;
		return -1;
	}

	//	Read in/Config
	//
	double x,y,z;
	int n, n1, n2, n3;
	int iMat;
	int emOnOff, fhOnOff, bhOnOff;

	//	Input the HGCal On/Off Settings
	//
	hgConfigFile >> emOnOff >> fhOnOff >> bhOnOff;
	hgcParameters.emOnOff =emOnOff;
	hgcParameters.fhOnOff = fhOnOff;
	hgcParameters.bhOnOff = bhOnOff;
	
	//	EM Input Part
	//
	hgConfigFile >> n >> n1 >> n2 >> n3 >> iMat;
	hgcParameters.em.nLayers_Total = n;
	hgcParameters.em.nLayers_1 = n1;
	hgcParameters.em.nLayers_2 = n2;
	hgcParameters.em.nLayers_3 = n3;
	hgcParameters.em.iAbsMaterial = iMat;

	hgConfigFile >> x >> y >> z;
	hgcParameters.em.fullEMAbsXYZ[0] = x*mm;
	hgcParameters.em.fullEMAbsXYZ[1] = y*mm;
	hgcParameters.em.fullEMAbsXYZ[2] = z*mm;

	hgConfigFile >> z;
	hgcParameters.em.fullEMAbsXYZ[3] = z*mm;
	hgConfigFile >> z;
	hgcParameters.em.fullEMAbsXYZ[4] = z*mm;

	hgConfigFile >> x >> y >> z;
	hgcParameters.em.fullEMPadXYZ[0] = x*mm;
	hgcParameters.em.fullEMPadXYZ[1] = y*mm;
	hgcParameters.em.fullEMPadXYZ[2] = z*mm;

	hgConfigFile >> x >> y;
	hgcParameters.em.fullEMPadXYZ[3] = x*mm;
	hgcParameters.em.fullEMPadXYZ[4] = y*mm;

	hgConfigFile >> x >> y >> z;
	hgcParameters.em.fullEMPadReadoutXYZ[0] = x*mm;
	hgcParameters.em.fullEMPadReadoutXYZ[1] = y*mm;
	hgcParameters.em.fullEMPadReadoutXYZ[2] = z*mm;

	//	FH Input Part
	//
	hgConfigFile >> n >> iMat;
	hgcParameters.fh.nLayers_Total = n;
	hgcParameters.fh.iAbsMaterial = iMat;

	hgConfigFile >> x >> y >> z;
	hgcParameters.fh.fullFHAbsXYZ[0] = x*mm;
	hgcParameters.fh.fullFHAbsXYZ[1] = y*mm;
	hgcParameters.fh.fullFHAbsXYZ[2] = z*mm;

	hgConfigFile >> x >> y >> z;
	hgcParameters.fh.fullFHPadXYZ[0] = x*mm;
	hgcParameters.fh.fullFHPadXYZ[1] = y*mm;
	hgcParameters.fh.fullFHPadXYZ[2] = z*mm;

	hgConfigFile >> x >> y >> z;
	hgcParameters.fh.fullFHPadReadoutXYZ[0] = x*mm;
	hgcParameters.fh.fullFHPadReadoutXYZ[1] = y*mm;
	hgcParameters.fh.fullFHPadReadoutXYZ[2] = z*mm;

	//	BH Input Part
	//
	hgConfigFile >> n >> iMat;
	hgcParameters.bh.nLayers_Total = n;
	hgcParameters.bh.iAbsMaterial = iMat;

	hgConfigFile >> x >> y >> z;
	hgcParameters.bh.fullBHAbsXYZ[0] = x*mm;
	hgcParameters.bh.fullBHAbsXYZ[1] = y*mm;
	hgcParameters.bh.fullBHAbsXYZ[2] = z*mm;

	hgConfigFile >> x >> y >> z;
	hgcParameters.bh.fullBHPadXYZ[0] = x*mm;
	hgcParameters.bh.fullBHPadXYZ[1] = y*mm;
	hgcParameters.bh.fullBHPadXYZ[2] = z*mm;

	hgConfigFile >> x >> y >> z;
	hgcParameters.bh.fullBHPadReadoutXYZ[0] = x*mm;
	hgcParameters.bh.fullBHPadReadoutXYZ[1] = y*mm;
	hgcParameters.bh.fullBHPadReadoutXYZ[2] = z*mm;

	return 1;
}//	end of Read HGCAL Config Data

//
//	Build HE
//
void SHDetectorConstruction::BuildHE()
{
	//	SD
	//
	G4SDManager *SDManager = G4SDManager::GetSDMpointer();
	SHSDCounter *SHSD = new SHSDCounter("SHSD", runParams, _heTree);
	SDManager->AddNewDetector(SHSD);


	//	Read in Config File
	//
	ReadHEConfigFile();

	//	Additional dimensions for HE: layer's dimensions
	//
	G4double fullHELayerX = heParameters.fullHEAbsXYZ[0];
	G4double fullHELayerY = heParameters.fullHEAbsXYZ[1];
	G4double fullHELayerZ = heParameters.fullHEAbsXYZ[2] +
		heParameters.fullHEActXYZ[2];

	//	HE Module itself
	//
	G4double fullHEModuleX = fullHELayerX;
	G4double fullHEModuleY = fullHELayerY;
	G4double fullHEModuleZ = heParameters.nLayers_Total*fullHELayerZ;

	//	Fot debugging...
	//
	G4cout << "### HE Module Dimensions: " << fullHEModuleX/cm << "  "
		<< fullHEModuleY/cm << "  " << fullHEModuleZ/cm
		<< G4endl
		<< "### HE Layer: " << fullHELayerX/cm << "  " << fullHELayerY/cm
		<< "  " << fullHELayerZ/cm
		<< G4endl
		<< "### HE Abs: " << heParameters.fullHEAbsXYZ[0]/cm
		<< "  " << heParameters.fullHEAbsXYZ[1]/cm << "  "
		<< heParameters.fullHEAbsXYZ[2]/cm
		<< G4endl;

	//	Place the HE Module itself first
	//
	solidHE_Module = new G4Box("solidHE_Module", 
			fullHEModuleX/2., fullHEModuleY/2., fullHEModuleZ/2.);
	logicHE_Module = new G4LogicalVolume(solidHE_Module,
			mVacuum, "logicHE_Module");
	physHE_Module = new G4PVPlacement(0, G4ThreeVector(), 
			logicHE_Module, "physHE_Module", logicWorld, 0, 0, true);

	//	Construct the Layer(Place it in the very end)
	//
	solidHE_Layer = new G4Box("solidHE_Layer", 
			fullHELayerX/2., fullHELayerY/2., fullHELayerZ/2.);
	logicHE_Layer = new G4LogicalVolume(solidHE_Layer, mVacuum,
			"logicHE_Layer");

	//	Abs
	//
	solidHE_Abs = new G4Box("solidHE_Abs", heParameters.fullHEAbsXYZ[0]/2,
			heParameters.fullHEAbsXYZ[1]/2., heParameters.fullHEAbsXYZ[2]/2.);
	logicHE_Abs = new G4LogicalVolume(solidHE_Abs, mHEAbsMat, "logicHE_Abs");
	G4double zpos = -fullHELayerZ/2. + heParameters.fullHEAbsXYZ[2]/2.;
	physHE_Abs = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
			logicHE_Abs, "physHE_Abs", logicHE_Layer, 0, 0, true);

	//	Scintillator
	//
	solidHE_Act = new G4Box("solidHE_Act", heParameters.fullHEActXYZ[0]/2.,
			heParameters.fullHEActXYZ[1]/2., heParameters.fullHEActXYZ[2]/2.);
	logicHE_Act = new G4LogicalVolume(solidHE_Act, mHEScintMain, 
			"logicHE_Act", 0, SHSD, 0);
	zpos = -fullHELayerZ/2. + heParameters.fullHEAbsXYZ[2] +
		heParameters.fullHEActXYZ[2]/2.;
	physHE_Act = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
		logicHE_Act, "physHE_Act", logicHE_Layer, 0, 0, true);

	//	Place all the Layeers
	//
	for (int iLayer=0; iLayer<heParameters.nLayers_Total; iLayer++)
	{
		zpos = -fullHEModuleZ/2. + (iLayer + 0.5)*fullHELayerZ;
		physHE_Layer = new G4PVPlacement(0, G4ThreeVector(0, 0, zpos),
				logicHE_Layer, "physHE_Layer", logicHE_Module,
				0, iLayer, true);
	}
}

//
//	Read in HE Parameters
//
int SHDetectorConstruction::ReadHEConfigFile()
{
	cout << "### Reading in HE Configuraiton File...." << endl;

	//	Init/Open/Check File
	//
	ifstream heConfigFile(runParams.heInputFileName);
	if (!heConfigFile)
	{
		cout << "### ERROR: File " << runParams.heInputFileName
			<< " hasn't been found!!!" << endl;
		return -1;
	}

	//	Read in/Config
	//
	double x,y,z;
	int n;

	heConfigFile >> n;
	heParameters.nLayers_Total = n;

	heConfigFile >> x >> y >> z;
	heParameters.fullHEAbsXYZ[0] = x*mm;
	heParameters.fullHEAbsXYZ[1] = y*mm;
	heParameters.fullHEAbsXYZ[2] = z*mm;

	heConfigFile >> x >> y >> z;
	heParameters.fullHEActXYZ[0] = x*mm;
	heParameters.fullHEActXYZ[1] = y*mm;
	heParameters.fullHEActXYZ[2] = z*mm;

	return 1;
}














