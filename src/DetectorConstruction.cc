#include "DetectorConstruction.hh"
#include "INDRA.hh"
#include "FAZIA.hh"

#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iostream>
#include "TMath.h"

using std::stringstream;
using namespace std;
using namespace CLHEP;
  
  DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
	logicWorld = nullptr;
	physWorld = nullptr;
	solidWorld = nullptr;
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //----------------Material --------------------
  G4NistManager* nist = G4NistManager::Instance();

  const G4double labTemp = CLHEP::STP_Temperature + 20.*kelvin;

  // -----------------------------------------------------
  // World

  G4Material* world_mat = nist -> FindOrBuildMaterial("G4_Galactic");
  //G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");		
  G4double world_size = 5*m;

  solidWorld =
	new G4Box("World",                       // its name
		3*world_size,                // half x
		3*world_size,                // half y
		3*world_size);               // half zNDDetectorConstruction::LANDDetectorConstruction()


  logicWorld =
	new G4LogicalVolume(solidWorld,          //its solid
		world_mat,           //its material
		"World");            //its name

  physWorld =
	new G4PVPlacement(0,                     //no rotation
		G4ThreeVector(),       //at (0,0,0)
		logicWorld,            //its logical volume
		"World",               //its name
		0,                     //its mother  volume
		false,                 //no boolean operation
		0,                     //copy number
		true);                 //overlaps checking

  logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());

  FAZIA forward;
  INDRA backward;

  if(1)
  {	
 	forward.faziaGeo(logicWorld);
  }	
  if(1)
  {	
 	backward.indraGeo(logicWorld);
  }	

  return physWorld;
}

