#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "INDRA.hh"
#include "FAZIA.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class FAZIA;
class INDRA;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
  
  	G4LogicalVolume* logicWorld;

  private:
  	G4Box* solidWorld;
	G4VPhysicalVolume* physWorld;
};

#endif
