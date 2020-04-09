#ifndef FAZIA_HH
#define FAZIA_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class FAZIA
{
  public:
    
	FAZIA();    
	virtual ~FAZIA();

	void faziaGeo(G4LogicalVolume* logicWorld);
};

#endif
