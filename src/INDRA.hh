#ifndef INDRA_HH
#define INDRA_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"

//using std::stringstream;
//using namespace std;
//using namespace CLHEP;

class G4VPhysicalVolume;
class G4LogicalVolume;

class INDRA 
{
  public:
    INDRA();
    virtual ~INDRA();
	
	void indraGeo(G4LogicalVolume* logicWorld);
};

#endif
