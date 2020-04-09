#ifndef FAZIAEVENTACTION_HH
#define FAZIAEVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "g4root.hh"

class FAZIAEventAction : public G4UserEventAction
{
  public:
	FAZIAEventAction();
	virtual ~FAZIAEventAction();

	// method from the base class
	virtual void BeginOfEventAction(const G4Event *);
	virtual void EndOfEventAction(const G4Event *);
	//G4int eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
	//G4int volumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();

	  void AddEnergyDeposit(G4int volumeID, G4double edep);
	  void FindPreStepKinetic(G4int volumeID, G4double kineticE);

  private:
	G4int eventID;
	//G4int volumeID;
	G4double edepTot[593];
	G4double edepTotModi[593];
	G4double PrekineticE[593];
	//G4int volI = 
};

#endif
