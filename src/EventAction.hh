#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "g4root.hh"

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    // method from the base class
    virtual void BeginOfEventAction(const G4Event *);
    virtual void EndOfEventAction(const G4Event *);
    //G4int eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
    //G4int volumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();

    void AddEnergyDeposit(G4int volume, G4double edep);
    void FindPreStepKinetic(G4int volume, G4double kineticE);

    void FindFAZIAMult(G4int partiID, G4int multp);
    void FindINDRAMult(G4int partiID, G4int mult);
    
    void FindPIDandVolumeID(G4int volume, G4int partiID);
    //void FindAnumber(G4int partiID);
  private:
    G4int eventID=0;
    G4int FAZIA_Mult=0;
    G4int INDRA_Mult=0;
    G4int partiID = 0;
    G4int PID = 0;
    //G4int Anum = 0;

    G4int mult=0;
    G4int multp=0;
    //G4int indramu;
    G4int VolumeID = 0;
    G4double edepTot[960];
    G4double edepTotModi[960];
    G4double PrekineticE[960];
    //G4int volI = 
};

#endif
