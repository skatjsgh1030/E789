#ifndef INDRAEVENTACTION_HH
#define INDRAEVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "g4root.hh"

class INDRAEventAction : public G4UserEventAction
{
  public:
    INDRAEventAction();
    virtual ~INDRAEventAction();

    // method from the base class
    virtual void BeginOfEventAction(const G4Event *);
    virtual void EndOfEventAction(const G4Event *);

    void AddEnergyDeposit(G4int volumeID, G4double edep);  
    void FindPreStepKinetic(G4int volumeID, G4double kineticE);

  private:
    G4double eventID;

    G4double INDRAedepTot[375];
    G4double INDRAedepTotModi[375];
    G4double INDRAPrekineticE[375];
};

#endif
