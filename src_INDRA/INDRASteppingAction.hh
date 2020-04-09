#ifndef INDRASTEPPINGACTION_HH
#define INDRASTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "g4root.hh"
//#include "g4cvs.hh"

class INDRASteppingAction : public G4UserSteppingAction
{
  public:
    INDRASteppingAction();
    virtual ~INDRASteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);
  
  private:
    //G4int fScintillationCounter;
    //G4int fCerenkovCounter;
    //G4int fEventNumber;
};

#endif
