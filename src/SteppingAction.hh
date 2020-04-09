#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "g4root.hh"
//#include "g4cvs.hh"

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
    virtual ~SteppingAction();
    G4int f_mult = 0;
    G4int i_mult = 0;
    G4int volume = 0;
    // method from the base class
    virtual void UserSteppingAction(const G4Step*);
  private:
    //G4int fScintillationCounter;
    //G4int fCerenkovCounter;
    //G4int fEventNumber;
};

#endif
