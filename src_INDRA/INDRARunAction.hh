#ifndef INDRARUNACTION_HH
#define INDRARUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "globals.hh"
#include "g4root.hh"

class INDRARunAction : public G4UserRunAction
{
  public:
    INDRARunAction();
    virtual ~INDRARunAction();

    // method from the base class
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
};

#endif
