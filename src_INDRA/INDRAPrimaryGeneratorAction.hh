#ifndef INDRAPRIMARYGENERATORACTION_HH
#define INDRAPRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "globals.hh"

class INDRAPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    INDRAPrimaryGeneratorAction();    
    virtual ~INDRAPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);    
  
    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  
  private:
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    G4int eventNum = 0;
};

#endif

