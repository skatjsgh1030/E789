#ifndef FAZIAPRIMARYGENERATORACTION_HH
#define FAZIAPRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "globals.hh"

class FAZIAPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    FAZIAPrimaryGeneratorAction();    
    virtual ~FAZIAPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  private:
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    G4double random;
    G4double phirand;
    G4double costheta;
    G4double Parti_ene;

    G4int eventNum = 0;

};

#endif
