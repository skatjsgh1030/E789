#include "INDRAPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

INDRAPrimaryGeneratorAction::INDRAPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
  G4ThreeVector gunposition = G4ThreeVector(0,0,-30*cm);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable -> FindParticle(particleName = "alpha");

  fParticleGun -> SetParticlePosition(gunposition);  
  fParticleGun -> SetParticleDefinition(particle);
  fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun -> SetParticleEnergy(170.*MeV);
}

INDRAPrimaryGeneratorAction::~INDRAPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void INDRAPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event

  fParticleGun -> GeneratePrimaryVertex(anEvent);
  std::cout << eventNum++ << G4endl;
}

