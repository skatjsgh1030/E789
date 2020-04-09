#include "PrimaryGeneratorAction.hh"
#include "MCEventGenerator.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

  PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
  fParticleGun = new G4ParticleGun();
  fEventGenerator = new MCEventGenerator("../../Si_Si_CsI/IQMD_data/IQMD_Soft250_1000.gen");    
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event
  G4int pdg;
  G4double vx, vy, vz, px, py, pz;

  fEventGenerator -> ReadNextEvent(vx, vy, vz);
  fParticleGun -> SetParticlePosition(G4ThreeVector(vx,vy,vz));

  while (fEventGenerator -> ReadNextTrack(pdg, px, py, pz)) {
	G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable() -> FindParticle(pdg);
	fParticleGun -> SetParticleDefinition(particle);

	G4ThreeVector momentum(px,py,pz);
	fParticleGun -> SetParticleMomentum(momentum.mag()*MeV);
	fParticleGun -> SetParticleMomentumDirection(momentum.unit());
	fParticleGun -> GeneratePrimaryVertex(anEvent);
  }
}
