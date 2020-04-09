#include "INDRASteppingAction.hh"
#include "INDRAEventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
//#include "G4RunManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#include "G4ThreeVector.hh"

#include "G4EventManager.hh"
//#include "G4cvs.hh"
//#include "G4root.hh"

  INDRASteppingAction::INDRASteppingAction()
: G4UserSteppingAction()
{
  //fScintillationCounter = 0;
  //fCerenkovCounter      = 0;
  //fEventNumber = 0;
}

INDRASteppingAction::~INDRASteppingAction()
{
}

void INDRASteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
  G4int volumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();
  G4int kineticE = step -> GetPreStepPoint() -> GetKineticEnergy();
  G4double edep = step -> GetTotalEnergyDeposit();
  G4double x = step -> GetDeltaPosition().x();  	
  G4double y = step -> GetDeltaPosition().y();  	
  G4double z = step -> GetDeltaPosition().z();

  //G4double moment = step -> GetDeltaMomentum();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> FillNtupleIColumn(0, eventID);
  analysisManager -> FillNtupleIColumn(1, volumeID);
  analysisManager -> FillNtupleDColumn(2, edep);
  analysisManager -> FillNtupleDColumn(3, x);
  analysisManager -> FillNtupleDColumn(4, y);
  analysisManager -> FillNtupleDColumn(5, z);
  //analysisManager -> FillNtupleDColumn(6, moment);
  analysisManager -> AddNtupleRow();

  INDRAEventAction *eventAction; 
  eventAction = (INDRAEventAction *) G4EventManager::GetEventManager() -> GetUserEventAction();
  if(step -> GetPostStepPoint() -> GetPhysicalVolume() != nullptr && step -> GetPreStepPoint() -> GetPhysicalVolume() != nullptr)
  {
    G4int prevolumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();
    G4int postvolumeID = step -> GetPostStepPoint() -> GetPhysicalVolume() -> GetCopyNo();

    for(int volI = 0 ; volI < 375 ; volI++)
    { 
      if(volumeID == volI)
      {
        eventAction -> AddEnergyDeposit(volI, edep);
      }
      if(prevolumeID == 0 && postvolumeID == volI) 
      {
        eventAction -> FindPreStepKinetic(volI, kineticE);
      }
    }
  }
}
