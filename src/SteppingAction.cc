#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#include "G4ThreeVector.hh"

#include "G4EventManager.hh"
//#include "G4cvs.hh"
//#include "G4root.hh"
  SteppingAction::SteppingAction()
: G4UserSteppingAction()
{
}

SteppingAction::~SteppingAction()
{
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
  G4int volumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();
  G4int kineticE = step -> GetPreStepPoint() -> GetKineticEnergy();
  G4int partiID = step -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
  //G4int trackID = step -> GetTrack() -> GetTrackID();
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
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=  
  // G4cout << "vol" << volumeID << G4endl;
  EventAction *eventAction;
  eventAction = (EventAction *) G4EventManager::GetEventManager() -> GetUserEventAction();
  if(step -> GetPostStepPoint() -> GetPhysicalVolume() != nullptr && step -> GetPreStepPoint() -> GetPhysicalVolume() != nullptr )// for avoid out of world information
  {
    G4int prevolumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();
    G4int postvolumeID = step -> GetPostStepPoint() -> GetPhysicalVolume() -> GetCopyNo();
    //-------------------------------------------------------------------------------------
    //  G4cout << "pre" << prevolumeID << G4endl;
    //if(volumeID == volI > 1  && volumeID == volI >= 592) //FAZIA Multi
    //G4cout << "pre if Mult :" << f_mult << G4endl;
    //G4cout << "Track ID : " << trackID << G4endl;
    if(partiID != 22 && partiID != 2112)
    {
      if(prevolumeID == 0 && 1 <= postvolumeID && postvolumeID <= 192 ) //FAZIA Multi
      {
        //G4cout << "FAZIA pre ++ Mult :" << f_mult << G4endl;
        f_mult = 1;
        //G4cout << "FAZIA post++ Mult :" << f_mult << G4endl;
        //G4cout << "FAZIA Vacuum :" << prevolumeID << G4endl;
        //G4cout << "FAZIA Tele :" << postvolumeID << G4endl;
        eventAction -> FindFAZIAMult(partiID,f_mult);
      }    
      if(prevolumeID == 0 && 701 <= postvolumeID && postvolumeID < 960) //INDRA Multi
      {
        //G4cout << "INDRA pre ++ Mult :" << i_mult << G4endl;
        i_mult = 1;
        //G4cout << "INDRA post ++ Mult :" << i_mult << G4endl;
        //G4cout << "INDRA Vacuum :" << prevolumeID << G4endl;
        //G4cout << "INDRA Tele :" << postvolumeID << G4endl;
        eventAction -> FindINDRAMult(partiID,i_mult);
      }
      //------------------------------------------------------------------------------------
      if(prevolumeID == 0 && 1<= postvolumeID && postvolumeID < 960)
      {
        //G4cout<< "volumestep :"<<postvolumeID <<G4endl;
        volume = postvolumeID;
        eventAction -> FindPIDandVolumeID(volume, partiID);
      }
      //------------------------------------------------------------------------------------
    }
    else if(partiID == 22 && partiID == 2112)
    {
      f_mult = 0;
      i_mult = 0;
      eventAction -> FindFAZIAMult(partiID, f_mult);
      eventAction -> FindINDRAMult(partiID, i_mult);
    }
    //-------------------------------------------------------------------------------------
    for(int volI = 0 ; volI < 960 ; volI++)
    {
      if (volumeID == volI)
      {
        eventAction -> AddEnergyDeposit(volI, edep);
      }
      if(prevolumeID == 0 && postvolumeID == volI)
      {
        eventAction -> FindPreStepKinetic(volI, kineticE);
      }
    }
    i_mult = 0; 
    f_mult = 0; 
  }
}
