#include "EventAction.hh"
#include "SteppingAction.hh"
#include "G4RunManager.hh"
#include "TMath.h"
#include "Randomize.hh"

using namespace CLHEP;
using namespace std;

  EventAction::EventAction()
: G4UserEventAction()
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event*)
{
  //eventID = 0;
  //G4int partiID = 0;
  eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
  G4cout << "event start : " <<eventID << G4endl;
  for(G4int volI =0 ; volI < 960 ; volI++){
    edepTot[volI] = 0.;
    edepTotModi[volI] = 0.;}
}

void EventAction::FindPreStepKinetic(G4int volume, G4double kineticE)
{
  for(G4int volI = 0 ; volI < 960 ; volI++)
  {	
    if(volI == volume){	
      PrekineticE[volI] = kineticE;}
  }
}

void EventAction::AddEnergyDeposit(G4int volume, G4double edep){
  for(G4int volI = 0 ; volI < 960 ; volI++){
    if(volI == volume){
      edepTot[volI] += edep;
      break;}}
}
//----------------------------------------------------------------------------------
void EventAction::FindFAZIAMult(G4int partiID, G4int multp){
  if(multp !=0 && partiID != 22 && partiID != 2112)
  {
    //multp++;
    FAZIA_Mult += multp;
    //cout<<" FAZIA particle ="<< partiID <<endl;
    //cout<<" eventMult of FAZIA "<< FAZIA_Mult <<endl;
  }
}
void EventAction::FindINDRAMult(G4int partiID, G4int mult){
  if(mult !=0 && partiID != 22 && partiID != 2112)
  {
    //mult++;
    INDRA_Mult += mult;
    //cout<<"  INDRA particle ="<< partiID <<endl;
    //cout<<" eventMult of INDRA "<< INDRA_Mult <<endl;
  }
}
//-----------------------------------------------------------------------------------
void EventAction::FindPIDandVolumeID(G4int volume, G4int partiID){

  VolumeID = volume;
  PID = partiID;
  //cout << "volumeID : "<<volume <<endl;
  //cout << "PID : "<<partiID <<endl;
}  
//-----------------------------------------------------------------------------------
void EventAction::EndOfEventAction(const G4Event*)
{
  G4double std_Sinoi = 0.05;
  G4double std_CsInoi = 0.1;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();

  //FAZIA_Mult = multp;
  //INDRA_Mult = mult;

  analysisManager -> FillNtupleDColumn(1, 0, eventID);
  analysisManager -> FillNtupleDColumn(1, 1, VolumeID);
  analysisManager -> FillNtupleDColumn(1, 2, PID);
  analysisManager -> FillNtupleDColumn(1, 3, FAZIA_Mult);
  analysisManager -> FillNtupleDColumn(1, 4, INDRA_Mult);

  // Manage Energy Deposite in Crystals
  for(G4int volI = 0 ; volI < 960 ; volI++){    
    if(volI >= 0 && volI <= 200){//FAZIA Si1 response
      edepTotModi[volI] = RandGauss::shoot(edepTot[volI],TMath::Sqrt((0.115*PrekineticE[volI]*0.0000036)+(std_Sinoi*std_Sinoi)));
      analysisManager -> FillNtupleDColumn(1, 5+volI, edepTotModi[volI]);
    }
    if(volI >= 201 && volI <= 400){//FAZIA Si2 response
      edepTotModi[volI] = RandGauss::shoot(edepTot[volI],TMath::Sqrt((0.115*PrekineticE[volI]*0.0000036)+(std_Sinoi*std_Sinoi)));
      analysisManager -> FillNtupleDColumn(1, 5+volI, edepTotModi[volI]);
    }
    if(volI >= 401 && volI <= 600){//FAZIA CsI response
      edepTotModi[volI] = RandGauss::shoot(edepTot[volI],TMath::Sqrt((0.115*PrekineticE[volI]*0.0000036)+(std_CsInoi*std_CsInoi)));
      analysisManager -> FillNtupleDColumn(1, 5+volI, edepTotModi[volI]);
    }
    if(601 <= volI && volI <= 860){//INDRA CsI response
      edepTotModi[volI] = RandGauss::shoot(edepTot[volI],TMath::Sqrt((0.115*PrekineticE[volI]*0.0000036)+(std_Sinoi*std_Sinoi)));
      analysisManager -> FillNtupleDColumn(1, 5+volI, edepTotModi[volI]);
    }
    if(861 <= volI && volI < 960){//INDRA Si response
      edepTotModi[volI] = RandGauss::shoot(edepTot[volI],TMath::Sqrt((0.115*PrekineticE[volI]*0.0000036)+(std_CsInoi*std_CsInoi)));
      analysisManager -> FillNtupleDColumn(1, 5+volI, edepTotModi[volI]);
    }
  }
  G4cout << "event end : " <<eventID << G4endl;
  analysisManager -> AddNtupleRow(1);
  FAZIA_Mult = 0;
  INDRA_Mult = 0;
  VolumeID = 0;
  PID = 0;
  //multp = 0;
  //mult = 0;
}
