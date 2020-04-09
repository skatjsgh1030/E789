#include "FAZIAEventAction.hh"
#include "FAZIASteppingAction.hh"
#include "G4RunManager.hh"
#include "TMath.h"
#include "Randomize.hh"

using namespace CLHEP;
using namespace std;

  FAZIAEventAction::FAZIAEventAction()
: G4UserEventAction()
{
}

FAZIAEventAction::~FAZIAEventAction()
{
}

void FAZIAEventAction::BeginOfEventAction(const G4Event*)
{
  eventID = 0;
  //volumeID = 0;

  for(G4int volI =0 ; volI < 593 ; volI++)
  {
    edepTot[volI] = 0.;
    edepTotModi[volI] = 0.;
  }
}

void FAZIAEventAction::FindPreStepKinetic(G4int volume, G4double kineticE)
{
  for(G4int volI = 0 ; volI < 593 ; volI++)
  {	
    if(volI == volume)
    {	
      PrekineticE[volI] = kineticE;
    }
  }	
}

void FAZIAEventAction::AddEnergyDeposit(G4int volume, G4double edep)
{
  for(G4int volI = 0 ; volI < 593 ; volI++)
  {
    if(volI == volume)
    {
      edepTot[volI] += edep/2;
      break;
    }
  }
}

void FAZIAEventAction::EndOfEventAction(const G4Event*)
{
  G4double std_Sinoi = 0.05;
  G4double std_CsInoi = 0.1;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();

  analysisManager -> FillNtupleDColumn(1, 0, eventID);

  for(G4int volI = 0 ; volI < 593 ; volI++)
  {
    if(volI >= 0 && volI <= 401)
    {
      edepTotModi[volI] = RandGauss::shoot(edepTot[volI],TMath::Sqrt((0.115*PrekineticE[volI]*0.0000036)+(std_Sinoi*std_Sinoi)));
      analysisManager -> FillNtupleDColumn(1, 1+volI, edepTotModi[volI]);
    }
    if(volI >= 402 && volI < 593)
    {
      edepTotModi[volI] = RandGauss::shoot(edepTot[volI],TMath::Sqrt((0.115*PrekineticE[volI]*0.0000036)+(std_CsInoi*std_CsInoi)));
      analysisManager -> FillNtupleDColumn(1, 1+volI, edepTotModi[volI]);
    }
  }

  analysisManager -> AddNtupleRow(1);
}
