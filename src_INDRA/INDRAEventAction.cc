#include "INDRAEventAction.hh"
#include "INDRASteppingAction.hh"
#include "G4RunManager.hh"
#include "TMath.h"
#include "Randomize.hh"

using namespace CLHEP;
using namespace std;

  INDRAEventAction::INDRAEventAction()
: G4UserEventAction()
{
}

INDRAEventAction::~INDRAEventAction()
{
}

void INDRAEventAction::BeginOfEventAction(const G4Event*)
{
  eventID = 0;

  for(G4int volI = 0 ; volI < 375 ; volI++)
  {
    INDRAedepTot[volI] = 0.;
    INDRAedepTotModi[volI] = 0.;
  }
}

void INDRAEventAction::FindPreStepKinetic(G4int volume, G4double kineticE)
{
  for(G4int volI = 0 ; volI < 375 ; volI++)
  {
    if(volI == volume)
    {  
      INDRAPrekineticE[volI] = kineticE;
    }
  }
}

void INDRAEventAction::AddEnergyDeposit(G4int volume, G4double edep)
{
  for(G4int volI = 0 ; volI < 375 ; volI++)
  {
    if(volI == volume)
    {  
      INDRAedepTot[volI] += edep/2;
      break;
    }
  }
}

void INDRAEventAction::EndOfEventAction(const G4Event*)
{
  G4double std_INDRASinoi = 0.05;
  G4double std_INDRACsInoi = 0.05;  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID(); 

  analysisManager -> FillNtupleDColumn(1, 0, eventID);

  for(G4int volI = 0 ; volI < 375 ; volI++)
  {
    if(volI >= 0 && volI < 275)
    {   
      INDRAedepTotModi[volI] = RandGauss::shoot(INDRAedepTot[volI],TMath::Sqrt((0.115*INDRAPrekineticE[volI]*0.0000036)+(std_INDRACsInoi*std_INDRACsInoi)));
      analysisManager -> FillNtupleDColumn(1, 1+volI, INDRAedepTotModi[volI]);
    }   
    if(volI >= 275 && volI < 375)
    {   
      INDRAedepTotModi[volI] = RandGauss::shoot(INDRAedepTot[volI],TMath::Sqrt((0.115*INDRAPrekineticE[volI]*0.0000036)+(std_INDRASinoi*std_INDRASinoi)));
      analysisManager -> FillNtupleDColumn(1, 1+volI, INDRAedepTotModi[volI]);
    }   
  }

  analysisManager -> AddNtupleRow(1);
}
