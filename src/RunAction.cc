#include "G4UIcommand.hh"
#include "RunAction.hh"
#include "g4root.hh"
//#include "g4cvs.hh"

  RunAction::RunAction()
: G4UserRunAction()
{
}

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> OpenFile("E789_Simdata.root");

  analysisManager -> CreateNtuple("step", "step");
  analysisManager -> CreateNtupleIColumn("eventID");
  analysisManager -> CreateNtupleIColumn("volumeID");
  analysisManager -> CreateNtupleDColumn("edep");
  analysisManager -> CreateNtupleDColumn("x");
  analysisManager -> CreateNtupleDColumn("y");
  analysisManager -> CreateNtupleDColumn("z");
  analysisManager -> FinishNtuple();

  analysisManager -> CreateNtuple("edep","edep");
  analysisManager -> CreateNtupleDColumn("eventID");
  analysisManager -> CreateNtupleDColumn("VolumeID");
  analysisManager -> CreateNtupleDColumn("PID");
  analysisManager -> CreateNtupleDColumn("FAZIA_Mult");
  analysisManager -> CreateNtupleDColumn("INDRA_Mult");
  //analysisManager -> CreateNtupleDColumn("edepTotModi[volumeID]");

  //energy deposit 
  G4int volID = 960;
  for(G4int volI = 0 ; volI < volID ; volI++)
  {
    G4String volIDInString = G4UIcommand::ConvertToString(volI);
    analysisManager -> CreateNtupleDColumn("edepTotModi"+volIDInString);
  }
  analysisManager -> FinishNtuple();

}

void RunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> Write();
  analysisManager -> CloseFile();
}
