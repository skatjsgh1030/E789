#include "G4UIcommand.hh"
#include "INDRARunAction.hh"
#include "g4root.hh"
//#include "g4cvs.hh"

INDRARunAction::INDRARunAction()
: G4UserRunAction()
{
}

INDRARunAction::~INDRARunAction()
{
  delete G4AnalysisManager::Instance();
}

void INDRARunAction::BeginOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> OpenFile("INDRASimdata.root");

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
  
  G4int volID = 375;
  for(G4int volI = 0 ; volI < volID ; volI++)
  {
    G4String volIDInString = G4UIcommand::ConvertToString(volI);
    analysisManager -> CreateNtupleDColumn("INDRAedepTotModi"+volIDInString);
  }  
  analysisManager -> FinishNtuple();  

}

void INDRARunAction::EndOfRunAction(const G4Run*)
{
  std::cout << "end and write" << std::endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> Write();
  analysisManager -> CloseFile();
}
