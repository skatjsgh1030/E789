#include "MCEventGenerator.hh"

using namespace std; 

MCEventGenerator::MCEventGenerator(G4String fileName) 
{
  fInputFile.open(fileName.data());
  fInputFile >> fNumEvents;
}
 
MCEventGenerator::~MCEventGenerator() 
{
  if(fInputFile.is_open()) fInputFile.close();
}
 
bool MCEventGenerator::ReadNextEvent(double &vx, double &vy, double &vz) 
{
  G4int eventID;
  if (!(fInputFile >> eventID >> fNumTracks >> vx >> vy >> vz))
    return false;
 
  fCurrentTrackID = 0;
  return true;
}
 
bool MCEventGenerator::ReadNextTrack(int &pdg, double &px, double &py, double &pz) 
{
  if (fCurrentTrackID >= fNumTracks)
    return false;
 
  fInputFile >> pdg >> px >> py >> pz;
  fCurrentTrackID++;
 
  return true;
}
