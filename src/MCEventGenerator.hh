#ifndef MCEVENTGENGENERATOR_HH
#define MCEVENTGENGENERATOR_HH
 
#include <fstream>
#include "globals.hh"

using namespace std;

class MCEventGenerator {
  public:
    MCEventGenerator(G4String fileName);
    virtual ~MCEventGenerator();
 
    bool ReadNextEvent(G4double &vx, G4double &vy, G4double &vz);
    bool ReadNextTrack(G4int &pdg, G4double &px, G4double &py, G4double &pz);
 
    G4int GetNumEvents() { return fNumEvents; };
 
  private:
    std::ifstream fInputFile;
    G4int fNumEvents;
    G4int fNumTracks;
    G4int fCurrentTrackID;
};
 
#endif
