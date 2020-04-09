#include "FAZIAPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "TMath.h"
#include <iostream>
//#define restM 105 //MeV

using namespace CLHEP;
using namespace std;

  FAZIAPrimaryGeneratorAction::FAZIAPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
}

FAZIAPrimaryGeneratorAction::~FAZIAPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void FAZIAPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  costheta = RandGauss::shoot(0,0.5); //mean , std 
  random = G4UniformRand();  
  phirand = G4UniformRand();
  Parti_ene = G4UniformRand();
  //cout <<"rand " <<random << endl;
  //cout <<"partene " <<Parti_ene << endl;
  //cout << "" << endl;
  ///////////////////////////////////////////////////////////////////////////
  G4int n_particle = 1;

  G4double Phi_rand = 2*pi*phirand;
  //G4double Theta_rand = (60*CLHEP::degree)*((TMath::Cos(90*CLHEP::degree*costheta))*(TMath::Cos(90*CLHEP::degree*costheta)) );
  G4double Theta_rand = (pi/2)*costheta;
  G4int r_value = 20; 
  G4double AMeV = 52; 
  G4double x = r_value*(TMath::Sin(Theta_rand))*(TMath::Cos(Phi_rand))*m;
  G4double y = r_value*(TMath::Sin(Theta_rand))*(TMath::Sin(Phi_rand))*m;
  G4double z = r_value*(TMath::Cos(Theta_rand))*m;

  fParticleGun = new G4ParticleGun(n_particle);
  G4ThreeVector gunposition = G4ThreeVector(0,0,0);

    if(0. <= random && random < 0.1)//select particle sort H
    {   
        G4double H1E = 0;
        G4int Z = 1;    
        G4int A = 1;    
         G4double mean = A*AMeV;
         G4double std = mean*0.3;
        H1E = RandGauss::shoot(mean,std); //mean , std 
    
        fParticleGun -> SetParticlePosition(gunposition);  
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(H1E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"1H"<< G4endl;
        std::cout <<"Energy" << H1E << G4endl;
    }   

    if(0.1 <= random && random < 0.2 )//d
    {   
        G4double H2E = 0;
        G4int Z = 1;    
        G4int A = 2;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        H2E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(H2E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"2H"<< G4endl;
        std::cout <<"Energy" << H2E << G4endl;
    }   

    if(0.2  <= random && random < 0.3)//t, 3He
    {   
        G4double H3E = 0;
        G4int Z = 1;    
        G4int A = 3;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        H3E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(H3E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"3H" <<  G4endl;
        std::cout <<"Energy" << H3E << G4endl;
    }   

    if(0.3 <= random && random <0.47 )//4He
    {   
        G4double He4E = 0;
        G4int Z = 2;    
        G4int A = 4;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        He4E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(He4E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"4He"<< G4endl;
        std::cout <<"Energy" << He4E << G4endl;
    }   
    if(0.47 <= random && random <0.50 )//6Li
    {   
        G4double Li6E = 0;
        G4int Z = 3;    
        G4int A = 6;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Li6E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Li6E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"6Li"<< G4endl;
        std::cout <<"Energy" << Li6E << G4endl;
    }   
    if(0.50 <= random && random <0.60 )//7Li
    {   
        G4double Li7E = 0;
        G4int Z = 3;    
        G4int A = 7;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Li7E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Li7E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"7Li"<< G4endl;
        std::cout <<"Energy" << Li7E << G4endl;
    }   
    if(0.60 <= random && random <0.74 )//8Li
    {   
        G4double Li8E = 0;
        G4int Z = 3;    
        G4int A = 8;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Li8E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Li8E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"8Li"<< G4endl;
        std::cout <<"Energy" << Li8E << G4endl;
    }   
    if(0.74 <= random && random <0.75 )//9Li
    {   
        G4double Li9E = 0;
        G4int Z = 3;    
        G4int A = 9;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Li9E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Li9E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"9Li"<< G4endl;
        std::cout <<"Energy" << Li9E << G4endl;
    }   
    if(0.75 <= random && random <0.80 )//7Be
    {   
        G4double Be7E = 0;
        G4int Z = 4;    
        G4int A = 7;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Be7E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Be7E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"7Be"<< G4endl;
        std::cout <<"Energy" << Be7E << G4endl;
    }   
    if(0.80 <= random && random <0.86 )//8Be
    {   
        G4double Be8E = 0;
        G4int Z = 4;    
        G4int A = 8;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Be8E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Be8E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"8Be"<< G4endl;
        std::cout <<"Energy" << Be8E << G4endl;
    }   
    if(0.86 <= random && random < 0.92)//9Be
    {   
        G4double Be9E = 0;
        G4int Z = 4;    
        G4int A = 9;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Be9E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Be9E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"9Be"  << G4endl;
        std::cout <<"Energy" << Be9E << G4endl;
    }   

    if(0.92 <= random && random <= 1.0 )//10Be
    {   
        G4double Be10E = 0;
        G4int Z = 4;    
        G4int A = 10;    
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Be10E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Be10E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"10Be"<< G4endl;
        std::cout <<"Energy" << Be10E << G4endl;
    }   

/*  if(0.6  <= random && random < 0.65)
    {
        G4double Ni64E = 0;
        Ni64E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(28,64,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Ni64E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"64Ni" <<  G4endl;
        std::cout <<"Energy" << Ni64E << G4endl;
    }

    if(0.6 <= random && random <0.7 )
    {
        G4double Ni58E = 0;
        G4int Z = 28;       
        G4int A = 58;       
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Ni58E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Ni58E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"58Ni"<< G4endl;
        std::cout <<"Energy" << Ni58E << G4endl;
    }
    if(0.7 <= random && random <0.75 )
    {
        G4double C11E = 0;
        G4int Z = 6;        
        G4int A = 11;       
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        C11E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(C11E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"11C"<< G4endl;
        std::cout <<"Energy" << C11E << G4endl;
    }
    if(0.75 <= random && random <0.8 )
    {
        G4double C12E = 0;
        G4int Z = 6;        
        G4int A = 12;       
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        C12E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(C12E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"12C"<< G4endl;
        std::cout <<"Energy" << C12E << G4endl;
    }
    if(0.8 <= random && random <0.85 )
    {
        G4double C14E = 0;
        G4int Z = 6;        
        G4int A = 14;       
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        C14E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(C14E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"14C"<< G4endl;
        std::cout <<"Energy" << C14E << G4endl;
    }
    if(0.85 <= random && random <0.9 )
    {
        G4double Mg24E = 0;
        G4int Z = 12;       
        G4int A = 24;       
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Mg24E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Mg24E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"24Mg"<< G4endl;
        std::cout <<"Energy" << Mg24E << G4endl;
    }
    if(0.9 <= random && random <0.95 )
    {
        G4double Mg26E = 0;
        G4int Z = 12;       
        G4int A = 26;       
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        Mg26E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(Mg26E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"26Mg"<< G4endl;
        std::cout <<"Energy" << Mg26E << G4endl;
    }
    if(0.95 <= random && random <=1.0 )
    {
        G4double S33E = 0;
        G4int Z = 16;       
        G4int A = 33;       
        G4double mean = A*AMeV;
         G4double std = mean*0.3;
        S33E = RandGauss::shoot(mean,std); //mean , std 

        fParticleGun -> SetParticlePosition(gunposition);  
        //fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
        fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
        fParticleGun -> SetParticleEnergy(S33E*MeV);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        std::cout <<"event" << eventNum++ << G4endl;
        std::cout <<"33S"<< G4endl;
        std::cout <<"Energy" << S33E << G4endl;
    }*/
  //this function is called at the begining of each event

  fParticleGun -> GeneratePrimaryVertex(anEvent);
  //std::cout << eventNum++ << G4endl;
}                                                                                                                                                                                                                                                                                                                                                                                                 
