#include "FAZIA.hh"
#include "DetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iostream>
#include "TMath.h"

#define csibackL 20.0
#define DistanceSi1 1000 //1000mm
#define DistanceSi2 1003.620 //1000.625mm
#define DistanceCsI 1061.5 //1051.058mm
#define halfL 11.35

using std::stringstream;
using namespace std;
using namespace CLHEP;

FAZIA::FAZIA()
{
}

FAZIA::~FAZIA()
{
}

void FAZIA::faziaGeo(G4LogicalVolume* logicWorld)
{
  G4int overlap =0;
  G4int RoDi = 1;
  G4int RoDiQ = 1;
  G4double gapTel = 1.0;//gap of each telescope 
  G4double gapBlk = 3.5;//gap of each Block
  G4double centX = 24.00;//22.477651;
  G4double centY = 26.00;//24.439468;

  //COR -> position of starting Core block(block number 0's 1,2,3,4 telescope)
  G4double C1x = halfL*3 + gapTel*(3/2)+centX; 
  G4double C1y = halfL*7+ gapTel*3 +centY;
  G4double C2x = halfL*1 + gapTel*(1/2)+centX;
  G4double C2y = halfL*7+ gapTel*3 +centY;
  G4double C3x = halfL*3 + gapTel*(3/2)+centX;
  G4double C3y = halfL*5 + gapTel*2 +centY;
  G4double C4x = halfL*1 + gapTel*(1/2)+centX;
  G4double C4y = halfL*5 + gapTel*2 +centY;

  G4double Blktheta1 = TMath::ATan( TMath::Sqrt(C1x*C1x + C1y*C1y)/1000.0 );  //0.10728646
  G4double Blktheta2 = TMath::ATan( TMath::Sqrt(C2x*C2x + C2y*C2y)/1000.0 ); //0.10162866
  G4double Blktheta3 = TMath::ATan( TMath::Sqrt(C3x*C3x + C3y*C3y)/1000.0 );  //0.08920
  G4double Blktheta4 = TMath::ATan( TMath::Sqrt(C4x*C4x + C4y*C4y)/1000.0 );  //0.0822759
  G4double Blkphi1 = TMath::ATan(C1y/C1x);  //1.1902899
  G4double Blkphi2 = TMath::ATan(C2y/C2x);  //1.3734008
  G4double Blkphi3 = TMath::ATan(C3y/C3x);  //1.1071487
  G4double Blkphi4 = TMath::ATan(C4y/C4x);  //1.3258177

  //MID -> position of starting Middle block(block number 4's 1,2,3,4 telescope)
  G4double M1x = halfL*11 + gapTel*(9/2) + gapBlk+centX; 
  G4double M1y = halfL*7 + gapTel*3 +centY;
  G4double M2x = halfL*9 + gapTel*(7/2) + gapBlk+centX;
  G4double M2y = halfL*7 + gapTel*3 +centY;
  G4double M3x = halfL*11 + gapTel*(9/2) + gapBlk+centX;
  G4double M3y = halfL*5 + gapTel*2 +centY;
  G4double M4x = halfL*9 + gapTel*(7/2) + gapBlk+centX;
  G4double M4y = halfL*5 + gapTel*2 +centY;

  G4double Blktheta41 = TMath::ATan( TMath::Sqrt(M1x*M1x + M1y*M1y)/1000.0 );  //0.15495
  G4double Blktheta42 = TMath::ATan( TMath::Sqrt(M2x*M2x + M2y*M2y)/1000.0 );  //0.140489
  G4double Blktheta43 = TMath::ATan( TMath::Sqrt(M3x*M3x + M3y*M3y)/1000.0 );  //0.143232
  G4double Blktheta44 = TMath::ATan( TMath::Sqrt(M4x*M4x + M4y*M4y)/1000.0 );  //0.1273692
  G4double Blkphi41 = TMath::ATan(M1y/M1x);  //0.6947382
  G4double Blkphi42 = TMath::ATan(M2y/M2x);  //0.78539816
  G4double Blkphi43 = TMath::ATan(M3y/M3x);  //0.58800260
  G4double Blkphi44 = TMath::ATan(M4y/M4x);  //0.67474094

  //OUTER -> position of starting Outer block(block number 5's 1,2,3,4 telescope)
  G4double O1x = halfL*3 + gapTel*(3/2)+centX; 
  G4double O1y = halfL*15 + gapTel*6 + gapBlk+centY; 
  G4double O2x = halfL*1 + gapTel*(1/2)+centX; 
  G4double O2y = halfL*15 + gapTel*6 + gapBlk+centY; 
  G4double O3x = halfL*3 + gapTel*(3/2)+centX; 
  G4double O3y = halfL*13 + gapTel*5 + gapBlk+centY; 
  G4double O4x = halfL*1 + gapTel*(1/2)+centX; 
  G4double O4y = halfL*13 + gapTel*5 + gapBlk+centY; 

  G4double Blktheta51 = TMath::ATan( TMath::Sqrt(O1x*O1x + O1y*O1y)/1000.0 );  //0.18234275
  G4double Blktheta52 = TMath::ATan( TMath::Sqrt(O2x*O2x + O2y*O2y)/1000.0 );  //0.17916567
  G4double Blktheta53 = TMath::ATan( TMath::Sqrt(O3x*O3x + O3y*O3y)/1000.0 );  //0.16345286
  G4double Blktheta54 = TMath::ATan( TMath::Sqrt(O4x*O4x + O4y*O4y)/1000.0 );  //0.15986910
  G4double Blkphi51 = TMath::ATan(O1y/O1x);  //1.3521274
  G4double Blkphi52 = TMath::ATan(O2y/O2x);  //1.4601391
  G4double Blkphi53 = TMath::ATan(O3y/O3x);  //1.3258177
  G4double Blkphi54 = TMath::ATan(O4y/O4x);  //1.4464413

  G4double CORx = gapTel*0.5+centX;
  G4double CORy = halfL*4 + gapTel*1.5 +centY;

  G4double CORtheta = TMath::ATan( TMath::Sqrt(CORx*CORx + CORy*CORy)/1000.0 ); 
  G4double CORphi = TMath::ATan(CORy/CORx); 

  G4double MIDx = halfL*8 + gapTel*3 + gapBlk+centX;
  G4double MIDy = halfL*4 + gapTel*1.5 +centY;

  G4double MIDtheta = TMath::ATan( TMath::Sqrt(MIDx*MIDx + MIDy*MIDy)/1000.0 ); 
  G4double MIDphi = TMath::ATan(MIDy/MIDx); 

  G4double OUTERx = gapTel*0.5+centX;
  G4double OUTERy = halfL*12 + gapTel*4.5 + gapBlk + centY;

  G4double OUTERtheta = TMath::ATan( TMath::Sqrt(OUTERx*OUTERx + OUTERy*OUTERy)/1000.0 ); 
  G4double OUTERphi = TMath::ATan(OUTERy/OUTERx); 
  //----------------------------------------------------------------------------------------------
  //----------------Material --------------------
  G4NistManager* nist = G4NistManager::Instance();
  const G4double labTemp = CLHEP::STP_Temperature + 20.*kelvin;

  G4Element*  elCs  = new G4Element("Caesium","Cs",  55,   132.90545*g/mole);
  G4Element*  elI   = new G4Element("Iodine", "I" , 53,   126.90447*g/mole);

  G4Material* CsI = new G4Material("CsI", 4.51*g/CLHEP::cm3, 2, kStateSolid, labTemp);
  CsI -> AddElement(elCs, 1);
  CsI -> AddElement(elI,  1);

  /*/ ------------ Generate & Add Material Properties Table ------------
  //
  G4double photonEnergy[] =
  {	1.90738*eV, 1.98368*eV, 2.06633*eV,
  2.15617*eV, 2.25418*eV, 2.36152*eV,
  2.47960*eV, 2.61011*eV, 2.75511*eV,
  2.91718*eV, 3.0995*eV};
  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  // CsI
  //fixed
  G4double refractiveIndex1[] =
  { 1.7789, 1.7820, 1.7856, 1.7897, 1.7945,
  1.8001, 1.8067, 1.8146, 1.8242, 1.8359,
  1.8506};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
  {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
  30.000*m, 28.500*m, 27.000*m,17.500*m, 14.500*m };

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
  { 1.00, 1.00, 1.00, 
  1.00, 1.00, 1.00, 
  1.00,1.00, 1.00,
  1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
  { 0.01, 1.00, 2.00,
  3.00, 4.00, 5.00,
  6.00,7.00, 6.00,
  5.00, 4.00 };

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
  ->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
  ->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
  ->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
  ->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  CsI->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator

  CsI->GetIonisation()->SetBirksConstant(0.126*mm/MeV);*/

  // Detector

  G4Material* scintillator_mat = nist -> FindOrBuildMaterial("CsI");
  G4Material* Silicon_mat = nist -> FindOrBuildMaterial("G4_Si");
  //G4Material* randomizer_mat = nist -> FindOrBuildMaterial("G4_Al");
  const int blkcor = 4;
  const int tel = 4;

  //////////////COR////////////////////////////////////////////////////////////////////detecor kind -> BlockNum -> Quad -> Tele
  G4int silicon_11Num[blkcor][tel] = {{1,6,11,16},// BLK 0
	{22,27,32,17},//+16 // add +16 and push a array elements // BLK 1
	{43,48,33,38},//+32 // add +32 and push two array elements
	{64,49,54,59}};//+48 // add +48 and push three array elements

  G4int silicon_12Num[blkcor][tel] = {{201,206,211,216},
	{222,227,232,217},
	{243,248,233,238},
	{264,249,254,259}};

  G4int csicrystal1Num[blkcor][tel] = {{401,406,411,416},
	{422,427,432,417},
	{443,448,433,438},
	{464,449,454,459}};
//------------------------------------------------------------------------------------
  G4int silicon_21Num[blkcor][tel] = {{2,7,12,13},//BLK 0
	{23,28,29,18},//+16// BLK 1
	{44,45,34,39},//+32
	{61,50,55,60}};//+48

  G4int silicon_22Num[blkcor][tel] = {{202,207,212,213},
	{223,228,229,218},
	{244,245,234,239},
	{261,250,255,260}};

  G4int csicrystal2Num[blkcor][tel] = {{402,407,412,413},
	{423,428,429,418},
	{444,445,434,439},
	{461,450,455,460}};
//------------------------------------------------------------------------------------
  G4int silicon_31Num[blkcor][tel] = {{4,5,10,15},//BLK 0
	{21,26,31,20},//+16// BLK 1
	{42,47,36,37},//+32
	{63,52,53,58}};//+48

  G4int silicon_32Num[blkcor][tel] = {{204,205,210,215},
	{221,226,231,220},
	{242,247,236,237},
	{263,252,253,258}};

  G4int csicrystal3Num[blkcor][tel] = {{404,405,410,415},
	{421,426,431,420},
	{442,447,436,437},
	{463,452,453,458}};
//------------------------------------------------------------------------------------
  G4int silicon_41Num[blkcor][tel] = {{3,8,9,14},//BLK 0
	{24,25,30,19},//+16// BLK 1
	{41,46,35,40},//+32
	{62,51,56,57}};//+48

  G4int silicon_42Num[blkcor][tel] = {{203,208,209,214},
	{224,225,230,219},
	{241,246,235,240},
	{262,251,256,257}};

  G4int csicrystal4Num[blkcor][tel] = {{403,408,409,414},
	{424,425,430,419},
	{441,446,435,440},
	{426,451,456,457}};
  //////////////////////////////////////////////////////////////////////////////
  string silicon_11[blkcor][tel] ={ {"S1B0Q1T1","S1B0Q2T2","S1B0Q3T3","S1B0Q4T4"},
	{"S1B1Q1T1","S1B1Q2T2","S1B1Q3T3","S1B1Q4T4"},
	{"S1B2Q1T1","S1B2Q2T2","S1B2Q3T3","S1B2Q4T4"},
	{"S1B3Q1T1","S1B3Q2T2","S1B3Q3T3","S1B3Q4T4"}};

  string silicon_12[blkcor][tel] ={ {"S2B0Q1T1","S2B0Q2T2","S2B0Q3T3","S2B0Q4T4"},
	{"S2B1Q1T1","S2B1Q2T2","S2B1Q3T3","S2B1Q4T4"},
	{"S2B2Q1T1","S2B2Q2T2","S2B2Q3T3","S2B2Q4T4"},
	{"S2B3Q1T1","S2B3Q2T2","S2B3Q3T3","S2B3Q4T4"}};

  string csicrystal1[blkcor][tel] ={ {"CsIB0Q1T1","CsIB0Q2T2","CsIB0Q3T3","CsIB0Q4T4"},
	{"CsIB1Q1T1","CsIB1Q2T2","CsIB1Q3T3","CsIB1Q4T4"},
	{"CsIB2Q1T1","CsIB2Q2T2","CsIB2Q3T3","CsIB2Q4T4"},
	{"CsIB3Q1T1","CsIB3Q2T2","CsIB3Q3T3","CsIB3Q4T4"}};

  string silicon_21[blkcor][tel] ={ {"S1B0Q1T2","S1B0Q2T3","S1B0Q3T4","S1B0Q4T1"},
	{"S1B1Q1T2","S1B1Q2T3","S1B1Q3T4","S1B1Q4T1"},
	{"S1B2Q1T2","S1B2Q2T3","S1B2Q3T4","S1B2Q4T1"},
	{"S1B3Q1T2","S1B3Q2T3","S1B3Q3T4","S1B3Q4T1"}};

  string silicon_22[blkcor][tel] ={ {"S2B0Q1T2","S2B0Q2T3","S2B0Q3T4","S2B0Q4T1"},
	{"S2B1Q1T2","S2B1Q2T3","S2B1Q3T4","S2B1Q4T1"},
	{"S2B2Q1T2","S2B2Q2T3","S2B2Q3T4","S2B2Q4T1"},
	{"S2B3Q1T2","S2B3Q2T3","S2B3Q3T4","S2B3Q4T1"}};

  string csicrystal2[blkcor][tel] = {{"CsIB0Q1T2","CsIB0Q2T3","CsIB0Q3T4","CsIB0Q4T1"},
	{"CsIB1Q1T2","CsIB1Q2T3","CsIB1Q3T4","CsIB1Q4T1"},
	{"CsIB2Q1T2","CsIB2Q2T3","CsIB2Q3T4","CsIB2Q4T1"},
	{"CsIB3Q1T2","CsIB3Q2T3","CsIB3Q3T4","CsIB3Q4T1"}};

  string silicon_31[blkcor][tel] ={ {"S1B0Q1T4","S1B0Q2T1","S1B0Q3T2","S1B0Q4T3"},
	{"S1B1Q1T4","S1B1Q2T1","S1B1Q3T2","S1B1Q4T3"},
	{"S1B2Q1T4","S1B2Q2T1","S1B2Q3T2","S1B2Q4T3"},
	{"S1B3Q1T4","S1B3Q2T1","S1B3Q3T2","S1B3Q4T3"}};

  string silicon_32[blkcor][tel] ={ {"S2B0Q1T4","S2B0Q2T1","S2B0Q3T2","S2B0Q4T3"},
	{"S2B1Q1T4","S2B1Q2T1","S2B1Q3T2","S2B1Q4T3"},
	{"S2B2Q1T4","S2B2Q2T1","S2B2Q3T2","S2B2Q4T3"},
	{"S2B3Q1T4","S2B3Q2T1","S2B3Q3T2","S2B3Q4T3"}};

  string csicrystal3[blkcor][tel] ={ {"CsIB0Q1T4","CsIB0Q2T1","CsIB0Q3T2","CsIB0Q4T3"},
	{"CsIB1Q1T4","CsIB1Q2T1","CsIB1Q3T2","CsIB1Q4T3"},
	{"CsIB2Q1T4","CsIB2Q2T1","CsIB2Q3T2","CsIB2Q4T3"},
	{"CsIB3Q1T4","CsIB3Q2T1","CsIB3Q3T2","CsIB3Q4T3"}};

  string silicon_41[blkcor][tel] ={ {"S1B0Q1T3","S1B0Q2T4","S1B0Q3T1","S1B0Q4T2"},
	{"S1B1Q1T3","S1B1Q2T4","S1B1Q3T1","S1B1Q4T2"},
	{"S1B2Q1T3","S1B2Q2T4","S1B2Q3T1","S1B2Q4T2"},
	{"S1B3Q1T3","S1B3Q2T4","S1B3Q3T1","S1B3Q4T2"}};

  string silicon_42[blkcor][tel] ={ {"S2B0Q1T3","S2B0Q2T4","S2B0Q3T1","S2B0Q4T2"},
	{"S2B1Q1T3","S2B1Q2T4","S2B1Q3T1","S2B1Q4T2"},
	{"S2B2Q1T3","S2B2Q2T4","S2B2Q3T1","S2B2Q4T2"},
	{"S2B3Q1T3","S2B3Q2T4","S2B3Q3T1","S2B3Q4T2"}};

  string csicrystal4[blkcor][tel] ={ {"CsIB0Q1T3","CsIB0Q2T4","CsIB0Q3T1","CsIB0Q4T2"},
	{"CsIB1Q1T3","CsIB1Q2T4","CsIB1Q3T1","CsIB1Q4T2"},
	{"CsIB2Q1T3","CsIB2Q2T4","CsIB2Q3T1","CsIB2Q4T2"},
	{"CsIB3Q1T3","CsIB3Q2T4","CsIB3Q3T1","CsIB3Q4T2"}};

  string TubeSurface[blkcor][tel] ={ {"CsI11","CsI12","CsI13","CsI14"},
	{"CsI111","CsI112","CsI113","CsI114"},
	{"CsI211","CsI212","CsI213","CsI214"},
	{"CsI311","CsI312","CsI313","CsI314"}};

  G4Box* Silicon_11[blkcor][tel];
  G4Box* Silicon_12[blkcor][tel];
  G4Trd* CsIcrystal1[blkcor][tel];

  G4Box* Silicon_21[blkcor][tel];
  G4Box* Silicon_22[blkcor][tel];
  G4Trd* CsIcrystal2[blkcor][tel];

  G4Box* Silicon_31[blkcor][tel];
  G4Box* Silicon_32[blkcor][tel];
  G4Trd* CsIcrystal3[blkcor][tel];

  G4Box* Silicon_41[blkcor][tel];
  G4Box* Silicon_42[blkcor][tel];
  G4Trd* CsIcrystal4[blkcor][tel];
  //G4Trd* CsIcrystal[1][qua][blkcor][tel];

  G4LogicalVolume* logicalDetector11[blkcor][tel];
  G4LogicalVolume* logicalDetector12[blkcor][tel];
  G4LogicalVolume* logicalDetector13[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector11[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector12[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector13[blkcor][tel];

  G4LogicalVolume* logicalDetector21[blkcor][tel];
  G4LogicalVolume* logicalDetector22[blkcor][tel];
  G4LogicalVolume* logicalDetector23[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector21[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector22[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector23[blkcor][tel];

  G4LogicalVolume* logicalDetector31[blkcor][tel];
  G4LogicalVolume* logicalDetector32[blkcor][tel];
  G4LogicalVolume* logicalDetector33[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector31[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector32[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector33[blkcor][tel];

  G4LogicalVolume* logicalDetector41[blkcor][tel];
  G4LogicalVolume* logicalDetector42[blkcor][tel];
  G4LogicalVolume* logicalDetector43[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector41[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector42[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector43[blkcor][tel];


  G4VisAttributes* Silicon_11VisAtt[4][4];
  G4VisAttributes* Silicon_12VisAtt[4][4];
  G4VisAttributes* CsIcrystal1VisAtt[4][4];
  G4VisAttributes* Silicon_21VisAtt[4][4];
  G4VisAttributes* Silicon_22VisAtt[4][4];
  G4VisAttributes* CsIcrystal2VisAtt[4][4];
  G4VisAttributes* Silicon_31VisAtt[4][4];
  G4VisAttributes* Silicon_32VisAtt[4][4];
  G4VisAttributes* CsIcrystal3VisAtt[4][4];
  G4VisAttributes* Silicon_41VisAtt[4][4];
  G4VisAttributes* Silicon_42VisAtt[4][4];
  G4VisAttributes* CsIcrystal4VisAtt[4][4];
  /////////////////////////////////////////////////////////////////////////////////////

  //////////////OUTER////////////////////////////////////////////////////////////////////5,7,9,11
  G4int silicon_511Num[blkcor][tel] = {{81,86,91,96},//+80 (1,6,11,16)
	{118,123,128,113},//+112 7
	{155,160,145,150},//+ 144 9
	{192,177,182,187}};//+176 11

  G4int silicon_512Num[blkcor][tel] = {{281,286,291,296},
	{318,323,328,313},
	{355,360,345,350},
	{392,377,382,387}};

  G4int csicrystal51Num[blkcor][tel] = {{481,486,491,496},
	{518,523,528,513},
	{555,560,545,550},
	{592,577,582,587}};
//------------------------------------------------------------------------------------
  G4int silicon_521Num[blkcor][tel] = {{82,87,92,93},//+80 (2,7,12,13)
	{119,124,125,114},//+112 7
	{156,157,146,151},//+144 9
	{189,178,183,188}};//176 11

  G4int silicon_522Num[blkcor][tel] = {{282,287,292,293},
	{319,324,325,314},
	{356,357,346,351},
	{389,378,383,388}};

  G4int csicrystal52Num[blkcor][tel] = {{482,487,492,493},
	{519,524,525,514},
	{556,557,546,551},
	{589,578,583,588}};
//------------------------------------------------------------------------------------
  G4int silicon_531Num[blkcor][tel] = {{84,85,90,95},//+80 (4,5,10,15)
	{117,122,127,116},//+112 7
	{154,159,148,149},//+144 9
	{191,180,181,186}};//176 11

  G4int silicon_532Num[blkcor][tel] = {{284,285,290,295},
	{317,322,327,316},
	{354,359,348,349},
	{391,380,381,386}};

  G4int csicrystal53Num[blkcor][tel] = {{484,485,490,495},
	{517,522,527,516},
	{554,559,548,549},
	{591,580,581,586}};
//------------------------------------------------------------------------------------
  G4int silicon_541Num[blkcor][tel] = {{83,88,89,94},//+80 (3,8,9,14)
	{120,121,126,115},//+112 7
	{153,158,147,152},//+144 9
	{190,179,184,185}};//+176 11

  G4int silicon_542Num[blkcor][tel] = {{283,288,289,294},
	{320,321,326,315},
	{353,358,347,352},
	{390,379,384,385}};

  G4int csicrystal54Num[blkcor][tel] = {{483,488,489,494},
	{520,521,526,515},
	{553,558,547,552},
	{590,579,584,585}};
  //////////////////////////////////////////////////////////////////////////////
  string silicon_511[blkcor][tel] ={ {"S1B5Q1T1","S1B5Q2T2","S1B5Q3T3","S1B5Q4T4"},
	{"S1B7Q1T1","S1B7Q2T2","S1B7Q3T3","S1B7Q4T4"},
	{"S1B9Q1T1","S1B9Q2T2","S1B9Q3T3","S1B9Q4T4"},
	{"S1B11Q1T1","S1B11Q2T2","S1B11Q3T3","S1B11Q4T4"}};

  string silicon_512[blkcor][tel] ={ {"S2B5Q1T1","S2B5Q2T2","S2B5Q3T3","S2B5Q4T4"},
	{"S2B7Q1T1","S2B7Q2T2","S2B7Q3T3","S2B7Q4T4"},
	{"S2B9Q1T1","S2B9Q2T2","S2B9Q3T3","S2B9Q4T4"},
	{"S2B11Q1T1","S2B11Q2T2","S2B11Q3T3","S2B11Q4T4"}};

  string csicrystal51[blkcor][tel] ={ {"CsIB5Q1T1","CsIB5Q2T2","CsIB5Q3T3","CsIB5Q4T4"},//-------------------
	{"CsIB7Q1T1","CsIB7Q2T2","CsIB7Q3T3","CsIB7Q4T4"},
	{"CsIB9Q1T1","CsIB9Q2T2","CsIB9Q3T3","CsIB9Q4T4"},
	{"CsIB11Q1T1","CsIB11Q2T2","CsIB11Q3T3","CsIB11Q4T4"}};

  string silicon_521[blkcor][tel] ={ {"S1B5Q1T2","S1B5Q2T3","S1B5Q3T4","S1B5Q4T1"},
	{"S1B7Q1T2","S1B7Q2T3","S1B7Q3T4","S1B7Q4T1"},
	{"S1B9Q1T2","S1B9Q2T3","S1B9Q3T4","S1B9Q4T1"},
	{"S1B11Q1T2","S1B11Q2T3","S1B11Q3T4","S1B11Q4T1"}};

  string silicon_522[blkcor][tel] ={ {"S2B5Q1T2","S2B5Q2T3","S2B5Q3T4","S2B5Q4T1"},
	{"S2B7Q1T2","S2B7Q2T3","S2B7Q3T4","S2B7Q4T1"},
	{"S2B9Q1T2","S2B9Q2T3","S2B9Q3T4","S2B9Q4T1"},
	{"S2B11Q1T2","S2B11Q2T3","S2B11Q3T4","S2B11Q4T1"}};

  string csicrystal52[blkcor][tel] = {{"CsIB5Q1T2","CsIB5Q2T3","CsIB5Q3T4","CsIB5Q4T1"},//-------------------
	{"CsIB7Q1T2","CsIB7Q2T3","CsIB7Q3T4","CsIB7Q4T1"},
	{"CsIB9Q1T2","CsIB9Q2T3","CsIB9Q3T4","CsIB9Q4T1"},
	{"CsIB11Q1T2","CsIB11Q2T3","CsIB11Q3T4","CsIB11Q4T1"}};

  string silicon_531[blkcor][tel] ={ {"S1B5Q1T4","S1B5Q2T1","S1B5Q3T2","S1B5Q4T3"},
	{"S1B7Q1T4","S1B7Q2T1","S1B7Q3T2","S1B7Q4T3"},
	{"S1B9Q1T4","S1B9Q2T1","S1B9Q3T2","S1B9Q4T3"},
	{"S1B11Q1T4","S1B11Q2T1","S1B11Q3T2","S1B11Q4T3"}};

  string silicon_532[blkcor][tel] ={ {"S2B5Q1T4","S2B5Q2T1","S2B5Q3T2","S2B5Q4T3"},
	{"S2B7Q1T4","S2B7Q2T1","S2B7Q3T2","S2B7Q4T3"},
	{"S2B9Q1T4","S2B9Q2T1","S2B9Q3T2","S2B9Q4T3"},
	{"S2B11Q1T4","S2B11Q2T1","S2B11Q3T2","S2B11Q4T3"}};

  string csicrystal53[blkcor][tel] ={ {"CsIB5Q1T4","CsIB5Q2T1","CsIB5Q3T2","CsIB5Q4T3"},//---------------------
	{"CsIB7Q1T4","CsIB7Q2T1","CsIB7Q3T2","CsIB7Q4T3"},
	{"CsIB9Q1T4","CsIB9Q2T1","CsIB9Q3T2","CsIB9Q4T3"},
	{"CsIB11Q1T4","CsIB11Q2T1","CsIB11Q3T2","CsIB11Q4T3"}};

  string silicon_541[blkcor][tel] ={ {"S1B5Q1T3","S1B5Q2T4","S1B5Q3T1","S1B5Q4T2"},
	{"S1B7Q1T3","S1B7Q2T4","S1B7Q3T1","S1B7Q4T2"},
	{"S1B9Q1T3","S1B9Q2T4","S1B9Q3T1","S1B9Q4T2"},
	{"S1B11Q1T3","S1B11Q2T4","S1B11Q3T1","S1B11Q4T2"}};

  string silicon_542[blkcor][tel] ={ {"S2B5Q1T3","S2B5Q2T4","S2B5Q3T1","S2B5Q4T2"},
	{"S2B7Q1T3","S2B7Q2T4","S2B7Q3T1","S2B7Q4T2"},
	{"S2B9Q1T3","S2B9Q2T4","S2B9Q3T1","S2B9Q4T2"},
	{"S2B11Q1T3","S2B11Q2T4","S2B11Q3T1","S2B11Q4T2"}};

  string csicrystal54[blkcor][tel] ={ {"CsIB5Q1T3","CsIB5Q2T4","CsIB5Q3T1","CsIB5Q4T2"},//------------------------
	{"CsIB7Q1T3","CsIB7Q2T4","CsIB7Q3T1","CsIB7Q4T2"},
	{"CsIB9Q1T3","CsIB9Q2T4","CsIB9Q3T1","CsIB9Q4T2"},
	{"CsIB11Q1T3","CsIB11Q2T4","CsIB11Q3T1","CsIB11Q4T2"}};

  string TubeSurface5[blkcor][tel] ={ {"CsI511","CsI512","CsI513","CsI514"},
	{"CsI711","CsI712","CsI713","CsI714"},
	{"CsI911","CsI912","CsI913","CsI914"},
	{"CsI1111","CsI1112","CsI1113","CsI1114"}};

  G4Box* Silicon_511[blkcor][tel];
  G4Box* Silicon_512[blkcor][tel];
  G4Trd* CsIcrystal51[blkcor][tel];

  G4Box* Silicon_521[blkcor][tel];
  G4Box* Silicon_522[blkcor][tel];
  G4Trd* CsIcrystal52[blkcor][tel];

  G4Box* Silicon_531[blkcor][tel];
  G4Box* Silicon_532[blkcor][tel];
  G4Trd* CsIcrystal53[blkcor][tel];

  G4Box* Silicon_541[blkcor][tel];
  G4Box* Silicon_542[blkcor][tel];
  G4Trd* CsIcrystal54[blkcor][tel];
  //G4Trd* CsIcrystal[1][qua][blkcor][tel];

  G4LogicalVolume* logicalDetector511[blkcor][tel];
  G4LogicalVolume* logicalDetector512[blkcor][tel];
  G4LogicalVolume* logicalDetector513[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector511[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector512[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector513[blkcor][tel];

  G4LogicalVolume* logicalDetector521[blkcor][tel];
  G4LogicalVolume* logicalDetector522[blkcor][tel];
  G4LogicalVolume* logicalDetector523[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector521[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector522[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector523[blkcor][tel];

  G4LogicalVolume* logicalDetector531[blkcor][tel];
  G4LogicalVolume* logicalDetector532[blkcor][tel];
  G4LogicalVolume* logicalDetector533[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector531[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector532[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector533[blkcor][tel];

  G4LogicalVolume* logicalDetector541[blkcor][tel];
  G4LogicalVolume* logicalDetector542[blkcor][tel];
  G4LogicalVolume* logicalDetector543[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector541[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector542[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector543[blkcor][tel];


  G4VisAttributes* Silicon_511VisAtt[4][4];
  G4VisAttributes* Silicon_512VisAtt[4][4];
  G4VisAttributes* CsIcrystal51VisAtt[4][4];
  G4VisAttributes* Silicon_521VisAtt[4][4];
  G4VisAttributes* Silicon_522VisAtt[4][4];
  G4VisAttributes* CsIcrystal52VisAtt[4][4];
  G4VisAttributes* Silicon_531VisAtt[4][4];
  G4VisAttributes* Silicon_532VisAtt[4][4];
  G4VisAttributes* CsIcrystal53VisAtt[4][4];
  G4VisAttributes* Silicon_541VisAtt[4][4];
  G4VisAttributes* Silicon_542VisAtt[4][4];
  G4VisAttributes* CsIcrystal54VisAtt[4][4];
  /////////////////////////////////////////////////////////////////////////////////////

  //////////////----MID----////////////////////////////////////////////////////////////////////4,6,8,10
  G4int silicon_411Num[blkcor][tel] = {{65,70,75,80},//+64 (1,6,11,16)
	{102,107,112,97},//+96 6
	{139,144,129,134},//+128 8
	{176,161,166,171}};//+160 10

  G4int silicon_412Num[blkcor][tel] = {{265,270,275,280},
	{302,307,312,297},
	{339,344,329,334},
	{376,361,366,371}};

  G4int csicrystal41Num[blkcor][tel] = {{465,470,475,480},
	{502,507,512,497},
	{539,544,529,534},
	{576,561,566,571}};
//------------------------------------------------------------------------------------
  G4int silicon_421Num[blkcor][tel] = {{66,71,76,77},//+64 (2,7,12,13)
	{103,108,109,98},//+96 6
	{140,141,130,135},//+128 8
	{173,162,167,172}};//+160 10

  G4int silicon_422Num[blkcor][tel] = {{266,271,276,277},
	{303,308,309,298},
	{340,341,330,335},
	{373,362,367,372}};

  G4int csicrystal42Num[blkcor][tel] = {{466,471,476,477},
	{503,508,509,498},
	{540,541,530,535},
	{573,562,567,572}};
//------------------------------------------------------------------------------------
  G4int silicon_431Num[blkcor][tel] = {{68,69,74,79},//+64 (4,5,10,15)
	{101,106,111,100},//+96 6
	{138,143,132,133},//+128 8
	{175,164,165,170}};//+160 10

  G4int silicon_432Num[blkcor][tel] = {{268,269,274,279},
	{301,306,311,300},
	{338,343,332,333},
	{375,364,365,370}};

  G4int csicrystal43Num[blkcor][tel] = {{468,469,474,479},
	{501,506,511,500},
	{538,543,532,533},
	{575,564,565,570}};
//------------------------------------------------------------------------------------
  G4int silicon_441Num[blkcor][tel] = {{67,72,73,78},//+64 (3,8,9,14)
	{104,105,110,99},//+96 6
	{137,142,131,136},//+128 8
	{174,163,168,169}};//+160 10

  G4int silicon_442Num[blkcor][tel] = {{267,272,273,278},
	{304,305,310,299},
	{337,342,331,336},
	{374,363,368,369}};

  G4int csicrystal44Num[blkcor][tel] = {{467,472,473,478},
	{504,505,510,499},
	{537,542,531,536},
	{574,563,568,569}};
  //////////////////////////////////////////////////////////////////////////////
  string silicon_411[blkcor][tel] ={ {"S1B4Q1T1","S1B4Q2T2","S1B4Q3T3","S1B4Q4T4"},
	{"S1B6Q1T1","S1B6Q2T2","S1B6Q3T3","S1B6Q4T4"},
	{"S1B8Q1T1","S1B8Q2T2","S1B8Q3T3","S1B8Q4T4"},
	{"S1B10Q1T1","S1B10Q2T2","S1B10Q3T3","S1B10Q4T4"}};

  string silicon_412[blkcor][tel] ={ {"S2B4Q1T1","S2B4Q2T2","S2B4Q3T3","S2B4Q4T4"},
	{"S2B6Q1T1","S2B6Q2T2","S2B6Q3T3","S2B6Q4T4"},
	{"S2B8Q1T1","S2B8Q2T2","S2B8Q3T3","S2B8Q4T4"},
	{"S2B10Q1T1","S2B10Q2T2","S2B10Q3T3","S2B10Q4T4"}};

  string csicrystal41[blkcor][tel] ={ {"CsIB4Q1T1","CsIB4Q2T2","CsIB4Q3T3","CsIB4Q4T4"},
	{"CsIB6Q1T1","CsIB6Q2T2","CsIB6Q3T3","CsIB6Q4T4"},
	{"CsIB8Q1T1","CsIB8Q2T2","CsIB8Q3T3","CsIB8Q4T4"},
	{"CsIB10Q1T1","CsIB10Q2T2","CsIB10Q3T3","CsIB10Q4T4"}};

  string silicon_421[blkcor][tel] ={ {"S1B4Q1T2","S1B4Q2T3","S1B4Q3T4","S1B4Q4T1"},
	{"S1B6Q1T2","S1B6Q2T3","S1B6Q3T4","S1B6Q4T1"},
	{"S1B8Q1T2","S1B8Q2T3","S1B8Q3T4","S1B8Q4T1"},
	{"S1B10Q1T2","S1B10Q2T3","S1B10Q3T4","S1B10Q4T1"}};

  string silicon_422[blkcor][tel] ={ {"S2B4Q1T2","S2B4Q2T3","S2B4Q3T4","S2B4Q4T1"},
	{"S2B6Q1T2","S2B6Q2T3","S2B6Q3T4","S2B6Q4T1"},
	{"S2B8Q1T2","S2B8Q2T3","S2B8Q3T4","S2B8Q4T1"},
	{"S2B10Q1T2","S2B10Q2T3","S2B10Q3T4","S2B10Q4T1"}};

  string csicrystal42[blkcor][tel] = {{"CsIB4Q1T2","CsIB4Q2T3","CsIB4Q3T4","CsIB4Q4T1"},
	{"CsIB6Q1T2","CsIB6Q2T3","CsIB6Q3T4","CsIB6Q4T1"},
	{"CsIB8Q1T2","CsIB8Q2T3","CsIB8Q3T4","CsIB8Q4T1"},
	{"CsIB10Q1T2","CsIB10Q2T3","CsIB10Q3T4","CsIB10Q4T1"}};

  string silicon_431[blkcor][tel] ={ {"S1B4Q1T4","S1B4Q2T1","S1B4Q3T2","S1B4Q4T3"},
	{"S1B6Q1T4","S1B6Q2T1","S1B6Q3T2","S1B6Q4T3"},
	{"S1B8Q1T4","S1B8Q2T1","S1B8Q3T2","S1B8Q4T3"},
	{"S1B10Q1T4","S1B10Q2T1","S1B10Q3T2","S1B10Q4T3"}};

  string silicon_432[blkcor][tel] ={ {"S2B4Q1T4","S2B4Q2T1","S2B4Q3T2","S2B4Q4T3"},
	{"S2B6Q1T4","S2B6Q2T1","S2B6Q3T2","S2B6Q4T3"},
	{"S2B8Q1T4","S2B8Q2T1","S2B8Q3T2","S2B8Q4T3"},
	{"S2B10Q1T4","S2B10Q2T1","S2B10Q3T2","S2B10Q4T3"}};

  string csicrystal43[blkcor][tel] ={ {"CsIB4Q1T4","CsIB4Q2T1","CsIB4Q3T2","CsIB4Q4T3"},
	{"CsIB6Q1T4","CsIB6Q2T1","CsIB6Q3T2","CsIB6Q4T3"},
	{"CsIB8Q1T4","CsIB8Q2T1","CsIB8Q3T2","CsIB8Q4T3"},
	{"CsIB10Q1T4","CsIB10Q2T1","CsIB10Q3T2","CsIB10Q4T3"}};

  string silicon_441[blkcor][tel] ={ {"S1B4Q1T3","S1B4Q2T4","S1B4Q3T1","S1B4Q4T2"},
	{"S1B6Q1T3","S1B6Q2T4","S1B6Q3T1","S1B6Q4T2"},
	{"S1B8Q1T3","S1B8Q2T4","S1B8Q3T1","S1B8Q4T2"},
	{"S1B10Q1T3","S1B10Q2T4","S1B10Q3T1","S1B10Q4T2"}};

  string silicon_442[blkcor][tel] ={ {"S2B4Q1T3","S2B4Q2T4","S2B4Q3T1","S2B4Q4T2"},
	{"S2B6Q1T3","S2B6Q2T4","S2B6Q3T1","S2B6Q4T2"},
	{"S2B8Q1T3","S2B8Q2T4","S2B8Q3T1","S2B8Q4T2"},
	{"S2B10Q1T3","S2B10Q2T4","S2B10Q3T1","S2B10Q4T2"}};

  string csicrystal44[blkcor][tel] ={ {"CsIB4Q1T3","CsIB4Q2T4","CsIB4Q3T1","CsIB4Q4T2"},
	{"CsIB6Q1T3","CsIB6Q2T4","CsIB6Q3T1","CsIB6Q4T2"},
	{"CsIB8Q1T3","CsIB8Q2T4","CsIB8Q3T1","CsIB8Q4T2"},
	{"CsIB10Q1T3","CsIB10Q2T4","CsIB10Q3T1","CsIB10Q4T2"}};

  string TubeSurface4[blkcor][tel] ={ {"CsI411","CsI412","CsI413","CsI414"},
	{"CsI611","CsI612","CsI613","CsI614"},
	{"CsI811","CsI812","CsI813","CsI814"},
	{"CsI1011","CsI1012","CsI1013","CsI1014"}};

  G4Box* Silicon_411[blkcor][tel];
  G4Box* Silicon_412[blkcor][tel];
  G4Trd* CsIcrystal41[blkcor][tel];

  G4Box* Silicon_421[blkcor][tel];
  G4Box* Silicon_422[blkcor][tel];
  G4Trd* CsIcrystal42[blkcor][tel];

  G4Box* Silicon_431[blkcor][tel];
  G4Box* Silicon_432[blkcor][tel];
  G4Trd* CsIcrystal43[blkcor][tel];

  G4Box* Silicon_441[blkcor][tel];
  G4Box* Silicon_442[blkcor][tel];
  G4Trd* CsIcrystal44[blkcor][tel];
  //G4Trd* CsIcrystal[1][qua][blkcor][tel];

  G4LogicalVolume* logicalDetector411[blkcor][tel];
  G4LogicalVolume* logicalDetector412[blkcor][tel];
  G4LogicalVolume* logicalDetector413[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector411[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector412[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector413[blkcor][tel];

  G4LogicalVolume* logicalDetector421[blkcor][tel];
  G4LogicalVolume* logicalDetector422[blkcor][tel];
  G4LogicalVolume* logicalDetector423[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector421[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector422[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector423[blkcor][tel];

  G4LogicalVolume* logicalDetector431[blkcor][tel];
  G4LogicalVolume* logicalDetector432[blkcor][tel];
  G4LogicalVolume* logicalDetector433[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector431[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector432[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector433[blkcor][tel];

  G4LogicalVolume* logicalDetector441[blkcor][tel];
  G4LogicalVolume* logicalDetector442[blkcor][tel];
  G4LogicalVolume* logicalDetector443[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector441[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector442[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector443[blkcor][tel];


  G4VisAttributes* Silicon_411VisAtt[4][4];
  G4VisAttributes* Silicon_412VisAtt[4][4];
  G4VisAttributes* CsIcrystal41VisAtt[4][4];
  G4VisAttributes* Silicon_421VisAtt[4][4];
  G4VisAttributes* Silicon_422VisAtt[4][4];
  G4VisAttributes* CsIcrystal42VisAtt[4][4];
  G4VisAttributes* Silicon_431VisAtt[4][4];
  G4VisAttributes* Silicon_432VisAtt[4][4];
  G4VisAttributes* CsIcrystal43VisAtt[4][4];
  G4VisAttributes* Silicon_441VisAtt[4][4];
  G4VisAttributes* Silicon_442VisAtt[4][4];
  G4VisAttributes* CsIcrystal44VisAtt[4][4];
  /////////////////////////////////////////////////////////////////////////////////////
  G4OpticalSurface* opTubeSurface[blkcor][tel];

  logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());

  G4RotationMatrix* Rot_1Si[blkcor][tel];
  G4RotationMatrix* Rot_2Si[blkcor][tel];
  G4RotationMatrix* Rot_3Si[blkcor][tel];
  G4RotationMatrix* Rot_4Si[blkcor][tel];

  G4RotationMatrix* Rot_41Si[blkcor][tel];
  G4RotationMatrix* Rot_42Si[blkcor][tel];
  G4RotationMatrix* Rot_43Si[blkcor][tel];
  G4RotationMatrix* Rot_44Si[blkcor][tel];

  G4RotationMatrix* Rot_51Si[blkcor][tel];
  G4RotationMatrix* Rot_52Si[blkcor][tel];
  G4RotationMatrix* Rot_53Si[blkcor][tel];
  G4RotationMatrix* Rot_54Si[blkcor][tel];
  //----------------------------------------------------------------------------------------------------------
  for(G4int Blkcor = 0 ; Blkcor < 4 ; Blkcor++)
  {
	for(G4int Tel = 0 ; Tel < 4 ; Tel++)
	{
	  //COR--------------------------------------------------------------------------
	  Rot_1Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_2Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_3Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_4Si[Blkcor][Tel] = new G4RotationMatrix;

	  if(Tel == 0){
		Rot_1Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta1)*TMath::Sin(Blkphi1+(pi*Blkcor/2)))/(TMath::Cos(Blktheta1))));
		Rot_1Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta1)*TMath::Cos(Blkphi1+(pi*Blkcor/2))));
		Rot_2Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta2)*TMath::Sin(Blkphi2+(pi*Blkcor/2)))/(TMath::Cos(Blktheta2))));
		Rot_2Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta2)*TMath::Cos(Blkphi2+(pi*Blkcor/2))));
		Rot_3Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta3)*TMath::Sin(Blkphi3+(pi*Blkcor/2)))/(TMath::Cos(Blktheta3))));
		Rot_3Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta3)*TMath::Cos(Blkphi3+(pi*Blkcor/2))));
		Rot_4Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta4)*TMath::Sin(Blkphi4+(pi*Blkcor/2)))/(TMath::Cos(Blktheta4))));
		Rot_4Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta4)*TMath::Cos(Blkphi4+(pi*Blkcor/2))));}
	  if(Tel == 1){
		Rot_1Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta1)*TMath::Sin(Blkphi1+(pi*Blkcor/2)))/(TMath::Cos(Blktheta1))));
		Rot_1Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta1)*TMath::Cos(Blkphi1+(pi*Blkcor/2))));
		Rot_2Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta2)*TMath::Sin(Blkphi2+(pi*Blkcor/2)))/(TMath::Cos(Blktheta2))));
		Rot_2Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta2)*TMath::Cos(Blkphi2+(pi*Blkcor/2))));
		Rot_3Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta3)*TMath::Sin(Blkphi3+(pi*Blkcor/2)))/(TMath::Cos(Blktheta3))));
		Rot_3Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta3)*TMath::Cos(Blkphi3+(pi*Blkcor/2))));
		Rot_4Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta4)*TMath::Sin(Blkphi4+(pi*Blkcor/2)))/(TMath::Cos(Blktheta4))));
		Rot_4Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta4)*TMath::Cos(Blkphi4+(pi*Blkcor/2))));}
	  if(Tel == 2){
		Rot_1Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta1)*TMath::Sin(Blkphi1+(pi*Blkcor/2)))/(TMath::Cos(Blktheta1))));
		Rot_1Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta1)*TMath::Cos(Blkphi1+(pi*Blkcor/2))));
		Rot_2Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta2)*TMath::Sin(Blkphi2+(pi*Blkcor/2)))/(TMath::Cos(Blktheta2))));
		Rot_2Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta2)*TMath::Cos(Blkphi2+(pi*Blkcor/2))));
		Rot_3Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta3)*TMath::Sin(Blkphi3+(pi*Blkcor/2)))/(TMath::Cos(Blktheta3))));
		Rot_3Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta3)*TMath::Cos(Blkphi3+(pi*Blkcor/2))));
		Rot_4Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta4)*TMath::Sin(Blkphi4+(pi*Blkcor/2)))/(TMath::Cos(Blktheta4))));
		Rot_4Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta4)*TMath::Cos(Blkphi4+(pi*Blkcor/2))));}
	  if(Tel == 3){
		Rot_1Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta1)*TMath::Sin(Blkphi1+(pi*Blkcor/2)))/(TMath::Cos(Blktheta1))));
		Rot_1Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta1)*TMath::Cos(Blkphi1+(pi*Blkcor/2))));
		Rot_2Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta2)*TMath::Sin(Blkphi2+(pi*Blkcor/2)))/(TMath::Cos(Blktheta2))));
		Rot_2Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta2)*TMath::Cos(Blkphi2+(pi*Blkcor/2))));
		Rot_3Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta3)*TMath::Sin(Blkphi3+(pi*Blkcor/2)))/(TMath::Cos(Blktheta3))));
		Rot_3Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta3)*TMath::Cos(Blkphi3+(pi*Blkcor/2))));
		Rot_4Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta4)*TMath::Sin(Blkphi4+(pi*Blkcor/2)))/(TMath::Cos(Blktheta4))));
		Rot_4Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan(TMath::Sin(Blktheta4)*TMath::Cos(Blkphi4+(pi*Blkcor/2))));}
	  //MID-------------------------------------------------------------------
	  Rot_41Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_42Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_43Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_44Si[Blkcor][Tel] = new G4RotationMatrix;

	  if(Tel == 0){
		Rot_41Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Sin(Blkphi41+(pi*Blkcor/2)))/(TMath::Cos(Blktheta41))));
		Rot_41Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Cos(Blkphi41+(pi*Blkcor/2)))));
		Rot_42Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Sin(Blkphi42+(pi*Blkcor/2)))/(TMath::Cos(Blktheta42))));
		Rot_42Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Cos(Blkphi42+(pi*Blkcor/2)))));
		Rot_43Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Sin(Blkphi43+(pi*Blkcor/2)))/(TMath::Cos(Blktheta43))));
		Rot_43Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Cos(Blkphi43+(pi*Blkcor/2)))));
		Rot_44Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Sin(Blkphi44+(pi*Blkcor/2)))/(TMath::Cos(Blktheta44))));
		Rot_44Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Cos(Blkphi44+(pi*Blkcor/2)))));}
	  if(Tel == 1){
		Rot_41Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Sin(Blkphi41+(pi*Blkcor/2)))/(TMath::Cos(Blktheta41))));
		Rot_41Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Cos(Blkphi41+(pi*Blkcor/2)))));
		Rot_42Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Sin(Blkphi42+(pi*Blkcor/2)))/(TMath::Cos(Blktheta42))));
		Rot_42Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Cos(Blkphi42+(pi*Blkcor/2)))));
		Rot_43Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Sin(Blkphi43+(pi*Blkcor/2)))/(TMath::Cos(Blktheta43))));
		Rot_43Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Cos(Blkphi43+(pi*Blkcor/2)))));
		Rot_44Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Sin(Blkphi44+(pi*Blkcor/2)))/(TMath::Cos(Blktheta44))));
		Rot_44Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Cos(Blkphi44+(pi*Blkcor/2)))));}
	  if(Tel == 2){
		Rot_41Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Sin(Blkphi41+(pi*Blkcor/2)))/(TMath::Cos(Blktheta41))));
		Rot_41Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Cos(Blkphi41+(pi*Blkcor/2)))));
		Rot_42Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Sin(Blkphi42+(pi*Blkcor/2)))/(TMath::Cos(Blktheta42))));
		Rot_42Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Cos(Blkphi42+(pi*Blkcor/2)))));
		Rot_43Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Sin(Blkphi43+(pi*Blkcor/2)))/(TMath::Cos(Blktheta43))));
		Rot_43Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Cos(Blkphi43+(pi*Blkcor/2)))));
		Rot_44Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Sin(Blkphi44+(pi*Blkcor/2)))/(TMath::Cos(Blktheta44))));
		Rot_44Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Cos(Blkphi44+(pi*Blkcor/2)))));}
	  if(Tel == 3){
		Rot_41Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Sin(Blkphi41+(pi*Blkcor/2)))/(TMath::Cos(Blktheta41))));
		Rot_41Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta41)*TMath::Cos(Blkphi41+(pi*Blkcor/2)))));
		Rot_42Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Sin(Blkphi42+(pi*Blkcor/2)))/(TMath::Cos(Blktheta42))));
		Rot_42Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta42)*TMath::Cos(Blkphi42+(pi*Blkcor/2)))));
		Rot_43Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Sin(Blkphi43+(pi*Blkcor/2)))/(TMath::Cos(Blktheta43))));
		Rot_43Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta43)*TMath::Cos(Blkphi43+(pi*Blkcor/2)))));
		Rot_44Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Sin(Blkphi44+(pi*Blkcor/2)))/(TMath::Cos(Blktheta44))));
		Rot_44Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta44)*TMath::Cos(Blkphi44+(pi*Blkcor/2)))));}
	  //OUTER-------------------------------------------------------
	  Rot_51Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_52Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_53Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_54Si[Blkcor][Tel] = new G4RotationMatrix;

	  if(Tel == 0){
		Rot_51Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Sin(Blkphi51+(pi*Blkcor/2)))/(TMath::Cos(Blktheta51))));
		Rot_51Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Cos(Blkphi51+(pi*Blkcor/2)))));
		Rot_52Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Sin(Blkphi52+(pi*Blkcor/2)))/(TMath::Cos(Blktheta52))));
		Rot_52Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Cos(Blkphi52+(pi*Blkcor/2)))));
		Rot_53Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Sin(Blkphi53+(pi*Blkcor/2)))/(TMath::Cos(Blktheta53))));
		Rot_53Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Cos(Blkphi53+(pi*Blkcor/2)))));
		Rot_54Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Sin(Blkphi54+(pi*Blkcor/2)))/(TMath::Cos(Blktheta54))));
		Rot_54Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Cos(Blkphi54+(pi*Blkcor/2)))));}
	  if(Tel == 1){
		Rot_51Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Sin(Blkphi51+(pi*Blkcor/2)))/(TMath::Cos(Blktheta51))));
		Rot_51Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Cos(Blkphi51+(pi*Blkcor/2)))));
		Rot_52Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Sin(Blkphi52+(pi*Blkcor/2)))/(TMath::Cos(Blktheta52))));
		Rot_52Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Cos(Blkphi52+(pi*Blkcor/2)))));
		Rot_53Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Sin(Blkphi53+(pi*Blkcor/2)))/(TMath::Cos(Blktheta53))));
		Rot_53Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Cos(Blkphi53+(pi*Blkcor/2)))));
		Rot_54Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Sin(Blkphi54+(pi*Blkcor/2)))/(TMath::Cos(Blktheta54))));
		Rot_54Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Cos(Blkphi54+(pi*Blkcor/2)))));}
	  if(Tel == 2){
		Rot_51Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Sin(Blkphi51+(pi*Blkcor/2)))/(TMath::Cos(Blktheta51))));
		Rot_51Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Cos(Blkphi51+(pi*Blkcor/2)))));
		Rot_52Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Sin(Blkphi52+(pi*Blkcor/2)))/(TMath::Cos(Blktheta52))));
		Rot_52Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Cos(Blkphi52+(pi*Blkcor/2)))));
		Rot_53Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Sin(Blkphi53+(pi*Blkcor/2)))/(TMath::Cos(Blktheta53))));
		Rot_53Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Cos(Blkphi53+(pi*Blkcor/2)))));
		Rot_54Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Sin(Blkphi54+(pi*Blkcor/2)))/(TMath::Cos(Blktheta54))));
		Rot_54Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Cos(Blkphi54+(pi*Blkcor/2)))));}
	  if(Tel == 3){
		Rot_51Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Sin(Blkphi51+(pi*Blkcor/2)))/(TMath::Cos(Blktheta51))));
		Rot_51Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta51)*TMath::Cos(Blkphi51+(pi*Blkcor/2)))));
		Rot_52Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Sin(Blkphi52+(pi*Blkcor/2)))/(TMath::Cos(Blktheta52))));
		Rot_52Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta52)*TMath::Cos(Blkphi52+(pi*Blkcor/2)))));
		Rot_53Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Sin(Blkphi53+(pi*Blkcor/2)))/(TMath::Cos(Blktheta53))));
		Rot_53Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta53)*TMath::Cos(Blkphi53+(pi*Blkcor/2)))));
		Rot_54Si[Blkcor][Tel] -> rotateX(1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Sin(Blkphi54+(pi*Blkcor/2)))/(TMath::Cos(Blktheta54))));
		Rot_54Si[Blkcor][Tel] -> rotateY(-1*TMath::ATan((TMath::Sin(Blktheta54)*TMath::Cos(Blkphi54+(pi*Blkcor/2)))));}
	  //one block's geometry-COR---------------------------------------------------------------------------------------------------------		
	  {//1

		//Si1
		Silicon_11[Blkcor][Tel] = new G4Box(silicon_11[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector11[Blkcor][Tel] = new G4LogicalVolume(Silicon_11[Blkcor][Tel],Silicon_mat,silicon_11[Blkcor][Tel]);
		PhysicalDetector11[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Cos(Blktheta1))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector11[Blkcor][Tel],//theta == 0.0144
			silicon_11[Blkcor][Tel],
			logicWorld,
			false,
			silicon_11Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_12[Blkcor][Tel] = new G4Box(silicon_12[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector12[Blkcor][Tel] = new G4LogicalVolume(Silicon_12[Blkcor][Tel] ,Silicon_mat ,silicon_12[Blkcor][Tel]);

		PhysicalDetector12[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Cos(Blktheta1))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector12[Blkcor][Tel],
			silicon_12[Blkcor][Tel],
			logicWorld,
			false,
			silicon_12Num[Blkcor][Tel],
			overlap);

		//CsI crystal		  
		CsIcrystal1[Blkcor][Tel] = new G4Trd(csicrystal1[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector13[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal1[Blkcor][Tel] ,scintillator_mat ,csicrystal1[Blkcor][Tel]);

		PhysicalDetector13[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Cos(Blktheta1))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector13[Blkcor][Tel],
			csicrystal1[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal1Num[Blkcor][Tel],
			overlap);
	  }
	  {//2
		Silicon_21[Blkcor][Tel] = new G4Box(silicon_21[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector21[Blkcor][Tel] = new G4LogicalVolume(Silicon_21[Blkcor][Tel],Silicon_mat,silicon_21[Blkcor][Tel]);
		PhysicalDetector21[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Cos(Blktheta2))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector21[Blkcor][Tel],//theta == 0.0316
			silicon_21[Blkcor][Tel],
			logicWorld,
			false,
			silicon_21Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_22[Blkcor][Tel] = new G4Box(silicon_22[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector22[Blkcor][Tel] = new G4LogicalVolume(Silicon_22[Blkcor][Tel] ,Silicon_mat ,silicon_22[Blkcor][Tel]);

		PhysicalDetector22[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Cos(Blktheta2))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector22[Blkcor][Tel],
			silicon_22[Blkcor][Tel],
			logicWorld,
			false,
			silicon_22Num[Blkcor][Tel],
			overlap);

		//CsI crystal		  
		CsIcrystal2[Blkcor][Tel] = new G4Trd(csicrystal2[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector23[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal2[Blkcor][Tel] ,scintillator_mat ,csicrystal2[Blkcor][Tel]);

		PhysicalDetector23[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Cos(Blktheta2))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector23[Blkcor][Tel],
			csicrystal2[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal2Num[Blkcor][Tel],
			overlap);
	  }
	  {//3
		Silicon_31[Blkcor][Tel] = new G4Box(silicon_31[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector31[Blkcor][Tel] = new G4LogicalVolume(Silicon_31[Blkcor][Tel],Silicon_mat,silicon_31[Blkcor][Tel]);
		PhysicalDetector31[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Cos(Blktheta3))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector31[Blkcor][Tel],//theta == 0.0147
			silicon_31[Blkcor][Tel],
			logicWorld,
			false,
			silicon_31Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_32[Blkcor][Tel] = new G4Box(silicon_32[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector32[Blkcor][Tel] = new G4LogicalVolume(Silicon_32[Blkcor][Tel] ,Silicon_mat ,silicon_32[Blkcor][Tel]);

		PhysicalDetector32[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Cos(Blktheta3))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector32[Blkcor][Tel],
			silicon_32[Blkcor][Tel],
			logicWorld,
			false,
			silicon_32Num[Blkcor][Tel],
			overlap);

		//CsI crystal		  
		CsIcrystal3[Blkcor][Tel] = new G4Trd(csicrystal3[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector33[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal3[Blkcor][Tel] ,scintillator_mat ,csicrystal3[Blkcor][Tel]);

		PhysicalDetector33[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Cos(Blktheta3))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector33[Blkcor][Tel],
			csicrystal3[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal3Num[Blkcor][Tel],
			overlap);
	  }
	  {//4
		Silicon_41[Blkcor][Tel] = new G4Box(silicon_41[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector41[Blkcor][Tel] = new G4LogicalVolume(Silicon_41[Blkcor][Tel],Silicon_mat,silicon_41[Blkcor][Tel]);
		PhysicalDetector41[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Cos(Blktheta4))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector41[Blkcor][Tel],//theta == 0.0147
			silicon_41[Blkcor][Tel],
			logicWorld,
			false,
			silicon_41Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_42[Blkcor][Tel] = new G4Box(silicon_42[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector42[Blkcor][Tel] = new G4LogicalVolume(Silicon_42[Blkcor][Tel] ,Silicon_mat ,silicon_42[Blkcor][Tel]);

		PhysicalDetector42[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Cos(Blktheta4))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector42[Blkcor][Tel],
			silicon_42[Blkcor][Tel],
			logicWorld,
			false,
			silicon_42Num[Blkcor][Tel],
			overlap);

		//CsI crystal		  
		CsIcrystal4[Blkcor][Tel] = new G4Trd(csicrystal4[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector43[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal4[Blkcor][Tel] ,scintillator_mat ,csicrystal4[Blkcor][Tel]);

		PhysicalDetector43[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Cos(Blktheta4))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(CORtheta)*TMath::Cos(CORphi+(RoDi*pi*Blkcor/2)),TMath::Sin(CORtheta)*TMath::Sin(CORphi+(RoDi*pi*Blkcor/2)),TMath::Cos(CORtheta))),
			logicalDetector43[Blkcor][Tel],
			csicrystal4[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal4Num[Blkcor][Tel],
			overlap);
	  }
	  //-------------------------------------------------------------------------------------------------------------------------------		
	  //one block's geometry-MID---------------------------------------------------------------------------------------------------------		
	  {//1
		//Si1
		Silicon_411[Blkcor][Tel] = new G4Box(silicon_411[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector411[Blkcor][Tel] = new G4LogicalVolume(Silicon_411[Blkcor][Tel],Silicon_mat,silicon_411[Blkcor][Tel]);
		PhysicalDetector411[Blkcor][Tel] = new G4PVPlacement(Rot_41Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta41))*TMath::Cos((Blkphi41+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Sin(Blktheta41))*TMath::Sin((Blkphi41+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Cos(Blktheta41))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector411[Blkcor][Tel],//theta4 == 0.0144
			silicon_411[Blkcor][Tel],
			logicWorld,
			false,
			silicon_411Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_412[Blkcor][Tel] = new G4Box(silicon_412[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector412[Blkcor][Tel] = new G4LogicalVolume(Silicon_412[Blkcor][Tel] ,Silicon_mat ,silicon_412[Blkcor][Tel]);

		PhysicalDetector412[Blkcor][Tel] = new G4PVPlacement(Rot_41Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta41))*TMath::Cos((Blkphi41+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Sin(Blktheta41))*TMath::Sin((Blkphi41+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Cos(Blktheta41))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector412[Blkcor][Tel],
			silicon_412[Blkcor][Tel],
			logicWorld,
			false,
			silicon_412Num[Blkcor][Tel],
			overlap);

		//CsI crystal4		  
		CsIcrystal41[Blkcor][Tel] = new G4Trd(csicrystal41[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector413[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal41[Blkcor][Tel] ,scintillator_mat ,csicrystal41[Blkcor][Tel]);

		PhysicalDetector413[Blkcor][Tel] = new G4PVPlacement(Rot_41Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta41))*TMath::Cos((Blkphi41+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Sin(Blktheta41))*TMath::Sin((Blkphi41+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Cos(Blktheta41))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector413[Blkcor][Tel],
			csicrystal41[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal41Num[Blkcor][Tel],
			overlap);
	  }
	  {//2
		Silicon_421[Blkcor][Tel] = new G4Box(silicon_421[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector421[Blkcor][Tel] = new G4LogicalVolume(Silicon_421[Blkcor][Tel],Silicon_mat,silicon_421[Blkcor][Tel]);
		PhysicalDetector421[Blkcor][Tel] = new G4PVPlacement(Rot_42Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta42))*TMath::Cos(Blkphi42+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Sin(Blktheta42))*TMath::Sin(Blkphi42+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Cos(Blktheta42))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector421[Blkcor][Tel],//theta4 == 0.0316
			silicon_421[Blkcor][Tel],
			logicWorld,
			false,
			silicon_421Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_422[Blkcor][Tel] = new G4Box(silicon_422[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector422[Blkcor][Tel] = new G4LogicalVolume(Silicon_422[Blkcor][Tel] ,Silicon_mat ,silicon_422[Blkcor][Tel]);

		PhysicalDetector422[Blkcor][Tel] = new G4PVPlacement(Rot_42Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta42))*TMath::Cos(Blkphi42+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Sin(Blktheta42))*TMath::Sin(Blkphi42+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Cos(Blktheta42))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector422[Blkcor][Tel],
			silicon_422[Blkcor][Tel],
			logicWorld,
			false,
			silicon_422Num[Blkcor][Tel],
			overlap);

		//CsI crystal4		  
		CsIcrystal42[Blkcor][Tel] = new G4Trd(csicrystal42[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector423[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal42[Blkcor][Tel] ,scintillator_mat ,csicrystal42[Blkcor][Tel]);

		PhysicalDetector423[Blkcor][Tel] = new G4PVPlacement(Rot_42Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta42))*TMath::Cos(Blkphi42+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Sin(Blktheta42))*TMath::Sin(Blkphi42+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Cos(Blktheta42))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector423[Blkcor][Tel],
			csicrystal42[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal42Num[Blkcor][Tel],
			overlap);
	  }
	  {//3
		Silicon_431[Blkcor][Tel] = new G4Box(silicon_431[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector431[Blkcor][Tel] = new G4LogicalVolume(Silicon_431[Blkcor][Tel],Silicon_mat,silicon_431[Blkcor][Tel]);
		PhysicalDetector431[Blkcor][Tel] = new G4PVPlacement(Rot_43Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta43))*TMath::Cos(Blkphi43+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Sin(Blktheta43))*TMath::Sin(Blkphi43+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Cos(Blktheta43))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector431[Blkcor][Tel],//theta4 == 0.0147
			silicon_431[Blkcor][Tel],
			logicWorld,
			false,
			silicon_431Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_432[Blkcor][Tel] = new G4Box(silicon_432[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector432[Blkcor][Tel] = new G4LogicalVolume(Silicon_432[Blkcor][Tel] ,Silicon_mat ,silicon_432[Blkcor][Tel]);

		PhysicalDetector432[Blkcor][Tel] = new G4PVPlacement(Rot_43Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta43))*TMath::Cos(Blkphi43+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Sin(Blktheta43))*TMath::Sin(Blkphi43+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Cos(Blktheta43))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector432[Blkcor][Tel],
			silicon_432[Blkcor][Tel],
			logicWorld,
			false,
			silicon_432Num[Blkcor][Tel],
			overlap);

		//CsI crystal4		  
		CsIcrystal43[Blkcor][Tel] = new G4Trd(csicrystal43[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector433[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal43[Blkcor][Tel] ,scintillator_mat ,csicrystal43[Blkcor][Tel]);

		PhysicalDetector433[Blkcor][Tel] = new G4PVPlacement(Rot_43Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta43))*TMath::Cos(Blkphi43+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Sin(Blktheta43))*TMath::Sin(Blkphi43+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Cos(Blktheta43))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector433[Blkcor][Tel],
			csicrystal43[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal43Num[Blkcor][Tel],
			overlap);
	  }
	  {//4
		Silicon_441[Blkcor][Tel] = new G4Box(silicon_441[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector441[Blkcor][Tel] = new G4LogicalVolume(Silicon_441[Blkcor][Tel],Silicon_mat,silicon_441[Blkcor][Tel]);
		PhysicalDetector441[Blkcor][Tel] = new G4PVPlacement(Rot_44Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta44))*TMath::Cos((Blkphi44+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Sin(Blktheta44))*TMath::Sin((Blkphi44+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Cos(Blktheta44))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector441[Blkcor][Tel],//theta4 == 0.0147
			silicon_441[Blkcor][Tel],
			logicWorld,
			false,
			silicon_441Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_442[Blkcor][Tel] = new G4Box(silicon_442[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector442[Blkcor][Tel] = new G4LogicalVolume(Silicon_442[Blkcor][Tel] ,Silicon_mat ,silicon_442[Blkcor][Tel]);

		PhysicalDetector442[Blkcor][Tel] = new G4PVPlacement(Rot_44Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta44))*TMath::Cos((Blkphi44+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Sin(Blktheta44))*TMath::Sin((Blkphi44+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Cos(Blktheta44))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector442[Blkcor][Tel],
			silicon_442[Blkcor][Tel],
			logicWorld,
			false,
			silicon_442Num[Blkcor][Tel],
			overlap);

		//CsI crystal4		  
		CsIcrystal44[Blkcor][Tel] = new G4Trd(csicrystal44[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector443[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal44[Blkcor][Tel] ,scintillator_mat ,csicrystal44[Blkcor][Tel]);

		PhysicalDetector443[Blkcor][Tel] = new G4PVPlacement(Rot_44Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta44))*TMath::Cos((Blkphi44+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Sin(Blktheta44))*TMath::Sin((Blkphi44+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Cos(Blktheta44))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(MIDtheta)*TMath::Cos(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Sin(MIDtheta)*TMath::Sin(MIDphi+(RoDi*pi*Blkcor/2)),TMath::Cos(MIDtheta))),
			logicalDetector443[Blkcor][Tel],
			csicrystal44[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal44Num[Blkcor][Tel],
			overlap);
	  }
	  //--------------------------------------------------------------------------------------------------------------------------------
	  //one block's geometry-OUTER---------------------------------------------------------------------------------------------------------		
	  {//1
		//Si1
		Silicon_511[Blkcor][Tel] = new G4Box(silicon_511[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector511[Blkcor][Tel] = new G4LogicalVolume(Silicon_511[Blkcor][Tel],Silicon_mat,silicon_511[Blkcor][Tel]);
		PhysicalDetector511[Blkcor][Tel] = new G4PVPlacement(Rot_51Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta51))*TMath::Cos((Blkphi51+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Sin(Blktheta51))*TMath::Sin((Blkphi51+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Cos(Blktheta51))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector511[Blkcor][Tel],//theta == 0.0144
			silicon_511[Blkcor][Tel],
			logicWorld,
			false,
			silicon_511Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_512[Blkcor][Tel] = new G4Box(silicon_512[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector512[Blkcor][Tel] = new G4LogicalVolume(Silicon_512[Blkcor][Tel] ,Silicon_mat ,silicon_512[Blkcor][Tel]);

		PhysicalDetector512[Blkcor][Tel] = new G4PVPlacement(Rot_51Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta51))*TMath::Cos((Blkphi51+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Sin(Blktheta51))*TMath::Sin((Blkphi51+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Cos(Blktheta51))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector512[Blkcor][Tel],
			silicon_512[Blkcor][Tel],
			logicWorld,
			false,
			silicon_512Num[Blkcor][Tel],
			overlap);

		//CsI crystal5		  
		CsIcrystal51[Blkcor][Tel] = new G4Trd(csicrystal51[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector513[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal51[Blkcor][Tel] ,scintillator_mat ,csicrystal51[Blkcor][Tel]);

		PhysicalDetector513[Blkcor][Tel] = new G4PVPlacement(Rot_51Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta51))*TMath::Cos((Blkphi51+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Sin(Blktheta51))*TMath::Sin((Blkphi51+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Cos(Blktheta51))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector513[Blkcor][Tel],
			csicrystal51[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal51Num[Blkcor][Tel],
			overlap);
	  }
	  {//2
		Silicon_521[Blkcor][Tel] = new G4Box(silicon_521[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector521[Blkcor][Tel] = new G4LogicalVolume(Silicon_521[Blkcor][Tel],Silicon_mat,silicon_521[Blkcor][Tel]);
		PhysicalDetector521[Blkcor][Tel] = new G4PVPlacement(Rot_52Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta52))*TMath::Cos(Blkphi52+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Sin(Blktheta52))*TMath::Sin(Blkphi52+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Cos(Blktheta52))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector521[Blkcor][Tel],//theta == 0.0316
			silicon_521[Blkcor][Tel],
			logicWorld,
			false,
			silicon_521Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_522[Blkcor][Tel] = new G4Box(silicon_522[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector522[Blkcor][Tel] = new G4LogicalVolume(Silicon_522[Blkcor][Tel] ,Silicon_mat ,silicon_522[Blkcor][Tel]);

		PhysicalDetector522[Blkcor][Tel] = new G4PVPlacement(Rot_52Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta52))*TMath::Cos(Blkphi52+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Sin(Blktheta52))*TMath::Sin(Blkphi52+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Cos(Blktheta52))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector522[Blkcor][Tel],
			silicon_522[Blkcor][Tel],
			logicWorld,
			false,
			silicon_522Num[Blkcor][Tel],
			overlap);

		//CsI crystal5		  
		CsIcrystal52[Blkcor][Tel] = new G4Trd(csicrystal52[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector523[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal52[Blkcor][Tel] ,scintillator_mat ,csicrystal52[Blkcor][Tel]);

		PhysicalDetector523[Blkcor][Tel] = new G4PVPlacement(Rot_52Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta52))*TMath::Cos(Blkphi52+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Sin(Blktheta52))*TMath::Sin(Blkphi52+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Cos(Blktheta52))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector523[Blkcor][Tel],
			csicrystal52[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal52Num[Blkcor][Tel],
			overlap);
	  }
	  {//3
		Silicon_531[Blkcor][Tel] = new G4Box(silicon_531[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector531[Blkcor][Tel] = new G4LogicalVolume(Silicon_531[Blkcor][Tel],Silicon_mat,silicon_531[Blkcor][Tel]);
		PhysicalDetector531[Blkcor][Tel] = new G4PVPlacement(Rot_53Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta53))*TMath::Cos(Blkphi53+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Sin(Blktheta53))*TMath::Sin(Blkphi53+(RoDi*pi*Blkcor/2))*mm, DistanceSi1*(TMath::Cos(Blktheta53))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector531[Blkcor][Tel],//theta == 0.0147
			silicon_531[Blkcor][Tel],
			logicWorld,
			false,
			silicon_531Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_532[Blkcor][Tel] = new G4Box(silicon_532[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector532[Blkcor][Tel] = new G4LogicalVolume(Silicon_532[Blkcor][Tel] ,Silicon_mat ,silicon_532[Blkcor][Tel]);

		PhysicalDetector532[Blkcor][Tel] = new G4PVPlacement(Rot_53Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta53))*TMath::Cos(Blkphi53+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Sin(Blktheta53))*TMath::Sin(Blkphi53+(RoDi*pi*Blkcor/2))*mm, DistanceSi2*(TMath::Cos(Blktheta53))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector532[Blkcor][Tel],
			silicon_532[Blkcor][Tel],
			logicWorld,
			false,
			silicon_532Num[Blkcor][Tel],
			overlap);

		//CsI crystal5		  
		CsIcrystal53[Blkcor][Tel] = new G4Trd(csicrystal53[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector533[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal53[Blkcor][Tel] ,scintillator_mat ,csicrystal53[Blkcor][Tel]);

		PhysicalDetector533[Blkcor][Tel] = new G4PVPlacement(Rot_53Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta53))*TMath::Cos(Blkphi53+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Sin(Blktheta53))*TMath::Sin(Blkphi53+(RoDi*pi*Blkcor/2))*mm, DistanceCsI*(TMath::Cos(Blktheta53))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector533[Blkcor][Tel],
			csicrystal53[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal53Num[Blkcor][Tel],
			overlap);
	  }
	  {//4
		Silicon_541[Blkcor][Tel] = new G4Box(silicon_541[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector541[Blkcor][Tel] = new G4LogicalVolume(Silicon_541[Blkcor][Tel],Silicon_mat,silicon_541[Blkcor][Tel]);
		PhysicalDetector541[Blkcor][Tel] = new G4PVPlacement(Rot_54Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi1*(TMath::Sin(Blktheta54))*TMath::Cos((Blkphi54+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Sin(Blktheta54))*TMath::Sin((Blkphi54+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi1*(TMath::Cos(Blktheta54))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector541[Blkcor][Tel],//theta == 0.0147
			silicon_541[Blkcor][Tel],
			logicWorld,
			false,
			silicon_541Num[Blkcor][Tel],
			overlap);

		//Si2		 
		Silicon_542[Blkcor][Tel] = new G4Box(silicon_542[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector542[Blkcor][Tel] = new G4LogicalVolume(Silicon_542[Blkcor][Tel] ,Silicon_mat ,silicon_542[Blkcor][Tel]);

		PhysicalDetector542[Blkcor][Tel] = new G4PVPlacement(Rot_54Si[Blkcor][Tel],
			G4ThreeVector(DistanceSi2*(TMath::Sin(Blktheta54))*TMath::Cos((Blkphi54+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Sin(Blktheta54))*TMath::Sin((Blkphi54+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceSi2*(TMath::Cos(Blktheta54))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector542[Blkcor][Tel],
			silicon_542[Blkcor][Tel],
			logicWorld,
			false,
			silicon_542Num[Blkcor][Tel],
			overlap);

		//CsI crystal5		  
		CsIcrystal54[Blkcor][Tel] = new G4Trd(csicrystal54[Blkcor][Tel],(20*mm)/2,(csibackL*mm)/2,(20*mm)/2,(csibackL*mm)/2,(100.*mm)/2);
		logicalDetector543[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal54[Blkcor][Tel] ,scintillator_mat ,csicrystal54[Blkcor][Tel]);

		PhysicalDetector543[Blkcor][Tel] = new G4PVPlacement(Rot_54Si[Blkcor][Tel],
			G4ThreeVector(DistanceCsI*(TMath::Sin(Blktheta54))*TMath::Cos((Blkphi54+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Sin(Blktheta54))*TMath::Sin((Blkphi54+(RoDi*pi*Blkcor/2)/1.0))*mm, DistanceCsI*(TMath::Cos(Blktheta54))*mm).rotate((RoDiQ*90*(Tel))*degree,G4ThreeVector(TMath::Sin(OUTERtheta)*TMath::Cos(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Sin(OUTERtheta)*TMath::Sin(OUTERphi+(RoDi*pi*Blkcor/2)),TMath::Cos(OUTERtheta))),
			logicalDetector543[Blkcor][Tel],
			csicrystal54[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal54Num[Blkcor][Tel],
			overlap);
	  }
	  //----------------------------------------------------------------------------------------------------------------------------------
	  /*opTubeSurface[Blkcor][Tel] = new G4OpticalSurface(TubeSurface[Blkcor][Tel]);
		opTubeSurface[Blkcor][Tel]->SetType(dielectric_metal);
		opTubeSurface[Blkcor][Tel]->SetFinish(polished);
		opTubeSurface[Blkcor][Tel]->SetModel(unified);

		new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector13[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector23[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector33[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector43[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);

		new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector413[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector423[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector433[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector443[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);

		new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector513[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector523[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector533[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
		new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector543[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);

		const G4int num = 2;
		G4double ephoton[num] = {2.034*eV, 4.136*eV};

	  //Optical alu rap Surface
	  G4double refractiveIndex[num] = {0.48, 1.55};
	  G4double specularLobe[num]    = {0.3, 0.3};
	  G4double specularSpike[num]   = {0.2, 0.2};
	  G4double backScatter[num]     = {0.2, 0.2};
	  G4double reflecitivity[num]   = {0.90, 0.92};

	  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

	  myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
	  myST1->AddProperty("REFLECITIVITY",         ephoton, reflecitivity,   num);
	  myST1->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
	  myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
	  myST1->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);

	  G4cout << "Tube Surface G4MaterialPropertiesTable" << G4endl;
	  myST1->DumpTable();

	  opTubeSurface[Blkcor][Tel]->SetMaterialPropertiesTable(myST1);*/

	  Silicon_11VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_11VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector11[Blkcor][Tel] -> SetVisAttributes (Silicon_11VisAtt[Blkcor][Tel]);

	  Silicon_12VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_12VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector12[Blkcor][Tel] -> SetVisAttributes (Silicon_12VisAtt[Blkcor][Tel]);

	  CsIcrystal1VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal1VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector13[Blkcor][Tel] -> SetVisAttributes (CsIcrystal1VisAtt[Blkcor][Tel]);
	  //
	  Silicon_21VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_21VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector21[Blkcor][Tel] -> SetVisAttributes (Silicon_21VisAtt[Blkcor][Tel]);

	  Silicon_22VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_22VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector22[Blkcor][Tel] -> SetVisAttributes (Silicon_22VisAtt[Blkcor][Tel]);

	  CsIcrystal2VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal2VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector23[Blkcor][Tel] -> SetVisAttributes (CsIcrystal2VisAtt[Blkcor][Tel]);
	  //
	  Silicon_31VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_31VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector31[Blkcor][Tel] -> SetVisAttributes (Silicon_31VisAtt[Blkcor][Tel]);

	  Silicon_32VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_32VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector32[Blkcor][Tel] -> SetVisAttributes (Silicon_32VisAtt[Blkcor][Tel]);

	  CsIcrystal3VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal3VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector33[Blkcor][Tel] -> SetVisAttributes (CsIcrystal3VisAtt[Blkcor][Tel]);
	  //		
	  Silicon_41VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_41VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector41[Blkcor][Tel] -> SetVisAttributes (Silicon_41VisAtt[Blkcor][Tel]);

	  Silicon_42VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_42VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector42[Blkcor][Tel] -> SetVisAttributes (Silicon_42VisAtt[Blkcor][Tel]);

	  CsIcrystal4VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal4VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector43[Blkcor][Tel] -> SetVisAttributes (CsIcrystal4VisAtt[Blkcor][Tel]);
	  //////////////////////////////////////
	  Silicon_411VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_411VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector411[Blkcor][Tel] -> SetVisAttributes (Silicon_411VisAtt[Blkcor][Tel]);

	  Silicon_412VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_412VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector412[Blkcor][Tel] -> SetVisAttributes (Silicon_412VisAtt[Blkcor][Tel]);

	  CsIcrystal41VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal41VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector413[Blkcor][Tel] -> SetVisAttributes (CsIcrystal41VisAtt[Blkcor][Tel]);
	  //
	  Silicon_421VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_421VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector421[Blkcor][Tel] -> SetVisAttributes (Silicon_421VisAtt[Blkcor][Tel]);

	  Silicon_422VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_422VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector422[Blkcor][Tel] -> SetVisAttributes (Silicon_422VisAtt[Blkcor][Tel]);

	  CsIcrystal42VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal42VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector423[Blkcor][Tel] -> SetVisAttributes (CsIcrystal42VisAtt[Blkcor][Tel]);
	  //
	  Silicon_431VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_431VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector431[Blkcor][Tel] -> SetVisAttributes (Silicon_431VisAtt[Blkcor][Tel]);

	  Silicon_432VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_432VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector432[Blkcor][Tel] -> SetVisAttributes (Silicon_432VisAtt[Blkcor][Tel]);

	  CsIcrystal43VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal43VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector433[Blkcor][Tel] -> SetVisAttributes (CsIcrystal43VisAtt[Blkcor][Tel]);
	  //		
	  Silicon_441VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_441VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector441[Blkcor][Tel] -> SetVisAttributes (Silicon_441VisAtt[Blkcor][Tel]);

	  Silicon_442VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_442VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector442[Blkcor][Tel] -> SetVisAttributes (Silicon_442VisAtt[Blkcor][Tel]);

	  CsIcrystal44VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal44VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector443[Blkcor][Tel] -> SetVisAttributes (CsIcrystal44VisAtt[Blkcor][Tel]);
	  //////////////////////////////
	  Silicon_511VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_511VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector511[Blkcor][Tel] -> SetVisAttributes (Silicon_511VisAtt[Blkcor][Tel]);

	  Silicon_512VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_512VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector512[Blkcor][Tel] -> SetVisAttributes (Silicon_512VisAtt[Blkcor][Tel]);

	  CsIcrystal51VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal51VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector513[Blkcor][Tel] -> SetVisAttributes (CsIcrystal51VisAtt[Blkcor][Tel]);
	  //
	  Silicon_521VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_521VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector521[Blkcor][Tel] -> SetVisAttributes (Silicon_521VisAtt[Blkcor][Tel]);

	  Silicon_522VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_522VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector522[Blkcor][Tel] -> SetVisAttributes (Silicon_522VisAtt[Blkcor][Tel]);

	  CsIcrystal52VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal52VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector523[Blkcor][Tel] -> SetVisAttributes (CsIcrystal52VisAtt[Blkcor][Tel]);
	  //
	  Silicon_531VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_531VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector531[Blkcor][Tel] -> SetVisAttributes (Silicon_531VisAtt[Blkcor][Tel]);

	  Silicon_532VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_532VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector532[Blkcor][Tel] -> SetVisAttributes (Silicon_532VisAtt[Blkcor][Tel]);

	  CsIcrystal53VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal53VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector533[Blkcor][Tel] -> SetVisAttributes (CsIcrystal53VisAtt[Blkcor][Tel]);
	  //		
	  Silicon_541VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_541VisAtt[Blkcor][Tel]-> SetForceWireframe(false);
	  logicalDetector541[Blkcor][Tel] -> SetVisAttributes (Silicon_541VisAtt[Blkcor][Tel]);

	  Silicon_542VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_542VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector542[Blkcor][Tel] -> SetVisAttributes (Silicon_542VisAtt[Blkcor][Tel]);

	  CsIcrystal54VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal54VisAtt[Blkcor][Tel] -> SetForceWireframe(false);
	  logicalDetector543[Blkcor][Tel] -> SetVisAttributes (CsIcrystal54VisAtt[Blkcor][Tel]);
	}
  }
  //G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
  //if (opticalSurface) opticalSurface->DumpInfo();
  // Generate & Add Material Properties Table attached to the optical surfaces
  //runManager -> SetSensitiveDetector(pvp);
  //return physWorld;
}
