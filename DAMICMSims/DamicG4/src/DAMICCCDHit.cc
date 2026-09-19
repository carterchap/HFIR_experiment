#include "DAMICCCDHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<DAMICCCDHit>* DAMICCCDHitAllocator;

DAMICCCDHit::DAMICCCDHit(G4int /*i */)
{
    fPos = G4ThreeVector(0,0,0);
    fgPos = G4ThreeVector(0,0,0);
    fTopBot = 0;
    fEnergyDeposit = 0;
    fPartIDElec = 0;
    fProcess = "NULL";
    //
    fPDGNumber = 0;
    fPartID = 0;
    fProdVolume = "NULL";
    fEnergy = 0;
    //Primary Nuc
    fTime = 0;
    fLocalTime=0;
    fPDGPrimaryNuc = 0;
    fEventID = 0;
    //Secondary Nuc
    fPDGSecondaryNuc = 0;
    //Part after Second
    fPDGPartASecond = 0;
    fEnergyPartASecond = 0;
    // CCD
    fCCDNum = 0;
}

DAMICCCDHit::~DAMICCCDHit()
{
}

void DAMICCCDHit::Draw()
{
}

void DAMICCCDHit::Print()
{
}
