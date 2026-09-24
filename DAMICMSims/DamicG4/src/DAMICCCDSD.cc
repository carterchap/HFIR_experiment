#include "DAMICCCDSD.hh"
#include "DAMICCCDHit.hh"
#include "DAMICRunAction.hh"
#include "DAMICTrackInformation.hh"
#include "DAMICAnalysisManager.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"

#include "Randomize.hh"
#include "TMath.h"

#include <cmath>

DAMICCCDSD::DAMICCCDSD(G4String name):
    G4VSensitiveDetector(name),
    _hitsCollection(nullptr),
    _hitID(-1)
{
    this->collectionName.insert("CCDColl");
}

DAMICCCDSD::~DAMICCCDSD()
{
}

void DAMICCCDSD::Initialize(G4HCofThisEvent* hce)
{
    _hitsCollection = new DAMICCCDHitsCollection(this->SensitiveDetectorName,this->collectionName[0]);

    if(_hitID<0)
    {
        _hitID = G4SDManager::GetSDMpointer()->GetCollectionID(this->collectionName[0]);
    }
    hce->AddHitsCollection(_hitID,_hitsCollection);
}

G4bool DAMICCCDSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    // FIXME MISSING DOCUMENTATION
    const G4double EnergyDep = step->GetTotalEnergyDeposit();
    // Nothing to do
    if(EnergyDep <=0 )
    {
        return true;
    }
    
    G4Track *track=step->GetTrack();
    const G4int partTrackID=track->GetTrackID();
    const G4int partStepID=track->GetCurrentStepNumber();
    const G4int partParentID=track->GetParentID();

    const G4StepPoint* thePrePoint  = step->GetPreStepPoint();
    const G4StepPoint* thePostPoint = step->GetPostStepPoint();

    G4String endproc = "UserLimit";
    G4String startproc = "UserGen";
    if(thePostPoint->GetProcessDefinedStep()!=0)
    {
        endproc = thePostPoint->GetProcessDefinedStep()->GetProcessName();
    }
    
    if(thePrePoint->GetProcessDefinedStep()!=0)
    {
        startproc =  thePrePoint->GetProcessDefinedStep()->GetProcessName();
    }
     
    const G4int partPDG = track->GetParticleDefinition()->GetPDGEncoding();
    const G4String nameProcess = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    const G4TouchableHandle theTouchablePre = step->GetPreStepPoint()->GetTouchableHandle();
    const G4TouchableHandle theTouchable = step->GetPostStepPoint()->GetTouchableHandle();

    const G4String volumeNname = theTouchable->GetVolume()->GetLogicalVolume()->GetName();

    const G4String volumeNnname = theTouchablePre->GetVolume()->GetLogicalVolume()->GetName();
   
    auto CCDSensorLV = theTouchablePre->GetVolume()->GetLogicalVolume();

    const G4ThreeVector preWorldPosition = step->GetPreStepPoint()->GetPosition();
    const G4ThreeVector postWorldPosition = step->GetPostStepPoint()->GetPosition();
    const G4ThreeVector localPrePosition = theTouchablePre->GetHistory()->GetTopTransform().TransformPoint(preWorldPosition);
    const G4ThreeVector localPostPosition = theTouchablePre->GetHistory()->GetTopTransform().TransformPoint(postWorldPosition);
    const G4double partStepLength = (postWorldPosition - preWorldPosition).mag();
    const G4TouchableHistory* touchable = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* motherPhysical = nullptr; // mother


    // Now the CCD Sensitive part (CCDSensor) is a daughter of a replicated volume which is called
    // CCDModulePair (if this name change on the GDML it should be propagated here!!!)
    // The copy number of the CCD is given by its mother, i.e. the CCDModulePair_PV
    auto CCDModulePVFound=false;
    G4int depthToGetCCDNum= 0;
    while(not CCDModulePVFound)
    {
        motherPhysical = touchable->GetVolume(depthToGetCCDNum);

        const auto motherName = motherPhysical->GetName();

        if(motherName=="CCDModulePair_PV" || motherName.find("CCDModule")==0 ||
                motherName=="ModulePV" || motherName=="ModulePVExt1")
        {
            CCDModulePVFound=true;
            break;
        }
        ++depthToGetCCDNum;
        // FIXME depthToGetCCDNum>10 should be eliminated, if the first condition works (i do not
        // see why can not work)
        if(motherName=="World_PV" || depthToGetCCDNum>10)
        {
            G4cout << "ERROR [DAMICCCDSD::ProcessHits] The physical volume 'CCDModulePair_PV' is "
                << "required to get the CCD ID number, but it is not found." << G4endl;
            exit(-1);
        }
    }
    // Start CCD id number from 1
    G4int copyNo = theTouchablePre->GetCopyNumber(depthToGetCCDNum);
    G4int copyNoCCD = theTouchablePre->GetCopyNumber(1);
   
    // if the volume used on the replica is made of two CCD, which are rotated
    // one wrt the other one 90 degree, the CCD number doesn't correspond to the real one
    // The replica volume contains two CCDs with the same ID number, however we can know which CCD
    // is due to the rotation property
    // In order to extract the correct CCD number id, we should access to the number ID of the
    // mother volume (replicate volume), as well as to the rotation parameter which allow us to know
    // if the CCD is rotated or not. Then, the number of the CCD is given by
    //    2*N_replica_mother * (n-1) 
    //    where n is 0 or 1 if the volume is rotated or not, resplectively
    //    and N_replica_mother the replica number given by the volume CCDModulePair

    const G4RotationMatrix* IsRotated = touchable->GetRotation(0);
   
    /*if(motherPhysical->GetName()=="CCDModulePair_PV")
    {
	if(copyNo ==40 || copyNo ==39 ){
	G4cout<< "is rotated: "<< (*IsRotated) << G4endl;
	G4cout<< "copy n: "<< copyNo << G4endl;
	}

        if((*IsRotated).isIdentity())
        {
            copyNo = 2*copyNo-1;
        }
        else
        {
            copyNo = 2*copyNo;
        }
    }*/

    if(motherPhysical->GetName()=="CCDModulePair_PV")
    {

	/*if(copyNo ==41 ){
	G4cout<< "copy n: "<< copyNo << G4endl;
	G4cout<< "copy n CCD: "<< copyNoCCD << G4endl;
	}*/

        if(copyNoCCD == 1)
        {
            copyNo = 4*copyNo-3;
        }
	if(copyNoCCD == 2)
        {
            copyNo = 4*copyNo-2;
        }
	if(copyNoCCD == 3)
        {
            copyNo = 4*copyNo-1;
        }
	if(copyNoCCD == 4)
        {
            copyNo = 4*copyNo;
        }
	
	//if(copyNo ==41*2 || copyNo == (41*2-1)  ){G4cout<< "FINAL copy n: "<< copyNo << G4endl;}
    }
  

    const G4int EventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    // Global time : track time wrt event, local time = wrt track time
    const G4double Time = step->GetTrack()->GetGlobalTime();
    const G4double local_time = thePostPoint->GetLocalTime();

    // ---------------------------------------------------
    // Get info to fill the hit
    //
    // GLOBAL POSITION 
    const G4ThreeVector Pos = thePostPoint->GetPosition();
    //G4cout << " copy num " << copyNo << " gzpos " <<  Pos << G4endl;

    //This gives the positon of the edep point inside the current volume
    const G4ThreeVector localPosition = theTouchablePre->GetHistory()->GetTopTransform().TransformPoint(Pos);
    
    //As the CCD silicon is now the only sensitive part, the translation (see DAMIC@SNOLAB code)
    //is not needed anymore
    //Move all coordinates to be positive
    auto shape = dynamic_cast<G4Box*>(CCDSensorLV->GetSolid());
    const G4double shapeHalfSizeX = shape->GetXHalfLength();
    const G4double shapeHalfSizeY = shape->GetYHalfLength();
    const G4double shapeHalfSizeZ = shape->GetZHalfLength();
    
    const G4double XCoord = localPosition.getX() + shapeHalfSizeX;
    const G4double YCoord = localPosition.getY() + shapeHalfSizeY;
    const G4double ZCoord = localPosition.getZ() + shapeHalfSizeZ;
    const G4ThreeVector PosRef(XCoord, YCoord, ZCoord);

    // TO BE DEPRECATED
    G4int _numMaxPix = TMath::Nint((2.0*shapeHalfSizeY)/CLHEP::mm/_pixSize);
    const G4int pixID = TMath::Nint(PosRef.x()/mm/_pixSize) + _numMaxPix*TMath::Nint(PosRef.y()/mm/_pixSize);
    
    DAMICCCDHit* hit = new DAMICCCDHit(copyNo);
    hit->SetPos(PosRef/CLHEP::mm);
    hit->SetGlobPos(Pos/CLHEP::mm);
    hit->SetPDGNumber(partPDG);
    hit->SetEnergyDeposit(EnergyDep);
    hit->SetMotherID(partParentID);
    hit->SetTrackID(partTrackID);
    hit->SetTime(Time);
    hit->SetLocalTime(local_time);
    hit->SetEventID(EventID);
    hit->SetCCDNum(copyNo);
    hit->SetPixID(pixID); // TO BE DEPRECATED
    _hitsCollection->insert(hit);

    // CCDOut tree is going to use this hit. DONE
    // ---------------------------------------------------

    // ---------------------------------------------------
    // Get info to fill the TRACK OUT branches
    //
    const G4int partVolID = GetVolumeID(volumeNnname);  //volumeNname; //to be fixed with the id of the volume or remove
    //if (partVolID<0) G4cout << volumeNnname << " not found " << G4endl;
    const G4double partDeltaE = (thePostPoint->GetKineticEnergy() - thePrePoint->GetKineticEnergy());
    const G4double partEnergy = thePrePoint->GetKineticEnergy();
    const G4double partStartProcess = GetProcessID(startproc);
    const G4double partEndProcess = GetProcessID(endproc);

    //if(partStartProcess == 1000 && (partEnergy/CLHEP::keV<1.8 && partEnergy/CLHEP::keV>1.5)) G4cout << "pdg: "<< partPDG << " energy[keV]: "<< partEnergy/CLHEP::keV <<" process: "<< track->GetCreatorModelName() << G4endl;
    //particle based tree - one entry per event, vector is per particle et step 
    DAMICAnalysisManager* man = DAMICAnalysisManager::GetInstance();
    man->FillTrackOut(
            partPDG, partTrackID, partStepID, partParentID,
            copyNo,
            partVolID, partStartProcess, partEndProcess, 
            partStepLength/CLHEP::mm,
            partEnergy/CLHEP::eV, partDeltaE/CLHEP::eV, EnergyDep/CLHEP::eV,
            //thePrePoint->GetPosition().x()/CLHEP::mm,
            //thePrePoint->GetPosition().y()/CLHEP::mm,
            //thePrePoint->GetPosition().z()/CLHEP::mm,
            XCoord/CLHEP::mm, YCoord/CLHEP::mm, ZCoord/CLHEP::mm,
            Time/CLHEP::s);
    
    // To store each Edep:
    /*man->FillCCDOut(
            partPDG,partTrackID,partParentID,copyNo,
            PosRef.x()/CLHEP::mm, PosRef.y()/CLHEP::mm,PosRef.z()/CLHEP::mm,
            Pos.x()/CLHEP::mm,Pos.y()/CLHEP::mm,Pos.z()/CLHEP::mm,
            EnergyDep/CLHEP::eV,Time/CLHEP::s,local_time/CLHEP::s);
    */
    return true;
}

