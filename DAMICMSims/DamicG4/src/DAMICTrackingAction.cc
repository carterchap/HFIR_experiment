#include "DAMICTrackingAction.hh"
#include "DAMICTrackInformation.hh"
#include "DAMICUtils.hh"
#include "DAMICAnalysisManager.hh"

#include "G4TrackingManager.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#include "G4SystemOfUnits.hh"


DAMICTrackingAction::DAMICTrackingAction():G4UserTrackingAction()
{
}

DAMICTrackingAction::~DAMICTrackingAction()
{
}

void DAMICTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    //Method only called once pr. track - no need to check if already added
    _trackid = aTrack->GetTrackID();

    _particlePDG = aTrack->GetParticleDefinition()->GetPDGEncoding();
    _partEKin = aTrack->GetKineticEnergy();
    _partPosition = aTrack->GetPosition();

    G4double vTime = aTrack->GetGlobalTime()/s;
    G4String VolumeName = aTrack->GetVolume()->GetName();
    _partCreVol = GetVolumeID(aTrack->GetVolume()->GetName());

    // fill maps: trackid--parentid 
    _parentMap[_trackid] = aTrack->GetParentID();
    // and trackid--pdg_code
    _pdgMap[_trackid] = _particlePDG;

    G4int last_replica = -1;
    
//    if(_trackid == 1)
    if (aTrack->GetParentID() == 0) 
    {
        fpTrackingManager->SetStoreTrajectory(true);        
        fpTrackingManager->SetUserTrackInformation(new DAMICTrackInformation(_particlePDG,_partEKin,_trackid,
                    VolumeName,_partPosition,vTime,last_replica));

    }
        
    
    //Access to Creator Process Name, and convert into an ID value
    G4String creatProc="";
    if(aTrack->GetCreatorProcess())
    {
         creatProc=aTrack->GetCreatorProcess()->GetProcessName();
    }
    else if(aTrack->GetParentID() == 0)
    {
        creatProc = "Generator";
    }
    else
    {
        creatProc = "unknown";
    }
    _creatProcID = GetProcessID(creatProc);
}


void DAMICTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if(secondaries)
    {
        DAMICTrackInformation* info = (DAMICTrackInformation*)(aTrack->GetUserInformation());
        size_t nSeco = secondaries->size();
        if(nSeco>0)
        {
            for(size_t i=0; i < nSeco; i++)
            {
                G4int ProcessCreatorSub = -1;
                G4int ProcessCreatorType = -1;
                G4int IDTrack = (*secondaries)[i]->GetTrackID();
                G4String VolumeName = "NULL";
                //G4double energykin = (*secondaries)[i]->GetKineticEnergy();
                if (IDTrack != 1)
                {
                    ProcessCreatorSub = (*secondaries)[i]->GetCreatorProcess()->GetProcessSubType();
                    ProcessCreatorType = (*secondaries)[i]->GetCreatorProcess()->GetProcessType();
                }
                
                //G4bool ElecIo = (ProcessCreatorSub == 2 && ProcessCreatorType == 2 && energykin < 10*keV);
                G4bool ElecIoSup = (ProcessCreatorSub == 2 && ProcessCreatorType == 2);
                if (!ElecIoSup)
                {
                    VolumeName = (*secondaries)[i]->GetVolume()->GetName();
                }
                
                G4int ParticlePDG = (*secondaries)[i]->GetParticleDefinition()->GetPDGEncoding();
                
                G4double Energy = (*secondaries)[i]->GetKineticEnergy();
                
                G4ThreeVector Position = (*secondaries)[i]->GetPosition();
                
                G4double vTime = (*secondaries)[i]->GetGlobalTime()/s;
                
               	G4int _last_replica_id = -1;

                DAMICTrackInformation* infoNew = new DAMICTrackInformation(info, ParticlePDG, Energy, IDTrack, VolumeName, Position, vTime,_last_replica_id);
                (*secondaries)[i]->SetUserInformation(infoNew);
            }
        }
    }

    G4StepPoint* thePostPoint = aTrack->GetStep()->GetPostStepPoint();
    G4String endproc("UserLimit");
    if(thePostPoint->GetProcessDefinedStep()!=0)
    {
        endproc = thePostPoint->GetProcessDefinedStep()->GetProcessName();    
    }
    
    G4int lastVolID = 0;
    G4int lastProcessID = 999;
    DAMICAnalysisManager * man = DAMICAnalysisManager::GetInstance();
    //if(man->GetPartInfoVectorSize() > 0)
    //{
        lastProcessID = GetProcessID(endproc);
        if(thePostPoint->GetPhysicalVolume()){
        	lastVolID = GetVolumeID(thePostPoint->GetPhysicalVolume()->GetName());
        }
    //}
    man->FillPartInfo(_particlePDG,_trackid,_parentMap[_trackid], _pdgMap[_parentMap[_trackid]],
            _partCreVol, lastVolID, _creatProcID, lastProcessID, 
            _partEKin/CLHEP::keV,_partPosition/CLHEP::mm);
 
}

void DAMICTrackingAction::ResetMaps()
{
    _parentMap.clear();
    _pdgMap.clear();
}

