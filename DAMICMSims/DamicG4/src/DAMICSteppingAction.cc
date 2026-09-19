#include "DAMICSteppingAction.hh"
#include "DAMICEventAction.hh"
#include "DAMICRunAction.hh"
#include "DAMICUtils.hh"
#include "DAMICAnalysisManager.hh"
#include "DAMICTrackInformation.hh"
#include <G4VUserTrackInformation.hh>
#include "G4Step.hh"
#include "DAMICAnalysisManager.hh"
#include <G4VProcess.hh>
#include "G4RunManager.hh"

DAMICSteppingAction::DAMICSteppingAction(DAMICRunAction *ra) : 
    G4UserSteppingAction(),
    _runAction(ra)
{
}

DAMICSteppingAction::~DAMICSteppingAction()
{
}

void DAMICSteppingAction::UserSteppingAction(const G4Step* aStep)
{

    DAMICAnalysisManager * man = DAMICAnalysisManager::GetInstance();    
    G4Track *track=aStep->GetTrack();
    G4int TrackID=track->GetTrackID();
    const G4int StepID= track->GetCurrentStepNumber();
    G4int counter = 0;
    G4int id_replica_last = -1;
    
    int volid = GetVolumeID(track->GetVolume()->GetName());

   //G4cout << "track def: "<< track->GetParticleDefinition()->GetParticleName() <<G4endl;
    

    if (!aStep->GetPostStepPoint()) { 
      return;
	  }

    if (!(aStep->GetPostStepPoint()->GetProcessDefinedStep())) { 
      return;
	  }

    G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent() -> GetEventID(); 
    
//G4cout <<aStep->GetTrack()->GetDefinition()->GetPDGEncoding()  << " pre point:"<< aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() <<"post point " << thePostPoint->GetPhysicalVolume()->GetName() << " proc: "<< thePostPoint->GetProcessDefinedStep()->GetProcessName()<< G4endl ;

    if (thePostPoint->GetProcessDefinedStep()->GetProcessName() == "nFission")
      G4cout << "fission ID" << TrackID << " kin " << thePostPoint->GetKineticEnergy()/CLHEP::eV <<  G4endl; 
    
    G4double time = (thePostPoint->GetGlobalTime() - _runAction->triggerTime);

    
    
    if(TrackID==1 && StepID==1)
    {
        _runAction->primary_volid = volid;
        _runAction->triggerTime = track->GetGlobalTime()/CLHEP::s;
       	counter = 0;
    }

    G4String Volume = "";
    G4String Material = "";

    if (aStep->GetPreStepPoint()) {
      Material = aStep->GetPreStepPoint()->GetMaterial()->GetName();
      Volume = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    }
    //G4cout << " Volume, material " << Volume << " " << Material << G4endl;
    
    G4String NextVolume = "";
    G4String NextMaterial = "";  

    if (!aStep->GetTrack()->GetNextVolume()) {
      
      return;}
    else {
       NextVolume = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();

      NextMaterial = aStep->GetPostStepPoint()->GetMaterial()->GetName();

    }

  
    //extract shieldingOut infos
    if (man->IsShieldOutON()) {
      //G4cout << "SHIELD OUT ON!!!"<< G4endl;
     int id_replica = -1;
     int id_nextreplica = -1;				      
    
    std::string delimiter = "_PV";
    std::string delimiter2 = "_";

    std::string token = "";
    token =  Volume.substr(0, Volume.find(delimiter));

    std::size_t pos = token.find(delimiter2);
    
    //G4cout << Volume << " " << NextVolume << G4endl; 

    if (pos != std::string::npos && (token.substr(0, pos)=="AncLeadShield" || token.substr(0, pos)=="PolyShield" || token.substr(0, pos)=="ColdCopper" ) ) {//||token.substr(0, pos)=="Cryo" 
      if (token.substr(pos+1)!="")  id_replica = stoi(token.substr(pos+1));
      if (token.substr(0, pos)=="AncLeadShield") id_replica+=500;
      if (token.substr(0, pos)=="PolyShield") id_replica+=900;
      if (token.substr(0, pos)=="ColdCopper") id_replica+=10;
      if (token.substr(0, pos)=="Cryo") id_replica+=100;
    
    } //end if npos
    else {
      id_replica = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();

      if (token=="PbLeftRightModule" || token=="PbTopBotModule" || token=="PbFrontBackModule")
	id_replica+=500; 
      else if (token=="PolyLeftRightModule" || token=="PolyTopBotModule"  || token=="PolyFrontBackModule" )
        id_replica+=900;

      if (token=="ExtLeadShield") id_replica+=500; 
      if (token=="ExtLeadShield_Ancient") id_replica+=550;
      else if (token=="ExtPolyShield" ) id_replica+=900;
      else if (token=="ExtBoratedPolyShield" ) id_replica+=950;	
      else if (token=="LabWall" ) id_replica+=2000;
      else if (token=="Lab" ) id_replica+=1000;
	
    }
 

    std::string token2 = "";
    token2 = NextVolume.substr(0, NextVolume.find(delimiter));
    std::size_t pos2 = token2.find(delimiter2);


    if (pos2 != std::string::npos && (token2.substr(0, pos2)=="AncLeadShield" || token2.substr(0, pos2)=="PolyShield" || token2.substr(0, pos2)=="ColdCopper"  ) ) {//|| token2.substr(0, pos2)=="Cryo"
      if (token2.substr(pos2+1)!="") id_nextreplica = stoi(token2.substr(pos2+1));
      if (token2.substr(0, pos2)=="AncLeadShield") id_nextreplica+=500;
      if (token2.substr(0, pos2)=="PolyShield") id_nextreplica+=900;
      if (token2.substr(0, pos2)=="ColdCopper") id_nextreplica+=10;
      if (token2.substr(0, pos2)=="Cryo") id_nextreplica+=100;
    
    }  //end npos2
    else {
      id_nextreplica = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
      if (token2=="PbLeftRightModule" || token2=="PbTopBotModule" || token2=="PolyFrontBackModule")
	id_nextreplica+=500;
      
      if (token2=="PolyLeftRightModule" || token2=="PolyTopBotModule"  || token2=="PolyFrontBackModule")
	id_nextreplica+=900;

      if (token2=="ExtLeadShield") id_nextreplica+=500; 
      if (token2=="ExtLeadShield_Ancient") id_nextreplica+=550;
      else if (token2=="ExtPolyShield" ) id_nextreplica+=900;
      else if (token2=="ExtBoratedPolyShield" ) id_nextreplica+=950;	
      else if (token2=="LabWall" ) id_nextreplica+=2000;
      else if (token2=="Lab" ) id_nextreplica+=1000;

    }

    DAMICTrackInformation* info = (DAMICTrackInformation*)(track->GetUserInformation());

    if (id_replica>0 ) {
      id_replica_last = info->GetReplicaID();
      info->SetReplicaID(id_replica);
    
     //G4cout  << "Track, Step " << TrackID << ", " <<  StepID  <<  " pdg " << aStep->GetTrack()->GetDefinition()->GetPDGEncoding() << " replica " <<  id_replica << " last " << id_replica_last << " next " << id_nextreplica << " size track info " << info->GetSize() << G4endl;
    }

 
     const G4double StepLength = (aStep->GetPostStepPoint()->GetPosition() - aStep->GetPreStepPoint()->GetPosition()).mag();
  
    G4String creatProc0="";
    if(aStep->GetTrack()->GetCreatorProcess())
    {
         creatProc0=aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
    }
    else if(aStep->GetTrack()->GetParentID() == 0)
    {
        creatProc0 = "Generator";
    }
    else
    {
        creatProc0 = "unknown";
    }
    int creatProcID0 = GetProcessID(creatProc0);

     if(man->IsLastVolumeON()){
	//G4cout << "Selected Last Volume & Stop Volume: " << man->GetShieldingLastVolume() <<" & " << man->GetShieldingStopVolume() << G4endl;
	
	if(Volume==man->GetShieldingLastVolume() && NextVolume == man->GetShieldingStopVolume()){


		if( (id_nextreplica < id_replica && id_replica <= id_replica_last && id_replica>500 && id_nextreplica>0 && id_replica_last>0 &&  aStep->GetPostStepPoint()->GetKineticEnergy()!=0)|| (id_nextreplica < id_replica && id_replica <= id_replica_last && id_replica==500 && id_nextreplica>=0 && id_replica_last>0 &&  aStep->GetPostStepPoint()->GetKineticEnergy()!=0) ){


    //G4cout <<"event id "<< eventID << " trackid" << aStep->GetTrack()->GetTrackID()<< " volume " << token << " id " << id_replica << " next volume " << token2 << " id "<< id_nextreplica  << " Last "<< id_replica_last << G4endl;
		    
		    man->FillShieldingOut(
		    aStep->GetTrack()->GetDefinition()->GetPDGEncoding(),
		    aStep->GetTrack()->GetTrackID(), 
		    aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::eV, aStep->GetPreStepPoint()->GetWeight(), 
		    aStep->GetPostStepPoint()->GetPosition().x(),
		    aStep->GetPostStepPoint()->GetPosition().y(), aStep->GetPostStepPoint()->GetPosition().z(),
		    aStep->GetPostStepPoint()->GetMomentum().x(), aStep->GetPostStepPoint()->GetMomentum().y(),
		    aStep->GetPostStepPoint()->GetMomentum().z(), time/CLHEP::s, id_replica,StepID,StepLength,creatProcID0);
	       
		    id_replica_last =  id_nextreplica;  
		    counter++;

		    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		    
		}
	}

     }
     else{
	     //G4cout << "You didn't select Last Volume & Stop Volume: particles information saved for each layer, particles killed when they pass from lead0 to lab "<< G4endl;
	     if( (id_nextreplica < id_replica && id_replica <= id_replica_last && id_replica>500 && id_nextreplica>0 && id_replica_last>0 &&  aStep->GetPostStepPoint()->GetKineticEnergy()!=0)|| (id_nextreplica < id_replica && id_replica <= id_replica_last && id_replica==500 && id_nextreplica>=0 && id_replica_last>0 &&  aStep->GetPostStepPoint()->GetKineticEnergy()!=0) ){
	     
		 //G4cout <<"event id "<< eventID << " trackid " << aStep->GetTrack()->GetTrackID()<< " step id: "<< StepID  <<" volume " << token << " id " << id_replica << " next volume " << token2 << " id "<< id_nextreplica  << " Last "<< id_replica_last << G4endl;
		 
		    man->FillShieldingOut(
		    aStep->GetTrack()->GetDefinition()->GetPDGEncoding(),
		    aStep->GetTrack()->GetTrackID(), 
		    aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::eV, aStep->GetPostStepPoint()->GetWeight(), 
		    aStep->GetPostStepPoint()->GetPosition().x(),
		    aStep->GetPostStepPoint()->GetPosition().y(), aStep->GetPostStepPoint()->GetPosition().z(),
		    aStep->GetPostStepPoint()->GetMomentum().x(), aStep->GetPostStepPoint()->GetMomentum().y(),
		    aStep->GetPostStepPoint()->GetMomentum().z(), time/CLHEP::s, id_replica,StepID,StepLength,creatProcID0);
	       
		    id_replica_last =  id_nextreplica;  
		    counter++;

		   //aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		   if(  (id_replica==500 && id_nextreplica==0) ){ aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
	     }
     }


    

    }//end if shieldout ON

  
      
}
