#include "DAMICTrackInformation.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


DAMICTrackInformation::DAMICTrackInformation(G4int VPDGParticle, G4double VEnergyParticle, G4int VParticleID, G4String VMaterialProd, G4ThreeVector VPosition, G4double  VTime, G4int vol_replica_id):G4VUserTrackInformation()
{
  PDGParticle = {VPDGParticle};
  EnergyParticle = {VEnergyParticle};
  ParticleID = {VParticleID};
  MaterialProd = {VMaterialProd};
  StartPosition = {VPosition};
  Time = {VTime};
  last_replica_id = {vol_replica_id};

}

DAMICTrackInformation::DAMICTrackInformation(const DAMICTrackInformation* InfoTrack): G4VUserTrackInformation()
{
  PDGParticle = InfoTrack->PDGParticle;
  EnergyParticle = InfoTrack->EnergyParticle;
  ParticleID = InfoTrack->ParticleID;
  MaterialProd = InfoTrack->MaterialProd;
  StartPosition = InfoTrack->StartPosition;
  Time = InfoTrack->Time;
  last_replica_id = InfoTrack->last_replica_id;

}

DAMICTrackInformation::DAMICTrackInformation(const DAMICTrackInformation* InfoTrack, G4int VPDGParticle, G4double VEnergyParticle, G4int VParticleID, G4String VMaterialProd, G4ThreeVector VPosition, G4double  VTime, G4int vol_replica_id): G4VUserTrackInformation()
{
  PDGParticle = InfoTrack->PDGParticle;
  EnergyParticle = InfoTrack->EnergyParticle;
  ParticleID = InfoTrack->ParticleID;
  MaterialProd = InfoTrack->MaterialProd;
  StartPosition = InfoTrack->StartPosition;
  Time = InfoTrack->Time;
  last_replica_id = InfoTrack->last_replica_id;
  AddPDGParticle(VPDGParticle);
  AddEnergyParticle(VEnergyParticle);
  AddParticleID(VParticleID);
  AddMaterialProd(VMaterialProd);
  AddStartPosition(VPosition);
  AddTime(VTime);
}

void DAMICTrackInformation::SetReplicaID(G4int val) {
  int j =  last_replica_id.size();
  last_replica_id.at(j-1) = val; 
}


G4int DAMICTrackInformation::GetReplicaID() {
  int j =  last_replica_id.size();
  return last_replica_id.at(j-1); 
}


G4int DAMICTrackInformation::GetSize() {
  int j =  last_replica_id.size();
  return j;
}

DAMICTrackInformation::~DAMICTrackInformation()
{}
