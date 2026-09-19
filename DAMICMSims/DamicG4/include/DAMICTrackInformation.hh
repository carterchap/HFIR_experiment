#ifndef DAMICTrackInformation_h
#define DAMICTrackInformation_h 1

#include "globals.hh"
#include "G4VUserTrackInformation.hh"
#include <vector>
#include "G4ThreeVector.hh"

class DAMICTrackInformation: public G4VUserTrackInformation
{
    public:
        DAMICTrackInformation(G4int VPDGParticle, G4double VEnergyParticle, G4int VParticleID, 
                G4String VMaterialProd, G4ThreeVector VPosition, G4double  VTime, G4int vol_replica_id);

        DAMICTrackInformation(const DAMICTrackInformation* InfoTrack);

        DAMICTrackInformation(const DAMICTrackInformation* InfoTrack, G4int VPDGParticle, G4double 
                VEnergyParticle, G4int VParticleID, G4String VMaterialProd, G4ThreeVector VPosition, G4double  VTime, G4int vol_replica_id);

        ~DAMICTrackInformation();
        
        virtual void Print() const {};
        
        void AddPDGParticle(G4int val) { PDGParticle.push_back(val); }
        std::vector<G4int> GetPDGParticle() const { return PDGParticle; }
        
        void AddEnergyParticle(G4double val) { EnergyParticle.push_back(val); }
        std::vector<G4double> GetEnergyParticle() const { return EnergyParticle; }
        
        void AddParticleID(G4int val) { ParticleID.push_back(val); }
        std::vector<G4int> GetParticleID() const { return ParticleID; }
        
        void AddMaterialProd(G4String val) { MaterialProd.push_back(val) ; }
        std::vector<G4String> GetMaterialProd() const { return MaterialProd; }
        
        void AddStartPosition(G4ThreeVector val) {StartPosition.push_back(val);}
        std::vector<G4ThreeVector> GetStartPosition() const {return StartPosition;}
        
        void AddTime(G4double val){ Time.push_back(val);}
        std::vector<G4double> GetTime() const {return Time;}
        
        void AddReplicaID(G4int val) { last_replica_id.push_back(val);}
        std::vector<G4int> GetReplicaID() const {return last_replica_id; }
    
        void SetReplicaID(G4int ) ; 
        G4int GetReplicaID(); 
  	G4int GetSize(); 
  	
    private:
        std::vector<G4int> PDGParticle;
        std::vector<G4double> EnergyParticle;
        std::vector<G4int> ParticleID;
        std::vector<G4String> MaterialProd;
        std::vector<G4ThreeVector> StartPosition;
        std::vector<G4double> Time;
        std::vector<G4int> last_replica_id;

};

#endif
