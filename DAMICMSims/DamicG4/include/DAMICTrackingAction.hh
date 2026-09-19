#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "DAMICUtils.hh"

#include "G4UserTrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"

#include "globals.hh"

class DetectorConstruction;

class DAMICTrackingAction : public G4UserTrackingAction 
{
    
    public:
        DAMICTrackingAction();
        ~DAMICTrackingAction();
        
        virtual void PreUserTrackingAction(const G4Track*);
        virtual void PostUserTrackingAction(const G4Track*);

        //Clear the maps for a new event
        void ResetMaps();

        inline G4int GetParentID(G4int child){ return _parentMap[child]; }
        inline G4int GetPDGbyTrackID(G4int ID){ return _pdgMap[ID]; }

    private:
        //What trackID <value> is parent to trackID <key> ?
        std::map<G4int,G4int> _parentMap;
        // trackID <key>, PDG<value>
        std::map<G4int,G4int> _pdgMap;

        //Store info to fill PartInfo at the PostUserTrackingAction
        //      -- last volume of the particle should also be informed
        //  FROM PreUserTrack:
        G4int _trackid; //ParticlePDG
        G4int _particlePDG;
        G4double _partEKin;
        G4ThreeVector _partPosition;
        G4int _partCreVol;
        G4int _creatProcID;

};

#endif
