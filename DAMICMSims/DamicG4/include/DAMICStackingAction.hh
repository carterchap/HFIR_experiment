#ifndef DAMICStackingAction_hh
#define DAMICStackingAction_hh

#include <stdio.h>

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "DAMICStackingActionMessenger.hh"
#include <vector>

class G4Track;

class DAMICStackingAction : public G4UserStackingAction
{
    public:
        DAMICStackingAction();
        ~DAMICStackingAction();
        
        void SetPartName(G4String );
        void SetPartEnergyKin(G4double);
        void SetPartEnergyTot(G4double);        
        G4double GetPartEnergyKin() const {return _partEnergyKin;}
        G4double GetPartEnergyTot() const {return _partEnergyTot;}
        G4String GetPartName () const {return _partName;}

        virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

        void SetLimitNucleus(G4int);
        G4int GetLimitNucleus() const {return _limitNucleus;}

        void AddParticlesKill(G4String);
        std::vector<G4int> GetParticlesKill() const {return _particlesKill;}
        
        void ResetNewRun();
        
        G4bool DoKill(G4int);
        G4bool DoKillNucleus(G4int);
        G4bool DoOneDecay(G4int, G4int);
        
        // Add function to unset the DoOneDecay process, to simulation the full decay chain
        void UnsetDoOneDecay(G4bool);
    
    private:
        G4String _partName;
        G4double _partEnergyKin;
        G4double _partEnergyTot;
        std::vector<G4int> _particlesKill;
        G4int _limitNucleus;
        DAMICStackingActionMessenger* _messenger;
        G4int _idpartPrimPartPDG;
        G4bool _oneDecay=true;
};

#endif
