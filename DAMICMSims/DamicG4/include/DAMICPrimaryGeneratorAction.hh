#ifndef DAMICPrimaryGeneratorAction_hh
#define DAMICPrimaryGeneratorAction_hh

#include <stdio.h>

#include "G4VUserPrimaryGeneratorAction.hh"
#include "DAMICParticleSource.hh"
#include "globals.hh"

using namespace std;

class G4Event;
class DAMICTrackingAction;

class DAMICPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
        DAMICPrimaryGeneratorAction();
        ~DAMICPrimaryGeneratorAction();
        
        void GeneratePrimaries(G4Event* anEvent);
        const long* GetEventSeeds() const  {return seeds;};
        G4double GetEnergyPrimary() const  {return energyPri;};

        void SetTrackingActionInstance(DAMICTrackingAction* trackingAction) {_trackingAction=trackingAction;};
 
    private:
        DAMICParticleSource* particleGun;
        long seeds[2];
        G4double energyPri;

        // needed to reset maps to PDG/ID parent
        DAMICTrackingAction* _trackingAction;
};

#endif
