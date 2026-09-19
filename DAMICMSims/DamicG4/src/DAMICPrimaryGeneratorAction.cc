#include "G4Event.hh"

#include "DAMICPrimaryGeneratorAction.hh"
#include "DAMICTrackingAction.hh"
#include "DAMICParticleSource.hh"

DAMICPrimaryGeneratorAction::DAMICPrimaryGeneratorAction():
    G4VUserPrimaryGeneratorAction(),
    _trackingAction(nullptr)
{

  particleGun = new DAMICParticleSource();
  energyPri=0;
  seeds[0] =-1;
  seeds[1] =-1;

}

DAMICPrimaryGeneratorAction::~DAMICPrimaryGeneratorAction()
{
    if(particleGun!=nullptr)
    { 
        delete particleGun;
        particleGun=nullptr;
    }
}

void DAMICPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{    
    _trackingAction->ResetMaps();
    particleGun->GeneratePrimaryVertex(anEvent);
}

