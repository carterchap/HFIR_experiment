#ifndef DAMICPhysicsListPenelope_h
#define DAMICPhysicsListPenelope_h

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

class G4LossTableManager;
//class DAMICMilliChargedPhysics;

class DAMICPhysicsListPenelope: public G4VUserPhysicsList
{
public:
  DAMICPhysicsListPenelope();
  ~DAMICPhysicsListPenelope();

public:
  virtual void SetCuts();


protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();

  // these methods Construct physics processes and register them
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructHad();
  virtual void ConstructOp();


  virtual void AddTransportation();

private:
  G4int VerboseLevel;
  G4int OpVerbLevel;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  // these methods Construct particles
  void ConstructMyBosons();
  void ConstructMyLeptons();
  void ConstructMyHadrons();
  void ConstructMyShortLiveds();

  //	DAMICMilliChargedPhysics *milliPhysics;
	
};


#endif
