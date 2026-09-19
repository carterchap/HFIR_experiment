#ifndef DAMICEventAction_h
#define DAMICEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class DAMICRunAction;

class DAMICEventAction : public G4UserEventAction
{
  public:
    DAMICEventAction(DAMICRunAction* runAction);
    virtual ~DAMICEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    void SetID(G4int);

    void SetNamePrimary(G4String);

    void SetEnergyPrimary(G4double nrj);

    void SetAtomicMass(G4int);
    void SetAtomicNumber(G4int);

    G4int GetID() const {return ID;}
    G4String GetNamePrimary() const {return NamePrimary;}
    G4double GetEnergyPrimary() const {return EnergyPrimary;}
    G4int GetAtomicNumber() const {return AtomicNumber;}
    G4int GetAtomicMass() const {return AtomicMass;}

  private:
    DAMICRunAction* fRunAction;
    G4int ID;
    G4String NamePrimary;
    G4double EnergyPrimary;
    G4int AtomicNumber;
    G4int AtomicMass;
    G4int S_hits;
    G4int CCDCollID;
};

#endif
