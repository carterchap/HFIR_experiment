#include <fstream>
#include <iomanip>

#include "DAMICStackingActionMessenger.hh"
#include "DAMICStackingAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"
#include <cmath>

DAMICStackingActionMessenger::DAMICStackingActionMessenger(DAMICStackingAction* fStackSource): 
    _stackAction(fStackSource)
{
    
    _particleTable = G4ParticleTable::GetParticleTable();
    
    _stackDirectory = new G4UIdirectory("/damic/stack/");
    _stackDirectory->SetGuidance("Stacking action control commands.");
    
    _listCmd = new G4UIcmdWithoutParameter("/damic/stack/list", this);
    _listCmd->SetGuidance("List available particles.");
    _listCmd->SetGuidance(" Invoke G4ParticleTable.");
    
    _resetCmd = new G4UIcmdWithoutParameter("/damic/stack/reset", this);
    _resetCmd->SetGuidance("To use in a macro file with several runs");
    
    IonLimitCmd = new G4UIcommand("/damic/stack/ionlimit",this);
    IonLimitCmd->SetGuidance("Set properties of ion to be generated.");
    IonLimitCmd->SetGuidance("[usage] /stack/ionlimit Z A");
    IonLimitCmd->SetGuidance("        Z:(int) AtomicNumber");
    IonLimitCmd->SetGuidance("        A:(int) AtomicMass");
    
    G4UIparameter* parami;
    parami = new G4UIparameter("Z",'i',true);
    parami->SetDefaultValue("1");
    IonLimitCmd->SetParameter(parami);
    parami = new G4UIparameter("A",'i',true);
    parami->SetDefaultValue("1");
    IonLimitCmd->SetParameter(parami);
    
    _addPartKillCmd = new G4UIcmdWithAString("/damic/stack/partkill", this);
    _addPartKillCmd->SetGuidance("Add particles to kill during the run");
    _addPartKillCmd->SetParameterName("Particule", true, true);
    
    _doFullDecayCmd = new G4UIcmdWithABool("/damic/stack/fulldecay",this);
    _doFullDecayCmd->SetGuidance("Activate full-decay option");
    _doFullDecayCmd->SetParameterName("Bool", false);
    _doFullDecayCmd->SetDefaultValue(false);

}


DAMICStackingActionMessenger::~DAMICStackingActionMessenger()
{
    // FIXME add deletes of all news and point to a nullptr FIXME
    delete _listCmd;
    delete IonLimitCmd;
    delete _addPartKillCmd;
    delete _resetCmd;
    delete _stackDirectory;
    delete _doFullDecayCmd;
}

void DAMICStackingActionMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
    //G4cout << "SetNew "  << G4endl;
    if (command == _listCmd)
    {
        _particleTable->DumpTable();
    }
    else if(command == IonLimitCmd)
    {
        G4Tokenizer next( newValues );
        // check argument
        G4int fAtomicNumber = StoI(next());
        G4int fAtomicMass = StoI(next());
        
        G4ParticleDefinition* ion;
        ion =  G4IonTable::GetIonTable()->GetIon(fAtomicNumber,fAtomicMass);
        if(ion==0) 
        {
            G4cout << "Ion with Z=" << fAtomicNumber;
            G4cout << " A=" << fAtomicMass << "is not well defined" << G4endl;
        }
        else 
        {
            _stackAction->SetLimitNucleus(ion->GetPDGEncoding());
        }
    }
    else if ( command == _addPartKillCmd)
    {
        _stackAction->AddParticlesKill(newValues);
    }
    else if (command == _resetCmd)
    {
        _stackAction->ResetNewRun();
    }
    else if (command == _doFullDecayCmd)
    {
        const G4bool boolValue = _doFullDecayCmd->GetNewBoolValue(newValues);
        _stackAction->UnsetDoOneDecay(boolValue);
    }
}
