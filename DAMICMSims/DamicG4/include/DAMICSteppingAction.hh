#ifndef DAMICSteppingAction_h
#define DAMICSteppingAction_h

#include "G4UserSteppingAction.hh"
#include "globals.hh"


class DAMICEventAction;
class G4LogicalVolume;
class DAMICRunAction;

class DAMICSteppingAction : public G4UserSteppingAction
{
    public:
        DAMICSteppingAction(DAMICRunAction *ra);
        virtual ~DAMICSteppingAction();
        
        virtual void UserSteppingAction(const G4Step*);
    
    private:
        DAMICRunAction*  _runAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
