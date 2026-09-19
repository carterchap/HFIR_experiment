#ifndef DAMICAnalysisManagerMessenger_h
#define DAMICAnalysisManagerMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//forward declaration
class DAMICAnalysisManager;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;


class DAMICAnalysisManagerMessenger: public G4UImessenger
{
    public:
        DAMICAnalysisManagerMessenger(DAMICAnalysisManager *_analManager);
        ~DAMICAnalysisManagerMessenger();
        
        void SetNewValue(G4UIcommand *command, G4String newValues);

    private:
        DAMICAnalysisManager *_analysisManager;
        
        G4UIdirectory* _analysisDirectory;
        G4UIdirectory* _storeShieldDir;
  
        /// Store TTree Track Command
        G4UIcmdWithABool *_storeTrackCmd;

        // Store TTree Part Command
        G4UIcmdWithABool *_storePartCmd;

        // Store TTree EventOut Light Command
        G4UIcmdWithABool * _storeEventLightCmd;

        ///Store particle information after shielding - volumeName provided
        G4UIcmdWithABool *_storeShieldCmd;
        //G4UIcmdWithAString *_StopShieldVolCmd;
        G4UIcommand *_lastShieldVolCmd;
  
};
#endif
