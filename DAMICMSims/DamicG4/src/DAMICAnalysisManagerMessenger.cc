#include "DAMICAnalysisManagerMessenger.hh"
#include "DAMICAnalysisManager.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

DAMICAnalysisManagerMessenger::DAMICAnalysisManagerMessenger(DAMICAnalysisManager* _analManager):
    _analysisManager(_analManager)
{
    // Analysis Directory
    _analysisDirectory = new G4UIdirectory("/damic/analysis/");
    _analysisDirectory->SetGuidance("Analysis control commands.");
   
    // To store TTree Track 
    _storeTrackCmd = new G4UIcmdWithABool("/damic/analysis/setstoretrack",this);
    _storeTrackCmd->SetGuidance("Sets whether ROOT TTree for Tracking must be stored into the output root file");
    _storeTrackCmd->SetParameterName("Bool", false);
    _storeTrackCmd->SetDefaultValue(false);

    // To store TTree Part
    _storePartCmd = new G4UIcmdWithABool("/damic/analysis/setstoreparticleinfo",this);
    _storePartCmd->SetGuidance("Sets whether PartInfo ROOT TTree for particle infor must be stored into the output root file");
    _storePartCmd->SetParameterName("Bool", false);
    _storePartCmd->SetDefaultValue(false);


    // To store TTree EventOut Light
    _storeEventLightCmd = new G4UIcmdWithABool("/damic/analysis/setstoreeventlight",this);
    _storeEventLightCmd->SetGuidance("Sets whether EventOut ROOT TTree for event info must be stored in a lighter version");
    _storeEventLightCmd->SetParameterName("Bool", false);
    _storeEventLightCmd->SetDefaultValue(false);


     // To store TTree Shielding
    
    _storeShieldDir = new G4UIdirectory("/damic/analysis/shielding/");
    _storeShieldDir->SetGuidance("shielding control commands.");

       // To store TTree Part
    _storeShieldCmd = new G4UIcmdWithABool("/damic/analysis/shielding/storeOutput",this);
    _storeShieldCmd->SetGuidance("Activate shielding output ");
    _storeShieldCmd->SetParameterName("Bool", true);
    _storeShieldCmd->SetDefaultValue(false);

    _lastShieldVolCmd =  new G4UIcommand("/damic/analysis/shielding/LastandStopVolumes",this);    
    _lastShieldVolCmd->SetGuidance("Sets last and stop volume name for shielding output, eg. AncLeadShield_0_PV Lab_PV");
     G4UIparameter* param;
     param = new G4UIparameter("ShieldOutVolName", 's', false);
     param->SetDefaultValue("AncLeadShield_0_PV");
     _lastShieldVolCmd->SetParameter(param);
     param = new G4UIparameter("ShieldStopVolName", 's', false);
     param->SetDefaultValue("Lab_PV");
     _lastShieldVolCmd->SetParameter(param);

}

DAMICAnalysisManagerMessenger::~DAMICAnalysisManagerMessenger()
{
    if(_storeTrackCmd!=nullptr)
    {        
        delete _storeTrackCmd;
        _storeTrackCmd=nullptr;    
    }
    if(_storePartCmd!=nullptr)
    {        
        delete _storePartCmd;
        _storePartCmd=nullptr;    
    }

    if(_storeEventLightCmd!=nullptr)
    {        
        delete _storeEventLightCmd;
        _storeEventLightCmd=nullptr;    
    }
    
    if(_lastShieldVolCmd!=nullptr)
    {        
        delete _lastShieldVolCmd;
        _lastShieldVolCmd=nullptr;
    }
    /*        
    if(_StopShieldVolCmd!=nullptr)
    {        
        delete _StopShieldVolCmd;
        _StopShieldVolCmd=nullptr;
    }*/

    if(_storeShieldCmd!=nullptr)
    {        
        delete _storeShieldCmd;
        _storeShieldCmd=nullptr;
    }

    
    if(_analysisDirectory!=nullptr)
    {
        delete _analysisDirectory;
        _analysisDirectory=nullptr;
    }

     if(_storeShieldDir!=nullptr)
    {
        delete _storeShieldDir;
        _storeShieldDir=nullptr;
    }

    
}

void DAMICAnalysisManagerMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
    if(command==_storeTrackCmd)
    {
        const G4bool boolValue = _storeTrackCmd->GetNewBoolValue(newValues);
        _analysisManager->SetStoreTrack(boolValue);
    }
    else if(command==_storePartCmd)
    {
        const G4bool boolValue = _storePartCmd->GetNewBoolValue(newValues);
        _analysisManager->SetStorePart(boolValue);
    }
    else if(command==_storeEventLightCmd)
    {
        const G4bool boolValue = _storeEventLightCmd->GetNewBoolValue(newValues);
        _analysisManager->SetEventLight(boolValue);
    }
    else if(command==_storeShieldCmd)
    {
        const G4bool boolValue = _storeShieldCmd->GetNewBoolValue(newValues);
        _analysisManager->SetStoreShielding(boolValue);
    }
    else if(command==_lastShieldVolCmd)
    {

	std::vector<G4String> tokens;
        G4Analysis::Tokenize(newValues, tokens);
	_analysisManager->LastVolumeON(true);
        _analysisManager->SetShieldingLastVolume(tokens[0],tokens[1]);
    }
    /*else if(command==_lastShieldVolCmd)
    {
        _analysisManager->SetShieldingLastVolume(newValues);
    }    

 else if(command==_StopShieldVolCmd)
    {
        _analysisManager->SetShieldingStopVolume(newValues);
    } */   

  }
