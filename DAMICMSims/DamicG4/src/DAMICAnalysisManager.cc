#include "DAMICAnalysisManager.hh"
#include "DAMICAnalysisManagerMessenger.hh"

#include "G4UnitsTable.hh"

#include <algorithm>


DAMICAnalysisManager* DAMICAnalysisManager::_instance = nullptr;

DAMICAnalysisManager * DAMICAnalysisManager::GetInstance()
{
    if( _instance == nullptr )
    {
        _instance = new DAMICAnalysisManager();
        _instance->_createAuxiliaryData();
        // The ROOT-G4 analysis manager, it contains
        // the root file, the trees, etc.. 
        _instance->_g4man = G4AnalysisManager::Instance();
        _instance->_g4man->SetFileName("test");
        _instance->_g4man->SetVerboseLevel(1);
        _instance->_g4man->SetActivation(true);
        // Create messenger (/damic/analysis)
        _instance->_messenger = new DAMICAnalysisManagerMessenger(_instance);
        // Book all the trees
        _instance->_treeTrackActivated=false;
	    _instance->_storeShielding = false;
	    _instance->_shieldVolumeLastValidName = "";
     	_instance->_shieldVolumeStopName = "";
	_instance->_storeEventLight = false;

	_instance->_bLastVolume = false;

	_instance->_storeEventLight = false;	
        _instance->_book();
	_instance->_bookEvent();
    }

    return _instance;
}

DAMICAnalysisManager::~DAMICAnalysisManager()
{
    if(_messenger!=nullptr)
    {
        delete _messenger;
        _messenger=nullptr;
    }
}

G4String Remove_empty_string(std::string str)
{
    unsigned int first_alfanum = 0;
    for(unsigned int i=0; i < str.size();++i)
    {
        if( str[i] != ' ')
        {
            first_alfanum = i;
            break;
        }
    }
    
    //    G4cout << first_alfanum << G4endl;

    if(first_alfanum>0)
    {
        str.erase(0, first_alfanum);
    }
    return str.substr( 0, str.find(" "));
}

void DAMICAnalysisManager::LastVolumeON(G4bool val) {
  
   _bLastVolume = val;
}

void DAMICAnalysisManager::ParseMacro(const G4String & macroName )
{
    // Initialize
    _pparticleName = "";
    _volumeName = "";
    _ion = "";

    // exit if the macro is not batch mode
    if( macroName=="" )
    {
        return;
    }

    // extract from macro (if exists) the name of the primary particle, and the volume where it is
    // simulated
    std::ifstream macrofile( macroName );
    std::string line;

    while( std::getline(macrofile, line) )
    {
        if(line.find("/damic/gun/particle") != std::string::npos)
        {
            std::string pname(line.substr( line.find(" ")+1 ));            
            pname.erase( std::remove(pname.begin(), pname.end(), ' '), pname.end() );
            _pparticleName = pname;
        }
        else if(line.find("/damic/gun/ion") != std::string::npos)
        {
            std::string ion_z(line.substr( line.find(" ")+1 ));             
            std::string ion_a( ion_z.substr( ion_z.find(" ")+1 ));
            ion_z = ion_z.substr(0,ion_z.find(" ") );
            ion_a = ion_a.substr(0,ion_a.find(" ") );
            _ion = ion_a+"a"+ion_z+"z";
        }
        else if(line.find("/damic/gun/position/addvolume") != std::string::npos)
        {
            std::string vname( line.substr( line.find(" ")+1 ) );
            _volumeName = Remove_empty_string( vname );
        }
    }
}

void DAMICAnalysisManager::_book()
{
    _runInfoId = _g4man->CreateNtuple("RunInfo", "RunInfo");
    _g4man->CreateNtupleIColumn(_runInfoId,"NEvts");
    _g4man->CreateNtupleIColumn(_runInfoId,"NCCDs");
    _g4man->CreateNtupleIColumn(_runInfoId,"Seed");
    _g4man->CreateNtupleIColumn(_runInfoId,"CCDVersion");         ///< CCD design version (thickness, size), 0 = no CCD
    _g4man->CreateNtupleIColumn(_runInfoId,"CCDFrameVersion");    ///< CCDFrame design version, 0 = no CCDFrame
    _g4man->CreateNtupleIColumn(_runInfoId,"CableVersion"); ///< Cable design version, 0 = no cable
    _g4man->CreateNtupleIColumn(_runInfoId,"VesselVersion");      ///< Vessel design version, 0 = no Vessel
    _g4man->CreateNtupleIColumn(_runInfoId,"CryoVersion");        ///< Cryostat design version, 0 = no cryostat
    _g4man->CreateNtupleIColumn(_runInfoId,"ShieldingVersion");   ///< Shielding design version, 0 = no cryostat
    _g4man->CreateNtupleSColumn(_runInfoId,"concatedVolumeNames");
    _g4man->CreateNtupleIColumn(_runInfoId,"volumeNameID", _runinfo_volName);
    _g4man->CreateNtupleDColumn(_runInfoId,"volumeMass", _runinfo_volMass);
    _g4man->CreateNtupleDColumn(_runInfoId,"volumeDensity", _runinfo_volDensity);
    _g4man->CreateNtupleDColumn(_runInfoId,"volumeVolume", _runinfo_volVolume);
    _g4man->CreateNtupleDColumn(_runInfoId,"volumeSurface", _runinfo_volSurface);
    _g4man->CreateNtupleSColumn(_runInfoId,"primaryParticle");
    _g4man->CreateNtupleSColumn(_runInfoId,"primaryIon");
    _g4man->CreateNtupleSColumn(_runInfoId,"simulatedVolume");
    _g4man->FinishNtuple(_runInfoId);


    //CCD based tree - one entry per event, vector is one per particle et pixel
    _ccdOutId = _g4man->CreateNtuple("CCDOut", "CCDOut");
					      
    _g4man->CreateNtupleIColumn(_ccdOutId, "EventID");      
    _g4man->CreateNtupleIColumn(_ccdOutId, "pdg", _ccdout_pdg );          
    _g4man->CreateNtupleIColumn(_ccdOutId, "trackid", _ccdout_trackID );      
    _g4man->CreateNtupleIColumn(_ccdOutId, "parentid", _ccdout_parentID );
    _g4man->CreateNtupleIColumn(_ccdOutId, "CCDid", _ccdout_ccdID ); 
    _g4man->CreateNtupleDColumn(_ccdOutId, "posx", _ccdout_posX ); 
    _g4man->CreateNtupleDColumn(_ccdOutId, "posy", _ccdout_posY );  
    _g4man->CreateNtupleDColumn(_ccdOutId, "posz", _ccdout_posZ );  
    _g4man->CreateNtupleDColumn(_ccdOutId, "gposx", _ccdout_posGlobX );
    _g4man->CreateNtupleDColumn(_ccdOutId, "gposy", _ccdout_posGlobY );
    _g4man->CreateNtupleDColumn(_ccdOutId, "gposz", _ccdout_posGlobZ );
    _g4man->CreateNtupleDColumn(_ccdOutId, "Edep", _ccdout_edep );
    _g4man->CreateNtupleDColumn(_ccdOutId, "time", _ccdout_time );
    _g4man->CreateNtupleDColumn(_ccdOutId, "local_time", _ccdout_local_time ); 
    _g4man->FinishNtuple(_ccdOutId);

}

void DAMICAnalysisManager::_bookEvent()
{
 //event  tree - one entry per event, vector is normally one if 
    _eventOutId = _g4man->CreateNtuple("EventOut", "EventOut");
    _g4man->CreateNtupleIColumn(_eventOutId, "EventID");       ///< event id

    if(_storeEventLight){
    _g4man->CreateNtupleDColumn(_eventOutId, "posx", _evtout_posx);          ///< primary posx
    _g4man->CreateNtupleDColumn(_eventOutId, "posy", _evtout_posy);          ///< primary posy
    _g4man->CreateNtupleDColumn(_eventOutId, "posz", _evtout_posz);          ///< primary posz
    }
    else{
    _g4man->CreateNtupleIColumn(_eventOutId, "pdg", _evtout_pdg);           ///< primary charge
    _g4man->CreateNtupleIColumn(_eventOutId, "charge", _evtout_charge);        ///< primary pdg
    _g4man->CreateNtupleIColumn(_eventOutId, "volid", _evtout_VolID);         ///< primary volumeID
    _g4man->CreateNtupleDColumn(_eventOutId, "energy", _evtout_ekin);        ///< primary energy eV
    _g4man->CreateNtupleDColumn(_eventOutId, "posx", _evtout_posx);          ///< primary posx
    _g4man->CreateNtupleDColumn(_eventOutId, "posy", _evtout_posy);          ///< primary posy
    _g4man->CreateNtupleDColumn(_eventOutId, "posz", _evtout_posz);          ///< primary posz
    _g4man->CreateNtupleDColumn(_eventOutId, "momx", _evtout_momx);          ///< primary momx :  for muons and prompt neutrons
    _g4man->CreateNtupleDColumn(_eventOutId, "momy", _evtout_momy);          ///< primary momy :  for muons and prompt neutrons
    _g4man->CreateNtupleDColumn(_eventOutId, "momz", _evtout_momz);          ///< primary momz :  for muons and prompt neutrons
    _g4man->CreateNtupleDColumn(_eventOutId, "triggerTime", _evtout_time);   ///< primary time :  for muons and prompt neutrons
    }
    _g4man->FinishNtuple(_eventOutId); //1 
}

void DAMICAnalysisManager::_bookPart()
{
    if(_storePart)
    {
        _partInfoId = _g4man->CreateNtuple("PartInfo","particle info");
        _g4man->CreateNtupleIColumn(_partInfoId,"EventID");
        _g4man->CreateNtupleIColumn(_partInfoId,"pdg",      _partinfo_pdg);
        _g4man->CreateNtupleIColumn(_partInfoId,"trackid",  _partinfo_trackID );
        _g4man->CreateNtupleIColumn(_partInfoId,"parentid", _partinfo_parentID) ;
        _g4man->CreateNtupleIColumn(_partInfoId,"parentpdg", _partinfo_parentpdg );
        _g4man->CreateNtupleIColumn(_partInfoId,"firstVolID", _partinfo_firstVolID);
        _g4man->CreateNtupleIColumn(_partInfoId,"lastVolID",   _partinfo_lastVolID);
        _g4man->CreateNtupleIColumn(_partInfoId,"CreatorProcessID",  _partinfo_creatorProcessID);
        _g4man->CreateNtupleIColumn(_partInfoId,"LastProcessID", _partinfo_lastProcessID);
        _g4man->CreateNtupleDColumn(_partInfoId,"initialEkin",  _partinfo_initialEkin);
        _g4man->CreateNtupleDColumn(_partInfoId,"posx", _partinfo_posX);
        _g4man->CreateNtupleDColumn(_partInfoId,"posy", _partinfo_posY);
        _g4man->CreateNtupleDColumn(_partInfoId,"posz", _partinfo_posZ);
        _g4man->FinishNtuple(_partInfoId);
    }
}

void DAMICAnalysisManager::_bookTrack()
{
    if(_storeTrack)
    {
        //particle based tree - one entry per event, vector is per particle et step
        _trackOutId = _g4man->CreateNtuple("TrackOut", "TrackOut");
        _g4man->CreateNtupleIColumn(_trackOutId, "EventID");     ///< primary pdg
        _g4man->CreateNtupleIColumn(_trackOutId, "pdg", _trackout_pdg );
        _g4man->CreateNtupleIColumn(_trackOutId, "trackid", _trackout_trackID );
        _g4man->CreateNtupleIColumn(_trackOutId, "stepid", _trackout_stepID );
        _g4man->CreateNtupleIColumn(_trackOutId, "parentid", _trackout_parentID );
        _g4man->CreateNtupleIColumn(_trackOutId, "CCDid", _trackout_ccdID );
        _g4man->CreateNtupleIColumn(_trackOutId, "volid", _trackout_volID );
        _g4man->CreateNtupleIColumn(_trackOutId, "startProcess", _trackout_startProcess );
        _g4man->CreateNtupleIColumn(_trackOutId, "endProcess", _trackout_endProcess );
        _g4man->CreateNtupleDColumn(_trackOutId, "StepLength", _trackout_stepLength );
        _g4man->CreateNtupleDColumn(_trackOutId, "energy", _trackout_energy );
        _g4man->CreateNtupleDColumn(_trackOutId, "deltaE", _trackout_deltaE );
        _g4man->CreateNtupleDColumn(_trackOutId, "Edep", _trackout_edep );
        _g4man->CreateNtupleDColumn(_trackOutId, "posx", _trackout_posX );
        _g4man->CreateNtupleDColumn(_trackOutId, "posy", _trackout_posY );
        _g4man->CreateNtupleDColumn(_trackOutId, "posz", _trackout_posZ );
        _g4man->CreateNtupleDColumn(_trackOutId, "time", _trackout_time );
        _g4man->FinishNtuple(_trackOutId);

       _treeTrackActivated=true;
    }
    else
    {
        if(_treeTrackActivated)
        {
            // deleteja el tree de alguna manare
            // in the GUI mode, the Track Tree can be booked by seting and unseting the command
            // several times, once you set the commant the track tree is booked, but if then is
            // unset, this should be removed, if i find the way to do it, buy now, i just keep it on
            // memory but not stored
        }
    }
}

void DAMICAnalysisManager::_bookShielding()
{
    // Tree with intermediate informations after shielding --
    // volume name of the shielding stage will be store in RunInfo
    // One entry per event, vector is one per particle
    if(_storeShielding)
    {
        _shieldOutId = _g4man->CreateNtuple("ShieldOut", "ShieldOut");
        _g4man->CreateNtupleIColumn(_shieldOutId, "EventID");
        _g4man->CreateNtupleIColumn(_shieldOutId, "pdg",    _shieldout_pdg);
        _g4man->CreateNtupleIColumn(_shieldOutId, "trackid",_shieldout_trackID);
        _g4man->CreateNtupleDColumn(_shieldOutId, "energy", _shieldout_energy);
        _g4man->CreateNtupleDColumn(_shieldOutId, "posx",   _shieldout_posX);
        _g4man->CreateNtupleDColumn(_shieldOutId, "posy",   _shieldout_posY);
        _g4man->CreateNtupleDColumn(_shieldOutId, "posz",   _shieldout_posZ);
        _g4man->CreateNtupleDColumn(_shieldOutId, "momx",   _shieldout_momX);
        _g4man->CreateNtupleDColumn(_shieldOutId, "momy",   _shieldout_momY);
        _g4man->CreateNtupleDColumn(_shieldOutId, "momz",   _shieldout_momZ);
        _g4man->CreateNtupleDColumn(_shieldOutId, "time",   _shieldout_time);
        _g4man->CreateNtupleIColumn(_shieldOutId, "layerID", _shieldout_layerID);
	_g4man->CreateNtupleIColumn(_shieldOutId, "step_id", _shieldout_stepid);
	_g4man->CreateNtupleDColumn(_shieldOutId, "step_lenght",   _shieldout_steplenght);
	_g4man->CreateNtupleIColumn(_shieldOutId, "creat_proc_id",   _shieldout_creat_proc_id);
        _g4man->FinishNtuple(_shieldOutId);
    }
}



void DAMICAnalysisManager::_createAuxiliaryData()
{
    // List of trees id in order to fill them per event
    // _eventwise_ids.push_back(_runInfoId); NO, because is per RUN
    _eventwise_ids.push_back(&_eventOutId);
    _eventwise_ids.push_back(&_ccdOutId);
    
    // Vectors related with Trees
    _aux_ivec.push_back(&_evtout_pdg);
    _aux_ivec.push_back(&_evtout_charge);
    _aux_ivec.push_back(&_evtout_VolID);

    _aux_ivec.push_back(&_ccdout_pdg);
    _aux_ivec.push_back(&_ccdout_trackID);
    _aux_ivec.push_back(&_ccdout_parentID);
    _aux_ivec.push_back(&_ccdout_ccdID);
        
    _aux_dvec.push_back(&_evtout_ekin);
    _aux_dvec.push_back(&_evtout_posx);
    _aux_dvec.push_back(&_evtout_posy);
    _aux_dvec.push_back(&_evtout_posz);
    _aux_dvec.push_back(&_evtout_momx);
    _aux_dvec.push_back(&_evtout_momy);
    _aux_dvec.push_back(&_evtout_momz);
    _aux_dvec.push_back(&_evtout_time);

    _aux_dvec.push_back(&_ccdout_posX);
    _aux_dvec.push_back(&_ccdout_posY);
    _aux_dvec.push_back(&_ccdout_posZ);
    _aux_dvec.push_back(&_ccdout_posGlobX);
    _aux_dvec.push_back(&_ccdout_posGlobY);
    _aux_dvec.push_back(&_ccdout_posGlobZ);
    _aux_dvec.push_back(&_ccdout_edep);
    _aux_dvec.push_back(&_ccdout_time);
    _aux_dvec.push_back(&_ccdout_local_time);
}

void DAMICAnalysisManager::_createAuxiliaryDataForPart()
{
    if(_storePart)
    {
        _eventwise_ids.push_back(&_partInfoId);

        // Vectors related with Trees
        _aux_ivec.push_back(&_partinfo_pdg);
        _aux_ivec.push_back(&_partinfo_trackID);
        _aux_ivec.push_back(&_partinfo_parentID);
        _aux_ivec.push_back(&_partinfo_parentpdg);
        _aux_ivec.push_back(&_partinfo_firstVolID);
        _aux_ivec.push_back(&_partinfo_lastVolID);
        _aux_ivec.push_back(&_partinfo_creatorProcessID);
        _aux_ivec.push_back(&_partinfo_lastProcessID);
        
        _aux_dvec.push_back(&_partinfo_initialEkin);
        _aux_dvec.push_back(&_partinfo_posX);
        _aux_dvec.push_back(&_partinfo_posY);
        _aux_dvec.push_back(&_partinfo_posZ);
    }
}

void DAMICAnalysisManager::_createAuxiliaryDataForTrack()
{
    if(_storeTrack)
    {
        _eventwise_ids.push_back(&_trackOutId);

        _aux_ivec.push_back(&_trackout_pdg);
        _aux_ivec.push_back(&_trackout_trackID);
        _aux_ivec.push_back(&_trackout_stepID);
        _aux_ivec.push_back(&_trackout_parentID);
        _aux_ivec.push_back(&_trackout_ccdID);
        _aux_ivec.push_back(&_trackout_volID);
        _aux_ivec.push_back(&_trackout_startProcess);
        _aux_ivec.push_back(&_trackout_endProcess);

        _aux_dvec.push_back(&_trackout_stepLength);
        _aux_dvec.push_back(&_trackout_energy);
        _aux_dvec.push_back(&_trackout_deltaE);
        _aux_dvec.push_back(&_trackout_edep);
        _aux_dvec.push_back(&_trackout_posX);
        _aux_dvec.push_back(&_trackout_posY);
        _aux_dvec.push_back(&_trackout_posZ);
        _aux_dvec.push_back(&_trackout_time);
    }
}

void DAMICAnalysisManager::_createAuxiliaryDataForShielding()
{
    if(_storeShielding)
    {
        _eventwise_ids.push_back(&_shieldOutId);

        _aux_ivec.push_back(&_shieldout_pdg);
	_aux_ivec.push_back(&_shieldout_layerID);
	
        _aux_ivec.push_back(&_shieldout_trackID);
	_aux_dvec.push_back(&_shieldout_weight);
        _aux_dvec.push_back(&_shieldout_energy);
        _aux_dvec.push_back(&_shieldout_posX);
        _aux_dvec.push_back(&_shieldout_posY);
        _aux_dvec.push_back(&_shieldout_posZ);
        _aux_dvec.push_back(&_shieldout_momX);
        _aux_dvec.push_back(&_shieldout_momY);
        _aux_dvec.push_back(&_shieldout_momZ);
        _aux_dvec.push_back(&_shieldout_time);
	_aux_ivec.push_back(&_shieldout_stepid);
	_aux_dvec.push_back(&_shieldout_steplenght);
	_aux_ivec.push_back(&_shieldout_creat_proc_id);
      }
}

void DAMICAnalysisManager::clear_vectors()
{
    for(const auto & vi: _aux_ivec)
    {
        vi->clear();
    }
    for(const auto & vd: _aux_dvec)
    {
        vd->clear();
    }
}

void DAMICAnalysisManager::FillTrees(G4int eventId)
{
    for(auto i: _eventwise_ids)
    {

        // Fill the event Id, in all trees MUST BE first element
        _g4man->FillNtupleIColumn(*i,0,eventId);
        // Fill the tree per event
        _g4man->AddNtupleRow(*i);
    }
}

void DAMICAnalysisManager::FillRunInfo(G4int nevents, G4int nccd, G4int seed,
        G4int ccdVersion, G4int frameVersion, 
        G4int cableVersion, G4int vesselVersion, 
        G4int cryoVersion, G4int shieldingVersion,
        std::string concatenatedVolNames,
        const std::vector<G4int> & volName, 
        const std::vector<G4float> & volMass_g, 
        const std::vector<G4float> & volDensity_gcm3, 
        const std::vector<G4float> & volVolume_cm3,
        const std::vector<G4float> & volSurface_cm2)
{
    for(unsigned int i=0; i<volMass_g.size(); ++i)
    {
        _runinfo_volName.push_back( volName[i] );
        _runinfo_volMass.push_back( volMass_g[i] );
        _runinfo_volDensity.push_back( volDensity_gcm3[i] );
        _runinfo_volVolume.push_back( volVolume_cm3[i] );
        _runinfo_volSurface.push_back( volSurface_cm2[i] );
    }
    _g4man->FillNtupleIColumn(_runInfoId,0,nevents);
    _g4man->FillNtupleIColumn(_runInfoId,1,nccd);
    _g4man->FillNtupleIColumn(_runInfoId,2,seed);
    _g4man->FillNtupleIColumn(_runInfoId,3,ccdVersion);
    _g4man->FillNtupleIColumn(_runInfoId,4,frameVersion);
    _g4man->FillNtupleIColumn(_runInfoId,5,cableVersion);
    _g4man->FillNtupleIColumn(_runInfoId,6,vesselVersion);
    _g4man->FillNtupleIColumn(_runInfoId,7,cryoVersion);
    _g4man->FillNtupleIColumn(_runInfoId,8,shieldingVersion);
    _g4man->FillNtupleSColumn(_runInfoId,9,concatenatedVolNames);
    _g4man->FillNtupleSColumn(_runInfoId,15,_pparticleName);
    _g4man->FillNtupleSColumn(_runInfoId,16,_ion);
    _g4man->FillNtupleSColumn(_runInfoId,17,_volumeName);

    _g4man->AddNtupleRow(_runInfoId);
}

void DAMICAnalysisManager::FillEventOut(G4int pdg, G4int charge,
        G4int volId, G4double kinEnergy,
        const G4ThreeVector & pos, 
        const G4ThreeVector & mom,
        G4double triggerTime)
{
    _evtout_pdg.push_back(pdg);
    _evtout_charge.push_back(charge);
    _evtout_VolID.push_back(volId);
    _evtout_posx.push_back(pos.getX());
    _evtout_posy.push_back(pos.getY());
    _evtout_posz.push_back(pos.getZ());
    _evtout_momx.push_back(mom.getX());
    _evtout_momy.push_back(mom.getY());
    _evtout_momz.push_back(mom.getZ());
    _evtout_ekin.push_back(kinEnergy);
    _evtout_time.push_back(triggerTime);

/*
    _g4man->FillNtupleIColumn(_eventOutId,1, pdg);     ///< primary pdg
    _g4man->FillNtupleIColumn(_eventOutId,2, charge);     ///< primary charge
    _g4man->FillNtupleIColumn(_eventOutId,3, volId);     ///< primary volumeID
    _g4man->FillNtupleDColumn(_eventOutId,4, kinEnergy);  ///< primary energy
    _g4man->FillNtupleDColumn(_eventOutId,5, pos.getX());   ///< primary posx
    _g4man->FillNtupleDColumn(_eventOutId,6, pos.getY());   ///< primary posy
    _g4man->FillNtupleDColumn(_eventOutId,7, pos.getZ());   ///< primary posz
    _g4man->FillNtupleDColumn(_eventOutId,8, momentum.getX());   ///< primary momx :  for muons and prompt neutrons
    _g4man->FillNtupleDColumn(_eventOutId,9, momentum.getY());   ///< primary momy :  for muons and prompt neutrons
    _g4man->FillNtupleDColumn(_eventOutId,10,momentum.getZ());   ///< primary momz :  for muons and prompt neutrons
    _g4man->FillNtupleDColumn(_eventOutId,11,triggerTime);   ///< primary momz :  for muons and prompt neutrons
    _g4man->AddNtupleRow(_eventOutId);  //keep this uncommented for events with multiple vertexes
*/
}

void DAMICAnalysisManager::FillCCDOut(G4int pdg, 
        G4int trackId, G4int parentId,G4int ccdId, 
        G4double posX, G4double posY, G4double posZ,
        G4double posGlobX, G4double posGlobY, G4double posGlobZ,
        G4double energy, G4double time, G4double ltime)
{

    _ccdout_pdg.push_back(pdg);
    _ccdout_trackID.push_back(trackId);
    _ccdout_parentID.push_back(parentId);
    _ccdout_ccdID.push_back(ccdId);
    _ccdout_posX.push_back(posX);
    _ccdout_posY.push_back(posY);
    _ccdout_posZ.push_back(posZ);
    _ccdout_posGlobX.push_back(posGlobX);
    _ccdout_posGlobY.push_back(posGlobY);
    _ccdout_posGlobZ.push_back(posGlobZ);
    _ccdout_edep.push_back(energy);
    _ccdout_time.push_back(time);
    _ccdout_local_time.push_back(ltime);
}

void DAMICAnalysisManager::FillTrackOut(G4int pdg, G4int trackId, G4int stepId, G4int parentId,
        G4int ccdId,
        G4int volId, G4int startProcess, G4int endProcess,
        G4double stepLength,
        G4double energy, G4double deltaE, G4double Edep,
        G4double posx, G4double posy, G4double posz,
        G4double time)
{
    if(_storeTrack)
    {
        _trackout_pdg.push_back(pdg);
        _trackout_trackID.push_back(trackId);
        _trackout_stepID.push_back(stepId);
        _trackout_parentID.push_back(parentId);
        _trackout_ccdID.push_back(ccdId);
        _trackout_volID.push_back(volId);
        _trackout_startProcess.push_back(startProcess);
        _trackout_endProcess.push_back(endProcess);
        _trackout_stepLength.push_back(stepLength);
        _trackout_energy.push_back(energy);
        _trackout_deltaE.push_back(deltaE);
        _trackout_edep.push_back(Edep);
        _trackout_posX.push_back(posx);
        _trackout_posY.push_back(posy);
        _trackout_posZ.push_back(posz);
       _trackout_time.push_back(time);
    }
}



void DAMICAnalysisManager::FillPartInfo(G4int pdg, G4int trackID, G4int parentID, G4int parentPDG,
        G4int firstVolID, G4int lastVolID,
        G4int creatorProcessID, G4int lastProcessID,
        G4double initialEkin,
        const G4ThreeVector & pos)
{
    if(_storePart)
    {
        _partinfo_pdg.push_back(pdg);
        _partinfo_trackID.push_back(trackID);
        _partinfo_parentID.push_back(parentID);
        _partinfo_parentpdg.push_back(parentPDG);
        _partinfo_firstVolID.push_back(firstVolID);
        _partinfo_lastVolID.push_back(lastVolID);
        _partinfo_creatorProcessID.push_back(creatorProcessID);
        _partinfo_lastProcessID.push_back(lastProcessID);
        _partinfo_initialEkin.push_back(initialEkin);
        _partinfo_posX.push_back(pos.getX());
        _partinfo_posY.push_back(pos.getY());
        _partinfo_posZ.push_back(pos.getZ());
    }
}

void DAMICAnalysisManager::FillShieldingOut(G4int pdg, G4int trackId, 
        G4double energy, G4double weight,  
        G4double posx, G4double posy, G4double posz,
        G4double momx, G4double momy, G4double momz,
        G4double time, G4int layerID, G4int stepid, G4double steplenght, G4int creat_proc_id)
{
    if(_storeShielding)
    {
        _shieldout_pdg.push_back(pdg);
        _shieldout_trackID.push_back(trackId);
        _shieldout_energy.push_back(energy);
	_shieldout_weight.push_back(weight);
        _shieldout_posX.push_back(posx);
        _shieldout_posY.push_back(posy);
        _shieldout_posZ.push_back(posz);
        _shieldout_momX.push_back(momx);
        _shieldout_momY.push_back(momy);
        _shieldout_momZ.push_back(momz);
        _shieldout_time.push_back(time);
	_shieldout_layerID.push_back(layerID);
	_shieldout_stepid.push_back(stepid);
	_shieldout_steplenght.push_back(steplenght);
	_shieldout_creat_proc_id.push_back(creat_proc_id);
    }
}

void DAMICAnalysisManager::SetStoreTrack(G4bool val)
{
    _storeTrack = val;
    _createAuxiliaryDataForTrack();
    _bookTrack();
}

void DAMICAnalysisManager::SetStorePart(G4bool val)
{
    _storePart = val;
    _createAuxiliaryDataForPart();
    _bookPart();
}

void DAMICAnalysisManager::SetEventLight(G4bool val)
{
    _storeEventLight = val;
    _bookEvent();
}

void DAMICAnalysisManager::SetStoreShielding(G4bool val)
{
   _storeShielding = val;
   _createAuxiliaryDataForShielding();
   _bookShielding();
}

void DAMICAnalysisManager::SetShieldingLastVolume(const G4String & Last, const G4String & Stop)
{
  _shieldVolumeLastValidName = Last;
   _shieldVolumeStopName = Stop;
}


