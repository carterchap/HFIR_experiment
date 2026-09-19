#ifndef DAMICAnalysisManager_hh
#define DAMICAnalysisManager_hh

#include "globals.hh"
#include "G4ThreeVector.hh"

//#include "g4root.hh" deprecated in v.11
#include "G4AnalysisManager.hh"

// system includes
#include<vector>

class DAMICAnalysisManagerMessenger;

class DAMICAnalysisManager
{
    public:
        static DAMICAnalysisManager * GetInstance();
        ~DAMICAnalysisManager();

        // Initialize the Manager
        void Initialize() { _g4man->OpenFile(); }
        void Finish() { _g4man->Write(); _g4man->CloseFile(); }
        
        // Clear all vectors
        void clear_vectors();
        // Fill all trees per event
        void FillTrees(G4int eventId);

        // Fill trees
        void FillRunInfo(G4int event, G4int nccd, G4int seed,
                G4int ccdVersion, G4int frameVersion, G4int cableVersion, 
                G4int vesselVersion, G4int cryoVersion, G4int shieldingVersion,
                std::string concatenatedVolNames,
                const std::vector<G4int> & volName, 
                const std::vector<G4float> & volMass_g, 
                const std::vector<G4float> & volDensity_gcm3, 
                const std::vector<G4float> & volVolume_cm3,
                const std::vector<G4float> & volSurface_cm2 );

        void FillEventOut(G4int pdg, G4int charge,
                G4int volId, G4double kinEnergy,
                const G4ThreeVector & pos, 
                const G4ThreeVector & momentum,
                G4double triggerTime);

        void FillCCDOut(G4int pdg, G4int trackId, G4int parentId,
                G4int ccdId, 
                G4double posX, G4double posY, G4double posZ,
                G4double posGlobX, G4double posGlobY, G4double posBlogZ,
                G4double energy, G4double time, G4double ltime);

        void FillTrackOut(G4int pdg, G4int trackId, G4int stepId, G4int parentId,
                G4int ccdId, 
                G4int volId, G4int startProcess, G4int endProcess, 
                G4double stepLength, 
                G4double energy, G4double deltaE, G4double Edep,                
                G4double Xcoord, G4double Ycoord, G4double Zcoord,
                G4double time);
        void FillPartInfo(G4int pdg, G4int trackID, G4int parentID, G4int parentpdg,
                G4int firstVolID, G4int lastVolID, 
                G4int creatorProcessID, G4int lastProcessID,
                G4double initialEkin, 
                const G4ThreeVector & pos );

       void FillShieldingOut(G4int pdg, G4int trackId, G4double energy, G4double weight, 
                G4double posX, G4double posY, G4double posZ,
                G4double momX, G4double momY, G4double momZ,
			       G4double time, G4int layerID, G4int stepid, G4double steplenght,G4int creat_proc_id);


        unsigned int GetPartInfoVectorSize() { return _partinfo_lastVolID.size(); }
        
        // Read from macro the name of the generated primary particle, and volume name
        void ParseMacro(const G4String & macroname);
        
        // To set/unset the store of tracking information into a ROOT TTree
        void SetStoreTrack(G4bool val);
        // To set/unset the store of particle information into a ROOT TTree
        void SetStorePart(G4bool val);

      // To set lighter version of EventOut tree
        void SetEventLight(G4bool val);

         // To set/unset the store of shielding intermediate information
        void SetStoreShielding(G4bool val);
        bool IsShieldOutON() {return _storeShielding;}; 

        // To set/unset the volume names for shielding intermediate information
        //void SetShieldingLastVolume(const G4String & val);
        //void SetShieldingStopVolume(const G4String & val);
        void SetShieldingLastVolume( const G4String & Last, const G4String & Stop );
  
      // return the 2 volume names (in/out) when ShieldOut tree is write out
      G4String GetShieldingLastVolume() { return _shieldVolumeLastValidName;}
      G4String GetShieldingStopVolume() { return _shieldVolumeStopName;}
      
      void LastVolumeON(G4bool val);
      
      bool IsLastVolumeON() {return _bLastVolume;};
        
    private:
        // Hidding constructor (singleton)
        DAMICAnalysisManager() { ; }
        // Be sure all TTree related vectors are included here
        void _createAuxiliaryData();
        // To include Track TTree if booked
        void _createAuxiliaryDataForTrack();
        // To include PartInfo TTree if booked
        void _createAuxiliaryDataForPart();
        // To include Shielding TTree if booked
        void _createAuxiliaryDataForShielding();

        // Initialization of the ROOT file and Trees        
        void _book();        
        // book independent for Track TTRee
        void _bookTrack();
        // book independent for PartInfo TTRee
        void _bookPart();
        // book independent for Shielding TTree
        void _bookShielding();
        // book independent for EventOut TTree
        void _bookEvent();

        static DAMICAnalysisManager * _instance;
        // The ROOT analysis manager (the actual manager with
        // root file and trees content and management
        G4AnalysisManager * _g4man;
        
        // XXX: BE SURE THAT all the TTree-related vectors
        // are listed in the _create_auxiliary_data !!
        // XXX

        // vector for the RunInfo tree: properties for each volume, mass, density and volume
        std::vector<G4int> _runinfo_volName;
        std::vector<G4double> _runinfo_volMass;
        std::vector<G4double> _runinfo_volDensity;
        std::vector<G4double> _runinfo_volVolume;
        std::vector<G4double> _runinfo_volSurface;


        // vector for the EventOut tree
        std::vector<G4int> _evtout_pdg;
        std::vector<G4int> _evtout_charge;
        std::vector<G4int> _evtout_VolID;
        std::vector<G4double> _evtout_posx;
        std::vector<G4double> _evtout_posy;
        std::vector<G4double> _evtout_posz;
        std::vector<G4double> _evtout_momx;
        std::vector<G4double> _evtout_momy;
        std::vector<G4double> _evtout_momz;
        std::vector<G4double> _evtout_ekin;
        std::vector<G4double> _evtout_time;


        // vector for the CCDOut tree
        std::vector<G4int> _ccdout_pdg;
        std::vector<G4int> _ccdout_trackID;
        std::vector<G4int> _ccdout_parentID;
        std::vector<G4int> _ccdout_ccdID;
        std::vector<G4double> _ccdout_posX;
        std::vector<G4double> _ccdout_posY;
        std::vector<G4double> _ccdout_posZ;
        std::vector<G4double> _ccdout_posGlobX;
        std::vector<G4double> _ccdout_posGlobY;
        std::vector<G4double> _ccdout_posGlobZ;
        std::vector<G4double> _ccdout_edep;
        std::vector<G4double> _ccdout_time;
        std::vector<G4double> _ccdout_local_time;

        // vector for the TrackOut tree
        std::vector<G4int> _trackout_pdg;
        std::vector<G4int> _trackout_trackID;
        std::vector<G4int> _trackout_stepID;
        std::vector<G4int> _trackout_parentID;
        std::vector<G4int> _trackout_ccdID;
        std::vector<G4int> _trackout_volID;
        std::vector<G4int> _trackout_startProcess;
        std::vector<G4int> _trackout_endProcess;
        std::vector<G4double> _trackout_stepLength;
        std::vector<G4double> _trackout_energy;
        std::vector<G4double> _trackout_deltaE;
        std::vector<G4double> _trackout_edep;
        std::vector<G4double> _trackout_posX;
        std::vector<G4double> _trackout_posY;
        std::vector<G4double> _trackout_posZ;
        std::vector<G4double> _trackout_time;
        
        // vectors for the PartInfo tree
        std::vector<G4int> _partinfo_pdg;
        std::vector<G4int> _partinfo_trackID;
        std::vector<G4int> _partinfo_parentID;
        std::vector<G4int> _partinfo_parentpdg;
        std::vector<G4int> _partinfo_firstVolID;
        std::vector<G4int> _partinfo_lastVolID;
        std::vector<G4int> _partinfo_creatorProcessID;
        std::vector<G4int> _partinfo_lastProcessID;
        std::vector<G4double> _partinfo_initialEkin;
        std::vector<G4double> _partinfo_posX;
        std::vector<G4double> _partinfo_posY;
        std::vector<G4double> _partinfo_posZ;

        // vector for the ShieldOut tree
        std::vector<G4int> _shieldout_pdg;
        std::vector<G4int> _shieldout_trackID;
        std::vector<G4double> _shieldout_energy;
        std::vector<G4double> _shieldout_weight;
        std::vector<G4double> _shieldout_posX;
        std::vector<G4double> _shieldout_posY;
        std::vector<G4double> _shieldout_posZ;
        std::vector<G4double> _shieldout_momX;
        std::vector<G4double> _shieldout_momY;
        std::vector<G4double> _shieldout_momZ;
        std::vector<G4double> _shieldout_time;
        std::vector<G4int> _shieldout_layerID;
	std::vector<G4int> _shieldout_stepid;
	std::vector<G4double> _shieldout_steplenght;
	std::vector<G4int> _shieldout_creat_proc_id;
	
        // auxiliary list of all vectors 
        std::vector<std::vector<G4int>* > _aux_ivec;
        std::vector<std::vector<G4double>* > _aux_dvec;

        // Identifiers of the each tree
        G4int _runInfoId;
        G4int _eventOutId;
        G4int _ccdOutId;
        G4int _trackOutId;
        G4int _partInfoId;
        G4int _shieldOutId;

        std::vector<G4int*> _eventwise_ids;

        // Name of the primary particle generator, and the volume/surface
        std::string _pparticleName;
        std::string _volumeName;
        std::string _ion;

        // Set STore Track
        G4bool _storeTrack;
        G4bool _treeTrackActivated;
        // to store partinfo
        G4bool _storePart;
	//to store light eventout
	G4bool _storeEventLight;
        // to store shielding 
        G4bool _storeShielding;
        G4String _shieldVolumeLastValidName;
        G4String _shieldVolumeStopName;
        G4bool _bLastVolume;
        // Messenger
        DAMICAnalysisManagerMessenger* _messenger;
};
#endif
