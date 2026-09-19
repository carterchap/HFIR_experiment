#ifndef DAMICRunAction_h
#define DAMICRunAction_h

#include "G4UserRunAction.hh"

#ifndef G4VERSION_NUMBER
#include "G4Version.hh"
#endif

#if G4VERSION_NUMBER < 1030
#include "G4Parameter.hh"
#elseif G4VERSION_NUMBER >= 1030
#include "G4Accumulable.hh"
#endif

#include "globals.hh"
#include <vector>

#include "DAMICVolumeCoordinateStore.hh"

class G4Run;


class DAMICRunAction : public G4UserRunAction
{

    public:
        DAMICRunAction();
        virtual ~DAMICRunAction();
        
        virtual void BeginOfRunAction(const G4Run* run);
        virtual void EndOfRunAction(const G4Run*);
        
        void AddEdep (G4double edep);

        // Return the volume unique id present in the RunInfo tree
        // or -1 if does not exist
        int GetVolumeID(const G4String & volname);
       
        // Navigate through all geometry tree to list all PV and its instances
        void VolumeStructure();
       
        G4double EventID;
        G4double triggerTime;
        G4int primary_volid;

    private:
        
        // For VolumeStructure
        // map between the PV_name with an unique ID value
        std::map<G4String, int> _mapVolumeTag;
        // map between the PV_name and its Logical Volume
        std::map<G4String,G4LogicalVolume*> _mapLVobject;
        // map between the PV_name and the number of instances (i.e. how many times appear on the
        // geometry)
        std::map<G4String,int> _mapPVNinstances;
};

#endif
