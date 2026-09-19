#include "DAMICRunAction.hh"
#include "DAMICPrimaryGeneratorAction.hh"
#include "DAMICDetectorConstruction.hh"
#include "DAMICAnalysisManager.hh"
#include "DAMICParticleSource.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

#include "G4AccumulableManager.hh"
//#include "G4ParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
//#include "G4LogicalVolume.hh"
//#include "G4UnitsTable.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4EmCalculator.hh"

#include "G4ios.hh"
#include "Randomize.hh"

#include <fstream>
#include <iostream>

using namespace std;

DAMICRunAction::DAMICRunAction(): G4UserRunAction()
{
}

DAMICRunAction::~DAMICRunAction()
{
    delete G4AnalysisManager::Instance();
}

void DAMICRunAction::BeginOfRunAction(const G4Run* run)
{
    G4int numberEvent = run->GetNumberOfEventToBeProcessed();
    G4cout << "Events to be processed: " << numberEvent << G4endl;
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    G4cout << "Random generator seed value=" << G4Random::getTheSeed() << G4endl;
    
     
    G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
    parameterManager->Reset();
    
    // XXX: PROVISIONAL. IT SHOULD BE TAKEN FROM THE MACRO?
    G4int CCDVersion = -1;
    G4int CCDFrameVersion = -1;
    G4int VesselVersion = -1;
    G4int CableVersion = -1;
    G4int CryoVersion = -1;
    G4int ShieldingVersion = -1; 
    //MS : XXXX: ADD HERE LAST_SHIELD_VOLUME

    // Read geomtry tree and get the number of instances for each PV
    VolumeStructure();
    
    // Using the number of instances for each PV, get some properties
    std::vector<G4int> volTag;
    std::vector<G4float> volMass_kg;
    std::vector<G4float> volDensity_gcm3;
    std::vector<G4float> volVolume_cm3;
    std::vector<G4float> volSurface_cm2;
    std::string volConcName;

    // Print volume properties and record them into the RunInfo tree
    G4cout << "--- List of PV:" << G4endl;
    G4cout << "Tag      PV name    N_instances     Density     Mass        Volume      Surface " << G4endl;
    G4cout << "[units]  -   -   -   g/cm^3   kg  cm^3   cm^2" << G4endl;
    for(auto& i: _mapPVNinstances)
    {
        G4cout << " "<< _mapVolumeTag[i.first] <<" ,  "<< i.first <<" ,  "<< i.second <<" , " 
            << _mapLVobject[i.first]->GetMaterial()->GetDensity()*(CLHEP::cm3/CLHEP::g) <<" , "
            << _mapLVobject[i.first]->GetMass()/CLHEP::kg * i.second <<" , "
            << _mapLVobject[i.first]->GetSolid()->GetCubicVolume()/(CLHEP::cm3) * i.second <<" , "
            << _mapLVobject[i.first]->GetSolid()->GetSurfaceArea()/(CLHEP::cm2) * i.second << G4endl;
        
        if(!volConcName.empty())
        {
            volConcName += ";";
        }
        volConcName+=i.first;
        volTag.push_back(_mapVolumeTag[i.first]);
        volDensity_gcm3.push_back(_mapLVobject[i.first]->GetMaterial()->GetDensity()*(CLHEP::cm3/CLHEP::g));
        volMass_kg.push_back(_mapLVobject[i.first]->GetMass()/CLHEP::kg * i.second);
        volVolume_cm3.push_back(_mapLVobject[i.first]->GetSolid()->GetCubicVolume()/(CLHEP::cm3) * i.second);
        volSurface_cm2.push_back(_mapLVobject[i.first]->GetSolid()->GetSurfaceArea()/(CLHEP::cm2) * i.second);
    }
    G4cout << "" << G4endl;

    // Number of CCDs, assuming a PV name CCDSensor_PV
    int Nccds=0;
    if(_mapPVNinstances.count("CCDSensor_PV")>0)
    {
        Nccds=_mapPVNinstances["CCDSensor_PV"];
    }
    DAMICAnalysisManager * man = DAMICAnalysisManager::GetInstance();
    man->Initialize();
    man->FillRunInfo(numberEvent,Nccds,G4Random::getTheSeed(),
            CCDFrameVersion,CCDVersion,CableVersion,
            VesselVersion, CryoVersion, ShieldingVersion, 
            volConcName,volTag,volMass_kg,volDensity_gcm3,volVolume_cm3,volSurface_cm2); 


}

void DAMICRunAction::EndOfRunAction(const G4Run* /*run*/)
{
    // Merge parameters
    G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
    parameterManager->Merge();
    
    DAMICAnalysisManager * man = DAMICAnalysisManager::GetInstance();
    man->Finish();
    // Complete clean-up
}

int DAMICRunAction::GetVolumeID(const G4String & volname)
{
    if(_mapVolumeTag.find(volname) == _mapVolumeTag.end())
    {
        return -1;
    }
    return _mapVolumeTag[volname];
}

void DAMICRunAction::VolumeStructure()
{
    // A unique tag is associated to each LV present to the geomtry tree
    int volume_tag=0;

    // depth: starting from the top node (i.e. from the World)
    int depth=0;
    // for each depth (int), a list of PV
    std::map<int,std::vector<G4VPhysicalVolume*>> pvolumes_at_depth;
    // for each depth (int), the index of the current PV should be tracked
    std::map<int,unsigned int> pvIndex_at_depth;

    // Algorithm initialization: starting from the World
    // It ends, when we go back to World_PV
    G4VPhysicalVolume *pv = G4PhysicalVolumeStore::GetInstance()->GetVolume("World_PV");
    pvIndex_at_depth[depth]=0;
    pvolumes_at_depth[depth].push_back(pv);
    G4LogicalVolume *lv = pv->GetLogicalVolume();

    while(depth>=0)
    {
        // for each PV in depth, track daughters until to find a PV without daughters 
        bool searching_daughters=false;
        while(pvIndex_at_depth[depth] <= pvolumes_at_depth[depth].size()-1)
        {
            // load PV and its lv
            pv = pvolumes_at_depth[depth][pvIndex_at_depth[depth]];
            lv = pv->GetLogicalVolume();
            
            // if any, add daughters to next depth
            if(pv->GetLogicalVolume()->GetNoDaughters()>0)
            {
                // Has daughter
                // including another depth from this node (index)
                ++depth;
                // check if this level already exist (added from another node), 
                // and initialize the index map (if does exists)
                if(pvIndex_at_depth.find(depth)==pvIndex_at_depth.end())
                {
                    pvIndex_at_depth[depth]=0;
                }
                // append the daughters of LV to the next depth
                for(auto k=0;k<lv->GetNoDaughters();++k)
                {
                    pvolumes_at_depth[depth].push_back(lv->GetDaughter(k));
                }
                // Stop tracking more volumes from this depth
                // check if some of the daughters, have more daughters
                // needed to avoid to change depth again
                searching_daughters=true;

                // Continue to search for more daughter in the added level
                break;
            }
            // no daughters, go to next volume within the same level
            ++(pvIndex_at_depth[depth]);

            // Fill a map with the PV name and its instances
            // mother volumes are not included (only those volumes without daughters)
            G4String pv_name=pv->GetName();
            if(_mapPVNinstances.find(pv_name) == _mapPVNinstances.end())
            {
                _mapPVNinstances[pv_name]=0;
                _mapLVobject[pv_name]=lv;
                ++volume_tag;
                _mapVolumeTag[pv_name]=volume_tag;
            }
            ++(_mapPVNinstances[pv_name]);
        }
        // a deeper depth has been found
        if(searching_daughters)
        {
            continue;
        }
        // or, no more volumes to track in this depth 
        // continue naviagation: 
        // go to the next volume, in the previous depth
        --depth;
        ++(pvIndex_at_depth[depth]);
    }
}


