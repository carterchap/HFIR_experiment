#ifndef DAMICVolumeCoordinateStore_h
#define DAMICVolumeCoordinateStore_h
//TODO: cleanup includes, move to .cc file
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include <vector>


class DAMICPositionFileReader;

class DAMICVolumeCoordinateStore 
{
   public:
        DAMICVolumeCoordinateStore(G4PhysicalVolumeStore*);
        ~DAMICVolumeCoordinateStore();
        
        G4bool AddVolume(G4String Volume, G4double Concentration=1);
        G4bool RemoveVolume(G4String Volume);
        
        G4ThreeVector GetPointOnSurfaceByName(G4String VolumeName);
        G4ThreeVector GetPointOnSurfaceByID(std::vector<G4VSolid>::size_type VolumeID);
        
        G4ThreeVector GetPointInVolumeByName(G4String VolumeName);
        G4ThreeVector GetPointInVolumeByID(std::vector<G4VSolid>::size_type VolumeID);
        
        G4ThreeVector GetPointInList();
        
        G4int GetNoVolumes() {return NumberOfVolumes;}
        
        void Clear();
        //  G4ThreeVector GetGlobalCoordinates(G4VPhysicalVolume *Volume);
        G4ThreeVector GetPointOnSurface();
        G4ThreeVector GetPointInVolume();
        void BuildHistogram();
        
        void SetEmbeddedDistance(G4double min_distance=100*nm, G4double max_distance=100*nm);
        void SetMinEmbeddedDistance(G4double min_distance=100*nm);
        void SetMaxEmbeddedDistance(G4double max_distance=100*nm);
        void SetCCDFacingOnly(G4bool val=false);
        void SetAdvDiffModel(G4bool val=false);
        void SetPositionList(G4String fname);
        
        G4double GetMinEmbeddedDistance() const {return minEmbeddedDistance;}
        G4double GetMaxEmbeddedDistance() const {return maxEmbeddedDistance;}
        G4bool GetCCDFacingOnly() const {return CCDFacingOnly;}
    
   private:
        
        G4int NumberOfVolumes; 
        
        // Update volumes, concentration, material, etc.. stores
        void _updateStore(G4VPhysicalVolume * pv, G4LogicalVolume * lv,
                G4double concentration, std::vector<G4VPhysicalVolume*> cur_mothers);
        // Stores
        //Vector to store the volumes that we're using in
        std::vector<G4VPhysicalVolume*> _physicalVolumes;
        std::vector<G4LogicalVolume*> _logicalVolumes;
        //Holds pointers to the solids
        std::vector<G4VSolid*> _solids;
        //Vector to store the origin of the volumes, in the *global* coordinate system
        std::vector<G4ThreeVector> Translations;
        std::vector<G4RotationMatrix*> Rotations;
        //Vectors to hold the geometric surface area/volume of our volumes
        std::vector<G4double> Masses;
        std::vector<G4double> CubicVolumes;
        std::vector<G4double> SurfaceAreas;
        std::vector<G4double> Concentrations;
        //Holds the mothers
        std::vector<std::vector<G4VPhysicalVolume*> > Mothers;
        G4PhysicalVolumeStore *theStore;

        //List of positions to use for simulation
        //Needs to be static for multithreading purposes
        static DAMICPositionFileReader *fileReader;

        //Holds information about probabiility for particle generator
        std::vector<G4double> surfaceProbabilities;
        std::vector<G4double> volumeProbabilities;
        G4double totalSurfaceProbability;
        G4double totalVolumeProbability;
        G4bool histogramBuilt;
        G4int GetRandomSurfaceGeneratorID();
        G4int GetRandomVolumeGeneratorID();
        G4double minEmbeddedDistance;
        G4double maxEmbeddedDistance;
        G4bool CCDFacingOnly;
        G4bool advDiffModel;
        G4ThreeVector CCDOrigin;

};

#endif
