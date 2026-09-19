#ifndef DAMICParticleSource_h
#define DAMICParticleSource_h

#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "DAMICParticleSourceMessenger.hh"
#include "TTree.h"
#include "TFile.h"

class DAMICVolumeCoordinateStore;

class DAMICParticleSource: public G4VPrimaryGenerator
{
    public:
        DAMICParticleSource();
        ~DAMICParticleSource();

        void GeneratePrimaryVertex(G4Event* evt);
        
        //Input file
        void SetReadingShieldFile(G4bool);
  
        //Particle Sets
        void SetParticleEnergy(G4double);
        void SetParticleCharge(G4double);
        void SetParticleMomentumDirection(G4ParticleMomentum);
        void SetParticleTime(G4double);
        void SetParticleDefinition(G4ParticleDefinition*);
        void SetParticlePosition(G4ThreeVector);

        //Position Sets
        void SetMaterial(G4String);
        void SetVolumeSource(G4String);
        void SetMaterialSource(G4String);
        void SetMotherVolume(G4String);
        void SetSourceShape(G4String);
        void SetSurface(G4bool);
        void SetVolume(G4bool);
        void SetMaterial(G4bool);
        void SetSource(G4bool);
        void SetShape(G4bool);
        void SetList(G4bool);

        void SetPositionNumSize(G4int);
        void SetPositionNumValue(G4double, G4int);

        void SetEmbeddedDistance(G4double min_Dist, G4double max_Dist);
        void SetMinEmbeddedDistance(G4double min_Dist);
        void SetMaxEmbeddedDistance(G4double max_Dist);
        void SetCCDFacingOnly(G4bool val);
        void SetAdvDiffModel(G4bool val);
        void SetPositionList(G4String fname);

        // Position Reset
        void ResetPositionNum();

        //Energy Sets
        void SetDistriNRJ(G4bool, G4String);
        void SetMonoNRJ(G4bool);
        void SetEnergyNumSize(G4int);
        void SetEnergyNumValue(G4double, G4int);
        void SetHistoNRJ(G4String);
        void ReadHistoNRJ();
        //Energy Reset
        void ResetEnergyNum();

        // Direction Sets
        void SetDistriD(G4bool, G4String);
        void SetOneD(G4bool);
        void SetDirectionNumSize(G4int);
        void SetDirectionNumValue(G4double, G4int);

        // Direction Reset
        void ResetDirectionNum();
        
        // ReadParticleFile
        void ReadParticleFile(G4String);
        
        //GENERATORS for Shielding
        void SetCosmogenic(G4bool);
        void SetCosmoPart(G4String);
        void SetCosmoAngDistType(G4String);
        void SetCosmoSrcCenter(G4ThreeVector);
        void SetCosmoSrcSize(G4ThreeVector);
        void GenerateMuonAngleDist();
        void GenerateMuonAngleGaussDist();
        void GenerateIsotropicFlux(); 
        void GenerateNeutronEnergyDist();
        
        // void SetLabDepth(G4double);
  
        // energy distribution 
        void GenerateMuonEnergyDist(); 
        inline G4double GetParticleEnergy() {return _particleEnergy;}
        void GenerateMax(); 
        
        // VOLUME
        void AddVolume(G4String, G4double);
        void CalculProba();
        void ChooseVolume();
        //Energy
        void EnergyValue();
        
        // Direction
        void DirectionValue();
        //bools 
        G4bool IsSurface();
        G4bool IsVolume();
        G4bool IsMaterial();
        G4bool IsShape();
        G4bool IsSource();
        G4bool IsList();
        G4bool IsCosmogenic();
        G4bool IsReadingShieldFile();

        // bools Energy
        G4bool IsDistriNRJ();
        G4bool IsMonoNRJ();
  
        // bools Direction
        G4bool IsDistriD();
        G4bool IsOneD();
  
        //changes bools Position
        void DoVolume();
        void DoMaterial();
        void DoSource();
        void DoShape(G4String);
        void DoSurface();
        void DoPositionList();
        void DoCosmogenic();
        
        //Changes bools Energy
        void DoMonoNRJ();
        void DoDistriNRJ(G4String);
        
        // Changes bools Momentum
        void DoOneD();
        void DoDistriD(G4String);
  
        //coordinatesCreation
        void MaterialCoordinates();
        void VolumeCoordinates();
        void SourceCoordinates();
        void ShapeCoordinates();
        void CalculPosition(G4String, G4String, G4String);
        void SurfaceCoordinates();
        void ListCoordinates();
        //Mother Volume Hit and Miss
        void FindInMother(); // hit and miss
        
        G4double SpectrumBeta(G4double);
        
        // Random Vectors Process
        void RandomVecProf(std::vector<G4double>*);
        G4int PickRandom(std::vector<G4double>*);
  
    private:
        TFile *fin;
        TTree *treein;
  
        // Particle properties
        G4double _particleEnergy;
        G4double _particleCharge;
        G4ParticleMomentum _particleMomentumDirection;
        G4double _particleTime;
        G4ParticleDefinition* _particleDefinition;
        G4ThreeVector _particlePosition;
        
        //Position
        G4bool _bSurface; //True on surface
        G4bool _bVolume; //True if Volume
        G4bool _bMaterial;
        G4bool _bShape;
        G4bool _bSource;
        G4bool _bList;
        G4String _volumeSource; // IF VOLUME MODE
        G4String _materialSource; // IF MATERIAL MODE
        G4String _motherVolume;
        G4String _material;
        G4String _shape;
        
        G4bool _bInputFile;
        int  in_EventID;  
        std::vector<double> *in_posx = NULL;
        std::vector<double> *in_posy = NULL;
        std::vector<double> *in_posz = NULL;
        std::vector<double> *in_momx = NULL;
        std::vector<double> *in_momy = NULL;
        std::vector<double> *in_momz = NULL;
        std::vector<double> *in_energy = NULL;
        std::vector<double> *in_time = NULL;
        std::vector<double> *in_weight = NULL;
	std::vector<int> *in_layerID = NULL;
        std::vector<int> *in_pdg = nullptr;
        std::vector<int> *in_trackid=nullptr;

        // Energy
        G4bool _bDistriNRJ;
        G4bool _bMonoNRJ;
        G4String _distriNRJ;
        G4String _histoNRJ;
        
        //Tritium
        G4bool _bcalculatedMaxTri;
        G4double _maxTritium;
  
        //Direction
        G4bool _bOneD;
        G4bool _bDistriD;
        G4String _distriD;
        
        std::vector <G4double> _positionNum;
        std::vector <G4double> _energyNum;
        std::vector <G4double> _energyProba;
        std::vector <G4double> _directionNum;
        
        std::vector<G4String> _volumesUse;
        std::vector<G4double> _volumesUseMass;
        std::vector<G4double> _volumesConcentration;
        std::vector<G4double> _proba;
        
        // for volumes
        G4bool _bchanges;
        
        //cosmogenic - shielding
        G4bool    _bCosmogenic;
        
        G4String  _CosmoPart;
        G4String  _CosmoAngDistType;
        G4String  _CosmoShapeType; 
        G4double  _CosmoMinTheta = 0.;
        G4double  _CosmoMaxTheta = CLHEP::pi;
        G4double  _CosmoMinPhi = 0.;
        G4double  _CosmoMaxPhi = CLHEP::twopi;
        
        G4ThreeVector  _CosmoCentreCoords;
        G4ThreeVector  _CosmoSourceCubeSizes;
        
        //Navigator for Geometry
        G4Navigator* _gNavigator;
        
        // Messenger
        DAMICParticleSourceMessenger* _messenger;
        DAMICVolumeCoordinateStore* _coordinateStore;
  
};

#endif
