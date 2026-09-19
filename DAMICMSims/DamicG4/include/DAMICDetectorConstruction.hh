#ifndef DAMICDetectorConstruction_h
#define DAMICDetectorConstruction_h

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

#include "globals.hh"
#include<vector>

class G4GDMLParser;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;

class DAMICDetectorConstruction: public G4VUserDetectorConstruction
{
    public:
        DAMICDetectorConstruction(const G4GDMLParser& parser);
        virtual ~DAMICDetectorConstruction();

        virtual G4VPhysicalVolume* Construct();
        virtual void ConstructSDandField();

        //void RegionInitialization();

    private:
        // detector geometry GDML input file
        const G4GDMLParser& _gdmlParser;

        G4ThreeVector _GetColor(const std::string& s );

        // Read Region name and length cut from GDML attribute
        std::tuple<std::string, double, double, double, double> _GetRegionAndCut(const std::string& s);
};

#endif
