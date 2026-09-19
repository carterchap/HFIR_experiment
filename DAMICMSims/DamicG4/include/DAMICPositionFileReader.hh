#ifndef DAMICPositionFileReader_h
#define DAMICPositionFileReader_h

#include "G4ThreeVector.hh"

#include <vector>

class DAMICPositionFileReader 
{
    private:
        std::vector<G4ThreeVector> positionList;
        int positionListIdx;
        G4String fname;
    public:
        DAMICPositionFileReader(G4String filename);
        ~DAMICPositionFileReader() {};
        
        G4ThreeVector GetNextPosition();
        bool LoadFile(G4String filename);
};

#endif
