#include "DAMICPositionFileReader.hh"

#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iostream>

DAMICPositionFileReader::DAMICPositionFileReader(G4String filename) 
{
    positionListIdx=-1;
    LoadFile(filename);
}


bool DAMICPositionFileReader::LoadFile(G4String filename)
{
    G4cout << "Loading file: " << filename << G4endl;
    if (fname==filename)
    {
        return true;
    }
    fname=filename;
    positionList.clear();
    positionListIdx=-1;
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        G4cout << "Error, file for coordinate list not found.." << G4endl;
        return false;
    }
    
    G4double x,y,z;
    while (file >> x >> y >>z) 
    {
        positionList.push_back(G4ThreeVector(x*mm, y*mm, z*mm));
    }
    positionListIdx=0;
    G4cout << "Loaded list of positions from file: " << filename << G4endl;
}
        
G4ThreeVector DAMICPositionFileReader::GetNextPosition()
{
    G4ThreeVector Pos(-1,-1,-1);
    if (positionListIdx==-1)
    {
        G4cout << "Error, position list not loaded" << G4endl;
    }
    else
    {
        Pos=positionList[positionListIdx];
        positionListIdx++;
    }
    return Pos;
}
