/**
 * @file DAMICUtils.hh
 * @author: Mariangela Settimo
 * @date 2018 - Subatech
 */


#ifndef DAMICUtils_h
#define DAMICUtils_h 1

#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"


G4bool FindVolume(G4String Vname);
//G4String IsPointInsideVolume(G4ThreeVector& point);
//void ListOfLogVolumes();
//bool IsVolInsideDetector(G4String volName);
//bool IsVolInsideCubeMod(G4String volName);
//bool IsVolInsideFibre(G4String volName);
int GetVolID(G4String volName);    ///< using map
int GetVolumeID(G4String volName);  ///<using if conditions
//double CubesDistance(int cube1, int cube2);
G4int GetProcessID(G4String procname); ///< get the process code ID 


#endif

