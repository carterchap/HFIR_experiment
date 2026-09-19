#ifndef dimensions_H
#define dimensions_H

#include "G4SystemOfUnits.hh"
#include "defines.hh"
static const G4double worldSize=2.0*m;
//thickness without deadlayer
static const G4double CCD_dead_thickness=.5*um;
static const G4double CCD_thickness_raw=.5*mm;
static const G4double CCD_thickness=CCD_thickness_raw-CCD_dead_thickness*2;
static const G4double CCDx=4.5*cm;
static const G4double CCDy=4.5*cm;
static const G4double CCD_separation=.5*mm;
