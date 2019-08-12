//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef APPARATUSXRAYABSORBERS_HH
#define APPARATUSXRAYABSORBERS_HH

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units
#include "G4IntersectionSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;

///////////////////////////////////////////////////////////////////////
// ApparatusXRayAbsorbers
///////////////////////////////////////////////////////////////////////
class ApparatusXRayAbsorbers
{
public:
    ApparatusXRayAbsorbers(G4int suppSwitch);
    ~ApparatusXRayAbsorbers();

    G4int Build();
    G4int Place(G4LogicalVolume* expHallLog, G4int selector);

private:
    G4LogicalVolume* fAbsorber;
    G4LogicalVolume* fAbsorberLayer1Log;
    G4LogicalVolume* fAbsorberLayer2Log;
    G4LogicalVolume* fAbsorberLayer3Log;

    G4AssemblyVolume* fAssemblyAbsorber;

private:
    // Materials
    G4Material* fAbsorberLayer1Mat;
    G4Material* fAbsorberLayer2Mat;
    G4Material* fAbsorberLayer3Mat;
    // Methods
    G4int BuildLayer1();
    G4int BuildLayer2();
    G4int BuildLayer3();
    G4Box* AbsorberLayer(G4int);
    G4double TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi);
    G4double TransZ(G4double x, G4double z, G4double theta);

    G4double halfLengthZ;

    // from suppressed
    // Dimensions
    G4double fInnerAbsorberThickness; // thickness of Cu and Sn layer
    G4double fOuterAbsorberThickness; // thickness of Ta layer
    G4double fSideLength; // x, y dimension of absorber
    G4double fRadialDistance; // distance from origin to detector facing side of abs

    G4double fGriffinCoords[16][5];

    G4bool fSurfCheck;

};

#endif
