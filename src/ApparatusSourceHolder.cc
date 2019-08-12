#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "Global.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ApparatusSourceHolder.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


ApparatusSourceHolder::ApparatusSourceHolder() :
    // LogicalVolumes
    fCeramicPelletLog(0),
    fSourceSphereLog(0),
    fSourceSupportLog(0)
{

/////////////////////////////////////////////////////////////////////
// Defining materials 
/////////////////////////////////////////////////////////////////////
    
    // instancing G4NistManager
    nistMan = G4NistManager::Instance();
    nistMan->SetVerbose(0);

    G4double density;            
    std::vector<G4int>  natoms;
    std::vector<G4String> elements;
    
    // Defining new material (Tungsten Carbide) for pellet
    elements.push_back("W");    natoms.push_back(1);
    elements.push_back("C");    natoms.push_back(1);

    density = 15.63 * g/cm3;
    G4Material* WC = nistMan->ConstructNewMaterial("WC", elements, natoms, density);

    elements.clear();
    natoms.clear();

    // Defining Peek for sphere
    elements.push_back("C");    natoms.push_back(19);
    elements.push_back("H");    natoms.push_back(12);
    elements.push_back("O");    natoms.push_back(3);

    density = 1.31 * g/cm3;
    G4Material* Peek = nistMan->ConstructNewMaterial("Peek", elements, natoms, density);

    elements.clear();
    natoms.clear();

    // defining Source plastic for support structure
    G4Material* delrin = nistMan->FindOrBuildMaterial("G4_POLYOXYMETHYLENE");

/////////////////////////////////////////////////////////////////////
// Defining physical parameters for geometries 
/////////////////////////////////////////////////////////////////////

    fullPhi                     = 360.0*deg;
    fullTheta                   = 180.0*deg;

    fPelletLength               = .292*2.54*cm;
    fPelletInnerRadius         = 0.*cm; 
    fPelletOuterRadius         = .25*2.54*cm; 

    fSphereInnerRadius         = 0.*cm; 
    fSphereOuterRadius         = .75*2.54*cm; 

    fLargeCyclinderLength       = 23.2*2.54*cm;
    fLargeCyclinderInnerRadius  = 0.*cm;
    fLargeCyclinderOuterRadius  = 1.25*2.54*cm;

    fMediumCyclinderLength      = 2.0*2.54*cm;
    fMediumCyclinderInnerRadius = 0.*cm;
    fMediumCyclinderOuterRadius = 0.5*2.54*cm;

    fSmallCyclinderLength       = 4.5*2.54*cm;
    fSmallCyclinderInnerRadius  = 0.*cm;
    fSmallCyclinderOuterRadius  = 0.315*2.54*cm;

    fPelletMaterial             = WC;
    fSphereMaterial             = Peek;
    fSupportMaterial            = delrin;
    
    // Left over lengths from LaBr Detector file, 
    // leaving in because some are used for placing the source holder
    
    fCrystalLengthZ            = 2.0*2.54*cm;
    fCrystalOuterRadius         = 1.0*2.54*cm;

    fPackingLengthZ            = fCrystalLengthZ;
    fPackingInnerRadius         = fCrystalOuterRadius;
    fPackingOuterRadius         = 0.0625*2.54*cm + fCrystalOuterRadius;

    fPackingLidInnerRadius    = 0.0*cm;
    fPackingLidOuterRadius      = fPackingOuterRadius;
    fPackingFrontLidThickness = 0.0625*2.54*cm;

    fDiscLidInnerRadius       = 0.0*cm;
    fDiscLidOuterRadius       = fPackingLidOuterRadius;
    fDiscFrontLidThickness = 0.006*2.54*cm;

    fSealLidInnerRadius       = 0.0*cm;
    fSealLidOuterRadius       = fPackingLidOuterRadius;
    fSealFrontLidThickness   = 0.04*2.54*cm;

    fCanLengthZ                = fCrystalLengthZ + fPackingFrontLidThickness + fDiscFrontLidThickness + fSealFrontLidThickness;
    fCanInnerRadius            = fPackingOuterRadius;
    fCanOuterRadius            = 0.019*2.54*cm + fPackingOuterRadius;

    fCanLidInnerRadius        = 0.0*cm;
    fCanLidOuterRadius        = fCanOuterRadius;
    fCanFrontLidThickness     = 0.019*2.54*cm;
    fCanBackLidThickness      = 0.019*2.54*cm;


    fDetectorLengthZ           = fCrystalLengthZ +
            fPackingFrontLidThickness +
            fDiscFrontLidThickness +
            fSealFrontLidThickness +
            fCanFrontLidThickness +
            fCanBackLidThickness;
     //G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))) )+atan((sqrt(2))/(1+sqrt(2))))*deg;
     G4double triangleThetaAngle = 54.735610317245360*deg;

     // theta
     fDetectorAngles[0][0]    = triangleThetaAngle;
     fDetectorAngles[1][0]    = triangleThetaAngle;
     fDetectorAngles[2][0]    = triangleThetaAngle;
     fDetectorAngles[3][0]    = triangleThetaAngle;
     fDetectorAngles[4][0]    = 180.0*deg - triangleThetaAngle;
     fDetectorAngles[5][0]    = 180.0*deg - triangleThetaAngle;
     fDetectorAngles[6][0]    = 180.0*deg - triangleThetaAngle;
     fDetectorAngles[7][0]    = 180.0*deg - triangleThetaAngle;
     // phi
     fDetectorAngles[0][1]    = 22.5*deg;
     fDetectorAngles[1][1]    = 112.5*deg;
     fDetectorAngles[2][1]    = 202.5*deg;
     fDetectorAngles[3][1]    = 292.5*deg;
     fDetectorAngles[4][1]    = 22.5*deg;
     fDetectorAngles[5][1]    = 112.5*deg;
     fDetectorAngles[6][1]    = 202.5*deg;
     fDetectorAngles[7][1]    = 292.5*deg;
     // yaw (alpha)
     fDetectorAngles[0][2]    = 0.0*deg;
     fDetectorAngles[1][2]    = 0.0*deg;
     fDetectorAngles[2][2]    = 0.0*deg;
     fDetectorAngles[3][2]    = 0.0*deg;
     fDetectorAngles[4][2]    = 0.0*deg;
     fDetectorAngles[5][2]    = 0.0*deg;
     fDetectorAngles[6][2]    = 0.0*deg;
     fDetectorAngles[7][2]    = 0.0*deg;
     // pitch (beta)
     fDetectorAngles[0][3]    = triangleThetaAngle;
     fDetectorAngles[1][3]    = triangleThetaAngle;
     fDetectorAngles[2][3]    = triangleThetaAngle;
     fDetectorAngles[3][3]    = triangleThetaAngle;
     fDetectorAngles[4][3]    = 180.0*deg - triangleThetaAngle;
     fDetectorAngles[5][3]    = 180.0*deg - triangleThetaAngle;
     fDetectorAngles[6][3]    = 180.0*deg - triangleThetaAngle;
     fDetectorAngles[7][3]    = 180.0*deg - triangleThetaAngle;
     // roll (gamma)
     fDetectorAngles[0][4]    = 22.5*deg;
     fDetectorAngles[1][4]    = 112.5*deg;
     fDetectorAngles[2][4]    = 202.5*deg;
     fDetectorAngles[3][4]    = 292.5*deg;
     fDetectorAngles[4][4]    = 22.5*deg;
     fDetectorAngles[5][4]    = 112.5*deg;
     fDetectorAngles[6][4]    = 202.5*deg;
     fDetectorAngles[7][4]    = 292.5*deg;
}

ApparatusSourceHolder::~ApparatusSourceHolder() {
    // LogicalVolumes
    delete fCeramicPelletLog;
    delete fSourceSphereLog;
    delete fSourceSupportLog;
}

G4int ApparatusSourceHolder::Build() {

    // Build assembly volume
    fAssembly = new G4AssemblyVolume();

    // Check surfaces to determine any problematic overlaps. Turn this on to have Geant4 check the surfaces.
    // Do not leave this on, it will slow the DetectorConstruction process!
    // This was last check on July 26, 2017. - GOOD!
    fSurfCheck = false;

    // Use this when simulating first source (pre July 2019)
    //G4cout << "BuildPelletVolume" << G4endl;
    //BuildCeramicPelletVolume();
    
    G4cout << "BuildSourceSphereVolume" << G4endl;
    BuildSourceSphereVolume();
    G4cout << "BuildSourceSupportVolume" << G4endl;
    BuildSourceSupportVolume();

    return 1;
}

G4int ApparatusSourceHolder::Place(G4LogicalVolume* expHallLog, G4int detectorNumber) {
    G4int detectorCopyID = 0;

    G4cout << "SourceHolder Position Number = " << detectorNumber << G4endl;

    fCopyNumber = detectorCopyID + detectorNumber;

    G4double position = 0.;

    G4double theta  = fDetectorAngles[detectorNumber][0];
    G4double phi    = fDetectorAngles[detectorNumber][1];
    //G4double alpha  = fDetectorAngles[detectorNumber][2]; // yaw
    G4double beta   = fDetectorAngles[detectorNumber][3]; // pitch
    G4double gamma  = fDetectorAngles[detectorNumber][4]; // roll

    G4double x = 0;
    G4double y = 0;
    G4double z = position;

    G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
    rotate->rotateY(M_PI);
    rotate->rotateY(M_PI+beta);
    rotate->rotateZ(gamma);

    G4ThreeVector move(TransX(x,y,z,theta,phi), TransY(x,y,z,theta,phi), TransZ(x,y,z,theta));

    fAssembly->MakeImprint(expHallLog, move, rotate, fCopyNumber, fSurfCheck);
    return 1;
}

G4int ApparatusSourceHolder::BuildCeramicPelletVolume() {
    G4Material* material = fPelletMaterial; 
    if( !material ) {
        G4cout << " ----> Material " << fPelletMaterial << " not found, cannot build the ceramic pellet! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
    visAtt->SetVisibility(true);
//    visAtt->SetForceAuxEdgeVisible(true);

    G4Tubs* pellet = BuildCeramicPellet();

    // Define rotation and movement objects
    G4ThreeVector direction   = G4ThreeVector(0,0,1);
    G4double zPosition         = .051*2.54*cm; 
    G4ThreeVector move         = zPosition * direction;
    G4RotationMatrix* rotate    = new G4RotationMatrix;

    //logical volume
    if( fCeramicPelletLog == NULL ) {
        fCeramicPelletLog = new G4LogicalVolume(pellet, material, "CeramicPelletLog", 0, 0, 0);
        fCeramicPelletLog->SetVisAttributes(visAtt);
    }

    fAssembly->AddPlacedVolume(fCeramicPelletLog, move, rotate);

    return 1;
}

G4int ApparatusSourceHolder::BuildSourceSphereVolume() {
    G4Material* material = fSphereMaterial;
    if( !material ) {
        G4cout << " ----> Material " << fSphereMaterial << " not found, cannot build the source holder sphere! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(52., 61., 87.));
    visAtt->SetVisibility(true);
//    visAtt->SetForceAuxEdgeVisible(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double zPosition;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Subtract Geometries
    /////////////////////////////////////////////////////////////////////
    sourceSphere = BuildSourceSphere();
    cyc          = BuildCeramicPellet();

    G4ThreeVector fTrans(0.,0.,0.051*2.54*cm);

    //logical volume
    if( fSourceSphereLog == NULL ) {

        // First source (pre July 2019) with ceramic backing
        //delSphere = new G4SubtractionSolid("delSphere",sphere, cyc, fRot, fTrans);
        //fSourceSphereLog = new G4LogicalVolume(delSphere, material, "SourceSphereLog", 0, 0, 0);
        
        // Second source (post July 2019) without ceramic backing
        fSourceSphereLog = new G4LogicalVolume(sourceSphere, material, "SourceSphereLog", 0, 0, 0);
        fSourceSphereLog->SetVisAttributes(visAtt);
    }

    // place front canLid
    zPosition  = 0;
    move       = zPosition * direction;

    //add physical cylinder
    fAssembly->AddPlacedVolume(fSourceSphereLog, move, rotate);

    return 1;
}

G4int ApparatusSourceHolder::BuildSourceSupportVolume() {
    G4Material* material = fSupportMaterial;
    if( !material ) {
        G4cout << " ----> Material " << fSupportMaterial << " not found, cannot build the Source Support Structure! " << G4endl;
        return 0;
    }

    // Set visualization attributes
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1., 1., 1.));
    visAtt->SetVisibility(true);

    G4ThreeVector direction =  G4ThreeVector(0,0,1);
    G4double zPosition;
    G4ThreeVector move;
    G4RotationMatrix* rotate = new G4RotationMatrix;

    /////////////////////////////////////////////////////////////////////
    // Build and Subtract Geometries
    /////////////////////////////////////////////////////////////////////
    G4Tubs* lrgCyc = BuildLargeSourceCyclinder();
    G4Tubs* medCyc = BuildMediumSourceCyclinder();
    G4Tubs* smlCyc = BuildSmallSourceCyclinder();

    G4ThreeVector fLrgMedTrans(0.,0.,-11.6*2.54*cm);
    G4RotationMatrix* fRot = new G4RotationMatrix;

    G4UnionSolid* lrgMedSupport = new G4UnionSolid("LrgMedSupport",lrgCyc, medCyc, fRot, fLrgMedTrans);

    G4ThreeVector fSmlMedTrans(0.,0.,-13.6*2.54*cm);
    G4RotationMatrix* fRotZ = new G4RotationMatrix;
//    fRotZ->rotateX(M_PI*rad);

    G4UnionSolid* totalSupport = new G4UnionSolid("TotalSupport",lrgMedSupport, smlCyc, fRotZ, fSmlMedTrans);

    //logical volume
    if( fSourceSupportLog == NULL ) {
        fSourceSupportLog = new G4LogicalVolume(totalSupport, material, "SourceSupportLog", 0, 0, 0);
        fSourceSupportLog->SetVisAttributes(visAtt);
    }

    // place front canLid
    zPosition  = 18.85*2.54*cm;
    move       = zPosition * direction;

    //add physical cylinder
    fAssembly->AddPlacedVolume(fSourceSupportLog, move, rotate);

    return 1;
}
///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Tubs* ApparatusSourceHolder::BuildCeramicPellet() {
    G4double startPhi   = 0.0;
    G4double endPhi     = fullPhi;

    G4double innerRadius = fPelletInnerRadius;
    G4double outerRadius = fPelletOuterRadius;
    G4double halfLengthZ = fPelletLength/2.0;

    G4Tubs* pellet = new G4Tubs("pellet", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

    return pellet;
}//end ::BuildPellet


G4Sphere* ApparatusSourceHolder::BuildSourceSphere() {
    G4double startPhi   = 0.0;
    G4double endPhi     = fullPhi;
    G4double startTheta = 0.0;
    G4double endTheta   = fullTheta;
    

    G4double innerRadius   = fSphereInnerRadius;
    G4double outerRadius   = fSphereOuterRadius;

    G4Sphere* sphere = new G4Sphere("Sphere", innerRadius, outerRadius, startPhi, endPhi, startTheta, endTheta);

    return sphere;
}//end ::BuildAluminumCan

G4Tubs* ApparatusSourceHolder::BuildLargeSourceCyclinder() {
    G4double startPhi   = 0.0;
    G4double endPhi     = fullPhi;

    G4double innerRadius   = fLargeCyclinderInnerRadius;
    G4double outerRadius   = fLargeCyclinderOuterRadius;
    G4double halfLengthZ    = fLargeCyclinderLength/2.0;

    G4Tubs* largeCyclinder = new G4Tubs("LargeCyclinder", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

    return largeCyclinder;
}//end ::BuildLargeSourceCyclinder

G4Tubs* ApparatusSourceHolder::BuildMediumSourceCyclinder() {
    G4double startPhi   = 0.0;
    G4double endPhi     = fullPhi;

    G4double innerRadius   = fMediumCyclinderInnerRadius;
    G4double outerRadius   = fMediumCyclinderOuterRadius;
    G4double halfLengthZ    = fMediumCyclinderLength;

    G4Tubs* medCyclinder = new G4Tubs("MediumCyclinder", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

    return medCyclinder;
}//end ::BuildMediumSourceCyclinder

G4Tubs* ApparatusSourceHolder::BuildSmallSourceCyclinder() {
    G4double startPhi   = 0.0;
    G4double endPhi     = fullPhi;

    G4double innerRadius   = fSmallCyclinderInnerRadius;
    G4double outerRadius   = fSmallCyclinderOuterRadius;
    G4double halfLengthZ    = fSmallCyclinderLength;

    G4Tubs* smallCyclinder = new G4Tubs("SmallCyclinder", innerRadius, outerRadius, halfLengthZ, startPhi, endPhi);

    return smallCyclinder;
}//end ::BuildSmallSourceCyclinder

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector ApparatusSourceHolder::GetDirectionXYZ(G4double theta, G4double phi) {
    G4double x,y,z;
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);

    G4ThreeVector direction = G4ThreeVector(x,y,z);

    return direction;
}//end ::GetDirection
