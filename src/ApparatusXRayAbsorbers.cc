#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4UnionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ApparatusXRayAbsorbers.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

ApparatusXRayAbsorbers::~ApparatusXRayAbsorbers() {
   delete fAbsorber;
}

G4int ApparatusXRayAbsorbers::Build() { 
   // for defining Materials
   G4NistManager* nistMan = G4NistManager::Instance();
   nistMan->SetVerbose(0);

   // material closest to detector face
   fAbsorberLayer1Mat = nistMan->FindOrBuildMaterial("G4_Cu");    
   fAbsorberLayer2Mat = nistMan->FindOrBuildMaterial("G4_Sn");    
   fAbsorberLayer3Mat = nistMan->FindOrBuildMaterial("G4_Ta");    

   // Setup assembly volumes
   fAssemblyAbsorber = new G4AssemblyVolume();

   //G4cout<<"BuildLayer1"<<G4endl;
   BuildLayer1();

   //G4cout<<"BuildLayer2"<<G4endl;
   BuildLayer2();

   ////G4cout<<"BuildLayer3"<<G4endl;
   BuildLayer3();

   return 1;
}//end ::Build


G4int ApparatusXRayAbsorbers::Place(G4LogicalVolume* expHallLog, G4int selector) {
   if(expHallLog == nullptr) {
      std::cerr<<__PRETTY_FUNCTION__<<": expHallLog == nullptr!"<<std::endl;
      exit(1);
   }

   G4double theta;
   G4double phi;
   G4double alpha;
   G4double beta;
   G4double gamma;

   G4int iAbsorber1 = 0;
   G4int iAbsorber2 = 15;

   if(selector == 0) { // Full Absorbers
      iAbsorber1 = 0;
      iAbsorber2 = 15;
   } else if(selector == 1) { // kept in for subset of absorbers (if needed)
      // only build 1 
      iAbsorber1 = 0;
      iAbsorber2 = 0;
   } 

  // loop over all pieces
  for(G4int positionNumber=iAbsorber1; positionNumber<=iAbsorber2; positionNumber++) {

     // position parameters 
      theta  = fGriffinCoords[positionNumber][0];
      phi    = fGriffinCoords[positionNumber][1];
      alpha  = fGriffinCoords[positionNumber][2]; // yaw
      beta   = fGriffinCoords[positionNumber][3]; // pitch
      gamma  = fGriffinCoords[positionNumber][4]; // roll

      G4RotationMatrix* rotate = new G4RotationMatrix;
      rotate->rotateY(M_PI/2.0);
      rotate->rotateX(M_PI/2.0);
      rotate->rotateX(alpha);
      rotate->rotateY(beta);
      rotate->rotateZ(gamma);

      G4double distFromOriginDet = fAirBoxBackLengthDet/2.0 +fAirBoxFrontLengthDet + fNewRhombiRadiusDet;

      G4double x = 0;
      G4double y = 0;
      G4double z = distFromOriginDet;

      G4ThreeVector move(ApparatusXRayAbsorbers::TransX(x,y,z,theta,phi), ApparatusXRayAbsorbers::TransY(x,y,z,theta,phi), ApparatusXRayAbsorbers::TransZ(x,z,theta));
      G4ThreeVector moveFrontFace(fShift + fAppliedBackShift, 0, 0);

      fAssemblyAbsorber->MakeImprint(expHallLog, moveFrontFace, rotate, 0, fSurfCheck);

   }

   return 1;
}

G4int ApparatusXRayAbsorbers::BuildLayer1() {
   if(!fAbsorberLayer1Mat) {
      G4cout<<" ----> Material "<< fAbsorberLayer1Mat <<" not found, cannot build layer 1 of the X-Ray absorbers! "<<G4endl;
      return 0;
   }

   // Set visualization attributes
   G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
   visAtt->SetVisibility(true);

   G4Box* layer1 = AbsorberLayer(1);

   // Define rotation and movement objects
   G4ThreeVector direction = G4ThreeVector(0,0,1);
   // place holder value for testing
   G4double zPosition      = -fInnerAbsorberThickness/2.0;
   G4ThreeVector move      = zPosition * direction;

   G4RotationMatrix* rotate = new G4RotationMatrix;

   //logical volume
   if(fAbsorberLayer1Log == nullptr) {
      fAbsorberLayer1Log = new G4LogicalVolume(layer1, fAbsorberLayer1Mat, "layer1Log", 0, 0, 0);
      fAbsorberLayer1Log->SetVisAttributes(visAtt);
   }

   fAssemblyAbsorber->AddPlacedVolume(fAbsorberLayer1Log, move, rotate);

   return 1;
}

G4int ApparatusXRayAbsorbers::BuildLayer2() {
   if(!fAbsorberLayer2Mat) {
      G4cout<<" ----> Material "<< fAbsorberLayer2Mat <<" not found, cannot build layer 2 of the X-Ray absorbers! "<<G4endl;
      return 0;
   }

   // Set visualization attributes
   G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
   visAtt->SetVisibility(true);

   G4Box* layer = AbsorberLayer(2);

   // Define rotation and movement objects
   G4ThreeVector direction = G4ThreeVector(0,0,1);
   // place holder value for testing
   G4double zPosition      = -fInnerAbsorberThickness*3/2.;
   G4ThreeVector move      = zPosition * direction;

   G4RotationMatrix* rotate = new G4RotationMatrix;

   //logical volume
   if(fAbsorberLayer2Log == nullptr) {
      fAbsorberLayer2Log = new G4LogicalVolume(layer, fAbsorberLayer2Mat, "layer2Log", 0, 0, 0);
      fAbsorberLayer2Log->SetVisAttributes(visAtt);
   }

   fAssemblyAbsorber->AddPlacedVolume(fAbsorberLayer2Log, move, rotate);

   return 1;
}

G4int ApparatusXRayAbsorbers::BuildLayer3() {
   if(!fAbsorberLayer3Mat) {
      G4cout<<" ----> Material "<< fAbsorberLayer3Mat <<" not found, cannot build layer 3 of the X-Ray absorbers! "<<G4endl;
      return 0;
   }

   // Set visualization attributes
   G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
   visAtt->SetVisibility(true);

   G4Box* layer = AbsorberLayer(3);

   // Define rotation and movement objects
   G4ThreeVector direction = G4ThreeVector(0,0,1);
   // place holder value for testing
   G4double zPosition      = -fInnerAbsorberThickness*2. - fOuterAbsorberThickness/2.0;
   G4ThreeVector move      = zPosition * direction;

   G4RotationMatrix* rotate = new G4RotationMatrix;

   //logical volume
   if(fAbsorberLayer3Log == nullptr) {
      fAbsorberLayer3Log = new G4LogicalVolume(layer, fAbsorberLayer3Mat, "layer3Log", 0, 0, 0);
      fAbsorberLayer3Log->SetVisAttributes(visAtt);
   }

   fAssemblyAbsorber->AddPlacedVolume(fAbsorberLayer3Log, move, rotate);

   return 1;
}

G4Box* ApparatusXRayAbsorbers::AbsorberLayer(G4int selector) {

   G4cout << "selector " << selector << G4endl;

   G4double halfLengthX = fSideLength/2.0; // half length of square side
   G4double halfLengthY  = halfLengthX;

   // setting layer thickness
   if(selector == 1) { // first layers
      halfLengthZ  = fInnerAbsorberThickness/2.0;
   } else if(selector == 2) { 
      halfLengthZ  = fInnerAbsorberThickness/2.0;
   } else if(selector == 3) { 
      halfLengthZ  = fOuterAbsorberThickness/2.0;
   }

   G4Box* xRayAbsorberLayer = new G4Box("xRayAbsorberLayer", halfLengthX, halfLengthY, halfLengthZ);

   return xRayAbsorberLayer;
}

G4double ApparatusXRayAbsorbers::TransX(G4double x, G4double y, G4double z, G4double theta, G4double phi) {

   return (x*cos(theta)+z*sin(theta))*cos(phi)-y*sin(phi);
}

G4double ApparatusXRayAbsorbers::TransY(G4double x, G4double y, G4double z, G4double theta, G4double phi) {
   return (x*cos(theta)+z*sin(theta))*sin(phi)+y*cos(phi);
}

G4double ApparatusXRayAbsorbers::TransZ(G4double x, G4double z, G4double theta) {
   return -x*sin(theta)+z*cos(theta);
}
