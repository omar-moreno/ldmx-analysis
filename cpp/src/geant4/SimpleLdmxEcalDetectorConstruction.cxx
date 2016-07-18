
#include <SimpleLdmxEcalDetectorConstruction.h>

SimpleLdmxEcalDetectorConstruction::SimpleLdmxEcalDetectorConstruction() 
    : G4VUserDetectorConstruction() {
}

SimpleLdmxEcalDetectorConstruction::~SimpleLdmxEcalDetectorConstruction() {
}

G4VPhysicalVolume* SimpleLdmxEcalDetectorConstruction::Construct() { 

    // Define materials
    this->DefineMaterials();

    // Define volumes
    return this->DefineVolumes(); 

}

void SimpleLdmxEcalDetectorConstruction::ConstructSDandField() { 



}

//--------------------//
//  Private Methods   //
//--------------------//

void SimpleLdmxEcalDetectorConstruction::DefineMaterials() { 

    // Use the NIST Manager to define some basic materials
    G4NistManager* nist_mng = G4NistManager::Instance();
    nist_mng->FindOrBuildMaterial("G4_W");  // Tungsten
    nist_mng->FindOrBuildMaterial("G4_Si"); // Silicon
    
    G4double a;  // Mass of a mole;
    G4double z;  // Mean number of protons;  
    G4double density; // Density of material
    
    // Vacuum
    new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* SimpleLdmxEcalDetectorConstruction::DefineVolumes() { 
  
    G4Material* defaultMaterial  = G4Material::GetMaterial("Galactic"); 
    G4Material* absorberMaterial = G4Material::GetMaterial("G4_W");
    G4Material* gapMaterial      = G4Material::GetMaterial("G4_Si");

    // World and detector dimensions
  
    G4double W_X0            = 0.35*cm;
    G4int nofLayers          = 40;
    G4double calorSizeXY     = 10.*cm;
    G4double absoThickness   = 1.*W_X0;
    G4double gapThickness    = 0.05*cm;
    G4double layerThickness  = absoThickness + gapThickness;
    G4double calorThickness  = nofLayers * layerThickness;
    
    G4double worldSizeXY     = 1.5 * calorSizeXY;
    G4double worldSizeZ      = 2.2 * calorThickness; 
   
    //-----------//
    //   World   //
    //-----------//
    
    // Dimensions
  
    G4VSolid* worldS  = new G4Box("World", worldSizeXY, worldSizeXY, worldSizeZ);
                         
    G4LogicalVolume* worldLV = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
    G4VPhysicalVolume* worldPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 true);            // checking overlaps 

    //-----------------//
    //   Calorimeter   //
    //-----------------//
    G4ThreeVector posCal = G4ThreeVector(0, 0, 16*cm);
    

    G4VSolid* calorimeterS = new G4Box("Calorimeter", calorSizeXY/2, calorSizeXY/2, calorThickness/2); 
                         
    G4LogicalVolume* calorLV = new G4LogicalVolume(
                 calorimeterS,    // its solid
                 defaultMaterial, // its material
                 "Calorimeter");  // its name
                                   
    new G4PVPlacement(0,                // no rotation
                      posCal,           // its position
                      calorLV,          // its logical volume                         
                      "Calorimeter",    // its name
                      worldLV,          // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      true);  // checking overlaps 

    // Layer
    G4VSolid* layerS = new G4Box("Layer", calorSizeXY/2, calorSizeXY/2, layerThickness/2);
                         
    G4LogicalVolume* layerLV = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

     new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness);  // witdth of replica
  
    // Absorber
    G4VSolid* absorberS = new G4Box("Abso", calorSizeXY/2, calorSizeXY/2, absoThickness/2);
                         
    G4LogicalVolume* absorberLV = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name
                                   
    new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness/2), //  its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 true);  // checking overlaps 

    // Gap
    G4VSolid* gapS = new G4Box("Gap", calorSizeXY/2, calorSizeXY/2, gapThickness/2);
                         
    G4LogicalVolume* gapLV = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name
                                   
    new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness/2), //  its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 true);  // checking overlaps 

    // Always return the physical World
    return worldPV;
}
