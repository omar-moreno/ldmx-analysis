
#include <HydrogenEcalDetectorConstruction.h>

HydrogenEcalDetectorConstruction::HydrogenEcalDetectorConstruction() 
    : G4VUserDetectorConstruction() {
}

HydrogenEcalDetectorConstruction::~HydrogenEcalDetectorConstruction() {
}

G4VPhysicalVolume* HydrogenEcalDetectorConstruction::Construct() { 

    // Define materials
    this->DefineMaterials();

    // Define volumes
    return this->DefineVolumes(); 

}

void HydrogenEcalDetectorConstruction::ConstructSDandField() { 



}

//--------------------//
//  Private Methods   //
//--------------------//

void HydrogenEcalDetectorConstruction::DefineMaterials() { 

    // Use the NIST Manager to define some basic materials
    G4NistManager* nist_mng = G4NistManager::Instance();
    nist_mng->FindOrBuildMaterial("G4_lH2"); // Liquid hydrogen
    
    G4double a;  // Mass of a mole;
    G4double z;  // Mean number of protons;  
    G4double density; // Density of material
    
    // Vacuum
    new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* HydrogenEcalDetectorConstruction::DefineVolumes() { 
 
    //--------------------------// 
    //   Material definitions   //
    //--------------------------// 
    G4Material* vacuum  = G4Material::GetMaterial("Galactic");
    G4Material* absorber_material = G4Material::GetMaterial("G4_lH2");
    G4Material* gap_material = G4Material::GetMaterial("G4_lH2");
    

    // Radiation length of liquid hydrogen
    G4double h_x0 = 890.4*cm;

    // Number of Ecal layers
    G4int n_layers = 1;

    // Thickness of the Ecal absorbers
    G4double abs_thickness = 1.*h_x0;

    // Thickness of the gap between absorber layers
    G4double gap_thickness = 0.05*cm;

    // Thickness of a calorimeter layer
    G4double layer_thickness = abs_thickness + gap_thickness;
    
    // Thickness of calorimeter
    G4double ecal_thickness = n_layers * layer_thickness;

    // Size of calorimeter
    G4double ecal_size_xy = 10.*h_x0;

    //-----------//
    //   World   //
    //-----------//
    
    // Dimensions

    G4double world_size_xy     = 1.5 * ecal_size_xy;
    G4double world_size_z      = 2.2 * ecal_thickness; 

    G4VSolid* world_solid  
        = new G4Box("World", world_size_xy, world_size_xy, world_size_z);
                         
    G4LogicalVolume* world_logical_volume 
        = new G4LogicalVolume(world_solid, vacuum, "World");
                                   
    G4VPhysicalVolume* world_placement 
        = new G4PVPlacement(0, G4ThreeVector(), world_logical_volume, "World", 0, false, 0, true);

    //-----------------//
    //   Calorimeter   //
    //-----------------//
    
    G4ThreeVector ecal_pos = G4ThreeVector(0, 0, 1000*cm); 

    G4VSolid* ecal_solid 
        = new G4Box("Calorimeter", ecal_size_xy/2, ecal_size_xy/2, ecal_thickness/2); 
                         
    G4LogicalVolume* ecal_logical_volume 
        = new G4LogicalVolume(ecal_solid, vacuum, "Calorimeter");
                                   
    new G4PVPlacement(0, ecal_pos, ecal_logical_volume, "Calorimeter", world_logical_volume, false, 0, true); 

    G4VSolid* ecal_layer_solid 
        = new G4Box("Layer", ecal_size_xy/2, ecal_size_xy/2, layer_thickness/2);
                         
    G4LogicalVolume* ecal_layer_logical_volume 
        = new G4LogicalVolume(ecal_layer_solid, vacuum, "Layer");

     new G4PVReplica("Layer", ecal_layer_logical_volume, ecal_logical_volume, kZAxis, n_layers, layer_thickness); 
  
    // Absorber
    G4VSolid* abs_solid 
        = new G4Box("Abso", ecal_size_xy/2, ecal_size_xy/2, abs_thickness/2);
                         
    G4LogicalVolume* abs_logical_volume 
        = new G4LogicalVolume(abs_solid, absorber_material, "AbsoLV");       
                                   
    new G4PVPlacement(0, G4ThreeVector(0., 0., -gap_thickness/2), abs_logical_volume, "Abso", ecal_layer_logical_volume, false, 0, true);  

    // Gap
    G4VSolid* gap_solid 
        = new G4Box("Gap", ecal_size_xy/2, ecal_size_xy/2, gap_thickness/2);
                         
    G4LogicalVolume* gapLV 
        = new G4LogicalVolume(gap_solid, gap_material, "GapLV");     
                                   
    new G4PVPlacement(0, G4ThreeVector(0., 0., abs_thickness/2), gapLV, "Gap", ecal_layer_logical_volume, false, 0, true); 

    // Always return the physical World
    return world_placement;
}
