
#include <SimpleLdmxEcalPrimaryGeneratorAction.h>

SimpleLdmxEcalPrimaryGeneratorAction::SimpleLdmxEcalPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      particle_gun(new G4ParticleGun(1 /* Use only a single particle */)) { 
    
    // Define the kinematics of the particle
    G4ParticleDefinition* particle_definition 
        = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    particle_gun->SetParticleDefinition(particle_definition);
    particle_gun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.)); 
    particle_gun->SetParticleEnergy(2.8*GeV);
}

SimpleLdmxEcalPrimaryGeneratorAction::~SimpleLdmxEcalPrimaryGeneratorAction() { 
    delete particle_gun; 
}

void SimpleLdmxEcalPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) { 
  
    G4double world_z_half_length = 0;
    G4LogicalVolume* world_logical_vol = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    G4Box* world_box = 0;
    
    if ( world_logical_vol) world_box = dynamic_cast< G4Box*>(world_logical_vol->GetSolid()); 
    
    if ( world_box ) {
        world_z_half_length = world_box->GetZHalfLength();  
    }else  {
        G4ExceptionDescription msg;
        msg << "World volume of box not found." << G4endl;
        msg << "Perhaps you have changed geometry." << G4endl;
        msg << "The gun will be place in the center.";
        G4Exception("BasicCalPrimaryGeneratorAction::GeneratePrimaries()",
        "MyCode0002", JustWarning, msg);
    } 
  
  
    // Set gun position
    particle_gun->SetParticlePosition(G4ThreeVector(0., 0., (-world_z_half_length + 20)));

    particle_gun->GeneratePrimaryVertex(event);
}


