
#include <LeadCalorimeterSteppingAction.h>

LeadCalorimeterSteppingAction::LeadCalorimeterSteppingAction() 
    : G4UserSteppingAction() { 
}

LeadCalorimeterSteppingAction::~LeadCalorimeterSteppingAction() {

}

void LeadCalorimeterSteppingAction::UserSteppingAction(const G4Step* step) { 
   
    //std::cout << "************" << std::endl; 
    //std::cout << "*   Step   *" << std::endl;
    //std::cout << "************" << std::endl; 

    // Get the track associated with this step
    G4Track* track = step->GetTrack(); 

    // Get the particle type
    G4String particle_name = track->GetParticleDefinition()->GetParticleName();

    // Get the volume
    G4VPhysicalVolume* volume = track->GetVolume();
    G4String volume_name = volume->GetName(); 

    // Get the secondaries associated with this track
    const G4TrackVector* secondaries = step->GetSecondary();

    double incident_particle_energy = step->GetPreStepPoint()->GetTotalEnergy();
    //std::cout << "Total energy of " << particle_name << " : " 
    //          << incident_particle_energy << std::endl;
    //std::cout << "Particle currently in " << volume_name << std::endl; 

    // In the case that electrons are being fired on target, stop tracking the
    // electron once it interacts in the target regardless of whether a brem
    // photon was produced or not.
    if ((particle_name.compareTo("e-") + volume_name.compareTo("Target")) == 0) { 
        
        // If the initial interaction didn't result in any secondaries e.g. a
        // brem photon, don't continue.
        if (secondaries->size() == 0) { 
            track->SetTrackStatus(fKillTrackAndSecondaries); 
            return;
        }

        G4String process_name = secondaries->at(0)->GetCreatorProcess()->GetProcessName(); 
        //std::cout << "Incident electron produced " << secondaries->size() << " particle via " 
        //          << process_name << " process." << std::endl;

        // If secondaries were produced via a process other than brem, stop 
        // tracking all secondary tracks.
        if (process_name.compareTo("eBrem") != 0) {
            track->SetTrackStatus(fKillTrackAndSecondaries); 
            for (G4TrackVector::const_iterator secondary = secondaries->begin(); 
                    secondary != secondaries->end(); ++secondary) {  
                (*secondary)->SetTrackStatus(fKillTrackAndSecondaries);
            }
        }
    }
    
    // Keep tracking the photon (whether it's from a particle gun or produced
    // via brem), until it enters the absorber.
    if (volume_name.compareTo("Abso") == 0 && particle_name.compareTo("gamma") == 0) {

        //std::cout << "Evolution of Brem photon track postponed." << std::endl;

        // Get analysis manager
        G4AnalysisManager* analysis_manager = G4AnalysisManager::Instance();
        analysis_manager->FillNtupleDColumn(0, incident_particle_energy);
        analysis_manager->FillNtupleDColumn(1, track->GetMomentumDirection().x());
        analysis_manager->FillNtupleDColumn(2, track->GetMomentumDirection().y());
        analysis_manager->FillNtupleDColumn(3, track->GetMomentumDirection().z());
        analysis_manager->FillNtupleIColumn(4, secondaries->size()); 

        // If the initial interaction didn't result in any secondaries don't 
        // continue. 
        if (secondaries->size() == 0) {
            //analysis_manager->AddNtupleRow();
            return;
        }
        G4String process_name = secondaries->at(0)->GetCreatorProcess()->GetProcessName(); 

        // Only record photonuclear events
        if (process_name.compareTo("photonNuclear") == 0) { 
            std::cout << "Incident gamma produced " << secondaries->size() << " particle via " 
                      << process_name << " process." << std::endl;

            // Increase the count of photon nuclear events
            LeadCalorimeterRunAction::IncrementPhotoNuclearCount(); 
            
            G4Track* pi0_track = nullptr;
            G4Track* proton_track = nullptr;
            if (secondaries->size() == 2) {
                for (G4TrackVector::const_iterator secondary = secondaries->begin(); 
                        secondary != secondaries->end(); ++secondary) {  
                    G4String particle_name = (*secondary)->GetParticleDefinition()->GetParticleName();
                    std::cout << "Particle name: " << particle_name << std::endl;
                    if (particle_name.compareTo("pi0") == 0) { 
                        pi0_track = (*secondary); 
                    } else if (particle_name.compareTo("proton") == 0) { 
                        proton_track = (*secondary); 
                    }
                }
            
                if (pi0_track != nullptr && proton_track != nullptr) {
                    analysis_manager->FillNtupleDColumn(5,   pi0_track->GetKineticEnergy()); 
                    analysis_manager->FillNtupleDColumn(6,   pi0_track->GetMomentumDirection().x());
                    analysis_manager->FillNtupleDColumn(7,   pi0_track->GetMomentumDirection().y());
                    analysis_manager->FillNtupleDColumn(8,   pi0_track->GetMomentumDirection().z());
                    analysis_manager->FillNtupleDColumn(9,   proton_track->GetKineticEnergy()); 
                    analysis_manager->FillNtupleDColumn(10,  proton_track->GetMomentumDirection().x());
                    analysis_manager->FillNtupleDColumn(11,  proton_track->GetMomentumDirection().y());
                    analysis_manager->FillNtupleDColumn(12,  proton_track->GetMomentumDirection().z());
                
                    analysis_manager->AddNtupleRow();
                    std::cout << "Writing row to ntuple" << std::endl;
                }
            }
        }
        
        // Check if the track survived
        if (track->GetTrackStatus() == fAlive) { 
            //track->SetTrackStatus(fKillTrackAndSecondaries); 
            track->SetTrackStatus(fKillTrackAndSecondaries); 
        }

        // Postpone the evolution of all secondaries 
        for (G4TrackVector::const_iterator secondary = secondaries->begin(); 
                secondary != secondaries->end(); ++secondary) {  

            (*secondary)->SetTrackStatus(fKillTrackAndSecondaries);
        }
    } 
    return;
}
