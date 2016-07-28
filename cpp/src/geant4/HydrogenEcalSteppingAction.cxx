
#include <HydrogenEcalSteppingAction.h>

HydrogenEcalSteppingAction::HydrogenEcalSteppingAction() 
    : G4UserSteppingAction() { 
}

HydrogenEcalSteppingAction::~HydrogenEcalSteppingAction() {

}

void HydrogenEcalSteppingAction::UserSteppingAction(const G4Step* step) { 
   
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
        track->SetTrackStatus(fPostponeToNextEvent); 
        
        // If the initial interaction didn't result in any secondaries e.g. a
        // brem photon, don't continue.
        if (secondaries->size() == 0) return;

        G4String process_name = secondaries->at(0)->GetCreatorProcess()->GetProcessName(); 
        //std::cout << "Incident electron produced " << secondaries->size() << " particle via " 
        //          << process_name << " process." << std::endl;

        // If secondaries were produced via a process other than brem, stop 
        // tracking all secondary tracks.
        if (process_name.compareTo("eBrem") != 0) {
            for (G4TrackVector::const_iterator secondary = secondaries->begin(); 
                    secondary != secondaries->end(); ++secondary) {  
                (*secondary)->SetTrackStatus(fPostponeToNextEvent);
            }
        }
    }

    // Keep tracking the photon (whether it's from a particle gun or produced
    // via brem), until it enters the absorber.
    if (volume_name.compareTo("Abso") == 0 && particle_name.compareTo("gamma") == 0) {

        // Check if the track survived
        if (track->GetTrackStatus() == fAlive) { 
            track->SetTrackStatus(fPostponeToNextEvent); 
        }
        //std::cout << "Evolution of Brem photon track postponed." << std::endl;

        // Get analysis manager
        G4AnalysisManager* analysis_manager = G4AnalysisManager::Instance();
        analysis_manager->FillNtupleDColumn(0, incident_particle_energy);

        // If the initial interaction didn't result in any secondaries don't 
        // continue. 
        if (secondaries->size() == 0) {
            analysis_manager->AddNtupleRow();
            return;
        }
        G4String process_name = secondaries->at(0)->GetCreatorProcess()->GetProcessName(); 
        
        // Postpone the evolution of all secondaries 
        for (G4TrackVector::const_iterator secondary = secondaries->begin(); 
                secondary != secondaries->end(); ++secondary) {  
            (*secondary)->SetTrackStatus(fPostponeToNextEvent);
        }

        // Only record photonuclear events
        if (process_name.compareTo("photonNuclear") == 0) { 
            std::cout << "Incident gamma produced " << secondaries->size() << " particle via " 
                      << process_name << " process." << std::endl;
            analysis_manager->FillNtupleDColumn(1, incident_particle_energy);
        }
        analysis_manager->AddNtupleRow();
    } 

    return;
}
