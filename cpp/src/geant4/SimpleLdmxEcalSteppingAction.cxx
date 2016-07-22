
#include <SimpleLdmxEcalSteppingAction.h>

SimpleLdmxEcalSteppingAction::SimpleLdmxEcalSteppingAction() 
    : G4UserSteppingAction() { 

}

SimpleLdmxEcalSteppingAction::~SimpleLdmxEcalSteppingAction() {

}

void SimpleLdmxEcalSteppingAction::UserSteppingAction(const G4Step* step) { 
   
    std::cout << "*************" << std::endl; 
    std::cout << "*   Event   *" << std::endl;
    std::cout << "*************" << std::endl; 
    
    // 
    double incident_particle_energy = step->GetPreStepPoint()->GetTotalEnergy();
    std::cout << "Total energy of incident particle: " << incident_particle_energy << std::endl;

    // Get the track associated with this step
    //G4Track* track = step->GetTrack();
    const G4TrackVector* secondaries = step->GetSecondary();

    // If the initial interaction results in no secondaries, skip the rest pf
    // the event.
    if (secondaries->size() == 0) return;


    std::cout << "Total number of secondaries: " << secondaries->size() << std::endl;
   
    // Get the track associated with this step
    G4Track* track = step->GetTrack(); 

    // Check if the track survived
    bool trackIsAlive = (track->GetTrackStatus() == fAlive);

    // Postpone the tracking of all particles
    if (trackIsAlive) track->SetTrackStatus(fPostponeToNextEvent); 
    for (G4TrackVector::const_iterator secondary = secondaries->begin(); 
            secondary != secondaries->end(); ++secondary) { 
        (*secondary)->SetTrackStatus(fPostponeToNextEvent);
    }

    G4String process_name = secondaries->at(0)->GetCreatorProcess()->GetProcessName(); 
    std::cout << "Process name: " << process_name << std::endl;

    return;
}
