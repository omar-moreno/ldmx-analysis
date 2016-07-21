
#include <SimpleLdmxEcalSteppingAction.h>

SimpleLdmxEcalSteppingAction::SimpleLdmxEcalSteppingAction() 
    : G4UserSteppingAction() { 

}

SimpleLdmxEcalSteppingAction::~SimpleLdmxEcalSteppingAction() {

}

void SimpleLdmxEcalSteppingAction::UserSteppingAction(const G4Step* step) { 
   
    // Get the track associated with this step
    //G4Track* track = step->GetTrack();
    const std::vector<G4Track*>* secondaries = step->GetSecondary();

    // If the initial interaction results in no secondaries, skip the rest pf
    // the event.
    if (secondaries->size() == 0) return;

    std::cout << "*************" << std::endl; 
    std::cout << "*   Event   *" << std::endl;
    std::cout << "*************" << std::endl; 

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

    return;
}
