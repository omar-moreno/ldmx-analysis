
#include <HydrogenEcalRunAction.h>

HydrogenEcalRunAction::HydrogenEcalRunAction()
    : G4UserRunAction() {

    // Instantiate the analysis manager.  The type of analysis manager is 
    // chosen by selecting the namespace in the header GeantAnalysis.h
    G4AnalysisManager* analysis_manager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysis_manager->GetType() << G4endl;

    // Print event number 
    G4RunManager::GetRunManager()->SetPrintProgress(1000);     

    // Get the analysis manager 
    //analysis_manager->SetVerboseLevel(4);

    // Set the output file name
    analysis_manager->SetFileName("hydrogen_ecal_test");

    analysis_manager->CreateNtuple("gamma_p_to_pi_p", "kinematics");
    analysis_manager->CreateNtupleDColumn("e_gamma");
    analysis_manager->CreateNtupleDColumn("gamma_p_dir_x");
    analysis_manager->CreateNtupleDColumn("gamma_p_dir_y");
    analysis_manager->CreateNtupleDColumn("gamma_p_dir_z");
    analysis_manager->CreateNtupleIColumn("n_secondaries");
    analysis_manager->CreateNtupleDColumn("pi0_ke");
    analysis_manager->CreateNtupleDColumn("pi0_p_dir_x");
    analysis_manager->CreateNtupleDColumn("pi0_p_dir_y");
    analysis_manager->CreateNtupleDColumn("pi0_p_dir_z");
    analysis_manager->CreateNtupleDColumn("p_ke");
    analysis_manager->CreateNtupleDColumn("p_p_dir_x");
    analysis_manager->CreateNtupleDColumn("p_p_dir_y");
    analysis_manager->CreateNtupleDColumn("p_p_dir_z");
    analysis_manager->FinishNtuple(); 
}

HydrogenEcalRunAction::~HydrogenEcalRunAction() { 

    delete G4AnalysisManager::Instance();
}

void HydrogenEcalRunAction::BeginOfRunAction(const G4Run* run) { 

    // Get the analysis manager
    G4AnalysisManager* analysis_manager = G4AnalysisManager::Instance();

    // Open the output file
    analysis_manager->OpenFile();     
}

void HydrogenEcalRunAction::EndOfRunAction(const G4Run* run) { 

    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "%   Total number of events processed: " 
              << run->GetNumberOfEvent() << std::endl;
    std::cout << "%   Total number of photonuclear events: " 
              << photo_nuclear_count << std::endl;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    // Get the analysis manager
    G4AnalysisManager* analysis_manager = G4AnalysisManager::Instance();

    // Save the ntuple and close the file
    analysis_manager->Write(); 
    analysis_manager->CloseFile(); 
}

int HydrogenEcalRunAction::photo_nuclear_count = 0;
