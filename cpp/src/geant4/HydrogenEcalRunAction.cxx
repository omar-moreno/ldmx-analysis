
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

    analysis_manager->CreateNtuple("hydrogen_ecal", "kinematics");
    analysis_manager->CreateNtupleDColumn("Egamma");
    analysis_manager->CreateNtupleDColumn("Egamma_pn");
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

    // Get the analysis manager
    G4AnalysisManager* analysis_manager = G4AnalysisManager::Instance();

    // Save the ntuple and close the file
    analysis_manager->Write(); 
    analysis_manager->CloseFile(); 
}
