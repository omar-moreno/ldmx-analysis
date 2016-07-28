/**
 *
 */


//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>
#include <getopt.h>
#include <cstdlib>

//------------//
//   Geant4   //
//------------//
#include "G4RunManager.hh"
#include "G4VModularPhysicsList.hh"
#include "QGSP_BERT.hh"
#include "G4UImanager.hh"

#include <HydrogenEcalDetectorConstruction.h>
#include <HydrogenEcalActionInitialization.h>

using namespace std;

int main(int argc, char** argv) { 

    // Name of the macro to process
    string macro_name = ""; 

	// Parse all the command line arguments.  If there are no valid command 
	// line arguments passed, print the usage and exit.
	static struct option long_options[] = { 
		{ "macro_name", required_argument, 0, 'm' }, 
		{ 0, 0, 0, 0 } 
	};

	int option_index = 0; 
	int option_char; 
	while ((option_char = getopt_long(argc, argv, "m:", long_options, &option_index)) != -1) { 
		switch (option_char) { 
			case 'm': 
				macro_name = optarg;
				break;
			default: 
				return EXIT_FAILURE; 	
		}
	}

    // Check if a macro was specified.  If not, warn the user and exit the 
    // application.
    if (macro_name.empty()) { 
        cerr << "\n Please specify a macro to process." << endl;
        return EXIT_FAILURE;
    }

    // Choose the Random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);

    // Construct the default run manager
    G4RunManager* run_manager = new G4RunManager();

    run_manager->SetUserInitialization(new HydrogenEcalDetectorConstruction()); 

    // Define which physics list will be used
    G4VModularPhysicsList* physics_list = new QGSP_BERT;
    physics_list->SetVerboseLevel(1);
    run_manager->SetUserInitialization(physics_list); 

    run_manager->SetUserInitialization(new HydrogenEcalActionInitialization());

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Process the specified macro
    UImanager->ApplyCommand("/control/execute " + macro_name);
    
    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted 
    // in the main() program !
    delete run_manager;
}
