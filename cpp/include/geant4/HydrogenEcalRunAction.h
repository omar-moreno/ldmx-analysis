
#ifndef __HYDROGEN_ECAL_RUN_ACTION_H__
#define __HYDROGEN_ECAL_RUN_ACTION_H__

#include <GeantAnalysis.h>

//------------//
//   Geant4   //
//------------//
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

class HydrogenEcalRunAction : public G4UserRunAction { 

    public: 

        /** Default constructor */
        HydrogenEcalRunAction(); 

        /** Destructor */
        ~HydrogenEcalRunAction(); 

        void BeginOfRunAction(const G4Run* run); 

        void EndOfRunAction(const G4Run* run);

        static void IncrementPhotoNuclearCount() { ++photo_nuclear_count; };

    private: 

        static int photo_nuclear_count; 
}; // HydrogenEcalRunAction

#endif // __HYDROGEN_ECAL_RUN_ACTION_H__
