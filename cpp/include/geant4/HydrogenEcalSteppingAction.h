
#ifndef __HYDROGEN_ECAL_STEPPING_ACTION_H__
#define __HYDROGEN_ECAL_STEPPING_ACTION_H__

#include <GeantAnalysis.h>

//------------//
//   Geant4   //
//------------//
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4UserSteppingAction.hh"
#include "G4VProcess.hh"

class HydrogenEcalSteppingAction : public G4UserSteppingAction { 

    public: 

        /** Default constructor */
        HydrogenEcalSteppingAction();

        /** Destructor */
        ~HydrogenEcalSteppingAction();

        /**
         */
        void UserSteppingAction(const G4Step* step);
    
    private: 

};

#endif // __HYDROGEN_ECAL_STEPPING_ACTION_H__
