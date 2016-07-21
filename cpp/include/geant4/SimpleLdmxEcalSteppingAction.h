
#ifndef __SIMPLE_LDMX_ECAL_STEPPING_ACTION_H__
#define __SIMPLE_LDMX_ECAL_STEPPING_ACTION_H__

//------------//
//   Geant4   //
//------------//
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4UserSteppingAction.hh"


class SimpleLdmxEcalSteppingAction : public G4UserSteppingAction { 

    public: 

        /** Default constructor */
        SimpleLdmxEcalSteppingAction();

        /** Destructor */
        ~SimpleLdmxEcalSteppingAction();

        /**
         */
        void UserSteppingAction(const G4Step* step);

};

#endif // __SIMPLE_LDMX_ECAL_STEPPING_ACTION_H__
