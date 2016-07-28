
#ifndef __HYDROGEN_ECAL_ACTION_INITIALIZATION_H__
#define __HYDROGEN_ECAL_ACTION_INITIALIZATION_H__

#include <SimpleLdmxEcalPrimaryGeneratorAction.h>
#include <HydrogenEcalSteppingAction.h>
#include <HydrogenEcalRunAction.h>

//------------//
//   Geant4   //
//------------//
#include "G4VUserActionInitialization.hh"


class HydrogenEcalActionInitialization : public G4VUserActionInitialization { 

    public: 

        /** Default constructor */
        HydrogenEcalActionInitialization();

        /** Destructor */
        ~HydrogenEcalActionInitialization();

        /**
         */
        void BuildForMaster() const; 

        /**
         */
        void Build() const; 
};


#endif // __HYDROGEN_ECAL_ACTION_INITIALIZATION_H__
