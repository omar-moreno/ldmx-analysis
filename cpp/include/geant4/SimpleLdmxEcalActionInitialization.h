
#ifndef __SIMPLE_LDMX_ECAL_ACTION_INITIALIZATION_H__
#define __SIMPLE_LDMX_ECAL_ACTION_INITIALIZATION_H__

#include <SimpleLdmxEcalPrimaryGeneratorAction.h>

//------------//
//   Geant4   //
//------------//
#include "G4VUserActionInitialization.hh"


class SimpleLdmxEcalActionInitialization : public G4VUserActionInitialization { 

    public: 

        /** Default constructor */
        SimpleLdmxEcalActionInitialization();

        /** Destructor */
        ~SimpleLdmxEcalActionInitialization();

        /**
         */
        void BuildForMaster() const; 

        /**
         */
        void Build() const; 
};


#endif // __SIMPLE_LDMX_ECAL_ACTION_INITIALIZATION_H__
