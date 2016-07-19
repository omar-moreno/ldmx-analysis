
#ifndef __SIMPLE_LDMX_ECAL_PRIMARY_GENERATOR_ACTION_H__
#define __SIMPLE_LDMX_ECAL_PRIMARY_GENERATOR_ACTION_H__

//------------//
//   Geant4   //
//------------//
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class SimpleLdmxEcalPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction { 

    public: 

        /** Default constructor */
        SimpleLdmxEcalPrimaryGeneratorAction();

        /** Destructor */
        ~SimpleLdmxEcalPrimaryGeneratorAction();

        /**
         */
        void GeneratePrimaries(G4Event* event);

    private: 

        /** Geant4 particle gun */
        G4ParticleGun* particle_gun; 

}; // SimpleLdmxEcalPrimaryGeneratorAction

#endif // __SIMPLE_LDMX_ECAL_PRIMARY_GENERATOR_ACTION_H__

