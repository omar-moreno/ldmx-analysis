

#ifndef __SIMPLE_LDMX_ECAL_DETECTOR_CONSTRUCTION_H__
#define __SIMPLE_LDMX_ECAL_DETECTOR_CONSTRUCTION_H__

//------------//
//   Geant4   //
//------------//
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

// Forward declarations
class G4VPhysicalVolume;

class SimpleLdmxEcalDetectorConstruction : public G4VUserDetectorConstruction { 

    public: 

        /** Default constructor */
        SimpleLdmxEcalDetectorConstruction(); 
        
        /** Destructor */
        ~SimpleLdmxEcalDetectorConstruction();

        /**
         */
        G4VPhysicalVolume* Construct();

        /**
         */
        void ConstructSDandField(); 

    private: 

        void DefineMaterials(); 

        G4VPhysicalVolume* DefineVolumes(); 

}; // SimpleLdmxEcalDetectorConstruction

#endif // __SIMPLE_LDMX_ECAL_DETECTOR_CONSTRUCTION_H__
