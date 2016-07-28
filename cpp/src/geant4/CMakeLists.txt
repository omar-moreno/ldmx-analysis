#
# CMakeLists.txt
#   Omar Moreno 
#   SLAC National Accelrator Laboratory
#

add_library(SimpleLdmxEcalDetectorConstruction SimpleLdmxEcalDetectorConstruction.cxx)
target_link_libraries(SimpleLdmxEcalDetectorConstruction
    ${Geant4_Libraries}
)

add_library(HydrogenEcalDetectorConstruction HydrogenEcalDetectorConstruction.cxx)
target_link_libraries(HydrogenEcalDetectorConstruction
    ${Geant4_Libraries}
)

add_library(SimpleLdmxEcalPrimaryGeneratorAction SimpleLdmxEcalPrimaryGeneratorAction.cxx)
target_link_libraries(SimpleLdmxEcalPrimaryGeneratorAction
    ${Geant4_Libraries}
)

add_library(SimpleLdmxEcalSteppingAction SimpleLdmxEcalSteppingAction.cxx)
target_link_libraries(SimpleLdmxEcalSteppingAction
    ${Geant4_Libraries}
)

add_library(HydrogenEcalSteppingAction HydrogenEcalSteppingAction.cxx)
target_link_libraries(HydrogenEcalSteppingAction
    ${Geant4_Libraries}
)

add_library(HydrogenEcalRunAction HydrogenEcalRunAction.cxx)
target_link_libraries(HydrogenEcalRunAction
    ${Geant4_Libraries}
)


add_library(SimpleLdmxEcalActionInitialization SimpleLdmxEcalActionInitialization.cxx)
target_link_libraries(SimpleLdmxEcalActionInitialization
    ${Geant4_Libraries}
    SimpleLdmxEcalPrimaryGeneratorAction
    SimpleLdmxEcalSteppingAction
)

add_library(HydrogenEcalActionInitialization HydrogenEcalActionInitialization.cxx)
target_link_libraries(HydrogenEcalActionInitialization
    ${Geant4_Libraries}
    SimpleLdmxEcalPrimaryGeneratorAction
    HydrogenEcalRunAction
    HydrogenEcalSteppingAction
)

add_executable(simple_ldmx_ecal simple_ldmx_ecal.cxx)
target_link_libraries(simple_ldmx_ecal
    ${Geant4_LIBRARIES}
    SimpleLdmxEcalDetectorConstruction
    SimpleLdmxEcalActionInitialization
)

add_executable(hydrogen_ecal hydrogen_ecal.cxx)
target_link_libraries(hydrogen_ecal
    ${Geant4_LIBRARIES}
    HydrogenEcalDetectorConstruction
    HydrogenEcalActionInitialization
)
