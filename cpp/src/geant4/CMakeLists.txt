#
# CMakeLists.txt
#   Omar Moreno 
#   SLAC National Accelrator Laboratory
#

add_library(SimpleLdmxEcalDetectorConstruction SimpleLdmxEcalDetectorConstruction.cxx)
target_link_libraries(SimpleLdmxEcalDetectorConstruction
    ${Geant4_Libraries}
)

add_library(SimpleLdmxEcalPrimaryGeneratorAction SimpleLdmxEcalPrimaryGeneratorAction.cxx)
target_link_libraries(SimpleLdmxEcalPrimaryGeneratorAction
    ${Geant4_Libraries}
)

add_library(SimpleLdmxEcalActionInitialization SimpleLdmxEcalActionInitialization.cxx)
target_link_libraries(SimpleLdmxEcalActionInitialization
    ${Geant4_Libraries}
    SimpleLdmxEcalPrimaryGeneratorAction
)


add_executable(simple_ldmx_ecal simple_ldmx_ecal.cxx)
target_link_libraries(simple_ldmx_ecal
    ${Geant4_LIBRARIES}
    SimpleLdmxEcalDetectorConstruction
    SimpleLdmxEcalActionInitialization
)

