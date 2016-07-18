#
# CMakeLists.txt
#   Omar Moreno 
#   SLAC National Accelrator Laboratory
#

add_library(SimpleLdmxEcalDetectorConstruction SimpleLdmxEcalDetectorConstruction.cxx)
target_link_libraries(SimpleLdmxEcalDetectorConstruction
    ${Geant4_Libraries}
)

add_executable(simple_ldmx_ecal simple_ldmx_ecal.cxx)
target_link_libraries(simple_ldmx_ecal
    ${Geant4_LIBRARIES}
    SimpleLdmxEcalDetectorConstruction
)

