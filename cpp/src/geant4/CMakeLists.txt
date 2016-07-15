
find_library(G4_RUN_LIB G4run $ENV{G4_LIB_PATH}/lib64)

add_executable(simple_ldmx_ecal simple_ldmx_ecal.cxx)
target_link_libraries(simple_ldmx_ecal
    ${G4_RUN_LIB}
)

