#
# CMakeLists.txt
# Omar Moreno <omoreno@slac.stanford.edu>
# SLAC National Accelerator Facility
# June 2016 
#


# Find the required ROOT libraries
find_library(ROOT_CORE_LIB Core $ENV{ROOTSYS}/lib)
find_library(ROOT_RIO_LIB RIO $ENV{ROOTSYS}/lib)
find_library(ROOT_GPAD_LIB Gpad $ENV{ROOTSYS}/lib)
find_library(ROOT_TREE_LIB Tree $ENV{ROOTSYS}/lib)
find_library(ROOT_PHYSICS_LIB Physics $ENV{ROOTSYS}/lib)
find_library(ROOT_HIST_LIB Hist $ENV{ROOTSYS}/lib)

# Find the necessary LCIO libraries
find_library(LCIO_LIB lcio $ENV{LCIO}/lib)

add_library(Plotter Plotter.cxx) 
target_link_libraries(Plotter
    ${ROOT_CORE_LIB}
    ${ROOT_RIO_LIB}
    ${ROOT_HIST_LIB}
    ${ROOT_GPAD_LIB}
)

add_library(FlatTupleMaker FlatTupleMaker.cxx)
target_link_libraries(FlatTupleMaker
    ${ROOT_CORE_LIBRARY}
    ${ROOT_RIO_LIB}
    ${ROOT_TREE_LIB}
)

add_library(TrackUtils TrackUtils.cxx)
target_link_libraries(TrackUtils
    ${LCIO_LIB}
)

add_executable(lcio_analyzer lcio_analyzer.cxx)
target_link_libraries(lcio_analyzer 
    ${LCIO_LIB}
    LcioAbstractAnalysis
    HitAnalysis
    TrackAnalysis
    TaggerTrackerAnalysis
)

