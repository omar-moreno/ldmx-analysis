# Find the necessary LCIO libraries
find_library(LCIO_LIB lcio $ENV{LCIO}/lib)

add_library(LcioAbstractAnalysis LcioAbstractAnalysis.cxx)
target_link_libraries(LcioAbstractAnalysis 
    ${LCIO_LIB}    
)

add_library(HitAnalysis HitAnalysis.cxx)
target_link_libraries(HitAnalysis
    Plotter
    LcioAbstractAnalysis
)

add_library(TrackAnalysis TrackAnalysis.cxx)
target_link_libraries(TrackAnalysis
    Plotter
    TrackUtils
    LcioAbstractAnalysis
)

add_library(TaggerTrackerAnalysis TaggerTrackerAnalysis.cxx)
target_link_libraries(TaggerTrackerAnalysis
    FlatTupleMaker    
    LcioAbstractAnalysis
    TrackUtils
)
