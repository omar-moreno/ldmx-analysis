
#include <ClusterAnalysis.h>

ClusterAnalysis::ClusterAnalysis() { 
    LcioAbstractAnalysis::class_name = "ClusterAnalysis";
}

ClusterAnalysis::~ClusterAnalysis() { 
}

void ClusterAnalysis::initialize() { 
}

void ClusterAnalysis::processEvent(EVENT::LCEvent* event) { 

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* fitted_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerFittedRawTrackerHits"); 
}

void ClusterAnalysis::finalize() { 
}
