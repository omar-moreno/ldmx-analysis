
#include <HitAnalysis.h>

HitAnalysis::HitAnalysis()
    : plotter(new Plotter()) {  
    LcioAbstractAnalysis::class_name = "HitAnalysis";
}

HitAnalysis::~HitAnalysis() { 
}

void HitAnalysis::initialize() { 
    
    TH1* plot = nullptr; 
    plot = plotter->build1DHistogram("Recoil SimTrackerHit Layer 1 - x", 100, -100, 100);
    plot->GetXaxis()->SetTitle("x position (mm)");
}

void HitAnalysis::processEvent(EVENT::LCEvent* event) { 

    // Get the collection of sim hits from the event.  If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* recoil_hits = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");

    // Loop over all of the sim hits in the event
    for (int recoil_hit_n = 0; recoil_hit_n < recoil_hits->getNumberOfElements(); ++recoil_hit_n) { 
        
        // Get a sim hit from the list of hits
        EVENT::SimTrackerHit* recoil_hit = (EVENT::SimTrackerHit*) recoil_hits->getElementAt(recoil_hit_n);     
   
        const double* hit_position = recoil_hit->getPosition();

        
        //plotter->get1DHistogram("Recoil SimTrackerHit Layer 1 - x")->Fill(hit_position[0]); 
    }
}

void HitAnalysis::finalize() { 
    plotter->saveToRootFile("hit_analysis.root");
}
