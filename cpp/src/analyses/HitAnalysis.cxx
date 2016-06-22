
#include <HitAnalysis.h>

HitAnalysis::HitAnalysis()
    : plotter(new Plotter()) {  
    LcioAbstractAnalysis::class_name = "HitAnalysis";
}

HitAnalysis::~HitAnalysis() { 
}

void HitAnalysis::initialize() { 
    
    TH1* plot = nullptr; 

    // SimTrackerHits plots
    plot = plotter->build2DHistogram("Simulated Tagger Hits per Layer", 14, 1, 15, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build2DHistogram("Simulated Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build2DHistogram("Tagger Shaper Amplitude", 6, 0, 6, 1000, 2000, 5000);
    //plot->GetXaxis()->SetTitle("x position (mm)");
}

void HitAnalysis::processEvent(EVENT::LCEvent* event) { 

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* tagger_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerTrackerHits");

    // Create a cell ID decoder
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> tagger_sim_hit_decoder(tagger_sim_hits);
    
    // Loop over all SimTrackerHits in the event
    std::vector<int> hits(14, 0);
    for (int tagger_sim_hit_n = 0; 
            tagger_sim_hit_n < tagger_sim_hits->getNumberOfElements(); ++tagger_sim_hit_n) { 
    
        // Get a SimTrackerHit from the collection of hits
        EVENT::SimTrackerHit* tagger_sim_hit 
            = (EVENT::SimTrackerHit*) tagger_sim_hits->getElementAt(tagger_sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = tagger_sim_hit_decoder(tagger_sim_hit)["layer"];
        hits[layer - 1]++;
    }
    
    for (int layer_n = 0; layer_n < 14; layer_n++) { 
        plotter->get2DHistogram("Simulated Tagger Hits per Layer")->Fill(layer_n+1, hits[layer_n]);
    }

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* recoil_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");

    // Create a cell ID decoder
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> recoil_sim_hit_decoder(recoil_sim_hits);
    
    // Loop over all SimTrackerHits in the event
    hits.clear();
    for (int recoil_sim_hit_n = 0; 
            recoil_sim_hit_n < recoil_sim_hits->getNumberOfElements(); ++recoil_sim_hit_n) { 
    
        // Get a SimTrackerHit from the collection of hits
        EVENT::SimTrackerHit* recoil_sim_hit 
            = (EVENT::SimTrackerHit*) recoil_sim_hits->getElementAt(recoil_sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = recoil_sim_hit_decoder(recoil_sim_hit)["layer"];
        hits[layer - 1]++;
    }
    
    for (int layer_n = 0; layer_n < 10; layer_n++) { 
        plotter->get2DHistogram("Simulated Recoil Hits per Layer")->Fill(layer_n+1, hits[layer_n]);
    }


    // Get the collection of raw hits associated with the Tagger tracker from 
    // the event.  If no such collection exist, a DataNotAvailableException is
    // thrown.
    EVENT::LCCollection* tagger_raw_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerRawTrackerHits");

    // Loop over all raw hits in the event
    for (int tagger_raw_hit_n = 0; 
            tagger_raw_hit_n < tagger_raw_hits->getNumberOfElements(); ++tagger_raw_hit_n) { 
        
        // Get a sim hit from the list of hits
        EVENT::TrackerRawData* tagger_raw_hit 
            = (EVENT::TrackerRawData*) tagger_raw_hits->getElementAt(tagger_raw_hit_n);     

        // Get the ADC values associated with this hit
        std::vector<short> adc_values = tagger_raw_hit->getADCValues();
        
        for (int sample_n = 0; sample_n < 6; ++sample_n) { 
            plotter->get2DHistogram("Tagger Shaper Amplitude")->Fill(sample_n, adc_values[sample_n]);
        } 
    }
}

void HitAnalysis::finalize() { 
    plotter->saveToRootFile("hit_analysis.root");
}
