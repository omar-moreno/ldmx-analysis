
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

    plot = plotter->build1DHistogram("Total Tagger tracks", 5, 0, 5);
    plot = plotter->build1DHistogram("Total track hits", 8, 0, 8);
    plot = plotter->build1DHistogram("Tagger track momentum", 100, 3., 5.);
    plot = plotter->build2DHistogram("Tagger track momentum vs chi2", 100, 0, 5., 50, 0, 50);
    plot = plotter->build1DHistogram("Tagger track momentum - 6 hit", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track momentum - 7 hit", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track momentum - truth momentum", 100, -.5, .5); 
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

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
       = (EVENT::LCCollection*) event->getCollection("MCParticle");
    
    EVENT::MCParticle* incident_e = nullptr;
    // Loop over all of the MC particles 
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {
       
        incident_e = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
        if (incident_e->getParents().size() == 0) break; 
    } 

    // Get the collection of Tagger tracker tracks from the event.  If no such
    // collection exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* tagger_tracks 
        = (EVENT::LCCollection*) event->getCollection("Tracks");

    plotter->get1DHistogram("Total Tagger tracks")->Fill(tagger_tracks->getNumberOfElements());

    // Loop over all of the Tagger tracks in the event
    for (int tagger_track_n = 0; 
            tagger_track_n < tagger_tracks->getNumberOfElements(); ++tagger_track_n) {
        
        EVENT::Track* tagger_track 
            = (EVENT::Track*) tagger_tracks->getElementAt(tagger_track_n); 

        std::vector<TrackerHit*> track_hits = tagger_track->getTrackerHits();
        plotter->get1DHistogram("Total track hits")->Fill(track_hits.size());
         
        double p = TrackUtils::getMomentum(tagger_track, 1.5); 
        plotter->get1DHistogram("Tagger track momentum")->Fill(p);
        plotter->get2DHistogram("Tagger track momentum vs chi2")->Fill(p, tagger_track->getChi2());

        const double* mc_p_vec = incident_e->getMomentum();
        double mc_p = sqrt(mc_p_vec[0]*mc_p_vec[0] + mc_p_vec[1]*mc_p_vec[1] + mc_p_vec[2]*mc_p_vec[2]);
        plotter->get1DHistogram("Tagger track momentum - truth momentum")->Fill(p - mc_p);
    }

}

void HitAnalysis::finalize() { 
    plotter->saveToRootFile("hit_analysis.root");
}
