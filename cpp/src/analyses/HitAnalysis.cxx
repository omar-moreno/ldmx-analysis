
#include <HitAnalysis.h>

HitAnalysis::HitAnalysis()
    : plotter(new Plotter()), 
      findable_tracks(0), 
      found_tracks(0){  
    LcioAbstractAnalysis::class_name = "HitAnalysis";
}

HitAnalysis::~HitAnalysis() { 
}

void HitAnalysis::initialize() { 
    
    TH1* plot = nullptr; 

    // Truth plots
    plot = plotter->build1DHistogram("findable particle p", 100, 3, 5);
    plot = plotter->build1DHistogram("found particle p", 100, 3, 5);

    // Tagger plots 
    plot = plotter->build2DHistogram("Simulated Tagger Hits per Layer", 14, 1, 15, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");
    
    for (int layer_n = 1; layer_n <= 14; ++layer_n) { 
        plot = plotter->build2DHistogram("Simulated tagger hit position - sensor " + std::to_string(layer_n), 100, -100, 100, 20, -20, 20);
    
    }
    //plot->GetYa


    plot = plotter->build2DHistogram("Simulated Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build2DHistogram("Simulated Recoil Hit Time per Layer", 10, 1, 11, 50, 0, 50); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Hit times");

    plot = plotter->build2DHistogram("Tagger Shaper Amplitude", 6, 0, 6, 1000, 2000, 5000);

    plot = plotter->build2DHistogram("Raw Tagger Hits per Layer", 14, 1, 15, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build1DHistogram("Total Tagger tracks", 5, 0, 5);
    plot = plotter->build1DHistogram("Tagger - Hits per track", 8, 0, 8);
    plot = plotter->build1DHistogram("Tagger track p", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track p - 6 hit", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track p - 7 hit", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track p - truth momentum", 100, -.5, .5); 
    plotter->build1DHistogram("Tagger track charge", 3, -1, 2)->GetXaxis()->SetTitle("Track charge");
    plotter->build1DHistogram("Tagger track chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("Tagger doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("Tagger z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("Tagger sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("Tagger curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("Tagger tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    plotter->build1DHistogram("Tagger cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");
    
    plot = plotter->build2DHistogram("Tagger track p vs chi2", 100, 0, 5., 50, 0, 50);

    plot = plotter->build2DHistogram("Raw Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build1DHistogram("Total Recoil tracks", 5, 0, 5);
    plot = plotter->build1DHistogram("Recoil - Hits per track", 8, 0, 8);
    plot = plotter->build1DHistogram("Recoil track p", 100, 0, 5.);
    plot = plotter->build1DHistogram("Recoil track p - truth momentum", 100, -1, 1); 
    plotter->build1DHistogram("Recoil track charge", 3, -1, 2)->GetXaxis()->SetTitle("Track charge");
    plotter->build1DHistogram("Recoil track chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("Recoil doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("Recoil z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("Recoil sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("Recoil curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("Recoil tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    plotter->build1DHistogram("Recoil cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plot = plotter->build2DHistogram("Recoil track p vs chi2", 100, 0, 5., 50, 0, 50);
    plot = plotter->build2DHistogram("Recoil track p - truth momentum vs p truth", 100, -1, 1, 100, 0, 5); 
}

void HitAnalysis::processEvent(EVENT::LCEvent* event) { 

    std::cout << "Event: " << event->getEventNumber() << std::endl;

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
       = (EVENT::LCCollection*) event->getCollection("MCParticle");
   
    // Create a pointer to the incident electron i.e. the beam electron. 
    EVENT::MCParticle* incident_e = nullptr;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simple done by looking for an electron that has no parent. 
    int n_electrons = 0;
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
    
        if (particle->getPDG() == 11) n_electrons++; 

        if (particle->getParents().size() == 0) {
           incident_e = particle; 
        } 

    } 

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* tagger_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> tagger_sim_hit_decoder(tagger_sim_hits);
    
    std::vector<int> hits(14, 0);
    std::vector<int> incident_e_hits(14, 0); 
    // Loop over all Tagger SimTrackerHits in the event.
    for (int tagger_sim_hit_n = 0; 
            tagger_sim_hit_n < tagger_sim_hits->getNumberOfElements(); ++tagger_sim_hit_n) { 
    
        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* tagger_sim_hit 
            = (EVENT::SimTrackerHit*) tagger_sim_hits->getElementAt(tagger_sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = tagger_sim_hit_decoder(tagger_sim_hit)["layer"];
        hits[layer - 1]++;

        if (incident_e == tagger_sim_hit->getMCParticle()) { 
            incident_e_hits[layer - 1]++;
        }
    }
    
    bool findable_tagger_track = true;  
    for (int layer_n = 0; layer_n < 14; layer_n++) { 
        plotter->get2DHistogram("Simulated Tagger Hits per Layer")->Fill(layer_n+1, hits[layer_n]);
    
        if (incident_e_hits[layer_n] == 0 && layer_n < 10) {
            findable_tagger_track = false; 
        } 
    }
    
    if ((incident_e_hits[10]*incident_e_hits[11] == 0) && incident_e_hits[12]*incident_e_hits[13] == 0)
        findable_tagger_track = false;

    double* mc_p_vec = (double*) incident_e->getMomentum();
    double mc_p = sqrt(mc_p_vec[0]*mc_p_vec[0] + mc_p_vec[1]*mc_p_vec[1] + mc_p_vec[2]*mc_p_vec[2]);
    if (findable_tagger_track) {
        findable_tracks++;
        plotter->get1DHistogram("findable particle p")->Fill(mc_p);
    }
    //  Tagger Tracks  //
    //-----------------//

    // Get the collection of Tagger tracker tracks from the event.  If no such
    // collection exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* tagger_tracks 
        = (EVENT::LCCollection*) event->getCollection("TaggerTracks");

    plotter->get1DHistogram("Total Tagger tracks")->Fill(tagger_tracks->getNumberOfElements());
    
    if (tagger_tracks->getNumberOfElements() > 0) { 
        found_tracks++; 
        plotter->get1DHistogram("found particle p")->Fill(mc_p); 
    }

    // Loop over all of the Tagger tracks in the event
    for (int tagger_track_n = 0; 
            tagger_track_n < tagger_tracks->getNumberOfElements(); ++tagger_track_n) {
        
        EVENT::Track* tagger_track 
            = (EVENT::Track*) tagger_tracks->getElementAt(tagger_track_n); 

        std::vector<TrackerHit*> track_hits = tagger_track->getTrackerHits();
        plotter->get1DHistogram("Tagger - Hits per track")->Fill(track_hits.size());
         
        double p = TrackUtils::getMomentum(tagger_track, -1.5); 
        
        plotter->get1DHistogram("Tagger track p")->Fill(p);
        if (track_hits.size() == 6) plotter->get1DHistogram("Tagger track p - 6 hit")->Fill(p);  
        else if (track_hits.size() == 7) plotter->get1DHistogram("Tagger track p - 7 hit")->Fill(p);

        plotter->get1DHistogram("Tagger curvature")->Fill(tagger_track->getOmega()); 
        plotter->get1DHistogram("Tagger doca")->Fill(tagger_track->getD0());
        plotter->get1DHistogram("Tagger z0")->Fill(tagger_track->getZ0());
        plotter->get1DHistogram("Tagger sin(phi0)")->Fill(sin(tagger_track->getPhi()));
        plotter->get1DHistogram("Tagger tan_lambda")->Fill(tagger_track->getTanLambda());
        //plotter->get1DHistogram("Tagger cos(theta)")->Fill(TrackExtrapolator::getCosTheta(tagger_track));

        plotter->get2DHistogram("Tagger track p vs chi2")->Fill(p, tagger_track->getChi2());

        if (tagger_tracks->getNumberOfElements() == 1 && n_electrons == 1) { 

            std::vector<MCParticle*> daughters = incident_e->getDaughters();
            for (int index = 0; index < daughters.size(); index++) { 
             if (daughters[index]->getVertex()[2] < -4.0) { 
                    double* mc_d_p = (double*) daughters[index]->getMomentum();
                    mc_p_vec[0] = mc_p_vec[0] - mc_d_p[0];
                    mc_p_vec[1] -= mc_d_p[1];
                    mc_p_vec[2] -= mc_d_p[2];
                }
            }
            mc_p = sqrt(mc_p_vec[0]*mc_p_vec[0] + mc_p_vec[1]*mc_p_vec[1] + mc_p_vec[2]*mc_p_vec[2]);
            plotter->get1DHistogram("Tagger track p - truth momentum")->Fill(p - mc_p);
        }
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

    // Get the collection of Tagger tracker tracks from the event.  If no such
    // collection exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* recoil_tracks 
        = (EVENT::LCCollection*) event->getCollection("RecoilTracks");

    plotter->get1DHistogram("Total Recoil tracks")->Fill(recoil_tracks->getNumberOfElements());

    // Loop over all of the Recoil tracks in the event
    for (int recoil_track_n = 0; 
            recoil_track_n < recoil_tracks->getNumberOfElements(); ++recoil_track_n) {
        
        EVENT::Track* recoil_track 
            = (EVENT::Track*) recoil_tracks->getElementAt(recoil_track_n); 

        std::vector<TrackerHit*> track_hits = recoil_track->getTrackerHits();
        plotter->get1DHistogram("Recoil - Hits per track")->Fill(track_hits.size());
         
        double p = TrackUtils::getMomentum(recoil_track, -1.5);
        std::cout << "Recoil p: " << p << std::endl; 
        plotter->get1DHistogram("Recoil track p")->Fill(p);
        plotter->get2DHistogram("Recoil track p vs chi2")->Fill(p, recoil_track->getChi2());

        plotter->get1DHistogram("Recoil curvature")->Fill(recoil_track->getOmega()); 
        plotter->get1DHistogram("Recoil doca")->Fill(recoil_track->getD0());
        plotter->get1DHistogram("Recoil z0")->Fill(recoil_track->getZ0());
        plotter->get1DHistogram("Recoil sin(phi0)")->Fill(sin(recoil_track->getPhi()));
        plotter->get1DHistogram("Recoil tan_lambda")->Fill(recoil_track->getTanLambda());
        //plotter->get1DHistogram("Recoil cos(theta)")->Fill(TrackExtrapolator::getCosTheta(recoil_track));

        if (recoil_tracks->getNumberOfElements() == 1 && n_electrons == 1) { 

        mc_p_vec = (double*) incident_e->getMomentum();
        std::vector<MCParticle*> daughters = incident_e->getDaughters();
        for (int index = 0; index < daughters.size(); index++) { 
                double* mc_d_p = (double*) daughters[index]->getMomentum();
                mc_p_vec[0] = mc_p_vec[0] - mc_d_p[0];
                mc_p_vec[1] -= mc_d_p[1];
                mc_p_vec[2] -= mc_d_p[2];
        }
        mc_p = sqrt(mc_p_vec[0]*mc_p_vec[0] + mc_p_vec[1]*mc_p_vec[1] + mc_p_vec[2]*mc_p_vec[2]);
        //double* mc_p_end = (double*) incident_e->getMomentumAtEndpoint();
        plotter->get1DHistogram("Recoil track p - truth momentum")->Fill(p - mc_p);
        plotter->get2DHistogram("Recoil track p - truth momentum vs p truth")->Fill(p - mc_p, mc_p);
    } 
    }


}

void HitAnalysis::finalize() { 
    plotter->saveToRootFile("hit_analysis.root");

    std::cout << "[ ]: Simple Tracking Efficiency: " << (found_tracks/findable_tracks)*100 << "%" << std::endl;
}
