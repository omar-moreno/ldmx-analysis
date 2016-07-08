
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

    //------------------//
    //   Tagger plots   //
    //------------------//
    
    plot = plotter->build2DHistogram("Tagger SimTrackerHit per layer", 14, 1, 15, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build2DHistogram("Tagger SimTrackerHit PDG ID per layer", 14, 1, 15, 60, -30, 30); 
    plot->GetXaxis()->SetTitle("Layer"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build2DHistogram("Tagger SimTrackerHit dE/dx per layer", 14, 1, 15, 100, 0, .001); 
    plot->GetXaxis()->SetTitle("Layer"); 
    plot->GetYaxis()->SetTitle("dE/dx");

    plot = plotter->build2DHistogram("Tagger SimTrackerHit time per layer", 14, 1, 15, 400, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer"); 
    plot->GetYaxis()->SetTitle("time (ns)");

    for (int layer_n = 1; layer_n <= 14; ++layer_n) { 
        plot = plotter->build2DHistogram("Tagger SimTrackerHit position - Layer " + std::to_string(layer_n),
                600, -30, 30, 200, -2, 8);
        plot = plotter->build2DHistogram(
                "Tagger SimTrackerHit position - secondaries - Layer - " + std::to_string(layer_n),
                600, -30, 30, 200, -10, 10); 

        plot = plotter->build2DHistogram("Tagger Shaper Amplitude - Layer " + std::to_string(layer_n),
                6, 0, 6, 1000, 2000, 5000);
    }

    plot = plotter->build2DHistogram("Raw Tagger Hits per layer", 14, 1, 15, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build2DHistogram("Raw Tagger Hits time per layer", 14, 1, 15, 100, 0, 20); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("time (ns)");

    plot = plotter->build1DHistogram("Tagger Cluster size", 6, 0, 6);

    plot = plotter->build2DHistogram("Tagger Clusters per Layer", 14, 1, 15, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer"); 
    plot->GetYaxis()->SetTitle("Total clusters");

    plot = plotter->build2DHistogram("Tagger Cluster size per Layer", 14, 1, 15, 6, 0, 6); 
    plot->GetXaxis()->SetTitle("Layer"); 
    plot->GetYaxis()->SetTitle("Cluster size");

    //------------------//
    //   Recoil plots   //
    //------------------//
    
    plot = plotter->build1DHistogram("Recoil Cluster size", 6, 0, 6);

    plot = plotter->build2DHistogram("Simulated Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build2DHistogram("Simulated Recoil Hit Time per Layer", 10, 1, 11, 50, 0, 50); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Hit times");

    plot = plotter->build2DHistogram("Raw Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");
}

void HitAnalysis::processEvent(EVENT::LCEvent* event) { 


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
    double sim_px = 0; 
    double sim_py = 0;
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
            if (layer == 14) { 
                sim_px  = tagger_sim_hit->getMomentum()[0]; 
                sim_py  = tagger_sim_hit->getMomentum()[1];
            }
        } else { 
            plotter->get2DHistogram(
                "Tagger SimTrackerHit position - secondaries - Layer - " + std::to_string(layer))->Fill(
                tagger_sim_hit->getPosition()[0], tagger_sim_hit->getPosition()[1]);
        
        }

        plotter->get2DHistogram("Tagger SimTrackerHit position - Layer " + std::to_string(layer))->Fill(
                tagger_sim_hit->getPosition()[0], tagger_sim_hit->getPosition()[1]);
        plotter->get2DHistogram("Tagger SimTrackerHit PDG ID per layer")->Fill(layer, tagger_sim_hit->getMCParticle()->getPDG()); 
        plotter->get2DHistogram("Tagger SimTrackerHit dE/dx per layer")->Fill(layer, tagger_sim_hit->getEDep());
        plotter->get2DHistogram("Tagger SimTrackerHit time per layer")->Fill(layer, tagger_sim_hit->getTime());
    }


    // Get the collection of raw hits associated with the Tagger tracker from 
    // the event.  If no such collection exist, a DataNotAvailableException is
    // thrown.
    EVENT::LCCollection* tagger_raw_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerRawTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger raw hits.
    std::string encoder_string = "system:6,barrel:3,layer:4,module:12,sensor:1,side:32:-2,strip:1";
    UTIL::CellIDDecoder<EVENT::TrackerRawData> tagger_raw_hit_decoder(encoder_string);

    // Loop over all raw hits in the event
    std::vector<int> tagger_raw_hit_count(14, 0);
    for (int tagger_raw_hit_n = 0; 
            tagger_raw_hit_n < tagger_raw_hits->getNumberOfElements(); ++tagger_raw_hit_n) { 
        
        // Get a raw hit from the list of raw hits
        EVENT::TrackerRawData* tagger_raw_hit 
            = (EVENT::TrackerRawData*) tagger_raw_hits->getElementAt(tagger_raw_hit_n);     

        // Get the layer number associated with this hit.
        int layer = tagger_raw_hit_decoder(tagger_raw_hit)["layer"];
        tagger_raw_hit_count[layer - 1]++;

        // Get the ADC values associated with this hit
        std::vector<short> adc_values = tagger_raw_hit->getADCValues();
        
        for (int sample_n = 0; sample_n < 6; ++sample_n) { 
            plotter->get2DHistogram("Tagger Shaper Amplitude - Layer " + std::to_string(layer))->Fill(
                    sample_n, adc_values[sample_n]);
        } 
        
        plotter->get2DHistogram("Raw Tagger Hits time per layer")->Fill(layer, tagger_raw_hit->getTime());
    }

    for (int layer_n = 0; layer_n < 14; layer_n++) { 
        plotter->get2DHistogram("Raw Tagger Hits per layer")->Fill(layer_n+1, tagger_raw_hit_count[layer_n]);
    }


    EVENT::LCCollection* clusters 
        = (EVENT::LCCollection*) event->getCollection("TaggerClusters");

    UTIL::CellIDDecoder<EVENT::TrackerHit> tagger_cluster_decoder(encoder_string);

    std::vector<int> tagger_cluster_count(14, 0);
    for (int cluster_n = 0; cluster_n < clusters->getNumberOfElements(); ++cluster_n) { 
        
        EVENT::TrackerHit* cluster = (EVENT::TrackerHit*) clusters->getElementAt(cluster_n);
       
        // Get a raw hit from the list of raw hits
        EVENT::TrackerRawData* tagger_raw_hit 
            = (EVENT::TrackerRawData*) cluster->getRawHits()[0]; 

        // Get the layer number associated with this hit.
        int layer = tagger_raw_hit_decoder(tagger_raw_hit)["layer"];
        //std::cout << "Cluster layer: " << layer << std::endl;
        tagger_cluster_count[layer - 1]++;

        plotter->get1DHistogram("Tagger Cluster size")->Fill(cluster->getRawHits().size());
        plotter->get2DHistogram("Tagger Cluster size per Layer")->Fill(layer, cluster->getRawHits().size());
    }
    
    for (int layer_n = 0; layer_n < 14; layer_n++) { 
        plotter->get2DHistogram("Tagger Clusters per Layer")->Fill(layer_n+1, tagger_cluster_count[layer_n]);
    }

    clusters 
        = (EVENT::LCCollection*) event->getCollection("RecoilClusters");

    for (int cluster_n = 0; cluster_n < clusters->getNumberOfElements(); ++cluster_n) { 
        
        EVENT::TrackerHit* cluster = (EVENT::TrackerHit*) clusters->getElementAt(cluster_n);
         
        plotter->get1DHistogram("Recoil Cluster size")->Fill(cluster->getRawHits().size());

         
    }
 
    bool findable_tagger_track = true;  
    for (int layer_n = 0; layer_n < 14; layer_n++) { 
        plotter->get2DHistogram("Tagger SimTrackerHit per layer")->Fill(layer_n+1, hits[layer_n]);
    
        if (incident_e_hits[layer_n] == 0 && layer_n < 10) {
            findable_tagger_track = false; 
        } 
    }
    
    if ((incident_e_hits[10]*incident_e_hits[11] == 0) && incident_e_hits[12]*incident_e_hits[13] == 0)
        findable_tagger_track = false;

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* recoil_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");

    // Create a cell ID decoder
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> recoil_sim_hit_decoder(recoil_sim_hits);
    
    // Loop over all SimTrackerHits in the event
    std::vector<int> recoil_hits(10, 0);
    double sim_recoil_px = 0; 
    double sim_recoil_py = 0;
    double sim_recoil_pz = 0;
    for (int recoil_sim_hit_n = 0; 
            recoil_sim_hit_n < recoil_sim_hits->getNumberOfElements(); ++recoil_sim_hit_n) { 
    
        // Get a SimTrackerHit from the collection of hits
        EVENT::SimTrackerHit* recoil_sim_hit 
            = (EVENT::SimTrackerHit*) recoil_sim_hits->getElementAt(recoil_sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = recoil_sim_hit_decoder(recoil_sim_hit)["layer"];
        recoil_hits[layer - 1]++;
        if (incident_e == recoil_sim_hit->getMCParticle()) { 
            incident_e_hits[layer - 1]++;
            if (layer == 1) { 
                sim_recoil_px  = recoil_sim_hit->getMomentum()[0]; 
                sim_recoil_py  = recoil_sim_hit->getMomentum()[1];
                sim_recoil_pz  = recoil_sim_hit->getMomentum()[2];
            }
        }
    }
    
    for (int layer_n = 0; layer_n < 10; layer_n++) { 
        plotter->get2DHistogram("Simulated Recoil Hits per Layer")->Fill(layer_n+1, recoil_hits[layer_n]);
    }


}

void HitAnalysis::finalize() { 
    plotter->saveToRootFile("hit_analysis.root");
}
