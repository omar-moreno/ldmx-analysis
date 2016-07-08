

#include <TrackAnalysis.h>

TrackAnalysis::TrackAnalysis()
    : plotter(new Plotter()), 
      findable_tracks(0), 
      found_tracks(0){  
    LcioAbstractAnalysis::class_name = "TrackAnalysis";
}

TrackAnalysis::~TrackAnalysis() { 
}

void TrackAnalysis::initialize() { 
    
    TH1* plot = nullptr;

    plot = plotter->build1DHistogram("Beam electron momentum", 600, 0, 6); 

    // Truth plots
    plot = plotter->build1DHistogram("findable particle p", 100, 3, 5);
    plot = plotter->build1DHistogram("found particle p", 100, 3, 5);
    
    plot = plotter->build1DHistogram("Total Tagger tracks", 5, 0, 5);
    plot = plotter->build1DHistogram("Tagger - Hits per track", 8, 0, 8);
    plot = plotter->build1DHistogram("Tagger track p", 150, 0, 5.);
    plot = plotter->build1DHistogram("Tagger track p_{t}", 150, 0, .03);
    plot = plotter->build1DHistogram("Tagger track p - 6 hit", 150, 0., 5.);
    plot = plotter->build1DHistogram("Tagger track p - 7 hit", 150, 0., 5.);
    plot = plotter->build1DHistogram("Tagger track p - truth p", 100, -.5, .5); 
    plot = plotter->build1DHistogram("Tagger track p_{t} - truth p_{t}", 100, -.02, .02); 
    plotter->build1DHistogram("Tagger track charge", 3, -1, 2)->GetXaxis()->SetTitle("Track charge");
    plotter->build1DHistogram("Tagger track chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("Tagger doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("Tagger z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("Tagger sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("Tagger curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("Tagger tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    plotter->build1DHistogram("Tagger cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");
   
    plot = plotter->build1DHistogram("Tagger track p - d0 cut", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track p_x - d0 cut", 100, -.01, .02);
    plot = plotter->build1DHistogram("Tagger track p_y - d0 cut", 100, -.01, .02);
    plot = plotter->build1DHistogram("Tagger track p_{t} - d0 cut", 100, 0, .03);
    plot = plotter->build1DHistogram("Tagger track p_{t} - truth p_{t} - d0 cut", 100, -.02, .02); 
    plot = plotter->build1DHistogram("Tagger track p - 6 hit - d0 cut", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track p - 7 hit - d0 cut", 100, 3., 5.);
    plot = plotter->build1DHistogram("Tagger track p - truth p - d0 cut", 100, -.5, .5); 
    plot = plotter->build1DHistogram("Tagger track p_x - truth p_x - d0 cut", 100, -.1, .1); 
    plot = plotter->build1DHistogram("Tagger track p_y - truth p_y - d0 cut", 100, -.1, .1); 

    plot = plotter->build2DHistogram("Tagger track p vs chi2", 100, 0, 5., 50, 0, 50);
    plot = plotter->build2DHistogram("Tagger track momentum vs chi2", 100, 0, 5., 50, 0, 50);
    plot = plotter->build2DHistogram("Tagger track momentum vs z0", 500, 0, 5., 100, 0, 10);
    plot = plotter->build2DHistogram("Tagger track momentum vs d0", 500, 0, 5., 100, 0, 10);
    
    plot = plotter->build2DHistogram("Tagger track momentum vs chi2 - d0 cut", 100, 0, 5., 50, 0, 50);
    plot = plotter->build2DHistogram("Tagger track momentum vs z0 - d0 cut", 500, 0, 5., 100, 0, 10);
    plot = plotter->build2DHistogram("Tagger track momentum vs d0 - d0 cut", 500, 0, 5., 100, 0, 10);
    plot = plotter->build2DHistogram("Tagger track p_{t} - truth p_{t} vs p_{t} - d0 cut", 100, 0, 0.03, 100, -.02, .02); 

    plot = plotter->build2DHistogram("Raw Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    plot = plotter->build1DHistogram("Total Recoil tracks", 5, 0, 5);
    plot = plotter->build1DHistogram("Recoil - Hits per track", 8, 0, 8);
    plot = plotter->build1DHistogram("Recoil track p", 100, 0, 6.);
    plot = plotter->build1DHistogram("Recoil track p - truth p", 100, -1, 1); 
    plotter->build1DHistogram("Recoil track charge", 3, -1, 2)->GetXaxis()->SetTitle("Track charge");
    plotter->build1DHistogram("Recoil track chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("Recoil doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("Recoil z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("Recoil sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("Recoil curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("Recoil tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    plotter->build1DHistogram("Recoil cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");
    plot = plotter->build1DHistogram("Recoil track p_{t}", 100, 0, .03);
    plot = plotter->build1DHistogram("Recoil track p_{t} - truth p_{t}", 100, -.02, .02); 

    plot = plotter->build2DHistogram("Recoil track p vs chi2", 100, 0, 5., 50, 0, 50);
    plot = plotter->build2DHistogram("Recoil track p - truth p vs p truth", 100, -1, 1, 100, 0, 5); 
    plot = plotter->build2DHistogram("Recoil track p_{t} - truth p_{t} vs p_{t}", 100, 0, 0.03, 100, -.02, .02);

    plot = plotter->build1DHistogram("Recoil d0 - Tagger d0", 100, -0.1, 0.1);  
    plot = plotter->build1DHistogram("Recoil z0 - Tagger z0", 100, -0.1, 0.1);  
    plot = plotter->build1DHistogram("Recoil pt - Tagger pt", 100, -0.1, 0.1);  
}

void TrackAnalysis::processEvent(EVENT::LCEvent* event) { 

    std::cout << "Event: " << event->getEventNumber() << std::endl;

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
       = (EVENT::LCCollection*) event->getCollection("MCParticle");
   
    // Create a pointer to the incident electron i.e. the beam electron. 
    EVENT::MCParticle* incident_e = nullptr;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    //int n_electrons = 0;
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
    
        //if (particle->getPDG() == 11) n_electrons++; 

        if (particle->getParents().size() == 0) {
           incident_e = particle; 
        } 
    }

    // Calculate the momentum of the incident electron
    double* beam_e_pvec = (double*) incident_e->getMomentum();
    double beam_e_p = sqrt(pow(beam_e_pvec[0], 2) + pow(beam_e_pvec[1], 2) + pow(beam_e_pvec[2], 2));
    plotter->get1DHistogram("Beam electron momentum")->Fill(beam_e_p);

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* tagger_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> tagger_sim_hit_decoder(tagger_sim_hits);
    
    std::vector<int> incident_e_hits(14, 0); 
    // Loop over all Tagger SimTrackerHits in the event.
    for (int tagger_sim_hit_n = 0; 
            tagger_sim_hit_n < tagger_sim_hits->getNumberOfElements(); ++tagger_sim_hit_n) { 
    
        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* tagger_sim_hit 
            = (EVENT::SimTrackerHit*) tagger_sim_hits->getElementAt(tagger_sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = tagger_sim_hit_decoder(tagger_sim_hit)["layer"];

        if (incident_e == tagger_sim_hit->getMCParticle()) { 
            incident_e_hits[layer - 1]++;
        }
    }
    
    bool findable_tagger_track = true;  
    for (int layer_n = 0; layer_n < 14; layer_n++) { 
        if (incident_e_hits[layer_n] == 0) { 
            findable_tagger_track = false;
        } else if (incident_e_hits[layer_n] > 1) { 
            std::cout << "Incident electron has multiple SimTrackerHits per layer associated with it." << std::endl;
            std::cout << "Layer number: " << layer_n << " Number of hits: " << incident_e_hits[layer_n] << std::endl;
        }
    }
    
    if (findable_tagger_track) {
        findable_tracks++;
        plotter->get1DHistogram("findable particle p")->Fill(beam_e_p);
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
        plotter->get1DHistogram("found particle p")->Fill(beam_e_p); 
    }

    // Loop over all of the Tagger tracks in the event
    double tagger_d0 = 0;
    double tagger_z0 = 0;
    double tagger_pt = 0;  
    for (int tagger_track_n = 0; 
            tagger_track_n < tagger_tracks->getNumberOfElements(); ++tagger_track_n) {
        
        EVENT::Track* tagger_track 
            = (EVENT::Track*) tagger_tracks->getElementAt(tagger_track_n); 

        std::vector<TrackerHit*> track_hits = tagger_track->getTrackerHits();
        plotter->get1DHistogram("Tagger - Hits per track")->Fill(track_hits.size());
         
        double p = TrackUtils::getMomentum(tagger_track, -1.5);
        std::vector<double> p_vec = TrackUtils::getMomentumVector(tagger_track, -1.5);
        double pt = sqrt(p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2]); 
        
        plotter->get1DHistogram("Tagger track p")->Fill(p);
        plotter->get1DHistogram("Tagger track p_{t}")->Fill(pt);
        if (track_hits.size() == 6) plotter->get1DHistogram("Tagger track p - 6 hit")->Fill(p);  
        else if (track_hits.size() == 7) plotter->get1DHistogram("Tagger track p - 7 hit")->Fill(p);

        plotter->get1DHistogram("Tagger curvature")->Fill(tagger_track->getOmega()); 
        plotter->get1DHistogram("Tagger doca")->Fill(tagger_track->getD0());
        plotter->get1DHistogram("Tagger z0")->Fill(tagger_track->getZ0());
        plotter->get1DHistogram("Tagger sin(phi0)")->Fill(sin(tagger_track->getPhi()));
        plotter->get1DHistogram("Tagger tan_lambda")->Fill(tagger_track->getTanLambda());
        //plotter->get1DHistogram("Tagger cos(theta)")->Fill(TrackExtrapolator::getCosTheta(tagger_track));

        plotter->get2DHistogram("Tagger track p vs chi2")->Fill(p, tagger_track->getChi2());

        //double mc_pt = sqrt(sim_px*sim_px + sim_py*sim_py);
        /*if (tagger_tracks->getNumberOfElements() == 1 && n_electrons == 1) { 

            mc_p = sqrt(beam_e_pvec[0]*beam_e_pvec[0] + beam_e_pvec[1]*beam_e_pvec[1] + beam_e_pvec[2]*beam_e_pvec[2]);
            //std::cout << "mc pt: " << mc_pt << " pt: " << pt << std::endl;
            plotter->get1DHistogram("Tagger track p - truth p")->Fill(p - mc_p);
            plotter->get1DHistogram("Tagger track p_{t} - truth p_{t}")->Fill(pt - mc_pt);
        }*/

        plotter->get2DHistogram("Tagger track momentum vs chi2")->Fill(p, tagger_track->getChi2());
        plotter->get2DHistogram("Tagger track momentum vs d0")->Fill(p, tagger_track->getD0());
        plotter->get2DHistogram("Tagger track momentum vs z0")->Fill(p, tagger_track->getZ0());

        if (tagger_track->getD0() >  1.2) continue; 
        plotter->get1DHistogram("Tagger track p - truth p - d0 cut")->Fill(p - beam_e_p);
        //plotter->get1DHistogram("Tagger track p_x - truth p_x - d0 cut")->Fill(p_vec[1] - sim_py);
        //plotter->get1DHistogram("Tagger track p_y - truth p_y - d0 cut")->Fill(p_vec[2] - sim_px);
        plotter->get1DHistogram("Tagger track p - d0 cut")->Fill(p);
        plotter->get1DHistogram("Tagger track p_x - d0 cut")->Fill(p_vec[2]);
        plotter->get1DHistogram("Tagger track p_y - d0 cut")->Fill(p_vec[1]);
        plotter->get1DHistogram("Tagger track p_{t} - d0 cut")->Fill(pt);
        //plotter->get1DHistogram("Tagger track p_{t} - truth p_{t} - d0 cut")->Fill(pt - beam_e_pt);
        if (track_hits.size() == 6) plotter->get1DHistogram("Tagger track p - 6 hit - d0 cut")->Fill(p);  
        else if (track_hits.size() == 7) plotter->get1DHistogram("Tagger track p - 7 hit - d0 cut")->Fill(p);
        plotter->get2DHistogram("Tagger track momentum vs chi2 - d0 cut")->Fill(p, tagger_track->getChi2());
        plotter->get2DHistogram("Tagger track momentum vs d0 - d0 cut")->Fill(p, tagger_track->getD0());
        plotter->get2DHistogram("Tagger track momentum vs z0 - d0 cut")->Fill(p, tagger_track->getZ0());
        //plotter->get2DHistogram("Tagger track p_{t} - truth p_{t} vs p_{t} - d0 cut")->Fill(beam_e_pt, pt - beam_e_pt);
        tagger_d0 = tagger_track->getD0(); 
        tagger_z0 = tagger_track->getZ0(); 
        tagger_pt = pt;  

    }


    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    /*EVENT::LCCollection* recoil_sim_hits 
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
    }*/

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
        std::vector<double> p_vec = TrackUtils::getMomentumVector(recoil_track, -1.5);
        double pt = sqrt(p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2]); 
        //std::cout << "Recoil p: " << p << std::endl; 
        plotter->get1DHistogram("Recoil track p")->Fill(p);
        plotter->get2DHistogram("Recoil track p vs chi2")->Fill(p, recoil_track->getChi2());

        plotter->get1DHistogram("Recoil curvature")->Fill(recoil_track->getOmega()); 
        plotter->get1DHistogram("Recoil doca")->Fill(recoil_track->getD0());
        plotter->get1DHistogram("Recoil z0")->Fill(recoil_track->getZ0());
        plotter->get1DHistogram("Recoil sin(phi0)")->Fill(sin(recoil_track->getPhi()));
        plotter->get1DHistogram("Recoil tan_lambda")->Fill(recoil_track->getTanLambda());
        //plotter->get1DHistogram("Recoil cos(theta)")->Fill(TrackExtrapolator::getCosTheta(recoil_track));

        //double beam_e_pt = sqrt(sim_recoil_px*sim_recoil_px + sim_recoil_py*sim_recoil_py);
        //if (recoil_tracks->getNumberOfElements() == 1 && n_electrons == 1) { 

        //beam_e_p = sqrt(sim_recoil_px*sim_recoil_px +sim_recoil_py*sim_recoil_py +sim_recoil_pz*sim_recoil_pz );
        //double* beam_e_p_end = (double*) incident_e->getMomentumAtEndpoint();
        plotter->get1DHistogram("Recoil track p - truth p")->Fill(p - beam_e_p);
        plotter->get2DHistogram("Recoil track p - truth p vs p truth")->Fill(p - beam_e_p, beam_e_p);
        plotter->get1DHistogram("Recoil track p_{t}")->Fill(pt);
        //plotter->get1DHistogram("Recoil track p_{t} - truth p_{t}")->Fill(pt - beam_e_pt);
        //plotter->get2DHistogram("Recoil track p_{t} - truth p_{t} vs p_{t}")->Fill(beam_e_pt, pt - beam_e_pt); 
       
        if (tagger_d0 == 0) continue;

        plotter->get1DHistogram("Recoil d0 - Tagger d0")->Fill(recoil_track->getD0()*-1 - tagger_d0); 
        plotter->get1DHistogram("Recoil z0 - Tagger z0")->Fill(recoil_track->getZ0()*-1 - tagger_z0); 
        plotter->get1DHistogram("Recoil pt - Tagger pt")->Fill(pt - tagger_pt); 

        //} 
    }
}

void TrackAnalysis::finalize() { 
    plotter->saveToRootFile("track_analysis.root");

    std::cout << "[ ]: Simple Tracking Efficiency: " << (found_tracks/findable_tracks)*100 << "%" << std::endl;
}
