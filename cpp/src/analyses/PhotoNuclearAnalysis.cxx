/**
 *	@file PhotoNuclearAnalysis.cxx
 *	@brief Analysis used to study the performance of the Tagger tracker.
 *	@author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 27, 2016
 */

#include <PhotoNuclearAnalysis.h>

PhotoNuclearAnalysis::PhotoNuclearAnalysis() 
    : tuple(new FlatTupleMaker("photonuclear_tuple.root", "results")), 
      plotter(new Plotter()) {
    LcioAbstractAnalysis::class_name = "PhotoNuclearAnalysis";
}

PhotoNuclearAnalysis::~PhotoNuclearAnalysis() { 
}

void PhotoNuclearAnalysis::initialize() {

    TH1* plot = nullptr; 
    plot = plotter->build2DHistogram("Simulated Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
    plot->GetXaxis()->SetTitle("Layer number"); 
    plot->GetYaxis()->SetTitle("Total hits");

    tuple->addVector("sim_hit_layer");
    tuple->addVector("sim_hit_pos_x");
    tuple->addVector("sim_hit_pos_y");
    tuple->addVector("sim_hit_pos_z");
    tuple->addVector("sim_hit_time");
}

void PhotoNuclearAnalysis::processEvent(EVENT::LCEvent* event) { 

    tuple->setVariableValue("event", event->getEventNumber());

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
        = (EVENT::LCCollection*) event->getCollection("MCParticle");

    if (mc_particles->getNumberOfElements() == 1) return;

    // Create a pointer to the incident electron i.e. the beam electron. 
    EVENT::MCParticle* incident_e = nullptr;

    //std::cout << "Number of MC particles: " << mc_particles->getNumberOfElements() << std::endl;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    bool is_photonuclear = false;
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
        //std::cout << "Particle PDG ID: " << particle->getPDG() << std::endl;

        if (particle->getParents().size() == 0) {
            //std::cout << "Incident particle: " << particle->getPDG() 
            //          << " has " << particle->getDaughters().size() << " daughters." << std::endl;
            //incident_e = particle;
            for (auto daughter : particle->getDaughters()) { 
                if (daughter->getPDG() > 100 && daughter->getPDG() < 3000) { 
                    is_photonuclear = true;
                    std::cout << "Daughter: " << daughter->getPDG() << std::endl;
                }
            }
        }
    }
    
    if (!is_photonuclear) return; 


    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> sim_hit_decoder(sim_hits);

    std::vector<int> sim_hits_vec(10, 0);
    for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 

        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];
        sim_hits_vec[layer - 1]++;

        tuple->addToVector("sim_hit_layer", layer);
        tuple->addToVector("sim_hit_pos_x", sim_hit->getPosition()[0]); 
        tuple->addToVector("sim_hit_pos_y", sim_hit->getPosition()[1]); 
        tuple->addToVector("sim_hit_pos_z", sim_hit->getPosition()[2]); 
        tuple->addToVector("sim_hit_time", sim_hit->getTime());
    }

    for (int layer_n = 0; layer_n < 10; layer_n++) { 
        plotter->get2DHistogram("Simulated Recoil Hits per Layer")->Fill(layer_n+1, sim_hits_vec[layer_n]);
    }

    // Check if the track is findable
    /*tuple->setVariableValue("is_findable", 0);
      std::vector<EVENT::SimTrackerHit*> findable_sim_hits; 
      if (TrackUtils::isTrackFindable(14, incident_e, sim_hits, findable_sim_hits)) {
      ++findable_track;
      tuple->setVariableValue("is_findable", 1);  
      }

      for (EVENT::SimTrackerHit* sim_hit : findable_sim_hits) { 

    // Get the layer number associated with this hit.
    int layer = sim_hit_decoder(sim_hit)["layer"];

    tuple->addToVector("sim_hit_layer", layer);
    tuple->addToVector("sim_hit_pos_x", sim_hit->getPosition()[0]); 
    tuple->addToVector("sim_hit_pos_y", sim_hit->getPosition()[1]); 
    tuple->addToVector("sim_hit_pos_z", sim_hit->getPosition()[2]); 
    tuple->addToVector("sim_hit_time", sim_hit->getTime());
    }

    // Get the collection of Tagger tracker tracks from the event.  If no such
    // collection exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* tracks 
    = (EVENT::LCCollection*) event->getCollection("TaggerTracks");

    tuple->setVariableValue("n_tracks", tracks->getNumberOfElements());

    for (int track_n = 0; track_n < tracks->getNumberOfElements(); ++track_n) {

    EVENT::Track* track = (EVENT::Track*) tracks->getElementAt(track_n); 

    std::vector<TrackerHit*> track_hits = track->getTrackerHits();

    double p = TrackUtils::getMomentum(track, -1.5);
    std::vector<double> p_vec = TrackUtils::getMomentumVector(track, -1.5);
    double pt = sqrt(p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2]); 

    tuple->setVariableValue("n_track_hits", track_hits.size());
    tuple->setVariableValue("chi2", track->getChi2());
    tuple->setVariableValue("p", p);
    tuple->setVariableValue("px", p_vec[0]);
    tuple->setVariableValue("py", p_vec[1]);
    tuple->setVariableValue("pz", p_vec[2]);
    tuple->setVariableValue("d0", track->getD0());
    tuple->setVariableValue("z0", track->getZ0());
    tuple->setVariableValue("omega", track->getOmega());
    tuple->setVariableValue("phi0", track->getPhi());
    tuple->setVariableValue("tan_lambda", track->getTanLambda());
    }
    */

    tuple->fill();
}

void PhotoNuclearAnalysis::finalize() { 
    tuple->close();
}
