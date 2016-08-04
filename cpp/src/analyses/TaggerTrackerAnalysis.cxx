/**
 *	@file TaggerTrackerAnalysis.cxx
 *	@brief Analysis used to study the performance of the Tagger tracker.
 *	@author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 27, 2016
 */

#include <TaggerTrackerAnalysis.h>

TaggerTrackerAnalysis::TaggerTrackerAnalysis() 
    : tuple(new FlatTupleMaker("tagger_tracker_tuple.root", "results")),
      findable_track(0) {
    LcioAbstractAnalysis::class_name = "TaggerTrackerAnalysis";
}

TaggerTrackerAnalysis::~TaggerTrackerAnalysis() { 
}

void TaggerTrackerAnalysis::initialize() {
    
    tuple->addVariable("chi2"); 
    tuple->addVariable("d0");
    tuple->addVariable("event");
    tuple->addVariable("n_tracks");
    tuple->addVariable("n_track_hits");
    tuple->addVariable("omega");
    tuple->addVariable("p");
    tuple->addVariable("phi0");
    tuple->addVariable("px");
    tuple->addVariable("py");
    tuple->addVariable("pz");
    tuple->addVariable("is_findable");
    tuple->addVariable("tan_lambda");
    tuple->addVariable("z0");

    tuple->addVector("sim_hit_layer");
    tuple->addVector("sim_hit_pos_x");
    tuple->addVector("sim_hit_pos_y");
    tuple->addVector("sim_hit_pos_z");
    tuple->addVector("sim_hit_time");
}

void TaggerTrackerAnalysis::processEvent(EVENT::LCEvent* event) { 
   
    tuple->setVariableValue("event", event->getEventNumber());

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
       = (EVENT::LCCollection*) event->getCollection("MCParticle");
   
    // Create a pointer to the incident electron i.e. the beam electron. 
    EVENT::MCParticle* incident_e = nullptr;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    //int n_electrons = 0;
    //double largest_brem_energy = 0; 
    //double brem_position = -1; 
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
    
        //if (particle->getPDG() == 11) n_electrons++; 

        if (particle->getParents().size() == 0) {
           incident_e = particle; 
        }

       /*if (particle->getPDG() == 22 && particle->getVertex()[2] < -3) { 
            double* brem_energy_vec = (double*) particle->getMomentum(); 
            double brem_energy = sqrt(pow(brem_energy_vec[0], 2) + pow(brem_energy_vec[1], 2) + pow(brem_energy_vec[2], 2));
            if (brem_energy > largest_brem_energy) {
                largest_brem_energy = brem_energy;
                brem_position = particle->getVertex()[2];
            }
       }*/ 
    }


    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* sim_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> sim_hit_decoder(sim_hits);

    /*for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 
        
        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit 
            = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];

        tuple->addToVector("sim_hit_layer", layer);
        tuple->addToVector("sim_hit_pos_x", sim_hit->getPosition()[0]); 
        tuple->addToVector("sim_hit_pos_y", sim_hit->getPosition()[1]); 
        tuple->addToVector("sim_hit_pos_z", sim_hit->getPosition()[2]); 
        tuple->addToVector("sim_hit_time", sim_hit->getTime());
    }*/

    // Check if the track is findable
    tuple->setVariableValue("is_findable", 0);
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

    tuple->fill();
}

void TaggerTrackerAnalysis::finalize() { 
    tuple->close();
}
