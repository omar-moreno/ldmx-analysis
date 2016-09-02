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
      filter_pn(true), 
      findable_track(0) {
    LcioAbstractAnalysis::class_name = "TaggerTrackerAnalysis";
}

TaggerTrackerAnalysis::~TaggerTrackerAnalysis() { 
}

void TaggerTrackerAnalysis::initialize() {
   
    tuple->addVariable("beam_e_p");  
    tuple->addVariable("beam_e_px");  
    tuple->addVariable("beam_e_py");  
    tuple->addVariable("beam_e_pz");  
    tuple->addVariable("beam_e_p_last");  
    tuple->addVariable("beam_e_px_last");  
    tuple->addVariable("beam_e_py_last");  
    tuple->addVariable("beam_e_pz_last");  
    tuple->addVariable("chi2"); 
    tuple->addVariable("d0");
    tuple->addVariable("event");
    tuple->addVariable("n_tracks");
    tuple->addVariable("n_track_hits");
    tuple->addVariable("n_mishits"); 
    tuple->addVariable("omega");
    tuple->addVariable("p");
    tuple->addVariable("phi0");
    tuple->addVariable("px");
    tuple->addVariable("py");
    tuple->addVariable("pz");
    tuple->addVariable("is_findable");
    tuple->addVariable("tan_lambda");
    tuple->addVariable("x_target");
    tuple->addVariable("y_target");
    tuple->addVariable("z0");

    tuple->addVector("tagger_sim_hit_layer");
    tuple->addVector("tagger_sim_hit_dedx");
    tuple->addVector("tagger_sim_hit_backscattered");
    tuple->addVector("tagger_sim_hit_time");
    tuple->addVector("tagger_sim_hit_pos_x");
    tuple->addVector("tagger_sim_hit_pos_y");
    tuple->addVector("tagger_sim_hit_pos_z");

    tuple->addVector("findable_sim_hit_layer");
    tuple->addVector("findable_sim_hit_pos_x");
    tuple->addVector("findable_sim_hit_pos_y");
    tuple->addVector("findable_sim_hit_pos_z");
    tuple->addVector("findable_sim_hit_time");
}

void TaggerTrackerAnalysis::processEvent(EVENT::LCEvent* event) { 
   
    tuple->setVariableValue("event", event->getEventNumber());

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
       = (EVENT::LCCollection*) event->getCollection("MCParticle");
   
    // Create a pointer to the incident electron i.e. the beam electron. 
    EVENT::MCParticle* beam_e = nullptr;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    //int n_electrons = 0;
    //double largest_brem_energy = 0; 
    //double brem_position = -1; 
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
    
        //if (particle->getPDG() == 11) n_electrons++; 

        if (particle->getParents().size() == 0) {
           beam_e = particle; 
        }

        if (filter_pn && particle->getPDG() > 100 && createdWithinTarget(particle)) { 
            std::cout << "Skipping event with photo-nuclear reaction in target" << std::endl;
            return;
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

    // Calculate the momentum of the incident electron
    double* beam_e_pvec = (double*) beam_e->getMomentum();
    double beam_e_p = sqrt(pow(beam_e_pvec[0], 2) + pow(beam_e_pvec[1], 2) + pow(beam_e_pvec[2], 2));
    tuple->setVariableValue("beam_e_p", beam_e_p); 
    tuple->setVariableValue("beam_e_px", beam_e_pvec[0]); 
    tuple->setVariableValue("beam_e_py", beam_e_pvec[1]); 
    tuple->setVariableValue("beam_e_pz", beam_e_pvec[2]); 

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* sim_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> sim_hit_decoder(sim_hits);

    std::vector<int> sim_hits_vec(14, 0);
    for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 
        

        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit 
            = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n); 

        if (sim_hit->getMCParticle() == beam_e) continue;

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];
        sim_hits_vec[layer - 1]++;

        tuple->addToVector("tagger_sim_hit_backscattered", 0);
        if (sim_hit->getMCParticle()->isBackscatter()) 
            tuple->addToVector("tagger_sim_hit_backscattered", 1);
        tuple->addToVector("tagger_sim_hit_layer", layer);
        tuple->addToVector("tagger_sim_hit_dedx",  sim_hit->getEDep());
        tuple->addToVector("tagger_sim_hit_time",  sim_hit->getTime());
        tuple->addToVector("tagger_sim_hit_pos_x", sim_hit->getPosition()[0]); 
        tuple->addToVector("tagger_sim_hit_pos_y", sim_hit->getPosition()[1]); 
        tuple->addToVector("tagger_sim_hit_pos_z", sim_hit->getPosition()[2]); 

    }

    // Check if the track is findable
    tuple->setVariableValue("is_findable", 0);
    std::vector<EVENT::SimTrackerHit*> findable_sim_hits;
    int strategy[14] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
    if (TrackUtils::isTrackFindable(14, strategy, beam_e, sim_hits, findable_sim_hits)) {
        ++findable_track;
        tuple->setVariableValue("is_findable", 1);  
    }
    
    for (EVENT::SimTrackerHit* findable_sim_hit : findable_sim_hits) { 
        
        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(findable_sim_hit)["layer"];

        tuple->addToVector("findable_sim_hit_layer", layer);
        tuple->addToVector("findable_sim_hit_pos_x", findable_sim_hit->getPosition()[0]); 
        tuple->addToVector("findable_sim_hit_pos_y", findable_sim_hit->getPosition()[1]); 
        tuple->addToVector("findable_sim_hit_pos_z", findable_sim_hit->getPosition()[2]); 
        tuple->addToVector("findable_sim_hit_time", findable_sim_hit->getTime());
    
        if ((beam_e == findable_sim_hit->getMCParticle()) && (layer == 14)) { 
    
            // Calculate the momentum of the incident electron
            float* findable_sim_hit_pvec = (float*) findable_sim_hit->getMomentum();
            double findable_sim_hit_p = sqrt(pow(findable_sim_hit_pvec[0], 2) + pow(findable_sim_hit_pvec[1], 2) + pow(findable_sim_hit_pvec[2], 2));
            tuple->setVariableValue("beam_e_p_last", findable_sim_hit_p); 
            tuple->setVariableValue("beam_e_px_last", findable_sim_hit_pvec[0]); 
            tuple->setVariableValue("beam_e_py_last", findable_sim_hit_pvec[1]); 
            tuple->setVariableValue("beam_e_pz_last", findable_sim_hit_pvec[2]); 
        }
    }

    // Get the collection of LCRelations between a raw hit and a SimTrackerHit.
    EVENT::LCCollection* true_hit_relations
        = (EVENT::LCCollection*) event->getCollection("TaggerTrueHitRelations");

    // Instantiate an LCRelation navigator which allows faster access to either
    // the SimTrackerHit or raw hits.
    UTIL::LCRelationNavigator* true_hit_relations_nav
        = new UTIL::LCRelationNavigator(true_hit_relations);

    // Get the collection of Tagger tracker tracks from the event.  If no such
    // collection exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* tracks 
        = (EVENT::LCCollection*) event->getCollection("TaggerTracks");
   
    tuple->setVariableValue("n_tracks", tracks->getNumberOfElements());
   
    for (int track_n = 0; track_n < tracks->getNumberOfElements(); ++track_n) {
        
        EVENT::Track* track = (EVENT::Track*) tracks->getElementAt(track_n); 
    
        std::vector<TrackerHit*> track_hits = track->getTrackerHits();

        double n_mishit = 0;  
        for (TrackerHit* track_hit : track_hits) {
            bool mishit_found = false;
            for (int raw_hit_n = 0; raw_hit_n < track_hit->getRawHits().size(); ++raw_hit_n) { 
                
                EVENT::TrackerRawData* tagger_track_raw_hit 
                    = (EVENT::TrackerRawData*) track_hit->getRawHits()[raw_hit_n];
                
                EVENT::LCObjectVec sim_hit_list 
                    = true_hit_relations_nav->getRelatedToObjects(tagger_track_raw_hit);
                
                if (((EVENT::SimTrackerHit*) sim_hit_list[0])->getMCParticle() != beam_e) { 
                    mishit_found = true;
                }
            }
            if (mishit_found) n_mishit++; 
        } 
        tuple->setVariableValue("n_mishits", n_mishit); 

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
    tuple->clear();
}

void TaggerTrackerAnalysis::finalize() { 
    tuple->close();
}

bool TaggerTrackerAnalysis::createdWithinTarget(MCParticle* particle) {
    //std::cout << "Checking if particle was created within the target." << std::endl;
    if (particle->getVertex()[2] > -.1750 && particle->getVertex()[2] < .1750) {
        //std::cout << "Particle created within the target." << std::endl;
        return true;
    }
    //std::cout << "Particle created downstream of the target." << std::endl;
    return false; 
}

void TaggerTrackerAnalysis::filterPhotoNuclearEvents(bool filter_pn) { 
    this->filter_pn = filter_pn;
}
