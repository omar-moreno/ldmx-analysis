/**
 *	@file RecoilTrackerAnalysis.cxx
 *	@brief Analysis used to study the performance of the Recoil tracker.
 *	@author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 27, 2016
 */

#include <RecoilTrackerAnalysis.h>

RecoilTrackerAnalysis::RecoilTrackerAnalysis() 
    : tuple(new FlatTupleMaker("recoil_tracker_tuple.root", "results")),
      filter_pn(false) { 
      //findable_track(0) {
    LcioAbstractAnalysis::class_name = "RecoilTrackerAnalysis";
}

RecoilTrackerAnalysis::~RecoilTrackerAnalysis() { 
}

void RecoilTrackerAnalysis::initialize() {
 
    // Event information 
    tuple->addVariable("event");
    tuple->addVariable("n_tracks");
    tuple->addVariable("n_raw_hits");
    tuple->addVariable("n_3d_hits");

    tuple->addVariable("target_energy");
    tuple->addVariable("trigger_pad_energy");

    // Recoil track parameters and kinematics

    tuple->addVariable("n_recoil_mishits"); 
    tuple->addVariable("n_recoil_hits");
    tuple->addVariable("recoil_is_findable");
    tuple->addVariable("recoil_is_found");
    tuple->addVariable("recoil_chi2"); 
    tuple->addVariable("recoil_d0");
    tuple->addVariable("recoil_omega");
    tuple->addVariable("recoil_phi0");
    tuple->addVariable("recoil_tan_lambda");
    tuple->addVariable("recoil_z0");
    
    tuple->addVariable("recoil_p");
    tuple->addVariable("recoil_pt");
    tuple->addVariable("recoil_px");
    tuple->addVariable("recoil_py");
    tuple->addVariable("recoil_pz");

    tuple->addVariable("recoil_truth_p_first");
    tuple->addVariable("recoil_truth_pt_first");
    tuple->addVariable("recoil_truth_px_first");
    tuple->addVariable("recoil_truth_py_first");
    tuple->addVariable("recoil_truth_pz_first");

    tuple->addVariable("recoil_truth_p_last");
    tuple->addVariable("recoil_truth_pt_last");
    tuple->addVariable("recoil_truth_px_last");
    tuple->addVariable("recoil_truth_py_last");
    tuple->addVariable("recoil_truth_pz_last");

    tuple->addVariable("hardest_brem_energy");
    tuple->addVariable("hardest_brem_pos_z");

    tuple->addVariable("n_recoil_sim_hits");
    tuple->addVector("recoil_sim_hit_layer");
    tuple->addVector("recoil_sim_hit_dedx");
    tuple->addVector("recoil_sim_hit_backscattered");
    tuple->addVector("recoil_sim_hit_time");
    tuple->addVector("recoil_sim_hit_pos_x");
    tuple->addVector("recoil_sim_hit_pos_y");
    tuple->addVector("recoil_sim_hit_pos_z");

    tuple->addVariable("n_sim_hits");
    tuple->addVector("sim_hit_layer");
    tuple->addVector("sim_hit_dedx");
    tuple->addVector("sim_hit_backscattered");
    tuple->addVector("sim_hit_time");
    tuple->addVector("sim_hit_pos_x");
    tuple->addVector("sim_hit_pos_y");
    tuple->addVector("sim_hit_pos_z");

    tuple->addVariable("n_recoil_clusters");
    tuple->addVector("recoil_cluster_size");
    tuple->addVector("recoil_cluster_layer");
    tuple->addVector("recoil_cluster_dedx");
    tuple->addVector("recoil_cluster_time");
    tuple->addVector("recoil_cluster_pos_x");
    tuple->addVector("recoil_cluster_pos_y");
    tuple->addVector("recoil_cluster_pos_z");
    


    /*/reco
    tuple->addVariable("x_target");
    tuple->addVariable("y_target");

    */
}

void RecoilTrackerAnalysis::processEvent(EVENT::LCEvent* event) { 
   
    tuple->setVariableValue("event", event->getEventNumber());
    //std::cout << "Event: " << event->getEventNumber() << std::endl;
   
    // Get the collection of TargetHits from the event.  If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* target_hits
        = (EVENT::LCCollection*) event->getCollection("TargetHits");
    
    if (target_hits->getNumberOfElements() != 0) { 
        EVENT::CalorimeterHit* target_hit 
            = (EVENT::CalorimeterHit*) target_hits->getElementAt(0);
    
        // Get the energy deposited in the target
        tuple->setVariableValue("target_energy", target_hit->getEnergy()*1000);
    }

    EVENT::LCCollection* trigger_pad_hits
        = (EVENT::LCCollection*) event->getCollection("TriggerPadHits");

    if (trigger_pad_hits->getNumberOfElements() != 0) { 
        EVENT::CalorimeterHit* trigger_pad_hit 
            = (EVENT::CalorimeterHit*) trigger_pad_hits->getElementAt(0);
   
        // Get the energy deposited in the trigger pad
        tuple->setVariableValue("trigger_pad_energy", trigger_pad_hit->getEnergy()*1000); 
    }

    /*if (filter_pn) {
    
        // Get the collection of MC particles from the event. If no such collection
        // exist, a DataNotAvailableException is thrown.
        EVENT::LCCollection* mc_particles 
            = (EVENT::LCCollection*) event->getCollection("MCParticle");
        for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {
       
            EVENT::MCParticle* mc_particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);

            if (mc_particle->getPDG() > 100 && createdWithinTarget(mc_particle)) {
                std::cout << "Skipping event with photo-nuclear reaction in target" << std::endl;
                return;
            }
        }
    }*/

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
       = (EVENT::LCCollection*) event->getCollection("MCParticle");
   
    // Create a pointer that will be used to store the recoil electron. 
    EVENT::MCParticle* recoil_e = nullptr;

    // Loop over all of the MC particles and find the recoil electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    double hardest_brem_energy = -10000; 
    double hardest_brem_position = -10000; 
    for (int particle_n = 0; particle_n < mc_particles->getNumberOfElements(); ++particle_n) {

        EVENT::MCParticle* particle = (EVENT::MCParticle*) mc_particles->getElementAt(particle_n);
    
        if (particle->getPDG() == 11 && particle->getParents().size() == 0) {
           recoil_e = particle; 
        }

        if (particle->getPDG() == 22 
                && particle->getParents()[0]->getPDG() == 11
                && particle->getParents()[0]->getParents().size() == 0) { 
            
            double* brem_energy_vec = (double*) particle->getMomentum(); 
            double brem_energy 
                = sqrt(pow(brem_energy_vec[0], 2) + pow(brem_energy_vec[1], 2) + pow(brem_energy_vec[2], 2));

            if (brem_energy > hardest_brem_energy) {
                hardest_brem_energy = brem_energy;
                hardest_brem_position = particle->getVertex()[2];
            }
        }
    }

    tuple->setVariableValue("hardest_brem_energy", hardest_brem_energy);
    tuple->setVariableValue("hardest_brem_pos_z", hardest_brem_position);

    // Get the collection of SimTrackerHits associated with the Recoil tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");
    tuple->setVariableValue("n_sim_hits", sim_hits->getNumberOfElements());

    // Create a cell ID decoder used to get specific properties associated with
    // Recoil SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> sim_hit_decoder(sim_hits);

    std::vector<int> sim_hits_vec(10, 0);
    int n_recoil_sim_hits = 0;
    for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 
        
        // Get a Recoil SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit 
            = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];

        if (recoil_e == sim_hit->getMCParticle()) { 
            
            ++n_recoil_sim_hits;

            tuple->addToVector("recoil_sim_hit_layer", layer);
            tuple->addToVector("recoil_sim_hit_dedx",  sim_hit->getEDep());
            tuple->addToVector("recoil_sim_hit_time",  sim_hit->getTime());
            tuple->addToVector("recoil_sim_hit_pos_x", sim_hit->getPosition()[0]); 
            tuple->addToVector("recoil_sim_hit_pos_y", sim_hit->getPosition()[1]); 
            tuple->addToVector("recoil_sim_hit_pos_z", sim_hit->getPosition()[2]); 
        
            tuple->addToVector("recoil_sim_hit_backscattered", 0);
            if (sim_hit->getMCParticle()->isBackscatter()) 
                tuple->addToVector("recoil_sim_hit_backscattered", 1);
        }
            
        tuple->addToVector("sim_hit_layer", layer);
        tuple->addToVector("sim_hit_dedx",  sim_hit->getEDep());
        tuple->addToVector("sim_hit_time",  sim_hit->getTime());
        tuple->addToVector("sim_hit_pos_x", sim_hit->getPosition()[0]); 
        tuple->addToVector("sim_hit_pos_y", sim_hit->getPosition()[1]); 
        tuple->addToVector("sim_hit_pos_z", sim_hit->getPosition()[2]); 

        tuple->addToVector("sim_hit_backscattered", 0);
        if (sim_hit->getMCParticle()->isBackscatter()) 
            tuple->addToVector("sim_hit_backscattered", 1);
    
    }

    // Check if the recoil track is findable
    tuple->setVariableValue("recoil_is_findable", 0);
    int strategy[10] = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0};

    std::vector<EVENT::SimTrackerHit*> findable_sim_hits; 
    EVENT::SimTrackerHit* last_recoil_hit = nullptr;
    if (TrackUtils::isTrackFindable(10, strategy, recoil_e, sim_hits, findable_sim_hits)) {
        tuple->setVariableValue("recoil_is_findable", 1);  
        
    }
    
    int last_layer = 0;
    int first_layer = 15;
    EVENT::SimTrackerHit* last_recoil_sim_hit = nullptr;
    EVENT::SimTrackerHit* first_recoil_sim_hit = nullptr;
    for (EVENT::SimTrackerHit* sim_hit : findable_sim_hits) {
        
        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];
        
        if (layer > last_layer) { 
            last_layer = layer;
            last_recoil_sim_hit = sim_hit;
        }

        if (layer < first_layer) { 
            first_layer = layer; 
            first_recoil_sim_hit = sim_hit;
        }
    } 

    if (first_recoil_sim_hit != nullptr) { 
        float* recoil_truth_pvec = (float*) first_recoil_sim_hit->getMomentum();
        double recoil_truth_p 
            = sqrt(pow(recoil_truth_pvec[0], 2) + pow(recoil_truth_pvec[1], 2) + pow(recoil_truth_pvec[2], 2));
        double recoil_truth_pt 
            = sqrt(pow(recoil_truth_pvec[0], 2) + pow(recoil_truth_pvec[1], 2));
        tuple->setVariableValue("recoil_truth_p_first", recoil_truth_p); 
        tuple->setVariableValue("recoil_truth_pt_first", recoil_truth_pt); 
        tuple->setVariableValue("recoil_truth_px_first", recoil_truth_pvec[0]); 
        tuple->setVariableValue("recoil_truth_py_first", recoil_truth_pvec[1]); 
        tuple->setVariableValue("recoil_truth_pz_first", recoil_truth_pvec[2]); 
    }

    if (last_recoil_sim_hit != nullptr) { 
        float* recoil_truth_pvec = (float*) last_recoil_sim_hit->getMomentum();
        double recoil_truth_p 
            = sqrt(pow(recoil_truth_pvec[0], 2) + pow(recoil_truth_pvec[1], 2) + pow(recoil_truth_pvec[2], 2));
        double recoil_truth_pt 
            = sqrt(pow(recoil_truth_pvec[0], 2) + pow(recoil_truth_pvec[1], 2));
        tuple->setVariableValue("recoil_truth_p_last", recoil_truth_p); 
        tuple->setVariableValue("recoil_truth_pt_last", recoil_truth_pt); 
        tuple->setVariableValue("recoil_truth_px_last", recoil_truth_pvec[0]); 
        tuple->setVariableValue("recoil_truth_py_last", recoil_truth_pvec[1]); 
        tuple->setVariableValue("recoil_truth_pz_last", recoil_truth_pvec[2]); 
    }

    // Get the collection of raw hits from the event
    EVENT::LCCollection* raw_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilRawTrackerHits");
    tuple->setVariableValue("n_raw_hits", raw_hits->getNumberOfElements());

    // Get the collection of LCRelations between a raw hit and a SimTrackerHit.
    EVENT::LCCollection* true_hit_relations
        = (EVENT::LCCollection*) event->getCollection("RecoilTrueHitRelations");

    // Instantiate an LCRelation navigator which allows faster access to either
    // the SimTrackerHit or raw hits.
    UTIL::LCRelationNavigator* true_hit_relations_nav
        = new UTIL::LCRelationNavigator(true_hit_relations);

    // Get the collection of clusters from the event
    EVENT::LCCollection* clusters 
        = (EVENT::LCCollection*) event->getCollection("RecoilClusters");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger raw hits.
    std::string encoder_string = "system:6,barrel:3,layer:4,module:12,sensor:1,side:32:-2,strip:1";
    UTIL::CellIDDecoder<EVENT::TrackerRawData> raw_hit_decoder(encoder_string);

    int n_clusters = 0;
    for (int cluster_n = 0; cluster_n < clusters->getNumberOfElements(); ++cluster_n) { 
        
        EVENT::TrackerHit* cluster = (EVENT::TrackerHit*) clusters->getElementAt(cluster_n);

        // Get a raw hit from the list of raw hits
        EVENT::TrackerRawData* raw_hit 
            = (EVENT::TrackerRawData*) cluster->getRawHits()[0]; 

        // Get the layer number associated with this hit.
        int layer = raw_hit_decoder(raw_hit)["layer"];

        std::vector<LCObject*> sim_hit_list 
            = true_hit_relations_nav->getRelatedToObjects(raw_hit);
        
        bool recoil_cluster_found = false; 
        for (LCObject* hit : sim_hit_list) {
            if (recoil_e == ((EVENT::SimTrackerHit*) hit)->getMCParticle()) {
                recoil_cluster_found = true;
                n_clusters++; 
            }
        }

        if (recoil_cluster_found) {
            tuple->addToVector("recoil_cluster_size", cluster->getRawHits().size()); 
            tuple->addToVector("recoil_cluster_layer", layer);
            tuple->addToVector("recoil_cluster_dedx",  cluster->getEDep());
            tuple->addToVector("recoil_cluster_time",  cluster->getTime());
            tuple->addToVector("recoil_cluster_pos_x", cluster->getPosition()[0]); 
            tuple->addToVector("recoil_cluster_pos_y", cluster->getPosition()[1]); 
            tuple->addToVector("recoil_cluster_pos_z", cluster->getPosition()[2]); 
        }
    } 

    EVENT::LCCollection* stereo_hits
        = (EVENT::LCCollection*) event->getCollection("RecoilHelicalTrackHits"); 
    tuple->setVariableValue("n_3d_hits", stereo_hits->getNumberOfElements()); 

    // Get the collection of Recoil tracker tracks from the event.  If no such
    // collection exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* tracks 
        = (EVENT::LCCollection*) event->getCollection("RecoilTracks");
 
    tuple->setVariableValue("n_tracks", tracks->getNumberOfElements());
    if (tracks->getNumberOfElements() == 0) { 
        tuple->fill();
        tuple->clear();
        return;
    }
 
    tuple->setVariableValue("recoil_is_found", 0);

    std::vector<TrackerHit*> recoil_track_hits; 
    EVENT::Track* recoil_e_track = nullptr;
    double n_recoil_hits = 0;
    for (int track_n = 0; track_n < tracks->getNumberOfElements(); ++track_n) {

        bool is_recoil_track = false; 
        EVENT::Track* track = (EVENT::Track*) tracks->getElementAt(track_n); 
        
        std::vector<TrackerHit*> track_hits = track->getTrackerHits();
        
        double p = TrackUtils::getMomentum(track, -0.75);
        std::vector<double> p_vec = TrackUtils::getMomentumVector(track, -0.75);
        double pt = sqrt(p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2]); 
        
        //std::cout << "[ RecoilTrackerAnalysis ]: Track momentum " << p << " GeV" << std::endl;
        //std::cout << "[ RecoilTrackerAnalysis ]: Track hits " << track_hits.size() << std::endl;

        double current_n_recoil_hits = 0; 
        for (TrackerHit* track_hit : track_hits) {
            //std::cout << "[ RecoilTrackerAnalysis ]: Total number of raw hits " << track_hit->getRawHits().size() << std::endl;
            bool recoil_hit_found = true;
            for (int raw_hit_n = 0; raw_hit_n < track_hit->getRawHits().size(); ++raw_hit_n) { 
                
                EVENT::TrackerRawData* raw_hit 
                    = (EVENT::TrackerRawData*) track_hit->getRawHits()[raw_hit_n];
                
                std::vector<LCObject*> sim_hit_list 
                    = true_hit_relations_nav->getRelatedToObjects(raw_hit);
               
                for (LCObject* hit : sim_hit_list) {
                    //std::cout << "[ RecoilTrackerAnalysis ]: Sim particle " << ((EVENT::SimTrackerHit*) hit)->getMCParticle()->getPDG() << std::endl;
                    if (recoil_e != ((EVENT::SimTrackerHit*) hit)->getMCParticle()) {
                        //std::cout << "[ RecoilTrackerAnalysis ]: This is not a recoil electron hit." << std::endl;
                        recoil_hit_found = false;
                        break;
                    } 
                }
            }
            if (recoil_hit_found) current_n_recoil_hits++; 
        }

        if (current_n_recoil_hits > n_recoil_hits) { 
            n_recoil_hits = current_n_recoil_hits;
            recoil_e_track = track; 
        }

    } 

    if (recoil_e_track != nullptr) { 
        
        double p = TrackUtils::getMomentum(recoil_e_track, -0.75);
        std::vector<double> p_vec = TrackUtils::getMomentumVector(recoil_e_track, -0.75);
        double pt = sqrt(p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2]); 
        
        tuple->setVariableValue("n_recoil_hits", n_recoil_hits);
        tuple->setVariableValue("recoil_is_found", 1); 
        tuple->setVariableValue("recoil_chi2", recoil_e_track->getChi2());
        tuple->setVariableValue("recoil_p", p);
        tuple->setVariableValue("recoil_pt", pt);
        tuple->setVariableValue("recoil_px", p_vec[0]);
        tuple->setVariableValue("recoil_py", p_vec[1]);
        tuple->setVariableValue("recoil_pz", p_vec[2]);
        tuple->setVariableValue("recoil_d0", recoil_e_track->getD0());
        tuple->setVariableValue("recoil_z0", recoil_e_track->getZ0());
        tuple->setVariableValue("recoil_omega", recoil_e_track->getOmega());
        tuple->setVariableValue("recoil_phi0", recoil_e_track->getPhi());
        tuple->setVariableValue("recoil_tan_lambda", recoil_e_track->getTanLambda());
   
        double n_mishit = 0;  
        bool mishit_found = false;
        for (TrackerHit* track_hit : recoil_e_track->getTrackerHits()) {
            for (int raw_hit_n = 0; raw_hit_n < track_hit->getRawHits().size(); ++raw_hit_n) { 
                EVENT::TrackerRawData* raw_hit 
                    = (EVENT::TrackerRawData*) track_hit->getRawHits()[raw_hit_n];
                
                std::vector<LCObject*> sim_hit_list 
                    = true_hit_relations_nav->getRelatedToObjects(raw_hit);
               
                for (LCObject* hit : sim_hit_list) {
                    if (recoil_e != ((EVENT::SimTrackerHit*) hit)->getMCParticle()) {
                        mishit_found = true;
                        break;
                    } 
                }
            }
            if (mishit_found) n_mishit++; 
        }
        tuple->setVariableValue("n_recoil_mishits", n_mishit); 
    }

    tuple->fill();
    tuple->clear();

}

void RecoilTrackerAnalysis::finalize() { 
    tuple->close();
}

bool RecoilTrackerAnalysis::createdWithinTarget(MCParticle* particle) {
    //std::cout << "Checking if particle was created within the target." << std::endl;
    if (particle->getVertex()[2] > -.1750 && particle->getVertex()[2] < .1750) {
        //std::cout << "Particle created within the target." << std::endl;
        return true;
    }
    //std::cout << "Particle created downstream of the target." << std::endl;
    return false; 
}

void RecoilTrackerAnalysis::filterPhotoNuclearEvents(bool filter_pn) { 
    this->filter_pn = filter_pn;
}
