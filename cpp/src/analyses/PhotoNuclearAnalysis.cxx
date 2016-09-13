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

    tuple->addVariable("event");

    tuple->addVariable("target_energy");
    tuple->addVariable("trigger_pad_energy");

    tuple->addVariable("recoil_is_findable");
    tuple->addVariable("recoil_energy");

    tuple->addVariable("pn_gamma_energy");
    tuple->addVariable("pn_particle_mult");

    tuple->addVariable("n_recoil_hits");
    tuple->addVariable("n_recoil_pn_hits");
    tuple->addVariable("n_tracks");

    tuple->addVector("pn_pdg_id");
    tuple->addVector("pn_theta");
    tuple->addVector("pn_phi");
    tuple->addVector("is_findable");

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

    tuple->addVector("recoil_sim_hit_layer");
    tuple->addVector("recoil_sim_hit_dedx");
    tuple->addVector("recoil_sim_hit_backscattered");
    tuple->addVector("recoil_sim_hit_time");
    tuple->addVector("recoil_sim_hit_pos_x");
    tuple->addVector("recoil_sim_hit_pos_y");
    tuple->addVector("recoil_sim_hit_pos_z");

    tuple->addVector("tagger_sim_hit_layer");
    tuple->addVector("tagger_sim_hit_dedx");
    tuple->addVector("tagger_sim_hit_backscattered");
    tuple->addVector("tagger_sim_hit_time");
    tuple->addVector("tagger_sim_hit_pos_x");
    tuple->addVector("tagger_sim_hit_pos_y");
    tuple->addVector("tagger_sim_hit_pos_z");
}

void PhotoNuclearAnalysis::processEvent(EVENT::LCEvent* event) { 

    tuple->setVariableValue("event", event->getEventNumber());
    //std::cout << "Event: " << event->getEventNumber() << std::endl;

    // Get the collection of TargetHits from the event.  If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* target_hits
        = (EVENT::LCCollection*) event->getCollection("TargetHits");
    EVENT::CalorimeterHit* target_hit 
        = (EVENT::CalorimeterHit*) target_hits->getElementAt(0);

    // Get the energy deposited in the target
    tuple->setVariableValue("target_energy", target_hit->getEnergy()*1000);

    EVENT::LCCollection* trigger_pad_hits
        = (EVENT::LCCollection*) event->getCollection("TriggerPadHits");
    EVENT::CalorimeterHit* trigger_pad_hit 
        = (EVENT::CalorimeterHit*) trigger_pad_hits->getElementAt(0);
   
    // Get the energy deposited in the trigger pad
    tuple->setVariableValue("trigger_pad_energy", trigger_pad_hit->getEnergy()*1000); 

    // Get the collection of SimTrackerHits associated with the recoil tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* recoil_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");
    tuple->setVariableValue("n_recoil_hits", recoil_sim_hits->getNumberOfElements());

    // Create a cell ID decoder used to get specific properties associated with
    // Recoil SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> recoil_sim_hit_decoder(recoil_sim_hits);

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
        = (EVENT::LCCollection*) event->getCollection("MCParticle");

    //std::cout << "Number of MC particles: " << mc_particles->getNumberOfElements() << std::endl;

    int strategy[10] = { 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 }; 
    int n_tracks = 0; 
    
    // Loop over all of the MC particles in the event and find the photon leading
    // to the photonuclear reaction and the recoil electron.
    EVENT::MCParticle* pn_gamma = nullptr;
    EVENT::MCParticle* e_recoil = nullptr;
    std::vector<EVENT::SimTrackerHit*> findable_recoil_sim_hits;
    std::vector<EVENT::MCParticle*> pn_particles;  
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        // Get an MC particle from the collection of MC particles
        EVENT::MCParticle* mc_particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
       
        // If the MC particle is an electron and the electron has no parents, 
        // then the recoil electron has been found.
        if (mc_particle->getPDG() == 11 && mc_particle->getParents().size() == 0) { 

            e_recoil = mc_particle;

            // Check if the energy depositions of the recoil electron are 
            // consistent with a particle whose track should be found.
            tuple->setVariableValue("recoil_is_findable", 0); 
            if (TrackUtils::isTrackFindable(10, strategy, e_recoil,
                        recoil_sim_hits, findable_recoil_sim_hits)) { 
                //std::cout << "Particle (PDG ID = " << e_recoil->getPDG() << ") is findable." << std::endl;
                tuple->setVariableValue("recoil_is_findable", 1);
                n_tracks++; 
                continue;
            }
        }

        // If the MC particle is charged, check if its energy depositions are
        // consistent with a particle whose track should be found. 
        if (TrackUtils::isTrackFindable(10, strategy, mc_particle, recoil_sim_hits)) { 
            n_tracks++;
        }
        
        //std::cout << "MC Particle PDG ID: " 
        //          << mc_particle->getPDG() << " vertex_z : " 
        //          << mc_particle->getVertex()[2] << std::endl;
        
        // If the MC particle is a photon that was created in the target and has
        // the recoil electron as a parent, then the gamma that initiated the
        // photonuclear reaction has been found.
        if (mc_particle->getPDG() == 22 && createdWithinTarget(mc_particle)
                && mc_particle->getParents().size() == 1 
                && (mc_particle->getParents()[0]->getPDG() == 11 
                    && mc_particle->getParents()[0]->getParents().size() == 0)
                && mc_particle->getEnergy() >= .5) {
            
            //std::cout << "[ PhotoNuclearAnalysis ]: Photon that initiated "
            //          << " photo-nuclear reaction was found." << std::endl;

            pn_gamma = mc_particle;
            
            tuple->setVariableValue("pn_gamma_energy", mc_particle->getEnergy());
            tuple->setVariableValue("pn_particle_mult", mc_particle->getDaughters().size());

            // Loop over all of the daughters
            for (auto daughter : mc_particle->getDaughters()) {
                
                tuple->addToVector("pn_pdg_id", daughter->getPDG());
                pn_particles.push_back(daughter);
                
                // Calculate the scattering angle of the PN daughter
                double* pvec = (double*) daughter->getMomentum();
                double p = sqrt(pow(pvec[0], 2) + pow(pvec[1], 2) + pow(pvec[2], 2));
                double theta = acos(pvec[2]/p)*180/3.14159;
                tuple->addToVector("pn_theta", theta);
                    
                // If the particle is charged, check that it's findable
                /*int is_findable = 0;
                if (daughter->getCharge() != 0
                        && TrackUtils::isTrackFindable(10, strategy, daughter, recoil_sim_hits)) {    
                        //std::cout << "Particle (PDG ID = " 
                        //          << daughter->getPDG() << ") is findable."
                        //          << std::endl;
                    is_findable = 1;
                    n_tracks++;
                }
                tuple->addToVector("is_findable", is_findable);*/
            }
        }

        // If the MC particle is a hadron and it wasn't created within the 
        // target, skip the event.
        if (mc_particle->getPDG() > 100 && !createdWithinTarget(mc_particle)) {
            //std::cout << "[ PhotoNuclearAnalysis ]: Skipping event with "
            //          << "photo-nuclear reaction outside of the target." 
            //          << std::endl;
            tuple->clear();
            return;
        }
    }
    tuple->setVariableValue("n_tracks", n_tracks);

    int last_layer = 0;
    int first_layer = 15;
    EVENT::SimTrackerHit* last_recoil_sim_hit = nullptr;
    EVENT::SimTrackerHit* first_recoil_sim_hit = nullptr;
    for (EVENT::SimTrackerHit* recoil_sim_hit : findable_recoil_sim_hits) {
        
        // Get the layer number associated with this hit.
        int layer = recoil_sim_hit_decoder(recoil_sim_hit)["layer"];
        
        if (layer > last_layer) { 
            last_layer = layer;
            last_recoil_sim_hit = recoil_sim_hit;
        }

        if (layer < first_layer) { 
            first_layer = layer; 
            first_recoil_sim_hit = recoil_sim_hit;
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

    int n_recoil_pn_hits = 0;
    for (int sim_hit_n = 0; sim_hit_n < recoil_sim_hits->getNumberOfElements(); ++sim_hit_n) {
        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* recoil_sim_hit 
            = (EVENT::SimTrackerHit*) recoil_sim_hits->getElementAt(sim_hit_n); 
    
        // Get the layer number associated with this hit.
        int layer = recoil_sim_hit_decoder(recoil_sim_hit)["layer"];

        tuple->addToVector("recoil_sim_hit_layer", layer);
        tuple->addToVector("recoil_sim_hit_dedx",  recoil_sim_hit->getEDep());
        tuple->addToVector("recoil_sim_hit_time",  recoil_sim_hit->getTime());
        tuple->addToVector("recoil_sim_hit_pos_x", recoil_sim_hit->getPosition()[0]); 
        tuple->addToVector("recoil_sim_hit_pos_y", recoil_sim_hit->getPosition()[1]); 
        tuple->addToVector("recoil_sim_hit_pos_z", recoil_sim_hit->getPosition()[2]); 
    } 
    tuple->setVariableValue("n_recoil_pn_hits", n_recoil_pn_hits);
    /*



    
        if (recoil_sim_hit->getMCParticle() != e_recoil) n_recoil_pn_hits++;
    }
    */

    /*

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* tagger_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("TaggerTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> tagger_sim_hit_decoder(tagger_sim_hits);

    std::vector<int> tagger_sim_hits_vec(14, 0);
    for (int sim_hit_n = 0; sim_hit_n < tagger_sim_hits->getNumberOfElements(); ++sim_hit_n) { 
        
        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* tagger_sim_hit 
            = (EVENT::SimTrackerHit*) tagger_sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = tagger_sim_hit_decoder(tagger_sim_hit)["layer"];
        tagger_sim_hits_vec[layer - 1]++;

        tuple->addToVector("tagger_sim_hit_backscattered", 0);
        if (tagger_sim_hit->getMCParticle()->isBackscatter()) 
            tuple->addToVector("tagger_sim_hit_backscattered", 1);
        tuple->addToVector("tagger_sim_hit_layer", layer);
        tuple->addToVector("tagger_sim_hit_dedx",  tagger_sim_hit->getEDep());
        tuple->addToVector("tagger_sim_hit_time",  tagger_sim_hit->getTime());
        tuple->addToVector("tagger_sim_hit_pos_x", tagger_sim_hit->getPosition()[0]); 
        tuple->addToVector("tagger_sim_hit_pos_y", tagger_sim_hit->getPosition()[1]); 
        tuple->addToVector("tagger_sim_hit_pos_z", tagger_sim_hit->getPosition()[2]); 

    }

    
    */
    tuple->fill();
    tuple->clear(); 
}

void PhotoNuclearAnalysis::finalize() { 
    tuple->close();
}

bool PhotoNuclearAnalysis::createdWithinTarget(MCParticle* particle) {
    //std::cout << "Checking if particle was created within the target." << std::endl;
    if (particle->getVertex()[2] > -.1750 && particle->getVertex()[2] < .1750) {
        //std::cout << "Particle created within the target." << std::endl;
        return true;
    }
    //std::cout << "Particle created downstream of the target." << std::endl;
    return false; 
}
