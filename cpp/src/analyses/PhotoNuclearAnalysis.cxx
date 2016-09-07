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

    /*TH1* plot = nullptr; 
      plot = plotter->build2DHistogram("Simulated Recoil Hits per Layer", 10, 1, 11, 10, 0, 10); 
      plot->GetXaxis()->SetTitle("Layer number"); 
      plot->GetYaxis()->SetTitle("Total hits");*/

    tuple->addVariable("event");

    tuple->addVariable("recoil_is_findable");
    tuple->addVariable("recoil_energy");

    tuple->addVariable("pn_gamma_energy");
    tuple->addVariable("pn_particle_mult");

    tuple->addVariable("n_recoil_hits");
    tuple->addVariable("n_recoil_pn_hits");
    tuple->addVariable("n_recoil_tracks");

    tuple->addVector("pn_pdg_id");
    tuple->addVector("pn_theta");
    tuple->addVector("pn_phi");
    tuple->addVector("is_findable");

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
    std::cout << "Event: " << event->getEventNumber() << std::endl;

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* recoil_sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");


    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
        = (EVENT::LCCollection*) event->getCollection("MCParticle");

    //std::cout << "Number of MC particles: " << mc_particles->getNumberOfElements() << std::endl;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    EVENT::MCParticle* pn_gamma = nullptr;
    EVENT::MCParticle* e_recoil = nullptr;
    double n_neutron = 0;
    double n_pi0 = 0; 
    double n_pip = 0;
    double n_pim = 0;
    double n_proton = 0;
    double n_ion = 0;
    double n_gamma = 0;

    int first_strategy[10] = { 1, 1, 1, 1, 1, 1, 0, 0, 1, 0 }; 
    int second_strategy[10] = { 1, 1, 1, 1, 1, 1, 0, 0, 0, 1 };
    int total_recoil_tracks = 0; 
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* mc_particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
       
        if (mc_particle->getPDG() == 11 && mc_particle->getParents().size() == 0) { 
            
            e_recoil = mc_particle;

            tuple->setVariableValue("recoil_is_findable", 0); 
            if (TrackUtils::isTrackFindable(10, first_strategy, mc_particle, recoil_sim_hits)
                    || TrackUtils::isTrackFindable(10, second_strategy, mc_particle, recoil_sim_hits)) {
                //std::cout << "Particle (PDG ID = " << daughter->getPDG() << ") is findable." << std::endl;
                tuple->setVariableValue("recoil_is_findable", 1);
                total_recoil_tracks++; 
            } 
        }
        
        std::cout << "MC Particle PDG ID: " << mc_particle->getPDG() << " vertex_z : " << mc_particle->getVertex()[2] << std::endl;
        if (mc_particle->getPDG() == 22 && createdWithinTarget(mc_particle)) {
            
            std::cout << "Found particle created within target." << std::endl;
            if (mc_particle->getParents().size() == 0) std::cout << "null" << std::endl;
            if (mc_particle->getParents()[0]->getParents().size() == 0) { 
    
                std::cout << "Photonuclear gamma found." << std::endl;

                tuple->setVariableValue("pn_gamma_energy", mc_particle->getEnergy());
                tuple->setVariableValue("pn_particle_mult", mc_particle->getDaughters().size());

                pn_gamma = mc_particle;

                /*
                for (auto daughter : mc_particle->getDaughters()) {

                    // Check if the track is findable
                    tuple->addToVector("pn_pdg_id", daughter->getPDG());

                    double* pvec = (double*) daughter->getMomentum();
                    double p = sqrt(pow(pvec[0], 2) + pow(pvec[1], 2) + pow(pvec[2], 2));
                    double theta = acos(pvec[2]/p)*180/3.14159;
                    tuple->addToVector("pn_theta", theta);

                    int is_findable = 0;
                    if (daughter->getCharge() != 0) {
                        if (TrackUtils::isTrackFindable(10, first_strategy, daughter, recoil_sim_hits)
                               || TrackUtils::isTrackFindable(10, second_strategy, daughter, recoil_sim_hits)) {
                            //std::cout << "Particle (PDG ID = " << daughter->getPDG() << ") is findable." << std::endl;
                            is_findable = 1;
                            total_recoil_tracks++;
                        } 
                    } 
                    tuple->addToVector("is_findable", is_findable);
                }*/
            }
        }

        /*
        tuple->setVariableValue("n_neutron", n_neutron);
        tuple->setVariableValue("n_pi0", n_pi0);
        tuple->setVariableValue("n_pip", n_pip);
        tuple->setVariableValue("n_pim", n_pim);
        tuple->setVariableValue("n_proton", n_proton);
        tuple->setVariableValue("n_ion", n_ion);
        tuple->setVariableValue("n_gamma", n_gamma);

        if (mc_particle->getPDG() > 100 && !createdWithinTarget(mc_particle)) {
            //std::cout << "Skipping event with photo-nuclear event downstream." << std::endl;
            tuple->fill();
            return;
        }*/
    }
    
    /*
    tuple->setVariableValue("n_recoil_tracks", total_recoil_tracks);

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

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> recoil_sim_hit_decoder(recoil_sim_hits);
    tuple->setVariableValue("n_recoil_hits", recoil_sim_hits->getNumberOfElements());

    std::vector<int> recoil_sim_hits_vec(10, 0);
    int n_recoil_pn_hits = 0;
    for (int sim_hit_n = 0; sim_hit_n < recoil_sim_hits->getNumberOfElements(); ++sim_hit_n) { 

        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* recoil_sim_hit = (EVENT::SimTrackerHit*) recoil_sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = recoil_sim_hit_decoder(recoil_sim_hit)["layer"];
        recoil_sim_hits_vec[layer - 1]++;

        tuple->addToVector("recoil_sim_hit_backscattered", 0);
        if (recoil_sim_hit->getMCParticle()->isBackscatter()) 
            tuple->addToVector("recoil_sim_hit_backscattered", 1);
        tuple->addToVector("recoil_sim_hit_layer", layer);
        tuple->addToVector("recoil_sim_hit_dedx",  recoil_sim_hit->getEDep());
        tuple->addToVector("recoil_sim_hit_time",  recoil_sim_hit->getTime());
        tuple->addToVector("recoil_sim_hit_pos_x", recoil_sim_hit->getPosition()[0]); 
        tuple->addToVector("recoil_sim_hit_pos_y", recoil_sim_hit->getPosition()[1]); 
        tuple->addToVector("recoil_sim_hit_pos_z", recoil_sim_hit->getPosition()[2]); 
    
        if (recoil_sim_hit->getMCParticle() != e_recoil) n_recoil_pn_hits++;
    }
    tuple->setVariableValue("n_recoil_pn_hits", n_recoil_pn_hits);
    
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
