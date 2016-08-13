
#include <SignalAnalysis.h>

int SignalAnalysis::A_PRIME_PDG = 622; 

int SignalAnalysis::ELECTRON_PDG = 11;

int SignalAnalysis::FINAL_STATE = 1;

SignalAnalysis::SignalAnalysis()
    : tuple(new FlatTupleMaker("signal_analysis_ntuple.root", "results")) { 
    
        LcioAbstractAnalysis::class_name = "SignalAnalysis";     
}

SignalAnalysis::~SignalAnalysis() { 
}

void SignalAnalysis::initialize() { 

    tuple->addVariable("ap_mass");

    tuple->addVariable("event");

    tuple->addVariable("is_within_acceptance");

    tuple->addVariable("recoil_p");
    tuple->addVariable("recoil_pt");
    tuple->addVariable("recoil_px");
    tuple->addVariable("recoil_py");
    tuple->addVariable("recoil_pz");
    
    tuple->addVector("sim_hit_layer");
    tuple->addVector("sim_hit_pos_x");
    tuple->addVector("sim_hit_pos_y");
    tuple->addVector("sim_hit_pos_z");
    tuple->addVector("sim_hit_time");
}

void SignalAnalysis::processEvent(EVENT::LCEvent* event) { 

    tuple->setVariableValue("event", event->getEventNumber());
    //std::cout << "[ SignalAnalysis ]: Event number: " << event->getEventNumber() << std::endl;

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
       = (EVENT::LCCollection*) event->getCollection("MCParticle");

    if (mc_particles->getNumberOfElements() == 0) return;
    
    // Recoil electron
    EVENT::MCParticle* recoil_electron = nullptr;
    
    // Find the recoil electron
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) { 
        
        EVENT::MCParticle* mc_particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);

        if (mc_particle->getPDG() == SignalAnalysis::ELECTRON_PDG
                && mc_particle->getGeneratorStatus() == SignalAnalysis::FINAL_STATE) { 
            recoil_electron =  mc_particle;
        }

        if (mc_particle->getPDG() == SignalAnalysis::A_PRIME_PDG) { 
            tuple->setVariableValue("ap_mass", mc_particle->getMass()); 
        }
    } 

    // Calculate the momentum of the recoil electron
    double* recoil_e_pvec = (double*) recoil_electron->getMomentum();
    double recoil_e_p = sqrt(pow(recoil_e_pvec[0], 2) + pow(recoil_e_pvec[1], 2) + pow(recoil_e_pvec[2], 2));
    double recoil_e_pt = sqrt(pow(recoil_e_pvec[0], 2) + pow(recoil_e_pvec[1], 2));
    tuple->setVariableValue("recoil_p", recoil_e_p);
    tuple->setVariableValue("recoil_pt", recoil_e_pt);
    tuple->setVariableValue("recoil_px", recoil_e_pvec[0]);
    tuple->setVariableValue("recoil_py", recoil_e_pvec[1]);
    tuple->setVariableValue("recoil_pz", recoil_e_pvec[2]);

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> sim_hit_decoder(sim_hits);
   
    std::vector<int> recoil_e_hit_count(10, 0); 
    for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 
        
        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit 
            = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];
        //std::cout << "Layer: " << layer << std::endl;

        if (sim_hit->getMCParticle() == recoil_electron) {
            recoil_e_hit_count[layer - 1]++;
            tuple->addToVector("sim_hit_layer", layer);
            tuple->addToVector("sim_hit_pos_x", sim_hit->getPosition()[0]); 
            tuple->addToVector("sim_hit_pos_y", sim_hit->getPosition()[1]); 
            tuple->addToVector("sim_hit_pos_z", sim_hit->getPosition()[2]); 
            tuple->addToVector("sim_hit_time", sim_hit->getTime());
        }
    }
    
    tuple->setVariableValue("is_within_acceptance", 0);
    bool is_within_acceptance = true; 
    for (int layer_n = 0; layer_n < 6; ++layer_n) { 
        if (recoil_e_hit_count[layer_n] == 0) is_within_acceptance = false;
    }

    if (is_within_acceptance && 
            (recoil_e_hit_count[6]*recoil_e_hit_count[7] > 0 
             || recoil_e_hit_count[8] > 0 
             || recoil_e_hit_count[9] > 0)) {
        
        tuple->setVariableValue("is_within_acceptance", 1);
    } 
    tuple->fill();  
}

void SignalAnalysis::finalize() { 
     
    tuple->close(); 
}

