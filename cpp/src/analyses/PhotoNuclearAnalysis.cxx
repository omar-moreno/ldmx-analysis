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

    tuple->addVariable("pn_gamma_energy");
    tuple->addVariable("pn_particle_mult");
    tuple->addVariable("n_recoil_hits");
    tuple->addVariable("n_neutron");
    tuple->addVariable("n_pi0");
    tuple->addVariable("n_pip");
    tuple->addVariable("n_pim");
    tuple->addVariable("n_proton");
    tuple->addVariable("n_ion");
    tuple->addVariable("n_gamma");

    tuple->addVector("recoil_sim_hit_layer");
    tuple->addVector("recoil_sim_hit_dedx");
    tuple->addVector("recoil_sim_hit_backscattered");
    tuple->addVector("recoil_sim_hit_time");
    tuple->addVector("recoil_sim_hit_pos_x");
    tuple->addVector("recoil_sim_hit_pos_y");
    tuple->addVector("recoil_sim_hit_pos_z");
}

void PhotoNuclearAnalysis::processEvent(EVENT::LCEvent* event) { 

    tuple->setVariableValue("event", event->getEventNumber());
    //std::cout << "Event: " << event->getEventNumber() << std::endl;

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
        = (EVENT::LCCollection*) event->getCollection("MCParticle");

    //std::cout << "Number of MC particles: " << mc_particles->getNumberOfElements() << std::endl;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    EVENT::MCParticle* pn_gamma = nullptr;
    double n_neutron = 0;
    double n_pi0 = 0; 
    double n_pip = 0;
    double n_pim = 0;
    double n_proton = 0;
    double n_ion = 0;
    double n_gamma = 0;

    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* mc_particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
        //std::cout << "MC Particle PDG ID: " << mc_particle->getPDG() << " vertex_z : " << mc_particle->getVertex()[2] << std::endl;
        if (mc_particle->getPDG() == 22 && createdWithinTarget(mc_particle)) {
            if (mc_particle->getParents()[0]->getParents().size() == 0) { 
    
                //std::cout << "Photonuclear gamma found." << std::endl;

                tuple->setVariableValue("pn_gamma_energy", mc_particle->getEnergy());
                tuple->setVariableValue("pn_particle_mult", mc_particle->getDaughters().size());

                pn_gamma = mc_particle;

                for (auto daughter : mc_particle->getDaughters()) { 
                    switch (daughter->getPDG()) { 
                        case 2112 : 
                            n_neutron++;
                            break;
                        case 2212 : 
                            n_proton++;
                            break;
                        case 111: 
                            n_pi0++;
                            break;
                        case 211: 
                            n_pip++;
                            break;
                        case -211: 
                            n_pim++;
                            break;
                        case 22: 
                            n_gamma; 
                            break;
                        default : 
                            n_ion++;
                            break;
                    }
                }
            }
        }


        tuple->setVariableValue("n_neutron", n_neutron);
        tuple->setVariableValue("n_pi0", n_pi0);
        tuple->setVariableValue("n_pip", n_pip);
        tuple->setVariableValue("n_pim", n_pim);
        tuple->setVariableValue("n_proton", n_proton);
        tuple->setVariableValue("n_ion", n_ion);
        tuple->setVariableValue("n_gamma", n_gamma);

        if (mc_particle->getPDG() > 100 && !createdWithinTarget(mc_particle)) {
            //std::cout << "Skipping event with photo-nuclear event downstream." << std::endl;
            return;
        }
    }

    // Get the collection of SimTrackerHits associated with the Tagger tracker
    // from the event.  If no such collection exist, a DataNotAvailableException
    // is thrown.
    EVENT::LCCollection* sim_hits 
        = (EVENT::LCCollection*) event->getCollection("RecoilTrackerHits");

    // Create a cell ID decoder used to get specific properties associated with
    // Tagger SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> sim_hit_decoder(sim_hits);
    tuple->setVariableValue("n_recoil_hits", sim_hits->getNumberOfElements());

    std::vector<int> sim_hits_vec(10, 0);
    for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 

        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];
        sim_hits_vec[layer - 1]++;

        tuple->addToVector("recoil_sim_hit_backscattered", 0);
        if (sim_hit->getMCParticle()->isBackscatter()) 
            tuple->addToVector("recoil_sim_hit_backscattered", 1);
        tuple->addToVector("recoil_sim_hit_layer", layer);
        tuple->addToVector("recoil_sim_hit_time",  sim_hit->getEDep());
        tuple->addToVector("recoil_sim_hit_time",  sim_hit->getTime());
        tuple->addToVector("recoil_sim_hit_pos_x", sim_hit->getPosition()[0]); 
        tuple->addToVector("recoil_sim_hit_pos_y", sim_hit->getPosition()[1]); 
        tuple->addToVector("recoil_sim_hit_pos_z", sim_hit->getPosition()[2]); 
    }

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
