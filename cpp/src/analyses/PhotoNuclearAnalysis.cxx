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

    tuple->addVector("recoil_sim_hit_layer");
    tuple->addVector("recoil_sim_hit_pos_x");
    tuple->addVector("recoil_sim_hit_pos_y");
    tuple->addVector("recoil_sim_hit_pos_z");
    tuple->addVector("recoil_sim_hit_time");
}

void PhotoNuclearAnalysis::processEvent(EVENT::LCEvent* event) { 

    tuple->setVariableValue("event", event->getEventNumber());

    // Get the collection of MC particles from the event. If no such collection
    // exist, a DataNotAvailableException is thrown.
    EVENT::LCCollection* mc_particles 
        = (EVENT::LCCollection*) event->getCollection("MCParticle");

    //std::cout << "Number of MC particles: " << mc_particles->getNumberOfElements() << std::endl;

    // Loop over all of the MC particles and find the beam electron. For now, 
    // this is simply done by looking for an electron that has no parent. 
    EVENT::MCParticle* pn_gamma = nullptr;
    for (int mc_particle_n = 0; mc_particle_n < mc_particles->getNumberOfElements(); ++mc_particle_n) {

        EVENT::MCParticle* mc_particle = (EVENT::MCParticle*) mc_particles->getElementAt(mc_particle_n);
        if (mc_particle->getPDG() == 22 
                && mc_particle->getParents().size() == 1
                && mc_particle->getParents()[0]->getParents().size() == 0) { 
            
            tuple->setVariableValue("pn_gamma_energy", mc_particle->getEnergy());
            tuple->setVariableValue("pn_particle_mult", mc_particle->getDaughters().size());

            pn_gamma = mc_particle;
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

    std::vector<int> sim_hits_vec(10, 0);
    for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 

        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n); 

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];
        sim_hits_vec[layer - 1]++;

        tuple->addToVector("recoil_sim_hit_layer", layer);
        tuple->addToVector("recoil_sim_hit_pos_x", sim_hit->getPosition()[0]); 
        tuple->addToVector("recoil_sim_hit_pos_y", sim_hit->getPosition()[1]); 
        tuple->addToVector("recoil_sim_hit_pos_z", sim_hit->getPosition()[2]); 
        tuple->addToVector("recoil_sim_hit_time", sim_hit->getTime());
    }
    
    tuple->fill();
}

void PhotoNuclearAnalysis::finalize() { 
    tuple->close();
}
