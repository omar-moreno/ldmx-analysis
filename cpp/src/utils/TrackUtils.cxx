/**
 *	@section purpose:
 *  @author: Omar Moreno <omoreno1@ucsc.edu>
 *			 Santa Cruz Institute for Particle Physics
 *			 University of California, Santa Cruz
 *  @date: December 16, 2013
 *  @version: 1.0
 *
 */

#include <TrackUtils.h>

namespace { 
		const double param = 2.99792458e-04; 	
}

namespace TrackUtils { 

    
    double getX0(EVENT::Track* track){
        return -1*getDoca(track)*sin(getPhi0(track));  
    }

    double getY0(EVENT::Track* track){
        return getDoca(track)*cos(getPhi0(track)); 
    };

    double getR(EVENT::Track* track){
        return 1.0/track->getOmega(); 
    }; 
    
    double getDoca(EVENT::Track* track){
        return track->getD0();     
    };
    
    double getPhi0(EVENT::Track* track){
        return track->getPhi(); 
    };

    double getPhi(EVENT::Track* track, std::vector<double> position){ 
          double x = sin(getPhi0(track)) - (1/getR(track))*(position[0] - getX0(track)); 
          double y = cos(getPhi0(track)) + (1/getR(track))*(position[1] - getY0(track)); 
    
        return atan2(x, y); 
    }; 
    
    double getZ0(EVENT::Track* track){
        return track->getZ0(); 
    
    }; 
    
    double getTanLambda(EVENT::Track* track){
        return track->getTanLambda(); 
    }; 
    
    double getSinTheta(EVENT::Track* track){
       return 1/sqrt(1 + pow(getTanLambda(track), 2));  
    }; 
    
    double getCosTheta(EVENT::Track* track){
        return getTanLambda(track)/sqrt(1 + pow(getTanLambda(track), 2)); 
    }; 

    double getXc(EVENT::Track* track){ 
        return (getR(track) - getDoca(track))*sin(getPhi0(track)); 
    };

    double getYc(EVENT::Track* track){
        return -(getR(track) - getDoca(track))*cos(getPhi0(track));   
    };

	std::vector<double> getMomentumVector(EVENT::Track* track, double b_field){
		std::vector<double> p(3,0); 
		double pt = std::abs(getR(track)*b_field*param);
		
		p[0] = pt*cos(getPhi0(track)); 
		p[1] = pt*sin(getPhi0(track)); 
		p[2] = pt*getTanLambda(track); 		
		
		return p; 	
	};

	double getMomentum(EVENT::Track* track, double b_field){
	
		std::vector<double> p_vector = getMomentumVector(track, b_field); 
		double p_sum = 0;
        
        for (double p : p_vector) { 
            p_sum += p*p;  
        } 

		return sqrt(p_sum); 
	};

	int getCharge(EVENT::Track* track){
		int charge; 
		track->getOmega() > 0 ? charge = 1 : charge = -1; 
		return charge; 		
	};
}

bool TrackUtils::isTrackFindable(int layers, int strategy[],  EVENT::MCParticle* particle, EVENT::LCCollection* sim_hits) {
    std::vector<EVENT::SimTrackerHit*> empty_vec; 
    return TrackUtils::isTrackFindable(layers, strategy, particle, sim_hits, empty_vec); 
}

bool TrackUtils::isTrackFindable(int layers, int strategy[], EVENT::MCParticle* particle, EVENT::LCCollection* sim_hits, 
       std::vector<EVENT::SimTrackerHit*> &findable_sim_hits) {

    // Create a cell ID decoder used to get specific properties associated with
    // SimTrackerHits i.e. layer.
    UTIL::CellIDDecoder<EVENT::SimTrackerHit> sim_hit_decoder(sim_hits);

    // Loop over all Tagger SimTrackerHits in the event.
    std::vector<int> hits(layers, 0);
    for (int sim_hit_n = 0; sim_hit_n < sim_hits->getNumberOfElements(); ++sim_hit_n) { 
        
        // Get a Tagger SimTrackerHit from the collection of hits.
        EVENT::SimTrackerHit* sim_hit = (EVENT::SimTrackerHit*) sim_hits->getElementAt(sim_hit_n);

        // Get the layer number associated with this hit.
        int layer = sim_hit_decoder(sim_hit)["layer"];

        if (particle == sim_hit->getMCParticle())  {
            findable_sim_hits.push_back(sim_hit); 
            hits[layer - 1]++;
        } 
    }

    bool findable_track = true;
    for (int layer_n = 0; layer_n < layers; ++layer_n) { 
        if (hits[layer_n] == 0 && strategy[layer_n] == 1) {
            findable_sim_hits.clear(); 
            findable_track = false;
        }
    }
    
    return findable_track;

}

