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
