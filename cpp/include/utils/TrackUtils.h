/**
 *
 * @file TrackUtils.h
 * @author: <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *          SLAC National Accelerator Facility
 * @date: December 16, 2013
 *
 */

#ifndef __TRACK_UTILS_H__
#define __TRACK_UTILS_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <cmath>
#include <vector>

//----------//
//   LCIO   //
//----------//
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>

namespace TrackUtils { 

    /**
     *
     */
    double getX0(EVENT::Track* track);

    /**
     *
     */ 
    double getY0(EVENT::Track* track); 

    /**
     *
     */ 
    double getR(EVENT::Track* track); 

    /**
     *
     */ 
    double getDoca(EVENT::Track* track);

    /**
     *
     */ 
    double getPhi0(EVENT::Track* track);

    /**
     *
     */
   double getPhi(EVENT::Track*, std::vector<double>); 

    /**
     *
     */ 
    double getZ0(EVENT::Track* track); 

    /**
     *
     */ 
    double getTanLambda(EVENT::Track* track); 

    /**
     *
     */ 
    double getSinTheta(EVENT::Track* track); 

    /**
     *
     */ 
    double getCosTheta(EVENT::Track* track);

    /**
     *
     */
    double getXc(EVENT::Track* track); 

    /**
     *
     */
    double getYc(EVENT::Track* track); 

	/**
	 *
	 */
    template <typename T> int signum(T val){
        return (T(0) < val) - (val < T(0));
    }

	/**
	 *
	 */
	std::vector<double> getMomentumVector(EVENT::Track*, double);

	/**
	 *
	 */
	double getMomentum(EVENT::Track*, double);

	/**
	 *
	 */
	int getCharge(EVENT::Track*);

	/**
	 *
	 */
	int getLayer(EVENT::TrackerHit*);

}

#endif // __TRACK_UTILS_H__
