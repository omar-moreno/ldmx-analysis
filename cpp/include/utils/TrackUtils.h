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
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>

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


    /**
     *
     */
    bool isTrackFindable(int layers, int strategy[], EVENT::MCParticle* particle, EVENT::LCCollection* sim_hit);

    bool isTrackFindable(int layers, int strategy[], EVENT::MCParticle* particle, EVENT::LCCollection* sim_hit,
            std::vector<EVENT::SimTrackerHit*> &findable_sim_hits);
}

#endif // __TRACK_UTILS_H__
