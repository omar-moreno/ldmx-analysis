
#ifndef __TRACK_EXTRAPOLATOR_H__
#define __TRACK_EXTRAPOLATOR_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>
#include <cmath>

//----------//
//   LCIO   //
//----------//
#include <EVENT/Track.h>


#include <TrackUtils.h>

namespace TrackExtrapolator { 

        /**
         *
         */
        double getR(EVENT::Track* track);

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
        double getXc(EVENT::Track* track);

        /**
         *
         */
        double getYc(EVENT::Track* track);

        /**
         *
         */
        double getPathLength(EVENT::Track* track, double x1, double y1, double x2, double y2);

        /**
         *
         */
        double getPathToXPlane(EVENT::Track* track, double x);

        /**
         *
         */
        double getPhi(EVENT::Track* track, std::vector<double> position);

        /**
         *
         */
        // TODO: Move this function to a track utility namespace
        double getSinTheta(EVENT::Track* track); 

        /**
         *
         */
        // TODO: Move this function to a track utility namespace
        double getCosTheta(EVENT::Track* track); 

        /**
         *
         */
        std::vector<double> getPointOnHelix(EVENT::Track* track, double path_length);

        /**
         *
         */
        std::vector<double> extrapolateHelixToXPlane(EVENT::Track* track, double x);   

        /**
         *
         */
        std::vector<double> extrapolateTrack(EVENT::Track* track, double z);  

        /**
         * The end of the box dipole field during the engineering run.  Note 
         * that the dipole length was changed such that the integral of the 
         * box dipole field matches that of the field map integral.  The 
         * dipole edge is set to 45.72 cm (center of dipole) + 
         * 108/2 cm (dipole length).
         *
         */
        static const int DIPOLE_EDGE = 997.2;

}

#endif // __TRACK_EXTRAPOLATOR_H__
