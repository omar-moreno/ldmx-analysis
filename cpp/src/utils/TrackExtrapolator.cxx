
#include <TrackExtrapolator.h>


double TrackExtrapolator::getPathLength(EVENT::Track* track, double x1, double y1, double x2, double y2) { 
    double phi1 = atan2(y1 - TrackUtils::getYc(track), x1 - TrackUtils::getXc(track));
    double phi2 = atan2(y2 - TrackUtils::getYc(track), x2 - TrackUtils::getXc(track));
    double dphi = phi2 - phi1;

    if (dphi > 3.14159) dphi -= 2*3.14159;
    else if (dphi < -3.14159) dphi += 2*3.14159;

    return -TrackUtils::getR(track)*dphi;
}

double TrackExtrapolator::getPathToXPlane(EVENT::Track* track, double x) {

    double r = TrackUtils::getR(track);
    double y = TrackUtils::getYc(track) 
        + TrackUtils::signum<double>(r)*sqrt(r*r - pow(x - TrackUtils::getXc(track), 2));

    return getPathLength(track, TrackUtils::getX0(track), TrackUtils::getY0(track), x, y); 
}

std::vector<double> TrackExtrapolator::getPointOnHelix(EVENT::Track* track, double path_length) {

    double phi = track->getPhi() - (path_length/TrackUtils::getR(track));
    double x = TrackUtils::getXc(track) - TrackUtils::getR(track)*sin(phi); 
    double y = TrackUtils::getYc(track) + TrackUtils::getR(track)*cos(phi);
    double z = track->getZ0() + path_length*track->getTanLambda();

    std::vector<double> position(3, 0);
    position[0] = x; 
    position[1] = y;
    position[2] = z;

    return position;  
}

std::vector<double> TrackExtrapolator::extrapolateHelixToXPlane(EVENT::Track* track, double x) { 

    double path_length = getPathToXPlane(track, x);

    return getPointOnHelix(track, path_length);
}

std::vector<double> TrackExtrapolator::extrapolateTrack(EVENT::Track* track, double z) { 

    std::vector<double> position(3,0); 
    double dz = 0;

    /*if (z >= TrackExtrapolator::DIPOLE_EDGE) {

        // If the point of extrapolation is outside of the dipole edge, then 
        // extrapolate the helix to the edge and then use a straight line 
        // extrapolation beyond that

        // Extrapolate the helix to the edge of the dipole 
        position = extrapolateHelixToXPlane(track, TrackExtrapolator::DIPOLE_EDGE);

        // Get the difference between the dipole edge and the extrapolation
        // point. The track will be extrapolated assuming no field for this
        // distance i.e. straight line extrapolation
        dz = z - TrackExtrapolator::DIPOLE_EDGE;
    } else if (z <= 0) {

        // If the extrapolation point is upstream of the target, do something
        // similar as above

        position = extrapolateHelixToXPlane(track, 0);
        dz = z - position[0];  
    
    } else { */
    
        // If the extrapolation point is inside of the field region, 
        // analytically extrapolate the helix and return the position
        position = extrapolateHelixToXPlane(track, z); 
   
        // FIXME: This position should be in the JLab frame     
        return position;
    /*}

    // Calculate the value of Phi at the track position
    double phi = TrackUtils::getPhi(track, position);
    
    // Calcualte the distance to the extrapolation point
    double r = dz/TrackUtils::getSinTheta(track)*cos(phi);

    // Get the delta x and y values at the point of extrapolation 
    double dx = r*TrackUtils::getSinTheta(track)*sin(phi);
    double dy = r*TrackUtils::getCosTheta(track);

    // Calculate the position of the track at the extrapolation point
    std::vector<double> extrapolated_position(3, 0);
    extrapolated_position[0] = position[1] + dx;
    extrapolated_position[1] = position[2] + dy;
    extrapolated_position[2] = z; 

    return extrapolated_position; */
}
