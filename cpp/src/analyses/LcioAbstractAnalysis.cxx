/**
 *	@file LcioAbstractAnalysis.cxx
 *	@brief Abstract class describing an analysis used to process LCEvents.
 *	@author Omar Moreno <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 20, 2016
 *
 */

#include <LcioAbstractAnalysis.h>


std::string LcioAbstractAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}
