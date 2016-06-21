/**
 *	@file LcioAbstractAnalysis.h
 *	@brief Abstract class describing an analysis used to process LCEvents.
 *	@author Omar Moreno <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 20, 2016
 *
 */

#ifndef __LCIO_ANALYSIS_H__
#define __LCIO_ANALYSIS_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <string>

//----------//
//   LCIO   //
//----------//
#include <EVENT/LCEvent.h>

class LcioAbstractAnalysis {

	public: 

		virtual ~LcioAbstractAnalysis() { }; 

		/**
		 *
		 */
		virtual void initialize() = 0; 

		/**
		 *
		 */
		virtual void processEvent(EVENT::LCEvent*) = 0;

		/**
		 *
		 */
		virtual void finalize() = 0; 

		/**
		 *
		 */
		std::string toString();

    protected:

       /** Class name */
        std::string class_name;  

}; // LcioAnalysis 

#endif // __LCIO_ANALYSIS_H__
