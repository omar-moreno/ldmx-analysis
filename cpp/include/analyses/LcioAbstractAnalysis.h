/**
 *	@file LcioAbstractAnalysis.h
 *	@brief Abstract class describing an analysis used to process LCEvents.
 *	@author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
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
		 * Method used to initialize an analysis.  This method is called once
         * before any events are processed. 
		 */
		virtual void initialize() = 0; 

        /**
         *  Method used to process an event (LCEvent).  This method is called 
         *  once for every event. 
         *
         *  @param event The event to process.
         */
		virtual void processEvent(EVENT::LCEvent*) = 0;

        /** 
         * Method called at the end of an analysis.  This method is called once
         * per run.
         */
		virtual void finalize() = 0; 

        /**
         *  Provide a string representation of this analysis.  Most of the time
         *  this is just the name of the analysis class.
         *
         *  @return String representation of this analysis.
         */
		std::string toString();

    protected:

       /** Class name */
        std::string class_name;  

}; // LcioAnalysis 

#endif // __LCIO_ANALYSIS_H__
