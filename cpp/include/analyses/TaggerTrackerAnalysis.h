/**
 *	@file TaggerTrackerAnalysis.h
 *	@brief Analysis used to study the performance of the Tagger tracker.
 *	@author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 27, 2016
 */

#ifndef __TAGGER_TRACKER_ANALYSIS_H__
#define __TAGGER_TRACKER_ANALYSIS_H__

//--------------//
//   Analyses   //
//--------------//
#include <LcioAbstractAnalysis.h>
#include <FlatTupleMaker.h>
#include <TrackUtils.h>

//----------//
//   LCIO   //
//----------//
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

class TaggerTrackerAnalysis : public LcioAbstractAnalysis { 

    public: 

        /** Default constructor */
        TaggerTrackerAnalysis();

        /** Destructor */
        ~TaggerTrackerAnalysis();

		/**
		 * Method used to initialize an analysis.  This method is called once
         * before any events are processed. 
		 */
		void initialize(); 

        /**
         *  Method used to process an event (LCEvent).  This method is called 
         *  once for every event. 
         *
         *  @param event The event to process.
         */
        void processEvent(EVENT::LCEvent* event);

        /** 
         * Method called at the end of an analysis.  This method is called once
         * per run.
         */
        void finalize();

    private: 

        /** Allows the creation of a ROOT ntuple */
        FlatTupleMaker* tuple;

        /** Findable track count */
        int findable_track;

}; // TaggerTrackerAnalysis

#endif // __TAGGER_TRACKER_ANALYSIS_H__
