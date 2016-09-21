/**
 *	@file RecoilTrackerAnalysis.h
 *	@brief Analysis used to study the performance of the Recoil tracker.
 *	@author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 27, 2016
 */

#ifndef __RECOIL_TRACKER_ANALYSIS_H__
#define __RECOIL_TRACKER_ANALYSIS_H__

//--------------//
//   Analyses   //
//--------------//
#include <LcioAbstractAnalysis.h>
#include <FlatTupleMaker.h>
#include <TrackUtils.h>

//----------//
//   LCIO   //
//----------//
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerRawData.h>
#include <EVENT/MCParticle.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCRelationNavigator.h>

class RecoilTrackerAnalysis : public LcioAbstractAnalysis { 

    public: 

        /** Default constructor */
        RecoilTrackerAnalysis();

        /** Destructor */
        ~RecoilTrackerAnalysis();

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

        /**
         * Enable/disable the filtering of photonuclear events
         *
         * @param filter_pn True to enable filtering, false otherwise
         */
        void filterPhotoNuclearEvents(bool filter_pn);

    private: 

        bool createdWithinTarget(MCParticle* particle);

        /** Allows the creation of a ROOT ntuple */
        FlatTupleMaker* tuple;

        /** */
        bool filter_pn;

}; // RecoilTrackerAnalysis

#endif // __RECOIL_TRACKER_ANALYSIS_H__
