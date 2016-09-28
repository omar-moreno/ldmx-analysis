/**
 *	@file SignalAnalysis.h
 *	@brief Analysis used to study the performance of the Tagger tracker.
 *	@author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 *	        SLAC National Accelerator Facility
 *  @date June 27, 2016
 */

#ifndef __SIGNAL_ANALYSIS_H__
#define __SIGNAL_ANALYSIS_H__

//--------------//
//   Analyses   //
//--------------//
#include <FlatTupleMaker.h>
#include <LcioAbstractAnalysis.h>
#include <TrackUtils.h>
#include <TrackExtrapolator.h>

//----------//
//   LCIO   //
//----------//
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>

class SignalAnalysis : public LcioAbstractAnalysis { 

    public: 

        /** Defautl constructor */
        SignalAnalysis(); 

        /** Destructor */
        ~SignalAnalysis();

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

        /** A' PDG ID */
        static int A_PRIME_PDG;

        /** Electron PDG ID */
        static int ELECTRON_PDG; 

        /** Final state status value */
        static int FINAL_STATE; 

        int event_number; 
}; // SignalAnalysis

#endif // __SIGNAL_ANALYSIS_H__
