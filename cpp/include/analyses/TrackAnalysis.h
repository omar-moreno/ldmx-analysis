
#ifndef __TRACK_ANALYSIS_H__
#define __TRACK_ANALYSIS_H__

//--------------//
//   Analyses   //
//--------------//
#include <LcioAbstractAnalysis.h>

//-----------//
//   Utils   //
//-----------//
#include <Plotter.h>
#include <TrackUtils.h>

//----------//
//   LCIO   //
//----------//
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerRawData.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

class TrackAnalysis : public LcioAbstractAnalysis { 

    public: 

        /** Default Ctor */    
        TrackAnalysis(); 

        /** Destructor */
        ~TrackAnalysis(); 

        /** 
         * Method used to initialize an HPS analysis. This method is called
         * once at the beginning of an analysis. 
         */
        void initialize();

        /**
         *  Method used to process an HpsEvent.
         *
         *  @param event : HpsEvent that will be processed
         */
        void processEvent(EVENT::LCEvent* event);

        /** 
         * Method used to finalize an HPS analysis. This method is called once
         * at the end of an analysis.
         */
        void finalize();

        /**
         *  Provide a string representation of this analysis.
         *
         *  @return String representation of this analysis.
         */
        std::string toString();

    private:

        /** Name of the class */
        std::string class_name;
    
        /** Plotter used to make ROOT histograms */
        Plotter* plotter; 

        double findable_tracks; 
        double found_tracks;  
};

#endif
