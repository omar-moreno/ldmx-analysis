/**
 * 
 */

#ifndef __CLUSTER_ANALYSIS_H__
#define __CLUSTER_ANALYSIS_H__

//--------------//
//   Analyses   //
//--------------//
#include <LcioAbstractAnalysis.h>

class ClusterAnalysis : public LcioAbstractAnalysis { 

    public: 

        /** Default Ctor */    
        ClusterAnalysis(); 

        /** Destructor */
        ~ClusterAnalysis(); 

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


};

#endif
