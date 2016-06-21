
#ifndef __PLOTTER_H__
#define __PLOTTER_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <string>
#include <map>
#include <stdexcept>
#include <iostream>
#include <typeinfo>
#include <cstdlib>

//------------//
//--- ROOT ---//
//------------//
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TCanvas.h>

class Plotter { 

    public: 

        /**
         * Constructor
         */
        Plotter();

        /**
         *
         */
        ~Plotter(); 

        /**
         *
         */
        Plotter* setType(std::string hist_type) { this->hist_type = hist_type; return this; };
        
        /**
         *
         */
        Plotter* setGraphType(std::string graph_type) { this->graph_type = graph_type; return this; };

        /**
         *
         */
        void add1DHistogram(TH1* histogram); 

        /**
         *
         */
        Plotter* setLineColor(int color) { this->color = color; return this; };

        /**
         *
         */
        TH1* build1DHistogram(std::string name, int n_bins, double x_min, double x_max); 
       
        /**
         *
         */ 
        TH2* build2DHistogram(std::string name, int n_bins_x, double x_min, double x_max,
                int n_bins_y, double y_min, double y_max);
    
        /**
         *
         */
        TGraph* buildGraph(std::string name);

        /**
         *
         */
        TH1* get1DHistogram(std::string name);
        
        /**
         *
         */
        TH2* get2DHistogram(std::string name);

        /**
         *
         */
        TGraph* getGraph(std::string name); 

        /** */
        bool has1DHistogram(std::string name); 
        
        /** */
        bool has2DHistogram(std::string name); 
        
        /** */
        bool hasGraph(std::string name); 

        /**
         *
         */
        void saveToRootFile(std::string file_name); 
        
        /**
         *
         */
        void saveToPdf(std::string file_name); 



    private: 

        std::map<std::string, TH1*> histogram1D_map;
        std::map<std::string, TH2*> histogram2D_map;
        std::map<std::string, TGraph*> graph_map; 

        std::string hist_type;
        std::string graph_type; 

        int color; 

}; // Plotter

#endif // __PLOTTER_H__
