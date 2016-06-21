
#include <Plotter.h>

Plotter::Plotter()
    : hist_type("float"),
      graph_type("graph"), 
      color(kAzure+3) {
}

Plotter::~Plotter() {
}

TH1* Plotter::build1DHistogram(std::string name, int n_bins, double x_min, double x_max) { 
    
    if (histogram1D_map[name] != NULL) { 
        throw std::runtime_error("The histogram " + name + " already exist!");
    }
    
    std::string root_name = name + "_" + std::to_string(rand()%10000);
    
    if (hist_type.compare("float") == 0) {
        histogram1D_map[name] = new TH1F(root_name.c_str(), name.c_str(), n_bins, x_min, x_max);
    } else if (hist_type.compare("double") == 0) { 
        histogram1D_map[name] = new TH1D(root_name.c_str(), name.c_str(), n_bins, x_min, x_max);
    }

    histogram1D_map[name]->SetLineColor(color);
    histogram1D_map[name]->SetLineWidth(2);
    histogram1D_map[name]->SetFillStyle(3003);
    histogram1D_map[name]->SetFillColor(color - 1);
    
    return histogram1D_map[name];
}


TH2* Plotter::build2DHistogram(std::string name, int n_bins_x, double x_min, double x_max, 
        int n_bins_y, double y_min, double y_max) { 
    
    if (histogram2D_map[name] != NULL) { 
        throw std::runtime_error("The histogram " + name + " already exist!");
    }

    std::string root_name = name + "_" + std::to_string(rand()%10000);

    if (hist_type.compare("float") == 0) {
        histogram2D_map[name] = new TH2F(root_name.c_str(), name.c_str(), n_bins_x, x_min, x_max,
                n_bins_y, y_min, y_max);
    } else if (hist_type.compare("double") == 0) { 
        histogram2D_map[name] = new TH2D(root_name.c_str(), name.c_str(), n_bins_x, x_min, x_max,
                n_bins_y, y_min, y_max);
    }
    
    histogram2D_map[name]->SetLineColor(color);
    histogram2D_map[name]->SetStats(0); 

    return histogram2D_map[name];
}

TGraph* Plotter::buildGraph(std::string name) { 
    
    if (graph_map[name] != NULL) { 
        throw std::runtime_error("The graph " + name + " already exist!"); 
    }

    std::string root_name = name + "_" + std::to_string(rand()%10000);

    if (graph_type.compare("graph") == 0) {
        graph_map[name] = new TGraph();
    } else if (graph_type.compare("errors") == 0) {
       graph_map[name] = new TGraphErrors();  
    } else if (graph_type.compare("asymm") == 0) { 
        graph_map[name] = new TGraphAsymmErrors();
    } else { 
        throw std::runtime_error("Graph of type " + graph_type + " cannot be created.");
    }

    graph_map[name]->SetNameTitle(root_name.c_str(), name.c_str()); 
    return graph_map[name]; 
}


TH1* Plotter::get1DHistogram(std::string name) { 
    
    if (histogram1D_map[name] == NULL) { 
        throw std::runtime_error("Histogram " + name + " has not been created.");
    }

    return histogram1D_map[name];
}

TH2* Plotter::get2DHistogram(std::string name) { 

    if (histogram2D_map[name] == NULL) { 
        throw std::runtime_error("Histogram " + name + " has not been created.");
    }

    return histogram2D_map[name]; 
}

TGraph* Plotter::getGraph(std::string name) { 

    if (graph_map[name] == NULL) { 
        throw std::runtime_error("Graph " + name + " has not be created."); 
    }
    
    return graph_map[name]; 
}


bool Plotter::has1DHistogram(std::string name) { 
    if (histogram1D_map[name] == NULL) return false; 

    return true;
}

bool Plotter::has2DHistogram(std::string name) { 
    if (histogram2D_map[name] == NULL) return false; 
    return true;
}

bool Plotter::hasGraph(std::string name) { 
    if (graph_map[name] == NULL) return false; 
    return true;
}

void Plotter::add1DHistogram(TH1* histogram) { 
    histogram1D_map[histogram->GetName()] = histogram; 
}

void Plotter::saveToRootFile(std::string file_name) { 

    TFile* output_file = new TFile(file_name.c_str(), "RECREATE");
    std::map<std::string, TH1*>::iterator histogram1D_it =  histogram1D_map.begin();
    for (histogram1D_it; histogram1D_it != histogram1D_map.end(); ++histogram1D_it) {
        histogram1D_it->second->Draw(); 
        histogram1D_it->second->SetName(histogram1D_it->second->GetTitle());
        histogram1D_it->second->Write();
    }

    std::map<std::string, TH2*>::iterator histogram2D_it =  histogram2D_map.begin();
    for (histogram2D_it; histogram2D_it != histogram2D_map.end(); ++histogram2D_it) {
        histogram2D_it->second->Draw(); 
        histogram2D_it->second->SetName(histogram2D_it->second->GetTitle());
        histogram2D_it->second->Write();
    }

    std::map<std::string, TGraph*>::iterator graph_it =  graph_map.begin();
    for (graph_it; graph_it != graph_map.end(); ++graph_it) {
        graph_it->second->Draw("A*e"); 
        graph_it->second->SetName(graph_it->second->GetTitle());
        graph_it->second->Write();
    }

    delete output_file;
}

void Plotter::saveToPdf(std::string file_name) { 

    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    canvas->Print((file_name + "[").c_str());
    std::map<std::string, TH1*>::iterator histogram1D_it =  histogram1D_map.begin();
    for (histogram1D_it; histogram1D_it != histogram1D_map.end(); ++histogram1D_it) {
        histogram1D_it->second->Draw(); 
        canvas->Print((file_name + "(").c_str());
    }

    std::map<std::string, TH2*>::iterator histogram2D_it =  histogram2D_map.begin();
    for (histogram2D_it; histogram2D_it != histogram2D_map.end(); ++histogram2D_it) {
        histogram2D_it->second->Draw("colz"); 
        canvas->Print((file_name + "(").c_str());
    }

    std::map<std::string, TGraph*>::iterator graph_it =  graph_map.begin();
    for (graph_it; graph_it != graph_map.end(); ++graph_it) {
        graph_it->second->Draw("A*"); 
        canvas->Print((file_name + "(").c_str());
    }

    canvas->Print((file_name + "]").c_str());
    delete canvas;
}
