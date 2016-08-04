/**
 *	@file lcio_analyzer.cxx
 *	@brief Main application that loops through LCIO events and processes 
 *		   them using analyses of type LcioAbstractAnalysis.
 *	@author Omar Moreno <omoreno1@ucsc.edu>
 *	@date November 16, 2014
 *
 */

//----------------//
//   C++ StdLib   //
//----------------//
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <string>
#include <list>

//----------//
//   LCIO   //
//----------//
#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <EVENT/LCEvent.h>

//-------------------//
//   LCIO Analysis   //
//-------------------//
#include <LcioAbstractAnalysis.h>
#include <HitAnalysis.h>
#include <TrackAnalysis.h>
#include <TaggerTrackerAnalysis.h>

using namespace std;

void printUsage(); 

int main(int argc, char** argv) { 

	// Name of the LCIO file to be processed.
	string lcio_file_name; 
    string file_list_name;
    int event_count = -1; 

	// Parse all the command line arguments.  If there are no valid command 
	// line arguments passed, print the usage and exit.
	static struct option long_options[] = { 
		{ "file_name", required_argument, 0, 'i' }, 
        {"file_list",  required_argument, 0, 'l'},
        {"events",     required_argument, 0, 'n'},
        {"help",       no_argument,       0, 'h'},
		{ 0, 0, 0, 0 } 
	};

	int option_index = 0; 
	int option_char; 
	while ((option_char = getopt_long(argc, argv, "i:l:n:h", long_options, &option_index)) != -1) { 
		switch (option_char) { 
			case 'i': 
				lcio_file_name = optarg;
				break;
            case 'l':
                file_list_name = optarg;
                break; 
            case 'n':
                event_count = atoi(optarg);
                break;
            case 'h':
                printUsage();
                return EXIT_SUCCESS; 
			default: 
                printUsage(); 
				return EXIT_FAILURE; 	
		}
	}

	// If an LCIO file was not specified, warn the user and exit the application.
	if (lcio_file_name.empty() && file_list_name.empty()) { 
		cerr << "\n[ LCIO ANALYZER ]: Please specify a file to process." << endl;
        cerr << "[ LCIO ANALYZER ]: Use --help for usage.\n" << endl;
		return EXIT_FAILURE;
	} else if (!lcio_file_name.empty() && !file_list_name.empty()) { 
        cerr << "\n[ LCIO ANALYZER ]: Cannot specify both an LCIO file name an " 
             << " a list of files." << endl;
        cerr << "[ LCIO ANALYZER ]: Use --help for usage.\n" << endl;
		return EXIT_FAILURE;
    }

    // Create a list of files to process
    list<string> files; 
    string file;
    if (!lcio_file_name.empty()) { 
        files.push_back(lcio_file_name); 
    } else if (!file_list_name.empty()) { 
        
        ifstream file_list(file_list_name.c_str(), ifstream::in);
        if (!file_list.is_open()) { 
            cerr << "\n[ LCIO ANALYZER ]: Failed to open file " << file_list_name << endl;
            return EXIT_FAILURE;
        }
        
        while (file_list >> file) { 
            files.push_back(file); 
        }
        file_list.close();
    }

	// Container to hold all analyses
	list<LcioAbstractAnalysis*> analyses;

	// All all analyses that are to be run.
	analyses.push_back(new HitAnalysis());
	analyses.push_back(new TrackAnalysis());
    analyses.push_back(new TaggerTrackerAnalysis());

	// Create the LCIO reader and open the file.  If the file can't be opened, 
	// warn the user and exit the application.
	IO::LCReader *lc_reader = IOIMPL::LCFactory::getInstance()->createLCReader();

	// Initialize all analyses
	for (list<LcioAbstractAnalysis*>::iterator analysis = analyses.begin(); 
			analysis != analyses.end(); ++analysis) { 
		cout << "[ LCIO ANALYZER ]: Initializing analysis: " << (*analysis)->toString() << endl;
		(*analysis)->initialize();
	}	

    // Loop over all input files and process them
    for (list<string>::iterator files_it = files.begin(); files_it != files.end(); ++files_it) { 
        
	    try {
		    lc_reader->open((*files_it).c_str());
	    } catch (IO::IOException &e) { 
		    cout << "[ LCIO ANALYZER ]: File " << lcio_file_name << " cannot be opened!" << endl;
	        cout << "[ LCIO ANALYZER ]: File will not be processed." << endl;
        } 
        cout << "[ LCIO ANALYZER ]: Processing file: " << lcio_file_name << endl;

	    EVENT::LCEvent* event = NULL;
	    int event_number = 0;
	    // Loop over all events in the file
	    while ((event = lc_reader->readNextEvent())) { 
		    ++event_number;
		    if (event_number%500 == 0) 
			    cout << "[ LCIO ANALYZER ]: Processing event " << event_number << endl;

		    for (list<LcioAbstractAnalysis*>::iterator analysis = analyses.begin();
			       	analysis != analyses.end(); ++analysis) { 
			    (*analysis)->processEvent(event);
		    }	
	        if (event_number + 1 == event_count) break;
        }
        lc_reader->close(); 
    }

	for (list<LcioAbstractAnalysis*>::iterator analysis = analyses.begin();
		   	analysis != analyses.end(); ++analysis) { 
		(*analysis)->finalize();
		delete *analysis;
	}

	analyses.clear();	

	return EXIT_SUCCESS; 
}

void printUsage() {
    cout << "\nUsage: hps_analysis [OPTIONS]" << endl;
}
