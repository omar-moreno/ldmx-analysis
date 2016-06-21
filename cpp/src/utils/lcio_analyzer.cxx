/**
 *	@file lcio_analyzer.cxx
 *	@brief Main application that loops through LCIO events and processes 
 *		   them using analyses of type LcioAbstractAnalysis.
 *	@author Omar Moreno <omoreno1@ucsc.edu>
 *	@date November 16, 2014
 *
 */

#include <iostream>
#include <getopt.h>
#include <list>
#include <cstdlib>

#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <EVENT/LCEvent.h>

#include <LcioAbstractAnalysis.h>

#include <HitAnalysis.h>

using namespace std;

int main(int argc, char** argv) { 

	// Name of the LCIO file to be processed.
	string lcio_file_name; 

	// Parse all the command line arguments.  If there are no valid command 
	// line arguments passed, print the usage and exit.
	static struct option long_options[] = { 
		{ "file_name", required_argument, 0, 'i' }, 
		{ 0, 0, 0, 0 } 
	};

	int option_index = 0; 
	int option_char; 
	while ((option_char = getopt_long(argc, argv, "i:", long_options, &option_index)) != -1) { 
		switch (option_char) { 
			case 'i': 
				lcio_file_name = optarg;
				break;
			default: 
				return EXIT_FAILURE; 	
		}
	}

	// If an LCIO file was not specified, warn the user and exit the application.
	if (lcio_file_name.empty()) { 
		cerr << "\n[ LCIO ANALYZER ]: Please specify a file to process." << endl;
        cerr << "[ LCIO ANALYZER ]: Use --help for usage.\n" << endl;
		return EXIT_FAILURE;
	}

	// Create the LCIO reader and open the file.  If the file can't be opened, 
	// warn the user and exit the application.
	IO::LCReader *lc_reader = IOIMPL::LCFactory::getInstance()->createLCReader();
	try {
		lc_reader->open(lcio_file_name.c_str());
	} catch (IO::IOException &e) { 
		cout << "[ LCIO ANALYZER ]: File " << lcio_file_name << " cannot be opened!" << endl;
		return EXIT_FAILURE;
	}

	cout << "[ LCIO ANALYZER ]: Processing file: " << lcio_file_name << endl;

	// Container to hold all analyses
	list<LcioAbstractAnalysis*> analyses;

	// All all analyses that are to be run.
	analyses.push_back(new HitAnalysis());

	// Initialize all analyses
	for (list<LcioAbstractAnalysis*>::iterator analysis = analyses.begin(); 
			analysis != analyses.end(); ++analysis) { 
		cout << "[ LCIO ANALYZER ]: Initializing analysis: " << (*analysis)->toString() << endl;
		(*analysis)->initialize();
	}	

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
	}

	for (list<LcioAbstractAnalysis*>::iterator analysis = analyses.begin();
		   	analysis != analyses.end(); ++analysis) { 
		(*analysis)->finalize();
		delete *analysis;
	}

	analyses.clear();	

	return EXIT_SUCCESS; 
}
