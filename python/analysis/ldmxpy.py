#!/usr/bin/python

import argparse
import importlib
import ROOT as r
import os
import sys
import yaml

def parse_config(config_file) :

    print "Loading configuration from " + str(config_file)
    config = open(config_file, 'r')
    return yaml.load(config)

def main() : 
   
    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", action='store', dest='config',
                        help="Configuration file.")
    args = parser.parse_args()

    if not args.config :
        parser.error('A configuration file needs to be specified.')


    # Load the LDMX event model library.
    if not os.getenv('LDMX_SW_DIR') :
        print '[ ldmxpy ]: Error! Location of LDMX event library needs to be set.'
        sys.exit(2)
    ldmx_lib_path = os.environ['LDMX_SW_DIR'] + "/install/lib/libEvent.so"
    r.gSystem.Load(ldmx_lib_path)

    # Parse the configuration file
    config = parse_config(args.config)

    analyses = config["Analyses"]
    analyses_instances = []
    for analysis in analyses : 
        analysis_module_name, analysis_class_name = analysis.rsplit(".", 1)
        print "[ ldmxpy ]: Adding analysis ==> Module: %s Class: %s" % (analysis_module_name, analysis_class_name)
        analysis_class = getattr(importlib.import_module(analysis_module_name), analysis_class_name)
        analyses_instances.append(analysis_class())
        #analysis_class().initialize()

    for input_file in config["Files"] : 
        print 'Processing file ' + str(input_file)
        
        root_file = r.TFile(str(input_file))

        tree = root_file.Get("LDMX_Event")

        ldmx_event = r.event.SimEvent()
    
        b_ldmx_event = tree.GetBranch("LdmxEvent")
        b_ldmx_event.SetAddress(r.AddressOf(ldmx_event))

        for entry in xrange(0, tree.GetEntries()) : 

            if (entry + 1)%500 == 0 : print "Event " + str(entry + 1)

            tree.GetEntry(entry)
        
            for analyses in analyses_instances : 
                analyses.process(ldmx_event)

    for analyses in analyses_instances : 
        analyses.finalize()


if __name__ == "__main__":
    main()
