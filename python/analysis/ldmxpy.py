#!/usr/bin/python

import argparse
import importlib
import sys
import yaml

def parse_config(config_file) :

    print "Loading configuration from " + str(config_file)
    config = open(config_file, 'r')
    return yaml.load(config)

def main() : 
   
    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--config", 
                        help="Configuration file describing analysis paramters.")
    args = parser.parse_args()

    if args.config is None: 
        print 'A configuration file needs to be specified!'
        sys.exit(2)
    
    config = parse_config(args.config)

    analyses = config["Analyses"]
    analyses_instances = []
    for analysis in analyses : 
        analysis_module_name, analysis_class_name = analysis.rsplit(".", 1)
        print "[ ldmxpy ]: Adding analysis ==> Module: " + str(analysis_module_name) + " Class: " + str(analysis_class_name)
        analysis_class = getattr(importlib.import_module(analysis_module_name), analysis_class_name)
        analyses_instances.append(analysis_class())

    for input_file in config["Files"] : 
        print 'Processing file ' + str(input_file)

        for analyses in analyses_instances : 
            analyses.process(input_file)

    for analyses in analyses_instances : 
        analyses.finalize()


if __name__ == "__main__":
    main()
