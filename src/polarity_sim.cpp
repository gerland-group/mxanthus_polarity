/*
 * This file is part of the M. xanthus polarity simulator
 * Copyright (c) 2020 Filipe Tostevin
 *
 * This file contains the main routine that creates all simulator
 * objects and calls the propagation routines.
 */

#include <iostream>
#include <cstdlib>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/program_options.hpp>
#include "state.h"
#include "propagator.h"

using namespace std;
namespace po = boost::program_options;
namespace pt = boost::property_tree;

int main(
	int argc,
	char* argv[]
){
	// Command line options and help
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "print help message")
		("conf,c", po::value<string>(), "parameter configuration file")
		("mutants,m", "simulate deletion mutants")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if(vm.count("help")) {
		cout << desc << "\n";
		return(EXIT_SUCCESS);
	}

	// Read the input file given by the "--conf" option, or otherwise from
	// "params.in"
	pt::ptree ptree;
	string defaultconfig="params.in";
	if(vm.count("conf")) {
		cout << "Reading parameters from: " << vm["conf"].as<string>() << endl;
		pt::read_ini(vm["conf"].as<string>(), ptree);
	} else {
		cout << "No parameter file specified. Reading parameters from default file: " << defaultconfig << endl;
		pt::read_ini(defaultconfig, ptree);
	}

	// Get total simulation time from input file
	const double t_total=ptree.get<double>("t_total");
	const double sampling_interval=ptree.get<double>("sampling_interval");

	// Define the system state variable and initialize all elements to zero
	state conc(NSPEC,NPOOL,0);

	// Set up simulation infrastructure
	propagator P;
	// Read in parameter values from input file
	for(auto p: ptree.get_child("model_parameters")) {
		P.params.emplace(p.first, ptree.get<double>("model_parameters."+p.first));
	}
	//Run a test call of the propagator to check that all required parameters
	//exist
	try{ P(conc, conc, 0); }
	catch(...) {
		cout << "Error in model definition, check that all parameters are specified in the input file." << endl;
		exit(EXIT_FAILURE);
	}

	// Deletion mutants are defined by an integer that when represented as a
	// binary string defines the presence/absence of each protein,
	// i.e. if the number of proteins in NSPEC=3, then
	// mutant_number=1==b100 indicates only the first-listed protein is present,
	// mutant_number=7==b111 indicates that all proteins are present
	// Here we define the first simulation to run depending on whether we are
	// simulating all mutants, or just the wild-type
	int mutant_number;
	if(vm.count("mutants")) {
		cout << "Simulating deletion mutants" << endl;
		mutant_number=1;
	} else {
		cout << "Simulating only wild-type" << endl;
		mutant_number=(1<<NSPEC) -1;
	}

	// Loop over all mutants until we reach wild-type
	while(mutant_number < 1<<NSPEC) {
		cout << mutant_number << endl;
		//Each mutant will be output to a separate data file.
		//Decode mutant number bitwise into whether protein is present or not and
		//reconstruct readable output string from it, to use in the name of the
		//output file
		string mutant_string="";
		for(int s=0; s<NSPEC; ++s) {
			if( !( (mutant_number >> s) & 1) ) {
				mutant_string+=".d" + species_names[s];
			}
		}
		ofstream datafile("data.txt"+mutant_string);

		//Headings for output file
		datafile << "t ";
		for(int s=0; s<NSPEC; ++s) {
			datafile << species_names[s] << "1 "
			         << species_names[s] << "2 "
			         << species_names[s] << "t ";
		} datafile << endl;

		for(int s=0; s<NSPEC; ++s) {
			// Initialize total amounts of each protein to 1 or 0 according to the
			// deletion condition
			if( (mutant_number >> s) & 1 ) {
				conc(s,TOTAL)=1;
			} else {
				conc(s,TOTAL)=0;
			}
			// Initialize polar fractions
			conc(s,POLE1)=ptree.get<double>("initial_fractions_"+species_names[s]+".pole1")*conc(s,TOTAL);
			conc(s,POLE2)=ptree.get<double>("initial_fractions_"+species_names[s]+".pole2")*conc(s,TOTAL);
		}

		// Start at t=0
		double t=0;
		write_state(conc, t, datafile);

		// Run until the end of the simulation time t_total
		while(t<t_total) {
			// Call the ODE integrator to run until the next output time point
			boost::numeric::odeint::integrate(boost::ref(P), conc, t, t+sampling_interval, 0.2*sampling_interval);
			t+=sampling_interval;

			write_state(conc, t, datafile);
		}

		datafile.close();
		++mutant_number;
	}

	return EXIT_SUCCESS;
}
