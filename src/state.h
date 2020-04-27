/*
 * This file is part of the M. xanthus polarity simulator
 * Copyright (c) 2020 Filipe Tostevin
 *
 * Various arrays and functions associated with the state of the dynamical
 * system.
 */

#ifndef STATE_H
#define STATE_H

#include <fstream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>

// Enumerate the protein species in the system so we can reference them by name
enum species{
	A,
	B,
	R,
	NSPEC
};
// Name the protein species so we can get them from the input file and write to
// the output file by name
const std::string species_names[NSPEC] = {
	"A",
	"B",
	"R"
};

// Enumerate the possible pools to which proteins can contribute
enum pool{
	POLE1,
	POLE2,
	TOTAL,
	NPOOL
};

// Matrix object that will contain the state of the system
typedef boost::numeric::ublas::matrix<double> state;

// Utility functions
void write_state(const state&, const double, std::ofstream&);

#endif
