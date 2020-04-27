/*
 * This file is part of the M. xanthus polarity simulator
 * Copyright (c) 2020 Filipe Tostevin
 *
 * Definition of a propagator class that wraps the ODE integrator.
 */

#ifndef PROP_H
#define PROP_H

#include <unordered_map>
#include <boost/numeric/odeint.hpp>
#include "state.h"

// Define a map object that associates parameter values with a name
typedef std::unordered_map<std::string, double> params_array;

class propagator{
	public:
		// System parameters
		params_array params;

		// Default constructor does nothing
		propagator() {};
		// Function passed to the integrator that actually does the evaluation of
		// the ODEs
		void operator() (const state&, state&, const double);
};

#endif
