/*
 * This file is part of the M. xanthus polarity simulator
 * Copyright (c) 2020 Filipe Tostevin
 *
 * Definitions of utility functions related to the state of the system.
 */

#include "state.h"

// Write the current state (polar fractions and total of each protein) to the
// output stream
void write_state(
	const state& conc,
	const double t,
	std::ofstream& out
){
	out << t << " ";
	for(int s=0; s<NSPEC; ++s) {
		for(int p=0; p<NPOOL; ++p) {
			out << conc(s,p) << " ";
		}
	}
	out << std::endl;
}
