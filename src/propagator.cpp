/*
 * This file is part of the M. xanthus polarity simulator
 * Copyright (c) 2020 Filipe Tostevin
 *
 * Definitions of propagator function
 */

#include "propagator.h"


// Function passed to the integrator that actually does the evaluation of the
// dynamics of the ODE
void propagator::operator() (
	const state& conc,
	state& dconcdt,
	const double /*t*/
){
	for(int s=0; s<NSPEC; ++s) {
		dconcdt(s,TOTAL)=0.;
	}
	for(int p=0; p<TOTAL; ++p) {
		dconcdt(A,p)=( params.at("kA")
		              +params.at("kAR")*conc(R,p)
		             )*(conc(A,TOTAL)-(conc(A,POLE1)+conc(A,POLE2)))
		             -params.at("dA")*conc(A,p)
		             -params.at("dAB")*conc(A,p)*conc(B,p)*conc(B,p);
		dconcdt(B,p)=( params.at("kB")
		              +params.at("kBR")*conc(R,p)
		             )*(conc(B,TOTAL)-(conc(B,POLE1)+conc(B,POLE2)))
		             -params.at("dB")*conc(B,p)
		             -params.at("dBA")*conc(A,p)*conc(B,p)*conc(B,p);
		dconcdt(R,p)=( params.at("kR")
		              +params.at("kRR")*conc(R,p)
		              +params.at("kRB")*conc(B,p)/(1.+conc(A,p)/params.at("K"))
		             )*(conc(R,TOTAL)-(conc(R,POLE1)+conc(R,POLE2)))
		             -params.at("dR")*conc(R,p)/(1.+((p==POLE1 && (conc(A,TOTAL)<0.5 || conc(B,TOTAL)<0.5 || conc(R,TOTAL)<0.5))?params.at("Rbias"):0.));
	}
}

