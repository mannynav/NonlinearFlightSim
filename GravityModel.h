#pragma once

#include <functional>

#ifndef GRAVITY_MODEL_H
#define GRAVITY_MODEL_H

using GravityModel = std::function<double(double)>;

//ISA-based inverse-square gravity model
inline double isa_gravity_model(double h) {

	const double g_0 = 9.80665; //Standard gravity (m/s^2)
	const double R_e = 6371000.0; //Earth's mean radius (m)
	return g_0 * pow(R_e / (R_e + h), 2);

}


//Simplified 1976 US Standard Atmosphere gravity model
inline double us1976_gravity_model_simplified(double h) {

	const double g_0 = 9.80665; //m/s^2
	const double R_e = 6371000.0; //Earths's mean radius
	return g_0 * (1.0 - 2.0 * h / R_e);
}





#endif