#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

////////////////////////////////////////////////////////////////////////////////////////////////////
////// This header contains functions to set states of linear Civil Aircraft to certain conditions
///// to investigate flight modes
////////////////////////////////////////////////////////////////////////////////////////////////////


Eigen::VectorXd phugoid() {

	Eigen::VectorXd Xlong_o(4);
	Xlong_o[0] = -1.99;
	Xlong_o[1] = 0.183;
	Xlong_o[2] = -0.0038;
	Xlong_o[3] = 0.0038;

	return Xlong_o;
}

Eigen::VectorXd short_period() {

	Eigen::VectorXd Xlong_o(4);
	Xlong_o[0] = 0.03;
	Xlong_o[1] = 1.99;
	Xlong_o[2] = -1.0048;
	Xlong_o[3] = 0.0199;

	return Xlong_o;
}