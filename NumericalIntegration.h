#pragma once

#include <Eigen/Dense>
#include <numbers>
#include <iostream>
#include <fstream>

using Eigen::VectorXd;
using Eigen::MatrixXd;

//Simple forward euler integration scheme

Eigen::MatrixXd forward_euler_simulate_ss(
	const std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)>& state_space_model,
	const VectorXd& initial_x,
	const VectorXd& control_input,
	double start_time,
	double end_time,
	double steps,
	const MatrixXd& A_matrix,
	const MatrixXd& B_matrix)
{
	if (steps <= 0.0) {
		throw std::invalid_argument("Steps must be positive.");
	}

	int num_steps = static_cast<int>(steps) + 1;
	double timeIncrement = end_time / steps;

	Eigen::MatrixXd solution_matrix(initial_x.size(), num_steps);
	solution_matrix.col(0) = initial_x;

	VectorXd current_x = initial_x;
	double current_time = start_time;


	for (int i = 0; i < steps; i++)
	{
		solution_matrix.col(i + 1) = solution_matrix.col(i) + timeIncrement * state_space_model(solution_matrix.col(i), control_input, A_matrix, B_matrix);

	}

	return solution_matrix;

}



//4th order RK scheme

Eigen::VectorXd rk4_step_ss(
	const std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)>& state_space_model,
	const VectorXd& x_n,
	const VectorXd& u_n,
	double h,
	const MatrixXd& A,
	const MatrixXd& B) {

	VectorXd k1 = h * state_space_model(x_n, u_n, A, B);
	VectorXd k2 = h * state_space_model(x_n + 0.5 * k1, u_n, A, B);
	VectorXd k3 = h * state_space_model(x_n + 0.5 * k2, u_n, A, B);
	VectorXd k4 = h * state_space_model(x_n + k3, u_n, A, B);

	return x_n + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

Eigen::MatrixXd rk4_simulate_ss(
	const std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)>& state_space_model,
	const VectorXd& initial_x,
	const VectorXd& control_input,
	double start_time,
	double end_time,
	double steps,
	const MatrixXd& A_matrix,
	const MatrixXd& B_matrix) {

	if (steps <= 0.0) {
		throw std::invalid_argument("Steps must be positive.");
	}

	int num_steps = static_cast<int>(steps) + 1;
	double timeIncrement = end_time / steps;

	Eigen::MatrixXd solution_matrix(4, num_steps);
	solution_matrix.col(0) = initial_x;

	VectorXd current_x = initial_x;
	double current_time = start_time;

	for (int i = 1; i < num_steps; ++i) {

		current_x = rk4_step_ss(state_space_model, solution_matrix.col(i - 1), control_input, timeIncrement, A_matrix, B_matrix);
		solution_matrix.col(i) = current_x;
		current_time += timeIncrement;
	}

	return solution_matrix;
}




//4th order adaptive RK scheme

// Butcher tableau for embedded Runge-Kutta method (e.g., Dormand-Prince 4(5))
struct RKTableauSS {
	Eigen::MatrixXd a;
	Eigen::VectorXd b;   // Coefficients for the higher-order method
	Eigen::VectorXd b_alt; // Coefficients for the lower-order embedded method
	Eigen::VectorXd c;
};


// Dormand-Prince 4(5) tableau
RKTableauSS dp45_tableau_ss() {
	RKTableauSS tableau;
	tableau.a.resize(7, 6);
	tableau.b.resize(7);
	tableau.b_alt.resize(7);
	tableau.c.resize(7);

	tableau.c << 0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0;

	tableau.a <<
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.2, 0.0, 0.0, 0.0, 0.0, 0.0,
		3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0,
		44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0,
		19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0,
		9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0,
		35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 2565.0, 11.0 / 84.0;

	tableau.b << 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 2565.0, 11.0 / 84.0, 0.0;
	tableau.b_alt << 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0;

	return tableau;
}


std::pair<VectorXd, double> adaptive_rk_step_ss(
	const std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)>& state_space_model,
	const VectorXd& x_n,
	double t_n,
	double h,
	double tolerance,
	const RKTableauSS& tableau,
	const VectorXd& u_n, // Current control input
	const MatrixXd& A,
	const MatrixXd& B) {

	int s = tableau.a.rows();
	Eigen::MatrixXd k(x_n.size(), s);

	for (int i = 0; i < s; ++i) {
		VectorXd sum_ak = VectorXd::Zero(x_n.size());
		for (int j = 0; j < i; ++j) {
			sum_ak += tableau.a(i, j) * k.col(j);
		}
		k.col(i) = state_space_model(x_n + h * sum_ak, u_n, A, B);
	}

	VectorXd x_np1 = x_n;
	for (int i = 0; i < s; ++i) {
		x_np1 += h * tableau.b(i) * k.col(i);
	}

	VectorXd x_np1_alt = x_n;
	for (int i = 0; i < s; ++i) {
		x_np1_alt += h * tableau.b_alt(i) * k.col(i);
	}

	VectorXd error_estimate = (x_np1 - x_np1_alt).cwiseAbs();
	double error_norm = error_estimate.maxCoeff(); // Using max norm

	double safety_factor = 0.9;
	double p = tableau.b.size() - 1; // Order of the higher method

	double h_new;
	if (error_norm <= tolerance) {
		h_new = h * safety_factor * std::pow(tolerance / error_norm, 1.0 / (p + 1));
		return { x_np1, h_new };
	}
	else {
		h_new = h * safety_factor * std::pow(tolerance / error_norm, 1.0 / p);
		return { x_n, -std::abs(h_new) }; // Negative h_new indicates rejection
	}
}



Eigen::MatrixXd adaptive_rk_simulate_ss(
	const std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)>& state_space_model,
	VectorXd initial_x,
	const VectorXd& control_input,
	double start_time,
	double end_time,
	double initial_step_size,
	double tolerance,
	const RKTableauSS& tableau,
	const MatrixXd& A_matrix,
	const MatrixXd& B_matrix) {

	if (initial_step_size <= 0.0 || tolerance <= 0.0) {
		throw std::invalid_argument("Initial step size and tolerance must be positive.");
	}

	std::vector<VectorXd> solution_points;
	std::vector<double> time_points;

	solution_points.push_back(initial_x);
	time_points.push_back(start_time);

	VectorXd current_x = initial_x;
	double current_time = start_time;
	double h = initial_step_size;

	while (current_time < end_time) {
		if (current_time + h > end_time) {
			h = end_time - current_time;
		}

		VectorXd current_u = control_input;
		std::pair<VectorXd, double> result = adaptive_rk_step_ss(state_space_model, current_x, current_time, h, tolerance, tableau, current_u, A_matrix, B_matrix);
		VectorXd next_x = result.first;
		double h_new = result.second;

		if (h_new > 0) {
			solution_points.push_back(next_x);
			time_points.push_back(current_time + h);
			current_x = next_x;
			current_time += h;
			h = h_new;
		}
		else {
			h = -h_new; // Retry with the smaller step size
		}

		if (h < std::numeric_limits<double>::epsilon() * std::abs(current_time)) {
			std::cerr << "Warning: Step size became too small. Simulation might not converge." << std::endl;
			break;
		}
	}

	Eigen::MatrixXd solution_matrix(initial_x.size(), solution_points.size());
	for (size_t i = 0; i < solution_points.size(); ++i) {
		solution_matrix.col(i) = solution_points[i];
	}

	return solution_matrix;
}
