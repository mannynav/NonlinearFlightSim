

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> // Include Eigen's Eigenvalues module for eigenvalue and eigenvector computation
#include <numbers>
#include <iostream>
#include <fstream>

#include "Aircraft.h"
#include "NumericalIntegration.h"
#include "Output.h"
#include "FlightModes.h"

#include "map"

using Eigen::VectorXd;
using Eigen::MatrixXd;


Eigen::VectorXd Aircraft_Sim(CivilAircraft& aircraft, Eigen::VectorXd& X, Eigen::VectorXd& U)
{

	double pi = std::numbers::pi;

	double u = X[0]; //u
	double v = X[1]; //v
	double w = X[2]; //w
	double p = X[3]; //p
	double q = X[4]; //q
	double r = X[5]; //r
	double phi = X[6]; //phi
	const double theta = X[7]; //theta
	double psi = X[8]; //psi

	Eigen::VectorXd states = aircraft.initializeStates(X);
	u = states[0];
	v = states[1];
	w = states[2];
	p = states[3];
	q = states[4];
	r = states[5];
	phi = states[6];
	//theta = states[7];
	psi = states[8];



	Eigen::VectorXd controls = aircraft.intializeControlInputs(U);
	double aileron_inp = controls[0];
	double stabilizer_inp = controls[1];
	double rudder_inp = controls[2];
	double throttle1_inp = controls[3];
	double throttle2_inp = controls[4];

	////////////////////// Model Constants ///////////////////


	//Retrieve aircraft dictionary for constants
	std::map<std::string, double> aircraft_specs = aircraft.getAircraftSpecs();


	//Constants from aircraft specs
	double mass = aircraft_specs["mass"];
	double mean_aerodynamic_chord = aircraft_specs["mean_aerodynamic_chord"];
	double tail_body_AC_distance = aircraft_specs["lt"];
	double wing_planform_area = aircraft_specs["wing_planform_area"];
	double tail_planform_area = aircraft_specs["tail_planform_area"];

	double x_cg_pos_Fm = aircraft_specs["x_cg_pos_Fm"];
	double y_cg_pos_Fm = aircraft_specs["y_cg_pos_Fm"];
	double z_cg_pos_Fm = aircraft_specs["z_cg_pos_Fm"];

	double x_aero_pos_Fm = aircraft_specs["x_aero_pos_Fm"];
	double y_aero_pos_Fm = aircraft_specs["y_aero_pos_Fm"];
	double z_aero_pos_Fm = aircraft_specs["z_aero_pos_Fm"];


	//Engine constants
	double eng1_x_pos_Fm = aircraft_specs["eng1_x_pos_Fm"];
	double eng1_y_pos_Fm = aircraft_specs["eng1_y_pos_Fm"];
	double eng1_z_pos_Fm = aircraft_specs["eng1_z_pos_Fm"];

	double eng2_x_pos_Fm = aircraft_specs["eng2_x_pos_Fm"];
	double eng2_y_pos_Fm = aircraft_specs["eng2_y_pos_Fm"];
	double eng2_z_pos_Fm = aircraft_specs["eng2_z_pos_Fm"];




	//Other constants
	double sea_level_air_density = 1.225; //this will change when altitude model is implemented as a state
	double g = 9.81; //this will change when gravity model is implemented


	double depsda = 0.25;
	double alpha_initial = -11.5 * (pi / 180);
	double slope_coefficient_of_lift = 5.5; 
	double a3 = -768.5;
	double a2 = 609.2;
	double a1 = -155.2;
	double a0 = -15.212;
	double alpha_switch = 14.5 * (pi / 180); //lift slope goes from linear to non-linear


	//////////////// Saturation for control inputs - not used for now /////////////////
	
	double aileron_min = -25 * pi / 180;
	double aileron_max = 25 * pi / 180;

	double stabilizer_min = -25 * pi / 180;
	double stabilizer_max = 10 * pi / 180;

	double rudder_min = -30 * pi / 180;
	double rudder_max = 30 * pi / 180;

	double throttle1_min = 0.5 * pi / 180;
	double throttle1_max = 10 * pi / 180;

	double throttle2_min = 0.5 * pi / 180;
	double throttle2_max = 10 * pi / 180;


	//////////////// Other initializations using states /////////////////////
	double airSpeed = sqrt(u * u + v * v + w * w);
	double alpha = atan2(w, u);
	double beta = asin(v / airSpeed);
	double dynamicPressure = 0.5 * sea_level_air_density * pow(airSpeed, 2); //will change with airspeed and air density.

	Eigen::Vector3d angular_vel_be_body_frame({ p,q,r });
	Eigen::Vector3d translational_vel_body_frame({ u,v,w });



	////////////////// Aerodynamic Force Coefficients in stability axis /////////////////////////
	
	//Calculate coefficient of lift for wing body in stability axis
	double coefficient_of_lift_wb{};
	if (alpha <= alpha_switch)
	{
		coefficient_of_lift_wb = slope_coefficient_of_lift * (alpha - alpha_initial);
	}
	else {
		coefficient_of_lift_wb = a3 * pow(alpha, 3) + a2 * pow(alpha, 2) + a1 * alpha + a0;
	}

	////We need to rotate aerodynamic forces from the stability frame to the body frame via a rotation matrix
	Eigen::Matrix3d rot_stab_to_body(3, 3);
	rot_stab_to_body << cos(alpha), 0, -sin(alpha),
					0, 1, 0,
					sin(alpha), 0, cos(alpha);



	////Calculate gravitational forces in body frame. This does not cause a moment.
	Eigen::Vector3d gravity_body({ -g * sin(theta), g * cos(theta) * sin(phi), g * cos(theta) * cos(phi) });
	Eigen::Vector3d gravity_forces_body = aircraft.mass * gravity_body;
	


	//From F_b (all forces resolved in Fb) and calculate the udot, vdot, wdot
	Eigen::Vector3d forces_body = gravity_forces_body + aircraft.engine_forces_body_frame(g) + rot_stab_to_body*aircraft.aerodynamic_forces_stability_axis(alpha,beta,airSpeed,dynamicPressure);
	Eigen::Vector3d x1tox3dot = (1.0 / aircraft.mass) * forces_body - angular_vel_be_body_frame.cross(translational_vel_body_frame); 


	Eigen::Vector3d x4tox6dot = aircraft.cg_moments_body_frame(rot_stab_to_body, angular_vel_be_body_frame, airSpeed, dynamicPressure, alpha, beta,g);
	


	//Calculate phidot, thetadot, psidot
	Eigen::Matrix3d H_pi(3, 3);
	H_pi << 1, sin(phi)* tan(theta), cos(phi)* tan(theta),
		0, cos(phi), sin(phi),
		0, sin(phi) / cos(theta), cos(phi) / cos(theta);


	//Eigen::Vector3d x7tox9dot = H_pi * angular_vel_be_body_frame;
	Eigen::Vector3d x7tox9dot = aircraft.euler_angles_euler_kinematics(angular_vel_be_body_frame, phi, theta, psi);



	Eigen::VectorXd XDOT(9);
	XDOT[0] = x1tox3dot[0];
	XDOT[1] = x1tox3dot[1];
	XDOT[2] = x1tox3dot[2];
	XDOT[3] = x4tox6dot[0];
	XDOT[4] = x4tox6dot[1];
	XDOT[5] = x4tox6dot[2];
	XDOT[6] = x7tox9dot[0];
	XDOT[7] = x7tox9dot[1];
	XDOT[8] = x7tox9dot[2];

	return XDOT;

}

Eigen::VectorXd Implicit_Model(CivilAircraft& ac, Eigen::VectorXd& XDOT, Eigen::VectorXd& X, Eigen::VectorXd& U) {

	
	return Aircraft_Sim(ac, X, U) - XDOT;
}

void initialStatesControls(Eigen::VectorXd& states, Eigen::VectorXd& controls)
{
	states[0] = 85;
	states[1] = 0.0;
	states[2] = 0.0;
	states[3] = 0.0;
	states[4] = 0.0;
	states[5] = 0.0;
	states[6] = 0.0;
	states[7] = 0.1;
	states[8] = 0.0;

	controls[0] = 0.0;
	controls[1] = -0.1;
	controls[2] = 0.0;
	controls[3] = 0.08;
	controls[4] = 0.08;
}

void initializeTrimStatesInputsAndIncrements(Eigen::VectorXd& Xdoto, Eigen::VectorXd& Xo, Eigen::VectorXd& Uo, Eigen::MatrixXd& dxdot, Eigen::MatrixXd& dx, Eigen::MatrixXd& du)
{
	Xdoto.setZero();

	//Initialize states
	Xo[0] = 85;
	Xo[1] = 0.0;
	Xo[2] = 1.2713;
	Xo[3] = 0.0;
	Xo[4] = 0.0;
	Xo[5] = 0.0;
	Xo[6] = 0.0;
	Xo[7] = 0.0150;
	Xo[8] = 0.0;

	//Initialize inputs
	Uo[0] = 0.0;
	Uo[1] = -0.1780;
	Uo[2] = 0.0;
	Uo[3] = 0.0821;
	Uo[4] = 0.0821;

	dxdot.setConstant(0.00001);
	dx.setConstant(0.00001);
	du.setConstant(0.00001);
}

Eigen::MatrixXd LinearizeSystem_A(CivilAircraft& ac, Eigen::VectorXd& XDOTo, Eigen::VectorXd& Xo, Eigen::VectorXd& Uo, Eigen::MatrixXd& DXDOT, Eigen::MatrixXd& DX, Eigen::MatrixXd& DU)
{
	double n = XDOTo.size(); //number of states
	//double m = Uo.size(); //number of inputs

	Eigen::MatrixXd A(9, 9);

	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++) 
		{

			//get the perturbation or increment
			double dx = DX(i, j);

			//create the perturbation vector
			Eigen::VectorXd x_plus = Xo;
			Eigen::VectorXd x_minus = Xo;

			x_plus(j) = x_plus(j) + dx;
			x_minus(j) = x_minus(j) - dx;

			Eigen::VectorXd F = Implicit_Model(ac,XDOTo, x_plus, Uo);
			double F_plus_keep = F[i];

			F = Implicit_Model(ac,XDOTo, x_minus, Uo);
			double F_minus_keep = F[i];

			A(i, j) = (F_plus_keep - F_minus_keep) / (2 * dx);

		}
	}

	return A;
	
}

Eigen::MatrixXd LinearizeSystem_B(CivilAircraft& ac, Eigen::VectorXd& XDOTo, Eigen::VectorXd& Xo, Eigen::VectorXd& Uo, Eigen::MatrixXd& DXDOT, Eigen::MatrixXd& DX, Eigen::MatrixXd& DU)
{
	double n = XDOTo.size(); //number of states
	double m = Uo.size(); //number of inputs

	Eigen::MatrixXd B(9, 5);

	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < m; j++)
		{

			//get the perturbation or increment
			double du = DU(i, j);

			//create the perturbation vector
			Eigen::VectorXd u_plus = Uo;
			Eigen::VectorXd u_minus = Uo;

			u_plus(j) = u_plus(j) + du;
			u_minus(j) = u_minus(j) - du;

			Eigen::VectorXd F = Implicit_Model(ac,XDOTo, Xo, u_plus);

			double F_plus_keep = F[i];

			F = Implicit_Model(ac,XDOTo, Xo, u_minus);
	
			double F_minus_keep = F[i];
			
			B(i, j) = (F_plus_keep - F_minus_keep) / (2.0 * du);



		}
	}

	return B;

}

Eigen::VectorXd LinearStateSpace(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::MatrixXd& A, const Eigen::MatrixXd& B)
{
	return A * x + B * u;
}

Eigen::VectorXd readVectorFromCsv(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return Eigen::VectorXd(); // Return an empty vector
	}

	std::vector<double> data_vec;
	std::string line;

	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::string cell;

		// Assuming a single column, read the entire line as one value
		if (std::getline(ss, cell, ',')) { // Read until a comma (or end of line for single column)
			try {
				data_vec.push_back(std::stod(cell)); // Convert string to double
			}
			catch (const std::invalid_argument& e) {
				std::cerr << "Warning: Invalid number in CSV file '" << filename << "': " << cell << " (" << e.what() << ")" << std::endl;
				// You might choose to skip or handle this error differently
			}
			catch (const std::out_of_range& e) {
				std::cerr << "Warning: Number out of range in CSV file '" << filename << "': " << cell << " (" << e.what() << ")" << std::endl;
			}
		}
	}

	// Convert std::vector<double> to Eigen::VectorXd
	Eigen::VectorXd eigen_vector(data_vec.size());
	for (size_t i = 0; i < data_vec.size(); ++i) {
		eigen_vector(i) = data_vec[i];
	}

	return eigen_vector;
}

std::vector<double> readNumericCSVColumn(const std::string& filename, int columnIndex = 0)
{
	std::vector<double> data;
	std::ifstream file(filename);

	if (!file.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}

	

	std::string line;
	//Skip the header
	if (std::getline(file, line)) {}

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		std::vector<std::string> rowData;

		while (std::getline(lineStream, cell, ',')) {
			rowData.push_back(cell);
		}

		if (columnIndex >= rowData.size()) {
			throw std::runtime_error("Column index out of range in row: " + line);
		}

		try {
			double value = std::stod(rowData[columnIndex]);
			data.push_back(value);
		}
		catch (const std::invalid_argument& e) {
			throw std::runtime_error("Invalid numeric value in row: " + line + ", column: " + std::to_string(columnIndex));
		}
		catch (const std::out_of_range& e) {
			throw std::runtime_error("Numeric value out of range in row: " + line + ", column: " + std::to_string(columnIndex));
		}
	}

	return data;
}

int main()
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Initialization of non linear simulation
	////////////////////////////////////////////////////////////////////////////////////////////////////

	Eigen::VectorXd initialX(9); //Initial states for non linear simulation
	Eigen::VectorXd initialU(5); //Initial controls for nonlinear simlulation. Controls will be dependent on aircraft/object - 1 engine, 2 engine, etc
	initialStatesControls(initialX, initialU);
	
	CivilAircraft ac = CivilAircraft(initialX, initialU); 
	std::map<std::string, double> aircraft_map = ac.getAircraftSpecs();

	double simLength = 150;
	double steps = 200;
	double timeIncrement = simLength / steps;

	int rows = initialX.size();
	int time_steps = steps + 1;

	Eigen::MatrixXd nonLinearSolutionMatrix(rows, time_steps);
	nonLinearSolutionMatrix.setZero().col(0) = initialX;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Start non linear simulation
	////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i = 0; i < steps; i++)
	{
		Eigen::VectorXd nextk1 = nonLinearSolutionMatrix.col(i);
		Eigen::VectorXd k1 = Aircraft_Sim(ac,nextk1, initialU);

		Eigen::VectorXd nextk2 = nonLinearSolutionMatrix.col(i) + (timeIncrement / 2.0) * k1;
		Eigen::VectorXd k2 = Aircraft_Sim(ac,nextk2, initialU);

		Eigen::VectorXd nextk3 = nonLinearSolutionMatrix.col(i) + (timeIncrement / 2.0) * k2;
		Eigen::VectorXd k3 = Aircraft_Sim(ac,nextk3, initialU);

		Eigen::VectorXd nextk4 = nonLinearSolutionMatrix.col(i) + timeIncrement * k3;
		Eigen::VectorXd k4 = Aircraft_Sim(ac,nextk4, initialU);

		nonLinearSolutionMatrix.col(i + 1) = nonLinearSolutionMatrix.col(i) + (timeIncrement / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

	}

	

	/*for (int i = 0; i < steps; i++)
	{
		Eigen::VectorXd nextStep = nonLinearSolutionMatrix.col(i);
		nonLinearSolutionMatrix.col(i + 1) = nonLinearSolutionMatrix.col(i) + 0.100 * Aircraft_Sim(ac,nextStep, initialU);
	}*/

	////////////////////////////////////////////////////////////////////////////////////////////////////
	//// End of non linear simulation
	////////////////////////////////////////////////////////////////////////////////////////////////////

	std::string fileName = "Non Linear Solution.csv";
	outputToFile(nonLinearSolutionMatrix, fileName);


	std::cout << "----------------------------------------- End of non linear simulation -------------------------------------------------" << std::endl;




	std::cout << "----------------------------------------- Initializing trim point for Xo and Uo (Get from Matlab) -------------------------------------------------" << std::endl;

	//Read in state and control values from trim routine
	std::string states_file = "trim_states.csv";
	std::string controls_file = "trim_controls.csv";

	//Eigen::VectorXd states = readVectorFromCsv(states_file);
	//Eigen::VectorXd controls = readVectorFromCsv(controls_file);

	std::vector<double> states = readNumericCSVColumn(states_file);
	std::vector<double> controls = readNumericCSVColumn(controls_file);

	for (auto i : states) {
		std::cout << i << std::endl;
	}

	for (auto i : controls) {
		std::cout << i << std::endl;
	}


	//Verify trim conditions are satisfied using nonlinear simulation


	

	////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Initialize structures and matrices for linearization about straight and level trim
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//Vectors for states and inputs
	Eigen::VectorXd Xdoto(9);
	Eigen::VectorXd Xo(9);
	Eigen::VectorXd Uo(5);

	//Matrices for increments for the finite difference step
	Eigen::MatrixXd dxdot(9, 9);
	Eigen::MatrixXd dx(9, 9);
	Eigen::MatrixXd du(9, 5);

	initializeTrimStatesInputsAndIncrements(Xdoto, Xo, Uo, dxdot, dx, du);

	CivilAircraft acl = CivilAircraft(Xo, Uo);


	////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Linearize system to obtain A and B
	////////////////////////////////////////////////////////////////////////////////////////////////////

	std::cout << "----------------------------------------- Linearize system for A -------------------------------------------------" << std::endl;
	Eigen::MatrixXd Atest = LinearizeSystem_A(acl,Xdoto, Xo, Uo, dxdot, dx, du);
	std::cout << "Eigenvalues of A: " << std::endl;
	std::cout << Atest.eigenvalues() << std::endl;


	std::cout << "----------------------------------------- Linearize system for B -------------------------------------------------" << std::endl;
	Eigen::MatrixXd Btest = LinearizeSystem_B(acl,Xdoto, Xo, Uo, dxdot, dx, du);


	std::cout << "----------------------------------------- Start Longitudinal Model Process -------------------------------------------------" << std::endl;

	Eigen::MatrixXd S(9, 9); //Matrix for similarity transformation of states
	S.setZero();

	S << 1, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0,
		0, 1, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1;

	Eigen::MatrixXd T = S.inverse();

	Eigen::MatrixXd transformedA = T.inverse() * Atest * T;
	Eigen::MatrixXd transformedB = T.inverse() * Btest;

	std::cout << "transformed A: " << std::endl;
	std::cout << transformedA << std::endl;

	std::cout << "transformed B: " << std::endl;
	std::cout << transformedB << std::endl;

	Eigen::MatrixXd Alongitudinal = transformedA.block(0, 0, 4, 4);
	Eigen::MatrixXd Blongitudinal = transformedB.block(0, 0, 4, 5);

	std::cout << "A long: " << std::endl;
	std::cout << Alongitudinal << std::endl;

	std::cout << "Eigenvalues of A long: " << std::endl;
	std::cout << Alongitudinal.eigenvalues() << std::endl;

	Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(Alongitudinal); // Use EigenSolver to compute eigenvalues and eigenvectors
	std::cout << "Eigenvectors of A long: " << std::endl;
	std::cout << eigenSolver.eigenvectors() << std::endl;

	std::cout << "B long: " << std::endl;
	std::cout << Blongitudinal << std::endl;


	//Initialize simulation variables for longitudinal variables.
	double simTime = 150;
	double stepsLong = 150;
	double increment = simTime / stepsLong;

	//Eigen::VectorXd Xlong_o(4); //Phugoid.
	//Xlong_o = phugoid();

	//Eigen::VectorXd Xlong_o(4); //Short period
	//Xlong_o = short_period();

	Eigen::VectorXd Xlong_o(4); //Test Doublet
	Xlong_o[0] = 0.0;
	Xlong_o[1] = 0.0;
	Xlong_o[2] = 0.0;
	Xlong_o[3] = 0.0;

	int number_of_controls = 5;
	Eigen::MatrixXd longitudinal_control_inputs(number_of_controls, int(stepsLong) + 1);
	longitudinal_control_inputs.setZero();

	acl.initialize_deflections(0, longitudinal_control_inputs,0,0,0,0,stepsLong, increment); //initialize aileron
	acl.initialize_deflections(1, longitudinal_control_inputs,10,-10,10,10,stepsLong, increment); //initialize stabilizer
	acl.initialize_deflections(2, longitudinal_control_inputs,0,0,0,0,stepsLong, increment); //initialize rudder
	acl.initialize_thrusters(3, longitudinal_control_inputs, 0.5, 60, stepsLong, increment); //thruster 1
	acl.initialize_thrusters(4, longitudinal_control_inputs, 0.5, 60, stepsLong, increment); //thruster 2

	// Get the Dormand-Prince tableau
	//RKTableauSS tableau_ss = dp45_tableau_ss();
	std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)> ss_func = LinearStateSpace;


	//Perform simulation.
	Eigen::MatrixXd solution_linear_long = rk4_simulate_ss(ss_func, Xlong_o, longitudinal_control_inputs, 0, simTime, stepsLong, Alongitudinal, Blongitudinal);
	//Eigen::MatrixXd solution_linear_long = adaptive_rk_simulate_ss(ss_func, Xlong_o, Ulong_o, 0, 1000, 500, 0.03, tableau_ss,Alongitudinal, Blongitudinal);
	//Eigen::MatrixXd solution_linear_long = forward_euler_simulate_ss(ss_func, Xlong_o, Ulong_o, 0, 150, 300, Alongitudinal, Blongitudinal);

	fileName = "Linear Longitudinal Solution.csv";
	std::string file_longitudinal_control_inputs = "Linear Longitudinal Control Inputs.csv";

	outputToFile(solution_linear_long, fileName);
	outputToFile(longitudinal_control_inputs, file_longitudinal_control_inputs);




	std::cout << "----------------------------------------- Start Lateral Model Process -------------------------------------------------" << std::endl;
	
	//Initialize simulation variables.
	double lateral_sim_time = 200;
	double lateral_steps = 250;
	double lateral_increment = simTime / lateral_steps;

	//Initialize states
	Eigen::VectorXd Xlateral_o(4);
	Xlateral_o[0] = -1.99;
	Xlateral_o[1] = 0.183;
	Xlateral_o[2] = -0.0038;
	Xlateral_o[3] = 0.0038;

	Eigen::MatrixXd lateral_control_inputs(number_of_controls, int(lateral_steps) + 1);
	lateral_control_inputs.setZero();

	acl.initialize_deflections(0, lateral_control_inputs, 0, 0, 0, 0, lateral_steps, increment); //initialize aileron
	acl.initialize_deflections(1, lateral_control_inputs, 10, -10, 10, 10, lateral_steps, increment); //initialize stabilizer
	acl.initialize_deflections(2, lateral_control_inputs, 0, 0, 0, 0, lateral_steps, increment); //initialize rudder
	acl.initialize_thrusters(3, lateral_control_inputs, 0.0, 60, lateral_steps, increment); //thruster 1
	acl.initialize_thrusters(4, lateral_control_inputs, 0.0, 60, lateral_steps, increment); //thruster 2
	
	//Extract A matrix for lateral system.
	Eigen::MatrixXd Alateral = transformedA.block(4,4, 4,4);
	std::cout << "Alateral: " << std::endl;
	std::cout << Alateral << std::endl;
	std::cout << "A lateral eigenvalues: " << std::endl;
	std::cout << Alateral.eigenvalues() << std::endl;

	//Eigen::EigenSolver<Eigen::MatrixXd> eigenSolverLateral (Alateral); // Use EigenSolver to compute eigenvalues and eigenvectors
	//std::cout << "Eigenvectors of A lateral: " << std::endl;
	//std::cout << eigenSolverLateral.eigenvectors() << std::endl;

	//Extract B matrix for lateral system.
	Eigen::MatrixXd Blateral = transformedB.block(4, 0, 4, 5);
	std::cout << "Blateral: " << std::endl;
	std::cout << Blateral << std::endl;


	//tableau for adaptive RK.
	//RKTableauSS tableau_ss = dp45_tableau_ss();
	std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)> ss_func_lateral = LinearStateSpace;

	//Perform simulation.
	Eigen::MatrixXd solution_linear_lateral = rk4_simulate_ss(ss_func_lateral, Xlateral_o, lateral_control_inputs, 0, lateral_sim_time, lateral_steps, Alateral, Blateral);
	//Eigen::MatrixXd solution_linear_long = adaptive_rk_simulate_ss(ss_func, Xlateral_o, Ulateral_inputs, 0, 1000, 500, 0.03, tableau_ss,Alateral, Blateral);
	//Eigen::MatrixXd solution_linear_long = forward_euler_simulate_ss(ss_func, Xlateral_o, Ulateral_inputs, 0, 150, 300, Alateral, Blateral);

	fileName = "Linear Lateral Solution.csv";
	std::string file_lateral_control_inputs = "Linear Lateral Control Inputs.csv";

	outputToFile(solution_linear_lateral, fileName);
	outputToFile(lateral_control_inputs, file_lateral_control_inputs);
 
	return 0;

}