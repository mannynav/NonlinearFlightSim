#pragma once



/// <summary>
/// This class will be used a base class for aircraft or objects that will be simulated. 
/// </summary>

class Aircraft
{
public:


private:


};




/// <summary>
/// Non-linear Civil Aircraft Model with 5 controls
/// </summary>

class CivilAircraft : Aircraft
{
public:

	//Constructor for initialization of aircraft constants.
	CivilAircraft(Eigen::VectorXd& states, Eigen::VectorXd& controls) : pi(std::numbers::pi)
	{

		u = states[0];
		v = states[1];
		w = states[2];
		p = states[3];
		q = states[4];
		r = states[5];
		phi = states[6];
		theta = states[7];
		psi = states[8];

		aileron_inp = controls[0];
		stabilizer_inp = controls[1];
		rudder_inp = controls[2];
		throttle1_inp = controls[3];
		throttle2_inp = controls[4];

		//aircraft specs
		mass = 120000;
		mean_aerodynamic_chord = 6.6;
		lt = 24.8;
		wing_planform_area = 260;
		tail_planform_area = 64;


		x_cg_pos_Fm = 0.23 * mean_aerodynamic_chord;
		y_cg_pos_Fm = 0.0;
		z_cg_pos_Fm = 0.10 * mean_aerodynamic_chord;

		x_aero_pos_Fm = 0.12 * mean_aerodynamic_chord;
		y_aero_pos_Fm = 0;
		z_aero_pos_Fm = 0;

		InertiaMatrix << 40.07, 0, -2.0923,
			0, 64, 0,
			-2.0923, 0, 99.2;
		InertiaMatrix = mass * InertiaMatrix;

		InverseInertiaMatrix << 0.02498, 0, 0.000523,
			0, 0.015625, 0,
			0.0005232, 0, 0.010019;
		InverseInertiaMatrix = (1.0 / mass) * InverseInertiaMatrix;


		//Engine applications points of thrust
		eng1_x_pos_Fm = 0;
		eng1_y_pos_Fm = -7.94;
		eng1_z_pos_Fm = -1.9;

		eng2_x_pos_Fm = 0.0;
		eng2_y_pos_Fm = 7.94;
		eng2_z_pos_Fm = -1.9;


		//Aerodynamic force coefficients in stability axis
		depsda = 0.25;
		alpha_initial = -11.5 * (pi / 180);
		slope_coefficient_of_lift = 5.5;
		a3 = -768.5;
		a2 = 609.2;
		a1 = -155.2;
		a0 = -15.212;
		alpha_switch = 14.5 * (pi / 180); //lift slope goes from linear to non-linear



	}

	Eigen::VectorXd initializeStates(Eigen::VectorXd& stateVector)
	{
		//make sure to check size
		u = stateVector[0];
		v = stateVector[1];
		w = stateVector[2];
		p = stateVector[3];
		q = stateVector[4];
		r = stateVector[5];
		phi = stateVector[6];
		theta = stateVector[7];
		psi = stateVector[8];
		Eigen::VectorXd states(9);
		states[0] = u;
		states[1] = v;
		states[2] = w;
		states[3] = p;
		states[4] = q;
		states[5] = r;
		states[6] = phi;
		states[7] = theta;
		states[8] = psi;
		return states;
	}

	Eigen::VectorXd intializeControlInputs(Eigen::VectorXd& controlVector)
	{
		//make sure to check size
		aileron_inp = controlVector[0];
		stabilizer_inp = controlVector[1];
		rudder_inp = controlVector[2];
		throttle1_inp = controlVector[3];
		throttle2_inp = controlVector[4];

		Eigen::VectorXd controls(5);
		controls[0] = aileron_inp;
		controls[1] = stabilizer_inp;
		controls[2] = rudder_inp;
		controls[3] = throttle1_inp;
		controls[4] = throttle2_inp;
		return controls;
	}

	void initialize_deflections(int U_index, Eigen::MatrixXd& U_inputs, double pos_def_deg, double neg_def_deg, double phase_time1, double phase_time2, double steps, double increment)
	{
		// Convert degrees to radians
		double pos_rad = pos_def_deg * (3.14 / 180.0);
		double neg_rad = neg_def_deg * (3.14 / 180.0);

		// Apply the desired stabilizer profile (index 1 for stabilizer)
		for (int i = 0; i <= steps; ++i) {
			double current_time = i * increment;

			if (current_time >= 0 && current_time < phase_time1) {
				
				U_inputs(U_index, i) = pos_rad;
			}
			else if (current_time >= phase_time1 && current_time < 2 * phase_time2) {
				
				U_inputs(U_index, i) = neg_rad;
			}
		}
	}

	void initialize_thrusters(int U_index, Eigen::MatrixXd& U_inputs, double thrust_value, double phase_time, double steps, double increment)
	{
		if(U_index != 3 && U_index != 4)
		{
			throw std::out_of_range("U_index in initialize_thrusters must be either 3 or 4");
		}
	
		
		double thrust = thrust_value;

		for (int i = 0; i <= steps; ++i) {
			double current_time = i * increment;

			if (current_time >= 0 && current_time < phase_time) {

				U_inputs(U_index, i) = thrust;
			}
		}
	}


	Eigen::Vector3d engine_forces_body_frame(double g)
	{
		double engine1_thrust = throttle1_inp * mass * g;
		double engine2_thrust = throttle2_inp * mass * g;

		Eigen::Vector3d engine1_force_body({ engine1_thrust,0,0 });
		Eigen::Vector3d engine2_force_body({ engine2_thrust,0,0 });
		Eigen::Vector3d engine_forces_body = engine1_force_body + engine2_force_body;

		return engine_forces_body;
	}

	Eigen::Vector3d aerodynamic_forces_stability_axis(double alpha, double beta, double airSpeed, double dynamicPressure)
	{

		//Calculate coefficient of lift for wing body in stability axis
		double coefficient_of_lift_wb{};
		if (alpha <= alpha_switch)
		{
			coefficient_of_lift_wb = slope_coefficient_of_lift * (alpha - alpha_initial);
		}
		else {
			coefficient_of_lift_wb = a3 * pow(alpha, 3) + a2 * pow(alpha, 2) + a1 * alpha + a0;
		}


		//Calculate coefficient of lift for tail in stability axis
		double downwash_angle = depsda * (alpha - alpha_initial);
		double angle_of_attack_tail = alpha - downwash_angle + stabilizer_inp + 1.3 * q * lt / airSpeed;
		double coefficient_of_lift_tail = 3.1 * (tail_planform_area / wing_planform_area) * angle_of_attack_tail;

		double total_lift = coefficient_of_lift_wb + coefficient_of_lift_tail;


		//Total drag (not including tail)
		double total_drag = 0.13 + 0.07 * pow((5.5 * alpha + 0.654), 2);
		double side_force = -1.6 * beta + 0.24 * rudder_inp;


		Eigen::Vector3d aerodynamic_forces_stability({ -total_drag * dynamicPressure * wing_planform_area, side_force * dynamicPressure * wing_planform_area, -total_lift * dynamicPressure * wing_planform_area });

		return aerodynamic_forces_stability;

	}

	Eigen::Vector3d aerodynamic_forces_body_frame(double alpha, double beta, double airSpeed, double dynamicPressure, Eigen::Matrix3d& rot_stab_to_body)
	{
		//We need to rotate aerodynamic forces from the stability frame to the body frame via a rotation matrix
		Eigen::Vector3d aerodynamic_forces_body = rot_stab_to_body * aerodynamic_forces_stability_axis(alpha, beta, airSpeed, dynamicPressure);

		return aerodynamic_forces_body;

	}

	//Moments about aircraft center of gravity expressed in the body frame
	Eigen::Vector3d cg_moments_body_frame(Eigen::Matrix3d& rot_stab_to_body, Eigen::Vector3d& angular_vel_be_body_frame, double airSpeed, double dynamicPressure, double alpha, double beta, double g)
	{
		double downwash_angle = depsda * (alpha - alpha_initial);
		//Calculate moments in body frame about center of gravity.
		double eta11 = -1.4 * beta;
		//std::cout << "ac alpha: " << downwash_angle << std::endl;
		double eta21 = -0.59 - (3.1 * (tail_planform_area * lt) / (wing_planform_area * mean_aerodynamic_chord)) * (alpha - downwash_angle);
		double eta31 = (1 - alpha * (180 / (15 * pi))) * beta;
		Eigen::Vector3d eta({ eta11, eta21, eta31 });


		Eigen::MatrixXd correctionVec(3, 3);
		correctionVec << -11, 0, 5,
			0, (-4.03 * (tail_planform_area * lt * lt) / (wing_planform_area * mean_aerodynamic_chord * mean_aerodynamic_chord)), 0,
			1.7, 0, -11.5;

		Eigen::MatrixXd dCMdx = (mean_aerodynamic_chord / airSpeed) * correctionVec;


		Eigen::MatrixXd dCMdu(3, 3);
		dCMdu << -0.6, 0, 0.22,
			0, (-3.1 * (tail_planform_area * lt) / (wing_planform_area * mean_aerodynamic_chord)), 0,
			0, 0, -0.63;


		Eigen::Vector3d uVec({ aileron_inp,stabilizer_inp,rudder_inp });


		Eigen::Vector3d CMac_b = eta + dCMdx * angular_vel_be_body_frame + dCMdu * uVec;
		//std::cout << " ac eta" << eta[1]<< std::endl;


		//Normalize to an aerodynamic moments
		Eigen::Vector3d MAac_b = CMac_b * dynamicPressure * wing_planform_area * mean_aerodynamic_chord;



		////Transfer moment to CG
		Eigen::Vector3d rcg_b({ x_cg_pos_Fm,y_cg_pos_Fm,z_cg_pos_Fm });
		Eigen::Vector3d rac_b({ x_aero_pos_Fm,y_aero_pos_Fm,z_aero_pos_Fm });
		Eigen::Vector3d crossProduct = aerodynamic_forces_body_frame(alpha, beta, airSpeed, dynamicPressure, rot_stab_to_body).cross(rcg_b - rac_b);

		Eigen::Vector3d MAcg_b = MAac_b + crossProduct;

		Eigen::Vector3d cg_aerodynamic_moments_body = (MAac_b + crossProduct);



		Eigen::Vector3d momentArmEngine1({ x_cg_pos_Fm - eng1_x_pos_Fm, eng1_y_pos_Fm - y_cg_pos_Fm, z_cg_pos_Fm - eng1_z_pos_Fm });
		Eigen::Vector3d momentArmEngine2({ x_cg_pos_Fm - eng2_x_pos_Fm, eng2_y_pos_Fm - y_cg_pos_Fm, z_cg_pos_Fm - eng2_z_pos_Fm });


		//Get thrust of each engine
		double eng1_thrust = throttle1_inp * mass * g;
		double eng2_thrust = throttle2_inp * mass * g;

		Eigen::Vector3d engine1_force_body({ eng1_thrust,0,0 });
		Eigen::Vector3d engine2_force_body({ eng2_thrust,0,0 });

		Eigen::Vector3d momentEng1 = momentArmEngine1.cross(engine1_force_body);
		Eigen::Vector3d momentEng2 = momentArmEngine2.cross(engine2_force_body);
		Eigen::Vector3d cg_engine_moments_body = momentEng1 + momentEng2;

		Eigen::Vector3d total_cg_moments_body = cg_aerodynamic_moments_body + cg_engine_moments_body;

		Eigen::Vector3d x4tox6dot = InverseInertiaMatrix * (total_cg_moments_body - angular_vel_be_body_frame.cross(InertiaMatrix * angular_vel_be_body_frame));

		return x4tox6dot;

	}

	//Euler angles using traditional euler kinematics
	Eigen::Vector3d euler_angles_euler_kinematics(Eigen::Vector3d& angular_vel_be_body_frame, double phi, double theta, double psi)
	{

		////Calculate phidot, thetadot, psidot
		Eigen::Matrix3d H_pi(3, 3);
		H_pi << 1, sin(phi)* tan(theta), cos(phi)* tan(theta),
			0, cos(phi), sin(phi),
			0, sin(phi) / cos(theta), cos(phi) / cos(theta);
		Eigen::Vector3d x7tox9dot = H_pi * angular_vel_be_body_frame;

		return x7tox9dot;

	}

	//Euler angles using quaternions - to be implemented
	Eigen::Vector3d euler_angles_quaternion_kinematics() {}

	double pi{};

	double u{}; //u
	double v{}; //v
	double w{};//w
	double p{}; //p
	double q{}; //q
	double r{}; //r
	double phi{}; //phi
	double theta{}; //theta
	double psi{}; //psi


	//controls
	double aileron_inp{};
	double stabilizer_inp{};
	double rudder_inp{};
	double throttle1_inp{};
	double throttle2_inp{};


	//aircraft specs
	double number_of_engines = 2.0;
	double mass{};
	double mean_aerodynamic_chord{};
	double lt{};
	double wing_planform_area{};
	double tail_planform_area{};


	double x_cg_pos_Fm{};
	double y_cg_pos_Fm{};
	double z_cg_pos_Fm{};

	double x_aero_pos_Fm{};
	double y_aero_pos_Fm{};
	double z_aero_pos_Fm{};

	Eigen::Matrix3d InertiaMatrix{}; //inertia matrix needed for total_cg_moments_body
	Eigen::Matrix3d InverseInertiaMatrix{};


	//Engine constants
	double eng1_x_pos_Fm{};
	double eng1_y_pos_Fm{};
	double eng1_z_pos_Fm{};

	double eng2_x_pos_Fm{};
	double eng2_y_pos_Fm{};
	double eng2_z_pos_Fm{};


	//Aerodynamic constants - used for aerodynamic forces in stability axis
	double depsda{};
	double alpha_initial{};
	double slope_coefficient_of_lift{};
	double a3{};
	double a2{};
	double a1{};
	double a0{};
	double alpha_switch{}; //lift slope goes from linear to non-linear beyond this angle


	//Aerodynamic force coefficients
	double total_lift{};
	double total_drag{};
	double total_side_force{};

	double downwash_angle{};
	double angle_of_attack_tail{};
	double coefficient_of_lift_wb{};
	double coefficient_of_lift_tail{};


	//Vectors for forces in body axis
	Eigen::Vector3d aerodynamic_forces_body{};


	//Vectors for moments about center of gravity (cg)
	Eigen::Vector3d cg_aerodynamic_moments_body{};
	Eigen::Vector3d cg_engine_moments_body{};
	Eigen::Vector3d total_cg_moments_body{};


};


