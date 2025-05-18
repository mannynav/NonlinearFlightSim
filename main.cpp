

#include <Eigen/Dense>
#include <numbers>
#include <iostream>
#include <fstream>

using Eigen::VectorXd;
using Eigen::MatrixXd;

class Aircraft
{
public:


private:


};

class CivilAircraft : Aircraft
{
public:

	//Constructor for initialization of aircraft constants.
	CivilAircraft(Eigen::VectorXd&& states, Eigen::VectorXd& controls) : pi(std::numbers::pi)
	{

		u = states[0]; //u
		v = states[1]; //v
		w = states[2]; //w
		p = states[3]; //p
		q = states[4]; //q
		r = states[5]; //r
		phi = states[6]; //phi
		theta = states[7]; //theta
		psi = states[8]; //psi

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

	//Force vector
	Eigen::Vector3d forces_body_frame(double airSpeed, double alpha, double beta, double dynamicPressure, double g, double q, double theta, double phi)
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


		//Calculate total lift
		double total_lift = coefficient_of_lift_wb + coefficient_of_lift_tail;


		//Total drag (not including tail)
		double total_drag = 0.13 + 0.07 * pow((5.5 * alpha + 0.654), 2);


		//Side force
		double side_force = -1.6 * beta + 0.24 * rudder_inp;

		Eigen::Vector3d aerodynamic_forces_stability({ -total_drag * dynamicPressure * wing_planform_area, side_force * dynamicPressure * wing_planform_area, -total_lift * dynamicPressure * wing_planform_area });



		//We need to rotate aerodynamic forces from the stability frame to the body frame via a rotation matrix
		Eigen::MatrixXd rot_stab_to_body(3, 3);
		rot_stab_to_body << cos(alpha), 0, -sin(alpha),
			0, 1, 0,
			sin(alpha), 0, cos(alpha);

		aerodynamic_forces_body = rot_stab_to_body * aerodynamic_forces_stability;



		//Engine forces
		double engine1_thrust = throttle1_inp * mass * g;
		double engine2_thrust = throttle2_inp * mass * g;

		engine1_force_body[0] = engine1_thrust; //x
		engine1_force_body[1] = 0.0; //y
		engine1_force_body[2] = 0.0; //z


		engine2_force_body[0] = engine2_thrust; //x
		engine2_force_body[1] = 0.0; //y
		engine2_force_body[2] = 0.0; //z


		engine_forces_body = engine1_force_body + engine2_force_body;



		////Calculate gravitational forces in body frame. This does not cause a moment.
		Eigen::Vector3d gravity_body({ -g * sin(theta), g * cos(theta) * sin(phi), g * cos(theta) * cos(phi) });
		gravity_forces_body = mass * gravity_body;



		return (1/mass)*(aerodynamic_forces_body + engine_forces_body + gravity_forces_body);
	}

	//Moments about aircraft center of gravity expressed in the body frame
	Eigen::Vector3d cg_moments_body_frame(Eigen::Vector3d& angular_vel_be_body_frame, double airSpeed, double dynamicPressure, double alpha, double beta)
	{

		//Calculate moments in body frame about center of gravity.
		double eta11 = -1.4 * beta;
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


		//Normalize to an aerodynamic moments
		Eigen::Vector3d MAac_b = CMac_b * dynamicPressure * wing_planform_area * mean_aerodynamic_chord;


		////Transfer moment to CG
		Eigen::Vector3d rcg_b({ x_cg_pos_Fm,y_cg_pos_Fm,z_cg_pos_Fm });
		Eigen::Vector3d rac_b({ x_aero_pos_Fm,y_aero_pos_Fm,z_aero_pos_Fm });
		Eigen::Vector3d crossProduct = aerodynamic_forces_body.cross(rcg_b - rac_b);

		Eigen::Vector3d MAcg_b = MAac_b + crossProduct;

		cg_aerodynamic_moments_body = (MAac_b + crossProduct);



		//Engine moment due to offset from CG
		momentArmEngine1[0] = x_cg_pos_Fm - eng1_x_pos_Fm;
		momentArmEngine1[1] = eng1_y_pos_Fm - y_cg_pos_Fm;
		momentArmEngine1[2] = z_cg_pos_Fm - eng1_z_pos_Fm;

		momentArmEngine2[0] = x_cg_pos_Fm - eng2_x_pos_Fm;
		momentArmEngine2[1] = eng2_y_pos_Fm - y_cg_pos_Fm;
		momentArmEngine2[2] = z_cg_pos_Fm - eng2_z_pos_Fm;


		Eigen::Vector3d momentEng1 = momentArmEngine1.cross(engine1_force_body);
		Eigen::Vector3d momentEng2 = momentArmEngine2.cross(engine2_force_body);
		cg_engine_moments_body = momentEng1 + momentEng2;


		total_cg_moments_body = cg_aerodynamic_moments_body + cg_engine_moments_body;


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

	}

	//Euler angles using quaternions - to be implemented
	Eigen::Vector3d euler_angles_quaternion_kinematics() {}

	double pi{};

	//states
	double u{}; //u
	double v{}; //v
	double w{}; //w
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


	//Forces from engine1 and engine2
	Eigen::Vector3d engine1_force_body{};
	Eigen::Vector3d engine2_force_body{};

	//Vectors for forces in body axis
	Eigen::Vector3d gravity_forces_body{};
	Eigen::Vector3d engine_forces_body{};
	Eigen::Vector3d aerodynamic_forces_stability{};
	Eigen::Vector3d aerodynamic_forces_body{};
	Eigen::Vector3d total_forces_body{};

	
	//Moment arms
	Eigen::Vector3d momentArmEngine1{};
	Eigen::Vector3d momentArmEngine2{};


	//Vectors for moments about center of gravity (cg)
	Eigen::Vector3d cg_aerodynamic_moments_body{};
	Eigen::Vector3d cg_engine_moments_body{};
	Eigen::Vector3d total_cg_moments_body{};


};

Eigen::VectorXd Aircraft_Sim(CivilAircraft& aircraft, Eigen::VectorXd&& X, Eigen::VectorXd& U)
{

	//aircraft is not used yet.
	
	double pi = std::numbers::pi;

	double u = X[0]; //u
	double v = X[1]; //v
	double w = X[2]; //w
	double p = X[3]; //p
	double q = X[4]; //q
	double r = X[5]; //r
	double phi = X[6]; //phi
	double theta = X[7]; //theta
	double psi = X[8]; //psi

	double aileron_inp = U[0];
	double stabilizer_inp = U[1];
	double rudder_inp = U[2];
	double throttle1_inp = U[3];
	double throttle2_inp = U[4];

	////////////////////// Constants ///////////////////

	double mass = 120000;
	double mean_aerodynamic_chord = 6.6;
	double tail_body_AC_distance = 24.8;
	double wing_planform_area = 260;
	double tail_planform_area = 64;

	double x_cg_pos_Fm = 0.23 * mean_aerodynamic_chord;
	double y_cg_pos_Fm = 0.0;
	double z_cg_pos_Fm = 0.10 * mean_aerodynamic_chord;

	double x_aero_pos_Fm = 0.12 * mean_aerodynamic_chord;
	double y_aero_pos_Fm = 0;
	double z_aero_pos_Fm = 0;


	//Engine constants
	double eng1_x_pos_Fm = 0;
	double eng1_y_pos_Fm = -7.94;
	double eng1_z_pos_Fm = -1.9;

	double eng2_x_pos_Fm = 0.0;
	double eng2_y_pos_Fm = 7.94;
	double eng2_z_pos_Fm = -1.9;


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



	//////////////// Aerodynamic Force Coefficients in stability axis /////////////////////////
	
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
	double angle_of_attack_tail = alpha - downwash_angle + stabilizer_inp + 1.3 * q * tail_body_AC_distance / airSpeed;
	double coefficient_of_lift_tail = 3.1 * (tail_planform_area / wing_planform_area) * angle_of_attack_tail;


	//Calculate total lift
	double total_lift = coefficient_of_lift_wb + coefficient_of_lift_tail;
	

	//Total drag (not including tail)
	double total_drag = 0.13 + 0.07 * pow((5.5 * alpha + 0.654), 2);


	//Side force
	double side_force = -1.6 * beta + 0.24 * rudder_inp;
	
	Eigen::Vector3d aerodynamic_forces_stability({ -total_drag * dynamicPressure * wing_planform_area, side_force * dynamicPressure * wing_planform_area, -total_lift * dynamicPressure * wing_planform_area });
	

	//We need to rotate aerodynamic forces from the stability frame to the body frame via a rotation matrix
	Eigen::MatrixXd rot_stab_to_body(3, 3);

	rot_stab_to_body << cos(alpha), 0, -sin(alpha),
					0, 1, 0,
					sin(alpha), 0, cos(alpha);

	Eigen::Vector3d aerodynamic_forces_body = rot_stab_to_body * aerodynamic_forces_stability;


	//Calculate moments in body frame about center of gravity.
	double eta11 = -1.4 * beta;
	double eta21 = -0.59 - (3.1 * (tail_planform_area * tail_body_AC_distance) / (wing_planform_area * mean_aerodynamic_chord)) * (alpha - downwash_angle);
	double eta31 = (1 - alpha * (180 / (15 * pi))) * beta;
	Eigen::Vector3d eta({ eta11, eta21, eta31 });
	

	Eigen::MatrixXd correctionVec(3, 3);
	correctionVec << -11, 0, 5,
		0, (-4.03 * (tail_planform_area * tail_body_AC_distance * tail_body_AC_distance) / (wing_planform_area * mean_aerodynamic_chord * mean_aerodynamic_chord)), 0,
		1.7, 0, -11.5;

	Eigen::MatrixXd dCMdx = (mean_aerodynamic_chord / airSpeed) * correctionVec;



	Eigen::MatrixXd dCMdu(3, 3);
	dCMdu << -0.6, 0, 0.22,
		0, (-3.1 * (tail_planform_area * tail_body_AC_distance) / (wing_planform_area * mean_aerodynamic_chord)), 0,
		0, 0, -0.63;
	

	Eigen::Vector3d uVec({ aileron_inp,stabilizer_inp,rudder_inp });


	Eigen::Vector3d CMac_b = eta + dCMdx * angular_vel_be_body_frame + dCMdu * uVec;


	//Normalize to an aerodynamic moments
	Eigen::Vector3d MAac_b = CMac_b * dynamicPressure * wing_planform_area * mean_aerodynamic_chord;


	////Transfer moment to CG
	Eigen::Vector3d rcg_b({ x_cg_pos_Fm,y_cg_pos_Fm,z_cg_pos_Fm });
	Eigen::Vector3d rac_b({ x_aero_pos_Fm,y_aero_pos_Fm,z_aero_pos_Fm });
	Eigen::Vector3d crossProduct = aerodynamic_forces_body.cross(rcg_b - rac_b);

	Eigen::Vector3d MAcg_b = MAac_b + crossProduct;


	//Get thrust of each engine
	double eng1_thrust = throttle1_inp * mass * g;
	double eng2_thrust = throttle2_inp * mass * g;

	Eigen::Vector3d FE1_b({ eng1_thrust,0,0 });
	Eigen::Vector3d FE2_b({ eng2_thrust,0,0 });
	Eigen::Vector3d engine_forces_body = FE1_b + FE2_b;
	
	

	//Engine moment due to offset from CG
	Eigen::Vector3d momentArmEng1({ x_cg_pos_Fm - eng1_x_pos_Fm, eng1_y_pos_Fm - y_cg_pos_Fm, z_cg_pos_Fm - eng1_z_pos_Fm });
	Eigen::Vector3d momentArmEng2({ x_cg_pos_Fm - eng2_x_pos_Fm, eng2_y_pos_Fm - y_cg_pos_Fm, z_cg_pos_Fm - eng2_z_pos_Fm });

	Eigen::Vector3d momentEng1 = momentArmEng1.cross(FE1_b);
	Eigen::Vector3d momentEng2 = momentArmEng2.cross(FE2_b);

	Eigen::Vector3d momentsFromEngineCG_body = momentEng1 + momentEng2;



	////Calculate gravitational forces in body frame. This does not cause a moment.
	Eigen::Vector3d gravity_body({ -g * sin(theta), g * cos(theta) * sin(phi), g * cos(theta) * cos(phi) });
	Eigen::Vector3d gravity_forces_body = mass * gravity_body;


	//Inertia matrix
	Eigen::Matrix3d inertiaMatrix(3, 3);
	inertiaMatrix << 40.07, 0, -2.0923,
		             0, 64, 0,
		            -2.0923, 0, 99.2;
	inertiaMatrix = mass * inertiaMatrix;
	


	Eigen::Matrix3d invInertiaMatrix(3, 3);
	invInertiaMatrix << 0.02498, 0, 0.000523,
		0, 0.015625, 0,
		0.0005232, 0, 0.010019;
	invInertiaMatrix = (1.0 / mass) * invInertiaMatrix;



	////From F_b (all forces resolved in Fb) and calculate the udot, vdot, wdot
	Eigen::Vector3d forces_body = gravity_forces_body + engine_forces_body + aerodynamic_forces_body;
	Eigen::Vector3d x1tox3dot = (1 / mass) * forces_body - angular_vel_be_body_frame.cross(translational_vel_body_frame);
	


	//From Mcg_b (all moments about the cog in Fb) and calculate pdot, qdot,rdot
	Eigen::Vector3d cg_moments_body = MAcg_b + momentsFromEngineCG_body;

	Eigen::Vector3d x4tox6dot = invInertiaMatrix*(cg_moments_body - angular_vel_be_body_frame.cross(inertiaMatrix * angular_vel_be_body_frame));


	//Calculate phidot, thetadot, psidot
	Eigen::Matrix3d H_pi(3, 3);
	H_pi << 1, sin(phi)* tan(theta), cos(phi)* tan(theta),
		0, cos(phi), sin(phi),
		0, sin(phi) / cos(theta), cos(phi) / cos(theta);
	Eigen::Vector3d x7tox9dot = H_pi * angular_vel_be_body_frame;


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

Eigen::VectorXd RCAM_model_implicit(CivilAircraft& ac, Eigen::VectorXd& XDOT, Eigen::VectorXd& X, Eigen::VectorXd& U) {

	return Aircraft_Sim(ac, std::move(X), U) - XDOT;
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

	dxdot.setConstant(0.0000000001);


	dx.setConstant(0.0000000001);


	du.setConstant(0.0000000001);

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

			Eigen::VectorXd F = RCAM_model_implicit(ac,XDOTo, x_plus, Uo);
			double F_plus_keep = F[i];

			F = RCAM_model_implicit(ac,XDOTo, x_minus, Uo);
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

			Eigen::VectorXd F = RCAM_model_implicit(ac,XDOTo, Xo, u_plus);
			double F_plus_keep = F[i];

			F = RCAM_model_implicit(ac,XDOTo, Xo, u_minus);
			double F_minus_keep = F[i];

			B(i, j) = (F_plus_keep - F_minus_keep) / (2 * du);

		}
	}

	return B;

}

Eigen::VectorXd LinearStateSpace(const Eigen::VectorXd& x, const Eigen::VectorXd& u, const Eigen::MatrixXd& A, const Eigen::MatrixXd& B)
{
	return A * x + B * u;

}

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






int main()

{
	

	Eigen::VectorXd initialX(9);
	initialX[0] = 85;
	initialX[1] = 0.0;
	initialX[2] = 0.0;
	initialX[3] = 0.0;
	initialX[4] = 0.0;
	initialX[5] = 0.0;
	initialX[6] = 0.0;
	initialX[7] = 0.1;
	initialX[8] = 0.0;


	Eigen::VectorXd initialU(5);
	initialU[0] = 0.0;
	initialU[1] = -0.1;
	initialU[2] = 0.0;
	initialU[3] = 0.08;
	initialU[4] = 0.08;


	CivilAircraft ac = CivilAircraft(std::move(initialX), initialU);

	double simLength = 150;
	double steps = 200;
	double timeIncrement = simLength / steps;

	Eigen::MatrixXd nonLinearSolutionMatrix(9, 201);
	nonLinearSolutionMatrix.setZero();
	nonLinearSolutionMatrix.col(0) = initialX;


	for (int i = 0; i < steps; i++)
	{
		Eigen::VectorXd k1 = Aircraft_Sim(ac,nonLinearSolutionMatrix.col(i), initialU);

		Eigen::VectorXd nextk2 = nonLinearSolutionMatrix.col(i) + (timeIncrement / 2.0) * k1;
		Eigen::VectorXd k2 = Aircraft_Sim(ac,std::move(nextk2), initialU);

		Eigen::VectorXd nextk3 = nonLinearSolutionMatrix.col(i) + (timeIncrement / 2.0) * k2;
		Eigen::VectorXd k3 = Aircraft_Sim(ac,std::move(nextk3), initialU);

		Eigen::VectorXd nextk4 = nonLinearSolutionMatrix.col(i) + timeIncrement * k3;
		Eigen::VectorXd k4 = Aircraft_Sim(ac,std::move(nextk4), initialU);

		nonLinearSolutionMatrix.col(i + 1) = nonLinearSolutionMatrix.col(i) + (timeIncrement / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

	}


	/*for (int i = 0; i < 10; i++)
	{
		Eigen::VectorXd nextStep = solution.col(i);
		solution.col(i + 1) = solution.col(i) + 0.100 * RCAM_model(std::move(nextStep), initialU);
	}*/


	// --- Write the non linear solution matrix to a CSV file ---
	std::ofstream outputFile("non linear solution.csv");
	if (outputFile.is_open()) {
		for (int i = 0; i < nonLinearSolutionMatrix.rows(); ++i) {
			for (int j = 0; j < nonLinearSolutionMatrix.cols(); ++j) {
				outputFile << nonLinearSolutionMatrix(i, j);
				if (j < nonLinearSolutionMatrix.cols() - 1) {
					outputFile << ","; // Add a comma as a separator
				}
			}
			outputFile << "\n"; // Add a newline character for the next row
		}
		outputFile.close();
		std::cout << "non linear solution written to solution.csv" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}

	std::cout << "----------------------------------------- End of non linear simulation -------------------------------------------------" << std::endl;




	std::cout << "----------------------------------------- Initializing trim point for Xo and Uo (Get from Matlab) -------------------------------------------------" << std::endl;

	//Vectors for states and inputs
	Eigen::VectorXd Xdoto(9);
	Eigen::VectorXd Xo(9);
	Eigen::VectorXd Uo(5);

	//Matrices for increments for the finite difference step
	Eigen::MatrixXd dxdot(9, 9);
	Eigen::MatrixXd dx(9, 9);
	Eigen::MatrixXd du(9, 5);

	initializeTrimStatesInputsAndIncrements(Xdoto, Xo, Uo, dxdot, dx, du);

	std::cout << "----------------------------------------- Linearize system for A -------------------------------------------------" << std::endl;
	Eigen::MatrixXd Atest = LinearizeSystem_A(ac,Xdoto, Xo, Uo, dxdot, dx, du);
	std::cout << Atest << std::endl;
	std::cout << Atest.eigenvalues() << std::endl;



	std::cout << "----------------------------------------- Linearize system for B -------------------------------------------------" << std::endl;
	Eigen::MatrixXd Btest = LinearizeSystem_B(ac,Xdoto, Xo, Uo, dxdot, dx, du);
	std::cout << Btest << std::endl;



	std::cout << "----------------------------------------- Start Longitudinal Model Process -------------------------------------------------" << std::endl;

	std::cout << "----------------------------------------- Initialize S for similarity transformation, compute transformed A and B -------------------------------------------------" << std::endl;
	Eigen::MatrixXd S(9, 9);
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

	std::cout << S << std::endl;

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

	std::cout << "B long: " << std::endl;
	std::cout << Blongitudinal << std::endl;

	Eigen::MatrixXd Clongitudinal(4, 4);
	Clongitudinal.setIdentity();

	Eigen::MatrixXd Dlongitudinal(4, 5);
	Dlongitudinal.setZero();


	std::cout << "--------------------- Initialize Xlong and Ulong for Longitudinal Model ------------------------------ " << std::endl;

	Eigen::VectorXd Xlong_o(4); //Phugoid.
	Xlong_o[0] = -1.99;
	Xlong_o[1] = 0.183;
	Xlong_o[2] = -0.0038;
	Xlong_o[3] = 0.0038;

	Eigen::VectorXd Ulong_o(5); //Initialize inputs. Constant of 0 for now.
	Ulong_o.setZero();



	std::cout << "--------------------- State space simulation For Longitudinal Model ------------------------------ " << std::endl;

	double simTime = 200;
	double stepsLong = 250;
	double increment = simTime / stepsLong;


	// Get the Dormand-Prince tableau
	RKTableauSS tableau_ss = dp45_tableau_ss();
	std::function<VectorXd(const VectorXd&, const VectorXd&, const MatrixXd&, const MatrixXd&)> ss_func = LinearStateSpace;

	Eigen::MatrixXd solution_linear_long = rk4_simulate_ss(ss_func, Xlong_o, Ulong_o, 0, 150, 150, Alongitudinal, Blongitudinal);
	//Eigen::MatrixXd solution_linear_long = adaptive_rk_simulate_ss(ss_func, Xlong_o, Ulong_o, 0, 1000, 500, 0.03, tableau_ss,Alongitudinal, Blongitudinal);
	//Eigen::MatrixXd solution_linear_long = forward_euler_simulate_ss(ss_func, Xlong_o, Ulong_o, 0, 150, 300, Alongitudinal, Blongitudinal);
	std::cout << solution_linear_long.leftCols(5) << std::endl;


	std::cout << "----------------------------------------- End of simulation for Longitudinal Model -------------------------------------------------" << std::endl;


	// --- Write the non linear solution matrix to a CSV file ---
	std::ofstream outFile("Linear_Longitudnal_Solution.csv");
	if (outFile.is_open()) {
		for (int i = 0; i < solution_linear_long.rows(); ++i) {
			for (int j = 0; j < solution_linear_long.cols(); ++j) {
				outFile << solution_linear_long(i, j);
				if (j < solution_linear_long.cols() - 1) {
					outFile << ","; // Add a comma as a separator
				}
			}
			outFile << "\n"; // Add a newline character for the next row
		}
		outFile.close();
		std::cout << "Linear Longitudinal solution matrix written to csv" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}


	//CivilAircraft ac = CivilAircraft();
	//std::cout << ac.InertiaMatrix << std::endl;


	return 0;

}