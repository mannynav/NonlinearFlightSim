
#pragma once

#include "AircraftModel.h"
#include "EngineModel.h"
#include <Eigen/Dense>
#include <numbers>
#include <map>
#include <cmath>
#include <memory>


// ============================================================
//  CivilAircraft — large twin-engine transport (~120,000 kg)
//
//  Two GravityScaledEngine units (throttle indices 3 and 4),
//  managed by the base-class engines_ vector. Engine force and
//  engine moment are computed by AircraftModel; this class only
//  implements aerodynamics and the aerodynamic moment buildup.
// ============================================================

class CivilAircraft : public AircraftModel
{
public:

    CivilAircraft(Eigen::VectorXd& states, Eigen::VectorXd& controls)
    {
        initializeStates(states);
        intializeControlInputs(controls);

        // ---- Aircraft specs ----
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

        eng1_x_pos_Fm = 0;     eng1_y_pos_Fm = -7.94;  eng1_z_pos_Fm = -1.9;
        eng2_x_pos_Fm = 0.0;   eng2_y_pos_Fm = 7.94;  eng2_z_pos_Fm = -1.9;

        InertiaMatrix << 40.07, 0, -2.0923,
            0, 64, 0,
            -2.0923, 0, 99.2;
        InertiaMatrix = mass * InertiaMatrix;

        InverseInertiaMatrix << 0.02498, 0, 0.000523,
            0, 0.015625, 0,
            0.0005232, 0, 0.010019;
        InverseInertiaMatrix = (1.0 / mass) * InverseInertiaMatrix;

        depsda = 0.25;
        alpha_initial = -11.5 * (pi / 180);
        slope_coefficient_of_lift = 5.5;
        a3 = -768.5;  a2 = 609.2;  a1 = -155.2;  a0 = -15.212;
        alpha_switch = 14.5 * (pi / 180);

        // ---- Engines ----
        // position_cg = (eng_x - cg_x, eng_y - cg_y, cg_z - eng_z)
        // The y and z components reproduce the original engine-moment.
        engines_.push_back(std::make_unique<GravityScaledEngine>(
            Eigen::Vector3d(eng1_x_pos_Fm - x_cg_pos_Fm,
                eng1_y_pos_Fm - y_cg_pos_Fm,
                z_cg_pos_Fm - eng1_z_pos_Fm),
            /*throttle_index=*/3));
        engines_.push_back(std::make_unique<GravityScaledEngine>(
            Eigen::Vector3d(eng2_x_pos_Fm - x_cg_pos_Fm,
                eng2_y_pos_Fm - y_cg_pos_Fm,
                z_cg_pos_Fm - eng2_z_pos_Fm),
            /*throttle_index=*/4));
    }


    const char* name() const override { return "CivilAircraft"; }


    // ---- Aerodynamic forces in stability axis ----
    Eigen::Vector3d aerodynamic_forces_stability_axis(double alpha, double beta,
        double airSpeed, double dynamicPressure) override
    {
        double CL_wb = (alpha <= alpha_switch)
            ? slope_coefficient_of_lift * (alpha - alpha_initial)
            : a3 * std::pow(alpha, 3) + a2 * std::pow(alpha, 2) + a1 * alpha + a0;

        double eps_dn = depsda * (alpha - alpha_initial);
        double alpha_t = alpha - eps_dn + stabilizer_inp + 1.3 * q * lt / airSpeed;
        double CL_t = 3.1 * (tail_planform_area / wing_planform_area) * alpha_t;
        double CL = CL_wb + CL_t;

        double CD = 0.13 + 0.07 * std::pow(5.5 * alpha + 0.654, 2);
        double CY = -1.6 * beta + 0.24 * rudder_inp;

        return Eigen::Vector3d(
            -CD * dynamicPressure * wing_planform_area,
            CY * dynamicPressure * wing_planform_area,
            -CL * dynamicPressure * wing_planform_area);
    }


    // ---- Moments about CG in body frame ----
    Eigen::Vector3d cg_moments_body_frame(Eigen::Matrix3d& rot_stab_to_body,
        Eigen::Vector3d& omega_b,
        double airSpeed, double dynamicPressure,
        double alpha, double beta, double g) override
    {
        double eps_dn = depsda * (alpha - alpha_initial);
        double eta11 = -1.4 * beta;
        double eta21 = -0.59 - (3.1 * (tail_planform_area * lt) /
            (wing_planform_area * mean_aerodynamic_chord)) * (alpha - eps_dn);
        double eta31 = (1 - alpha * (180 / (15 * pi))) * beta;
        Eigen::Vector3d eta(eta11, eta21, eta31);

        Eigen::MatrixXd corr(3, 3);
        corr << -11, 0, 5,
            0, -4.03 * (tail_planform_area * lt * lt) /
            (wing_planform_area * mean_aerodynamic_chord * mean_aerodynamic_chord), 0,
            1.7, 0, -11.5;
        Eigen::MatrixXd dCMdx = (mean_aerodynamic_chord / airSpeed) * corr;

        Eigen::MatrixXd dCMdu(3, 3);
        dCMdu << -0.6, 0, 0.22,
            0, -3.1 * (tail_planform_area * lt) /
            (wing_planform_area * mean_aerodynamic_chord), 0,
            0, 0, -0.63;

        Eigen::Vector3d uVec(aileron_inp, stabilizer_inp, rudder_inp);
        Eigen::Vector3d CMac_b = eta + dCMdx * omega_b + dCMdu * uVec;
        Eigen::Vector3d MAac_b = CMac_b * dynamicPressure * wing_planform_area * mean_aerodynamic_chord;

        Eigen::Vector3d rcg(x_cg_pos_Fm, y_cg_pos_Fm, z_cg_pos_Fm);
        Eigen::Vector3d rac(x_aero_pos_Fm, y_aero_pos_Fm, z_aero_pos_Fm);
        Eigen::Vector3d MAcg_b = MAac_b
            + aerodynamic_forces_body_frame(alpha, beta, airSpeed, dynamicPressure,
                rot_stab_to_body).cross(rcg - rac);

        // Engine moments now come from the base-class engine loop
        Eigen::Vector3d M_eng = engine_moments_body_frame(g);

        Eigen::Vector3d M_total = MAcg_b + M_eng;
        return InverseInertiaMatrix * (M_total - omega_b.cross(InertiaMatrix * omega_b));
    }


    std::map<std::string, double> getAircraftSpecs() const override
    {
        std::map<std::string, double> s;
        s["mass"] = mass;
        s["mean_aerodynamic_chord"] = mean_aerodynamic_chord;
        s["lt"] = lt;
        s["wing_planform_area"] = wing_planform_area;
        s["tail_planform_area"] = tail_planform_area;
        s["x_cg_pos_Fm"] = x_cg_pos_Fm; s["y_cg_pos_Fm"] = y_cg_pos_Fm; s["z_cg_pos_Fm"] = z_cg_pos_Fm;
        s["x_aero_pos_Fm"] = x_aero_pos_Fm; s["y_aero_pos_Fm"] = y_aero_pos_Fm; s["z_aero_pos_Fm"] = z_aero_pos_Fm;
        s["eng1_x_pos_Fm"] = eng1_x_pos_Fm; s["eng1_y_pos_Fm"] = eng1_y_pos_Fm; s["eng1_z_pos_Fm"] = eng1_z_pos_Fm;
        s["eng2_x_pos_Fm"] = eng2_x_pos_Fm; s["eng2_y_pos_Fm"] = eng2_y_pos_Fm; s["eng2_z_pos_Fm"] = eng2_z_pos_Fm;
        return s;
    }


    // ---- Airframe-specific data ----
    double mean_aerodynamic_chord{}, lt{};
    double wing_planform_area{}, tail_planform_area{};
    double x_cg_pos_Fm{}, y_cg_pos_Fm{}, z_cg_pos_Fm{};
    double x_aero_pos_Fm{}, y_aero_pos_Fm{}, z_aero_pos_Fm{};
    double eng1_x_pos_Fm{}, eng1_y_pos_Fm{}, eng1_z_pos_Fm{};
    double eng2_x_pos_Fm{}, eng2_y_pos_Fm{}, eng2_z_pos_Fm{};

    double depsda{}, alpha_initial{}, slope_coefficient_of_lift{};
    double a3{}, a2{}, a1{}, a0{}, alpha_switch{};
};