
#pragma once

#include "AircraftModel.h"
#include "EngineModel.h"
#include <Eigen/Dense>
#include <numbers>
#include <map>
#include <cmath>
#include <memory>


// ============================================================
//  F16Aircraft — single-engine fighter, polynomial aero
//
//  Stability-derivative model in the style of CivilAircraft.
//  Stability axes (CL/CD), small-perturbation linear aero with
//  a parabolic drag polar. Coefficients are published F-16
//  values that appear in standard controls textbooks (Nelson,
//  Etkin & Reid, Roskam) for the early-block F-16 with a
//  stable CG (~0.30 cbar, i.e. not relaxed-stability).
//
//  Single ConstantMaxEngine on centerline (throttle index 3).
//  Engine thrust is along body-x through the CG → zero thrust
//  moment, which is the correct physics for a centerline jet.
//
//  Published mass/inertia/geometry (converted from S&L Appendix
//  A slug-ft^2 units to SI):
//    mass  = 9295.44 kg     (20500 lb / g)
//    Ixx   = 12874.8 kg·m^2
//    Iyy   = 75673.6 kg·m^2
//    Izz   = 85552.1 kg·m^2
//    Ixz   = 1331.4 kg·m^2
//    S     = 27.87 m^2      (300 ft^2)
//    cbar  = 3.45 m         (11.32 ft)
//    b     = 9.144 m        (30 ft)
//
//  Engine: T = throttle * T_max, T_max ≈ 76,310 N (~17,150 lbf,
//  F100 mil-power class). A faithful S&L thrust table with Mach,
//  altitude, and spool-up dynamics can replace ConstantMaxEngine
//  later without touching this class.
// ============================================================
class F16Aircraft : public AircraftModel
{
public:

    F16Aircraft(Eigen::VectorXd& states, Eigen::VectorXd& controls)
    {
        initializeStates(states);
        intializeControlInputs(controls);

        // ---- Mass & geometry (SI) ----
        mass = 9295.44;
        wing_planform_area = 27.87;
        mean_aerodynamic_chord = 3.45;
        wing_span = 9.144;

        x_cg_pos_Fm = 0.30 * mean_aerodynamic_chord;
        y_cg_pos_Fm = 0.0;
        z_cg_pos_Fm = 0.0;
        x_aero_pos_Fm = x_cg_pos_Fm;   // collocate; Cm coefficients absorb the offset
        y_aero_pos_Fm = 0.0;
        z_aero_pos_Fm = 0.0;

        // ---- Inertia tensor (kg·m²) ----
        const double Ixx = 12874.8;
        const double Iyy = 75673.6;
        const double Izz = 85552.1;
        const double Ixz = 1331.4;

        InertiaMatrix << Ixx, 0, -Ixz,
            0, Iyy, 0,
            -Ixz, 0, Izz;

        const double det = Ixx * Izz - Ixz * Ixz;
        InverseInertiaMatrix << Izz / det, 0, Ixz / det,
            0, 1.0 / Iyy, 0,
            Ixz / det, 0, Ixx / det;

        // ---- Longitudinal aero (stability axis, per radian) ----
        CL_0 = 0.0;
        CL_alpha = 4.6;
        CL_de = 0.5;     // stabilator (mapped to "stabilizer_inp")
        CD_0 = 0.02;
        k_drag = 0.10;    // parabolic polar: CD = CD_0 + k * CL^2

        Cm_0 = 0.0;
        Cm_alpha = -0.50;    // statically stable (textbook F-16 baseline)
        Cm_de = -1.5;     // strong all-moving tail
        Cm_q = -5.0;     // pitch damping

        // ---- Lateral / directional aero (stability axis, per radian) ----
        CY_beta = -1.0;
        CY_dr = 0.18;

        Cl_beta = -0.10;    // dihedral effect
        Cl_da = -0.10;    // sign convention matches CivilAircraft
        Cl_dr = 0.015;
        Cl_p = -0.45;    // roll damping
        Cl_r = 0.15;

        Cn_beta = 0.25;    // weathercock stability (positive)
        Cn_da = -0.02;    // adverse yaw
        Cn_dr = -0.10;
        Cn_p = -0.06;
        Cn_r = -0.30;    // yaw damping

        // ---- Single centerline engine ----
        const double T_max = 76310.0;  // N
        engines_.push_back(std::make_unique<ConstantMaxEngine>(
            Eigen::Vector3d(0.0, 0.0, 0.0),   // on CG: zero thrust moment
            /*throttle_index=*/3,
            T_max));
    }


    const char* name() const override { return "F-16"; }


    // ---- Aerodynamic forces in stability axis ----
    Eigen::Vector3d aerodynamic_forces_stability_axis(double alpha, double beta,
        double airSpeed, double dynamicPressure) override
    {
        // Lift
        double CL = CL_0 + CL_alpha * alpha + CL_de * stabilizer_inp;

        // Drag (parabolic polar)
        double CD = CD_0 + k_drag * CL * CL;

        // Side force
        double CY = CY_beta * beta + CY_dr * rudder_inp;

        return Eigen::Vector3d(
            -CD * dynamicPressure * wing_planform_area,
            CY * dynamicPressure * wing_planform_area,
            -CL * dynamicPressure * wing_planform_area);
    }


    // ---- Moments about CG in body frame ----
    Eigen::Vector3d cg_moments_body_frame(
        Eigen::Matrix3d& rot_stab_to_body,
        Eigen::Vector3d& omega_b,
        double airSpeed, double dynamicPressure,
        double alpha, double beta, double g) override
    {
        const double p_body = omega_b[0];
        const double q_body = omega_b[1];
        const double r_body = omega_b[2];

        const double V_safe = (airSpeed > 1e-3) ? airSpeed : 1e-3;
        const double pbar = p_body * wing_span / (2.0 * V_safe);
        const double qbar = q_body * mean_aerodynamic_chord / (2.0 * V_safe);
        const double rbar = r_body * wing_span / (2.0 * V_safe);

        // Roll
        const double Cl_total = Cl_beta * beta
            + Cl_da * aileron_inp
            + Cl_dr * rudder_inp
            + Cl_p * pbar
            + Cl_r * rbar;
        // Pitch
        const double Cm_total = Cm_0
            + Cm_alpha * alpha
            + Cm_de * stabilizer_inp
            + Cm_q * qbar;
        // Yaw
        const double Cn_total = Cn_beta * beta
            + Cn_da * aileron_inp
            + Cn_dr * rudder_inp
            + Cn_p * pbar
            + Cn_r * rbar;

        const double qS = dynamicPressure * wing_planform_area;
        Eigen::Vector3d MAac_b;
        MAac_b[0] = qS * wing_span * Cl_total;
        MAac_b[1] = qS * mean_aerodynamic_chord * Cm_total;
        MAac_b[2] = qS * wing_span * Cn_total;

        // Aero center collocated with CG → no aero-force-to-CG moment arm

        // Single centerline engine: moment is zero (arm is zero), but
        // route through the base loop for uniformity.
        Eigen::Vector3d M_eng = engine_moments_body_frame(g);

        Eigen::Vector3d M_total = MAac_b + M_eng;
        return InverseInertiaMatrix * (M_total - omega_b.cross(InertiaMatrix * omega_b));
    }


    std::map<std::string, double> getAircraftSpecs() const override
    {
        return {
            {"mass", mass},
            {"mean_aerodynamic_chord", mean_aerodynamic_chord},
            {"wing_span", wing_span},
            {"wing_planform_area", wing_planform_area}
        };
    }


    // ---- Geometry ----
    double mean_aerodynamic_chord{}, wing_span{};
    double wing_planform_area{};
    double x_cg_pos_Fm{}, y_cg_pos_Fm{}, z_cg_pos_Fm{};
    double x_aero_pos_Fm{}, y_aero_pos_Fm{}, z_aero_pos_Fm{};

    // ---- Stability derivatives (per radian) ----
    double CL_0{}, CL_alpha{}, CL_de{};
    double CD_0{}, k_drag{};
    double CY_beta{}, CY_dr{};
    double Cl_beta{}, Cl_da{}, Cl_dr{}, Cl_p{}, Cl_r{};
    double Cm_0{}, Cm_alpha{}, Cm_de{}, Cm_q{};
    double Cn_beta{}, Cn_da{}, Cn_dr{}, Cn_p{}, Cn_r{};
};