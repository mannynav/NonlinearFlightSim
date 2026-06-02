
#pragma once

#include "AircraftModel.h"
#include "EngineModel.h"
#include "LookupTable.h"
#include "GTMAeroTables.h"
#include <Eigen/Dense>
#include <numbers>
#include <map>
#include <cmath>
#include <memory>


// ============================================================
//  NasaGTMAircraft — NASA Generic Transport Model T2
//
//  Twin GravityScaledEngine units (throttle indices 3 and 4) via
//  the base-class engines_ vector. Aerodynamic coefficients come
//  from the real NASA T2 polynomial aerodatabase through an
//  8-table lookup moment buildup (3 static 2D + 5 rate 1D).
//
//  Published: mass, wing area, MAC, Iyy (arXiv:2106.08850).
//  Estimated: Ixx, Izz, Ixz, span, geometry offsets.
// ============================================================

class NasaGTMAircraft : public AircraftModel
{
public:

    NasaGTMAircraft(Eigen::VectorXd& states, Eigen::VectorXd& controls)
    {
        initializeStates(states);
        intializeControlInputs(controls);

        // ---- PUBLISHED geometry & mass ----
        mass = 26.190;
        wing_planform_area = 0.550;
        mean_aerodynamic_chord = 0.280;
        wing_span = 2.0;

        // ---- ESTIMATED geometry ----
        lt = 1.00;
        tail_planform_area = 0.14;

        x_cg_pos_Fm = 0.25 * mean_aerodynamic_chord;
        y_cg_pos_Fm = 0.0;
        z_cg_pos_Fm = 0.10 * mean_aerodynamic_chord;
        x_aero_pos_Fm = 0.13 * mean_aerodynamic_chord;
        y_aero_pos_Fm = 0.0;
        z_aero_pos_Fm = 0.0;

        eng1_x_pos_Fm = 0.0;  eng1_y_pos_Fm = -0.44;  eng1_z_pos_Fm = -0.10;
        eng2_x_pos_Fm = 0.0;  eng2_y_pos_Fm = 0.44;  eng2_z_pos_Fm = -0.10;

        // ---- Inertia tensor (kg·m²) ----
        const double Ixx = 1.32;
        const double Iyy = 5.768;
        const double Izz = 6.928;
        const double Ixz = 0.0784;

        InertiaMatrix << Ixx, 0, -Ixz,
            0, Iyy, 0,
            -Ixz, 0, Izz;

        const double det = Ixx * Izz - Ixz * Ixz;
        InverseInertiaMatrix << Izz / det, 0, Ixz / det,
            0, 1.0 / Iyy, 0,
            Ixz / det, 0, Ixx / det;

        // ---- Static aero (2D) ----
        cl_table_ = gtm_aero::build_cl_table();
        cd_table_ = gtm_aero::build_cd_table();
        cm_table_ = gtm_aero::build_cm_table();
        cy_table_ = gtm_aero::build_cy_table();
        cl_roll_table_ = gtm_aero::build_cl_roll_table();
        cn_table_ = gtm_aero::build_cn_table();

        // ---- Rate damping (1D in alpha) ----
        cm_q_table_ = gtm_aero::build_cm_q_table();
        cl_p_table_ = gtm_aero::build_cl_p_table();
        cl_r_table_ = gtm_aero::build_cl_r_table();
        cn_p_table_ = gtm_aero::build_cn_p_table();
        cn_r_table_ = gtm_aero::build_cn_r_table();

        // ---- Engines (two, gravity-scaled) ----
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


    const char* name() const override { return "NasaGTM-T2"; }


    // ---- Aerodynamic forces in stability axis ----
    Eigen::Vector3d aerodynamic_forces_stability_axis(
        double alpha, double beta, double airSpeed, double dynamicPressure) override
    {
        const double CL_static = cl_table_(alpha, stabilizer_inp);
        const double CD = cd_table_(alpha, stabilizer_inp);
        const double CY = cy_table_(beta, rudder_inp);

        const double tail_factor = 3.2 * (tail_planform_area / wing_planform_area);
        const double CL_q_damp = tail_factor * 1.3 * q * lt / airSpeed;
        const double CL = CL_static + CL_q_damp;

        return Eigen::Vector3d(
            -CD * dynamicPressure * wing_planform_area,
            CY * dynamicPressure * wing_planform_area,
            -CL * dynamicPressure * wing_planform_area);
    }


    // ---- Moments about CG in body frame (8-table buildup) ----
    Eigen::Vector3d cg_moments_body_frame(
        Eigen::Matrix3d& rot_stab_to_body,
        Eigen::Vector3d& omega_b,
        double airSpeed, double dynamicPressure,
        double alpha, double beta, double g) override
    {
        const double p_body = omega_b[0];
        const double q_body = omega_b[1];
        const double r_body = omega_b[2];

        const double c = mean_aerodynamic_chord;
        const double b = wing_span;
        const double S = wing_planform_area;

        const double V_safe = (airSpeed > 1e-3) ? airSpeed : 1e-3;
        const double pbar = p_body * b / (2.0 * V_safe);
        const double qbar = q_body * c / (2.0 * V_safe);
        const double rbar = r_body * b / (2.0 * V_safe);

        // Rolling
        const double Cl_total = cl_roll_table_(beta, aileron_inp)
            + cl_p_table_(alpha) * pbar
            + cl_r_table_(alpha) * rbar
            + 0.22 * rudder_inp;
        // Pitching
        const double Cm_total = cm_table_(alpha, stabilizer_inp)
            + cm_q_table_(alpha) * qbar;
        // Yawing
        const double Cn_total = cn_table_(beta, rudder_inp)
            + cn_p_table_(alpha) * pbar
            + cn_r_table_(alpha) * rbar;

        const double qS = dynamicPressure * S;
        Eigen::Vector3d MAac_b;
        MAac_b[0] = qS * b * Cl_total;
        MAac_b[1] = qS * c * Cm_total;
        MAac_b[2] = qS * b * Cn_total;

        Eigen::Vector3d rcg(x_cg_pos_Fm, y_cg_pos_Fm, z_cg_pos_Fm);
        Eigen::Vector3d rac(x_aero_pos_Fm, y_aero_pos_Fm, z_aero_pos_Fm);
        Eigen::Vector3d MAcg_b = MAac_b
            + aerodynamic_forces_body_frame(alpha, beta, airSpeed, dynamicPressure,
                rot_stab_to_body).cross(rcg - rac);

        // Engine moments from the base-class loop
        Eigen::Vector3d M_eng = engine_moments_body_frame(g);

        Eigen::Vector3d M_total = MAcg_b + M_eng;
        return InverseInertiaMatrix * (M_total - omega_b.cross(InertiaMatrix * omega_b));
    }


    std::map<std::string, double> getAircraftSpecs() const override
    {
        return {
            {"mass", mass},
            {"mean_aerodynamic_chord", mean_aerodynamic_chord},
            {"wing_span", wing_span},
            {"lt", lt},
            {"wing_planform_area", wing_planform_area},
            {"tail_planform_area", tail_planform_area},
            {"x_cg_pos_Fm", x_cg_pos_Fm}, {"y_cg_pos_Fm", y_cg_pos_Fm}, {"z_cg_pos_Fm", z_cg_pos_Fm},
            {"x_aero_pos_Fm", x_aero_pos_Fm}, {"y_aero_pos_Fm", y_aero_pos_Fm}, {"z_aero_pos_Fm", z_aero_pos_Fm},
            {"eng1_x_pos_Fm", eng1_x_pos_Fm}, {"eng1_y_pos_Fm", eng1_y_pos_Fm}, {"eng1_z_pos_Fm", eng1_z_pos_Fm},
            {"eng2_x_pos_Fm", eng2_x_pos_Fm}, {"eng2_y_pos_Fm", eng2_y_pos_Fm}, {"eng2_z_pos_Fm", eng2_z_pos_Fm}
        };
    }


    // ---- Geometry ----
    double mean_aerodynamic_chord{}, wing_span{}, lt{};
    double wing_planform_area{}, tail_planform_area{};
    double x_cg_pos_Fm{}, y_cg_pos_Fm{}, z_cg_pos_Fm{};
    double x_aero_pos_Fm{}, y_aero_pos_Fm{}, z_aero_pos_Fm{};
    double eng1_x_pos_Fm{}, eng1_y_pos_Fm{}, eng1_z_pos_Fm{};
    double eng2_x_pos_Fm{}, eng2_y_pos_Fm{}, eng2_z_pos_Fm{};

    // ---- Static aero (2D lookup tables) ----
    LookupTable<2> cl_table_;
    LookupTable<2> cd_table_;
    LookupTable<2> cm_table_;
    LookupTable<2> cy_table_;
    LookupTable<2> cl_roll_table_;
    LookupTable<2> cn_table_;

    // ---- Rate damping (1D lookup tables in alpha) ----
    LookupTable<1> cm_q_table_;
    LookupTable<1> cl_p_table_;
    LookupTable<1> cl_r_table_;
    LookupTable<1> cn_p_table_;
    LookupTable<1> cn_r_table_;
};