#pragma once

#include "AircraftModel.h"
#include "Aircraft.h"          // for CivilAircraft used by makeAircraft
#include "LookupTable.h"
#include "GTMAeroTables.h"
#include <Eigen/Dense>
#include <numbers>
#include <map>
#include <cmath>
#include <memory>
#include <stdexcept>


// ============================================================
//  NasaGTMAircraft — NASA Generic Transport Model T2
//
//  NASA GTM is a 5.5% dynamically-scaled twin-engine subscale
//  aircraft used by NASA Langley for upset / loss-of-control
//  research. The full Simulink simulation is open-source:
//
//      https://github.com/nasa/GTM_DesignSim
//
//  Mass, wing area, MAC, and Iyy below are PUBLISHED values
//  (Cunis et al., arXiv:1905.09794 Table 1; Lai/Cunis/Burlion,
//  arXiv:2106.08850 Table 2). Off-diagonal inertias are
//  engineering ESTIMATES.
//
//  Static aerodynamic coefficients (CL, CD, CY, Cl, Cm, Cn) come
//  from lookup tables in GTMAeroTables.h. Those table values are
//  currently hand-authored transport-class placeholders — replace
//  with real GTM data extracted from MATLAB once available.
//
//  Dynamic derivatives (rate damping) are still analytical
//  approximations in the `corr` matrix below; a higher-fidelity
//  implementation would store them as separate scalars or as
//  1-D tables of α.
//
//  IMPORTANT: trim_states.csv (built for the 120,000 kg
//  CivilAircraft at 85 m/s @ 10,000 ft) will NOT trim this
//  26 kg subscale aircraft. Use the in-C++ trim solver, or
//  generate a GTM-specific trim CSV.
// ============================================================
class NasaGTMAircraft : public AircraftModel
{
public:

    NasaGTMAircraft(Eigen::VectorXd& states, Eigen::VectorXd& controls)
    {
        initializeStates(states);
        intializeControlInputs(controls);

        // ---- PUBLISHED GTM-T2 geometry & mass ----
        mass = 26.190;   // kg     (arXiv:2106.08850)
        wing_planform_area = 0.550;    // m^2    (arXiv:2106.08850)
        mean_aerodynamic_chord = 0.280;    // m      (arXiv:2106.08850)

        // ---- ESTIMATED geometry ----
        lt = 1.00;         // tail-body AC distance (~ 0.4 · length)
        tail_planform_area = 0.14;         // ~25% of wing area

        x_cg_pos_Fm = 0.25 * mean_aerodynamic_chord;
        y_cg_pos_Fm = 0.0;
        z_cg_pos_Fm = 0.10 * mean_aerodynamic_chord;
        x_aero_pos_Fm = 0.13 * mean_aerodynamic_chord;
        y_aero_pos_Fm = 0.0;
        z_aero_pos_Fm = 0.0;

        // Engine offsets, scaled from full-scale CivilAircraft by ~5.5%
        eng1_x_pos_Fm = 0.0;  eng1_y_pos_Fm = -0.44;  eng1_z_pos_Fm = -0.10;
        eng2_x_pos_Fm = 0.0;  eng2_y_pos_Fm = 0.44;  eng2_z_pos_Fm = -0.10;

        // ---- Inertia tensor (kg·m²) ----
        // Iyy is published; Ixx, Izz, Ixz are reasonable estimates.
        const double Ixx = 1.32;     // ESTIMATE
        const double Iyy = 5.768;    // PUBLISHED (arXiv:2106.08850)
        const double Izz = 6.928;    // ESTIMATE
        const double Ixz = 0.0784;   // ESTIMATE

        InertiaMatrix << Ixx, 0, -Ixz,
            0, Iyy, 0,
            -Ixz, 0, Izz;

        // Analytical inverse (Ixy = Iyz = 0)
        const double det = Ixx * Izz - Ixz * Ixz;
        InverseInertiaMatrix << Izz / det, 0, Ixz / det,
            0, 1.0 / Iyy, 0,
            Ixz / det, 0, Ixx / det;

        // ---- Build all six static aero lookup tables ----
        cl_table_ = gtm_aero::build_cl_table();        // CL(α, η)
        cd_table_ = gtm_aero::build_cd_table();        // CD(α, η)
        cm_table_ = gtm_aero::build_cm_table();        // Cm(α, η)
        cy_table_ = gtm_aero::build_cy_table();        // CY(β, δr)
        cl_roll_table_ = gtm_aero::build_cl_roll_table();   // Cl(β, δa)
        cn_table_ = gtm_aero::build_cn_table();        // Cn(β, δr)
    }


    const char* name() const override { return "NasaGTM-T2"; }


    // ============================================================
    //  Aerodynamic forces in stability axis
    //  Static coefficients from tables; pitch-rate damping analytical.
    // ============================================================
    Eigen::Vector3d aerodynamic_forces_stability_axis(
        double alpha, double beta, double airSpeed, double dynamicPressure) override
    {
        const double CL_static = cl_table_(alpha, stabilizer_inp);
        const double CD = cd_table_(alpha, stabilizer_inp);
        const double CY = cy_table_(beta, rudder_inp);

        // Pitch-rate damping (dynamic derivative — kept analytical for now)
        const double tail_factor = 3.2 * (tail_planform_area / wing_planform_area);
        const double CL_q_damp = tail_factor * 1.3 * q * lt / airSpeed;

        const double CL = CL_static + CL_q_damp;

        return Eigen::Vector3d(
            -CD * dynamicPressure * wing_planform_area,
            CY * dynamicPressure * wing_planform_area,
            -CL * dynamicPressure * wing_planform_area);
    }


    // ============================================================
    //  Engine forces — thrust along body x-axis from both engines
    // ============================================================
    Eigen::Vector3d engine_forces_body_frame(double g) override
    {
        const double T1 = throttle1_inp * mass * g;
        const double T2 = throttle2_inp * mass * g;
        return Eigen::Vector3d(T1 + T2, 0.0, 0.0);
    }


    // ============================================================
    //  Moments about CG in body frame
    //  Static moment coefficients from tables; rate damping analytical.
    // ============================================================
    Eigen::Vector3d cg_moments_body_frame(
        Eigen::Matrix3d& rot_stab_to_body,
        Eigen::Vector3d& omega_b,
        double airSpeed, double dynamicPressure,
        double alpha, double beta, double g) override
    {
        // Static moment coefficients from tables. Tables absorb the primary
        // control surface for each axis. Small cross-control terms (e.g. rudder
        // contribution to roll) added analytically.
        Eigen::Vector3d eta;
        eta[0] = cl_roll_table_(beta, aileron_inp) + 0.22 * rudder_inp;
        eta[1] = cm_table_(alpha, stabilizer_inp);
        eta[2] = cn_table_(beta, rudder_inp);

        // Rate damping (analytical placeholder — replace with published
        // GTM rate derivatives Cl_p, Cm_q, Cn_r, etc. when available)
        Eigen::MatrixXd corr(3, 3);
        corr << -12.0, 0, 5.5,
            0, -4.20 * (tail_planform_area * lt * lt) /
            (wing_planform_area * mean_aerodynamic_chord * mean_aerodynamic_chord), 0,
            1.8, 0, -12.0;
        Eigen::MatrixXd dCMdx = (mean_aerodynamic_chord / airSpeed) * corr;

        Eigen::Vector3d CMac_b = eta + dCMdx * omega_b;

        // Dimensionalize and transfer to CG
        Eigen::Vector3d MAac_b = CMac_b * dynamicPressure * wing_planform_area * mean_aerodynamic_chord;

        Eigen::Vector3d rcg(x_cg_pos_Fm, y_cg_pos_Fm, z_cg_pos_Fm);
        Eigen::Vector3d rac(x_aero_pos_Fm, y_aero_pos_Fm, z_aero_pos_Fm);
        Eigen::Vector3d MAcg_b = MAac_b
            + aerodynamic_forces_body_frame(alpha, beta, airSpeed,
                dynamicPressure, rot_stab_to_body)
            .cross(rcg - rac);

        // Engine thrust moments
        Eigen::Vector3d arm1(x_cg_pos_Fm - eng1_x_pos_Fm, eng1_y_pos_Fm - y_cg_pos_Fm, z_cg_pos_Fm - eng1_z_pos_Fm);
        Eigen::Vector3d arm2(x_cg_pos_Fm - eng2_x_pos_Fm, eng2_y_pos_Fm - y_cg_pos_Fm, z_cg_pos_Fm - eng2_z_pos_Fm);
        Eigen::Vector3d F1(throttle1_inp * mass * g, 0, 0);
        Eigen::Vector3d F2(throttle2_inp * mass * g, 0, 0);
        Eigen::Vector3d M_eng = arm1.cross(F1) + arm2.cross(F2);

        Eigen::Vector3d M_total = MAcg_b + M_eng;
        return InverseInertiaMatrix * (M_total - omega_b.cross(InertiaMatrix * omega_b));
    }


    std::map<std::string, double> getAircraftSpecs() const override
    {
        return {
            {"mass",                   mass},
            {"mean_aerodynamic_chord", mean_aerodynamic_chord},
            {"lt",                     lt},
            {"wing_planform_area",     wing_planform_area},
            {"tail_planform_area",     tail_planform_area},
            {"x_cg_pos_Fm",            x_cg_pos_Fm},
            {"y_cg_pos_Fm",            y_cg_pos_Fm},
            {"z_cg_pos_Fm",            z_cg_pos_Fm},
            {"x_aero_pos_Fm",          x_aero_pos_Fm},
            {"y_aero_pos_Fm",          y_aero_pos_Fm},
            {"z_aero_pos_Fm",          z_aero_pos_Fm},
            {"eng1_x_pos_Fm",          eng1_x_pos_Fm},
            {"eng1_y_pos_Fm",          eng1_y_pos_Fm},
            {"eng1_z_pos_Fm",          eng1_z_pos_Fm},
            {"eng2_x_pos_Fm",          eng2_x_pos_Fm},
            {"eng2_y_pos_Fm",          eng2_y_pos_Fm},
            {"eng2_z_pos_Fm",          eng2_z_pos_Fm}
        };
    }


    // ---- Geometry ----
    double mean_aerodynamic_chord{}, lt{};
    double wing_planform_area{}, tail_planform_area{};
    double x_cg_pos_Fm{}, y_cg_pos_Fm{}, z_cg_pos_Fm{};
    double x_aero_pos_Fm{}, y_aero_pos_Fm{}, z_aero_pos_Fm{};
    double eng1_x_pos_Fm{}, eng1_y_pos_Fm{}, eng1_z_pos_Fm{};
    double eng2_x_pos_Fm{}, eng2_y_pos_Fm{}, eng2_z_pos_Fm{};

    // ---- Static aero (lookup tables) ----
    LookupTable<2> cl_table_;        // CL(α, η)
    LookupTable<2> cd_table_;        // CD(α, η)
    LookupTable<2> cm_table_;        // Cm(α, η)
    LookupTable<2> cy_table_;        // CY(β, δr)
    LookupTable<2> cl_roll_table_;   // Cl(β, δa)
    LookupTable<2> cn_table_;        // Cn(β, δr)
};


// ============================================================
//  Aircraft factory
// ============================================================
enum class AircraftType { CivilAircraft, NasaGTM };

inline std::unique_ptr<AircraftModel> makeAircraft(AircraftType type,
    Eigen::VectorXd& X,
    Eigen::VectorXd& U)
{
    switch (type) {
    case AircraftType::CivilAircraft: return std::make_unique<CivilAircraft>(X, U);
    case AircraftType::NasaGTM:       return std::make_unique<NasaGTMAircraft>(X, U);
    }
    throw std::invalid_argument("Unknown aircraft type");
}