
#pragma once

#include "LookupTable.h"
#include <Eigen/Dense>
#include <numbers>


// ============================================================
//  GTM Aerodynamic Coefficient Tables
//
//  Hand-authored placeholder data with transport-class
//  characteristics. Built so a NasaGTM simulation runs and
//  behaves plausibly through stall.
//
//  Replace with real NASA GTM data extracted from MATLAB tables
//  in github.com/nasa/GTM_DesignSim. Grid layout and units below
//  are designed to drop in cleanly once real data is available.
//
//  Units throughout: angles in RADIANS, coefficients dimensionless.
//  Coefficient sign conventions match the original CivilAircraft
//  polynomial model:
//    Cm_η < 0      — positive stabilizer pitches the nose down
//    Cl_δa < 0     — positive aileron rolls the body negatively
//    Cn_δr < 0     — positive rudder yaws the body negatively
// ============================================================
namespace gtm_aero {

    constexpr double D2R = std::numbers::pi / 180.0;


    // ============================================================
    //  Longitudinal — functions of α and η (stabilizer)
    //  α grid: 9 points, refined near stall at 14°
    //  η grid: 6 points across full stabilizer travel
    // ============================================================
    inline Eigen::VectorXd alpha_grid_long()
    {
        Eigen::VectorXd a(9);
        a << -5, 0, 4, 8, 12, 14, 16, 18, 20;
        return a * D2R;
    }

    inline Eigen::VectorXd eta_grid_long()
    {
        Eigen::VectorXd e(6);
        e << -25, -15, -5, 0, 5, 10;
        return e * D2R;
    }


    // ---------- CL(α, η) — total static lift ----------
    inline LookupTable<2> build_cl_table()
    {
        auto alpha_grid = alpha_grid_long();
        auto eta_grid = eta_grid_long();

        // CL_base(α) at η = 0 — transport lift curve with stall at ~14°
        const double cl_base[9] = {
            -0.30, 0.20, 0.70, 1.10, 1.40, 1.50, 1.40, 1.20, 1.05
        };
        const double dcl_deta = 0.8;   // stabilizer effectiveness (per radian)

        Eigen::VectorXd data(9 * 6);
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 6; ++j)
                data[i * 6 + j] = cl_base[i] + dcl_deta * eta_grid[j];

        return LookupTable<2>({ alpha_grid, eta_grid }, data);
    }


    // ---------- CD(α, η) — total drag ----------
    //   Base drag polar grows roughly with CL²; post-stall drag rises
    //   sharply. Stabilizer trim drag added as a quadratic in η.
    inline LookupTable<2> build_cd_table()
    {
        auto alpha_grid = alpha_grid_long();
        auto eta_grid = eta_grid_long();

        // CD_base(α) at η = 0 — transport drag polar with post-stall rise
        const double cd_base[9] = {
            0.045, 0.024, 0.045, 0.090, 0.150, 0.190, 0.270, 0.380, 0.510
        };

        Eigen::VectorXd data(9 * 6);
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 6; ++j) {
                const double eta_deg = eta_grid[j] / D2R;
                const double dcd_eta = 1.5e-4 * eta_deg * eta_deg;  // quadratic trim drag
                data[i * 6 + j] = cd_base[i] + dcd_eta;
            }
        }
        return LookupTable<2>({ alpha_grid, eta_grid }, data);
    }


    // ---------- Cm(α, η) — pitch moment about CG ----------
    //   Statically stable through stall; sign convention matches the
    //   original polynomial (Cm_η ≈ -2.9/rad).
    inline LookupTable<2> build_cm_table()
    {
        auto alpha_grid = alpha_grid_long();
        auto eta_grid = eta_grid_long();

        // Cm_base(α) at η = 0 — slight pitch-up past stall
        const double cm_base[9] = {
            +0.178, -0.467, -0.609, -0.751, -0.894, -0.965, -1.100, -1.200, -1.100
        };
        const double dcm_deta = -2.91;   // stabilizer effectiveness (per radian)

        Eigen::VectorXd data(9 * 6);
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 6; ++j)
                data[i * 6 + j] = cm_base[i] + dcm_deta * eta_grid[j];

        return LookupTable<2>({ alpha_grid, eta_grid }, data);
    }


    // ============================================================
    //  Lateral — functions of β and primary lateral control
    //  β grid: 7 points (-15° to +15°)
    //  Primary control varies by coefficient
    // ============================================================
    inline Eigen::VectorXd beta_grid_lat()
    {
        Eigen::VectorXd b(7);
        b << -15, -10, -5, 0, 5, 10, 15;
        return b * D2R;
    }


    // ---------- CY(β, δr) — side force ----------
    inline LookupTable<2> build_cy_table()
    {
        auto beta_grid = beta_grid_lat();

        Eigen::VectorXd dr_grid(5);
        dr_grid << -30, -15, 0, 15, 30;
        dr_grid *= D2R;

        // Linear in both variables: CY = -1.7·β + 0.25·δr
        const double Cy_beta = -1.7;
        const double Cy_dr = 0.25;

        Eigen::VectorXd data(7 * 5);
        for (int i = 0; i < 7; ++i)
            for (int j = 0; j < 5; ++j)
                data[i * 5 + j] = Cy_beta * beta_grid[i] + Cy_dr * dr_grid[j];

        return LookupTable<2>({ beta_grid, dr_grid }, data);
    }


    // ---------- Cl(β, δa) — roll moment ----------
    inline LookupTable<2> build_cl_roll_table()
    {
        auto beta_grid = beta_grid_lat();

        Eigen::VectorXd da_grid(5);
        da_grid << -25, -10, 0, 10, 25;
        da_grid *= D2R;

        // Cl_β (dihedral effect) and Cl_δa (aileron) matching original sign convention
        const double Cl_beta = -1.5;
        const double Cl_da = -0.6;

        Eigen::VectorXd data(7 * 5);
        for (int i = 0; i < 7; ++i)
            for (int j = 0; j < 5; ++j)
                data[i * 5 + j] = Cl_beta * beta_grid[i] + Cl_da * da_grid[j];

        return LookupTable<2>({ beta_grid, da_grid }, data);
    }


    // ---------- Cn(β, δr) — yaw moment ----------
    inline LookupTable<2> build_cn_table()
    {
        auto beta_grid = beta_grid_lat();

        Eigen::VectorXd dr_grid(5);
        dr_grid << -30, -15, 0, 15, 30;
        dr_grid *= D2R;

        // Cn_β (weathercock stability) and Cn_δr (rudder)
        const double Cn_beta = 1.0;
        const double Cn_dr = -0.66;

        Eigen::VectorXd data(7 * 5);
        for (int i = 0; i < 7; ++i)
            for (int j = 0; j < 5; ++j)
                data[i * 5 + j] = Cn_beta * beta_grid[i] + Cn_dr * dr_grid[j];

        return LookupTable<2>({ beta_grid, dr_grid }, data);
    }


} // namespace gtm_aero