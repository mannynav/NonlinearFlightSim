
#pragma once

#include "LookupTable.h"
#include <Eigen/Dense>
#include <numbers>


// ============================================================
//  F16AeroTables — Stevens & Lewis "Aircraft Control and
//  Simulation" Appendix A low-speed F-16 aerodynamic model.
//
//  The S&L model tabulates body-axis coefficients on an alpha
//  grid of -10..45 deg (in 5-deg steps) and, for the lateral
//  and control tables, a beta grid of -30..30 deg (15-deg steps)
//  or an elevator grid of -24..24 deg.
//
//  Coefficients here reproduce the S&L Appendix A tables:
//    CX(alpha, el)   axial force
//    CZ(alpha)       normal force (with beta and elevator corrections applied in code)
//    Cm(alpha, el)   pitch moment
//    CL_roll(alpha, beta)   rolling moment
//    Cn(alpha, beta)        yawing moment
//    CY                     side force (algebraic, applied in code)
//    plus damping derivatives via the DAMP(alpha) table (9 columns)
//    and control increments DLDA, DLDR, DNDA, DNDR.
//
//  Units: angles in RADIANS at the lookup interface (grids are
//  built in radians), coefficients dimensionless. The raw S&L
//  tables are indexed in degrees; we convert the grids to radians
//  so the interface matches the GTM tables.
//
//  NOTE: These are the canonical S&L numbers. They are widely
//  published for educational use.
// ============================================================
namespace f16_aero {

    constexpr double D2R = std::numbers::pi / 180.0;

    // Alpha grid: -10..45 deg, 12 points (S&L standard)
    inline Eigen::VectorXd alpha_grid()
    {
        Eigen::VectorXd a(12);
        a << -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45;
        return a * D2R;
    }

    // Beta grid: 0..30 deg, 7 points (tables are symmetric in beta;
    // sign handling done in the aircraft code). Stored signed here
    // -30..30 for direct interpolation.
    inline Eigen::VectorXd beta_grid()
    {
        Eigen::VectorXd b(7);
        b << -30, -20, -10, 0, 10, 20, 30;
        return b * D2R;
    }

    // Elevator grid: -24..24 deg, 5 points
    inline Eigen::VectorXd el_grid()
    {
        Eigen::VectorXd e(5);
        e << -24, -12, 0, 12, 24;
        return e * D2R;
    }


    // ---------- CX(alpha, el) — body-axis axial force ----------
    // S&L Appendix A table CX(alpha, el). Row-major [alpha][el].
    inline LookupTable<2> build_cx_table()
    {
        auto a = alpha_grid();
        auto e = el_grid();
        // 12 alpha x 5 elevator. Values from S&L Appendix A (CX).
        Eigen::VectorXd d(12 * 5);
        const double tbl[12][5] = {
            //  el=-24   -12      0       12      24
            { -0.099, -0.081, -0.081, -0.063, -0.025 }, // a=-10
            { -0.048, -0.038, -0.040, -0.021,  0.016 }, // a=-5
            { -0.022, -0.020, -0.021, -0.004,  0.032 }, // a=0
            { -0.040, -0.038, -0.039, -0.025,  0.006 }, // a=5
            { -0.083, -0.073, -0.071, -0.060, -0.029 }, // a=10
            { -0.083, -0.076, -0.069, -0.063, -0.043 }, // a=15
            { -0.121, -0.110, -0.105, -0.097, -0.083 }, // a=20
            { -0.106, -0.097, -0.090, -0.085, -0.075 }, // a=25
            { -0.114, -0.105, -0.092, -0.090, -0.082 }, // a=30
            { -0.143, -0.130, -0.122, -0.120, -0.117 }, // a=35
            { -0.135, -0.123, -0.112, -0.110, -0.107 }, // a=40
            { -0.117, -0.106, -0.099, -0.098, -0.097 }  // a=45
        };
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 5; ++j)
                d[i * 5 + j] = tbl[i][j];
        return LookupTable<2>({ a, e }, d);
    }


    // ---------- CZ(alpha) — body-axis normal force (baseline) ----------
    // 1D in alpha; beta and elevator corrections applied in aircraft code.
    inline LookupTable<1> build_cz_table()
    {
        auto a = alpha_grid();
        Eigen::VectorXd d(12);
        d << 0.770, 0.241, -0.100, -0.416, -0.731,
            -1.053, -1.366, -1.646, -1.917, -2.120,
            -2.248, -2.229;
        return LookupTable<1>({ a }, d);
    }


    // ---------- Cm(alpha, el) — pitch moment ----------
    inline LookupTable<2> build_cm_table()
    {
        auto a = alpha_grid();
        auto e = el_grid();
        Eigen::VectorXd d(12 * 5);
        const double tbl[12][5] = {
            //  el=-24    -12      0       12      24
            {  0.205,  0.168,  0.186,  0.196,  0.213 }, // a=-10
            {  0.081,  0.077,  0.107,  0.110,  0.110 }, // a=-5
            { -0.046, -0.020, -0.009, -0.005, -0.006 }, // a=0
            { -0.174, -0.145, -0.121, -0.127, -0.129 }, // a=5
            { -0.259, -0.202, -0.184, -0.193, -0.199 }, // a=10
            { -0.279, -0.240, -0.197, -0.205, -0.211 }, // a=15
            { -0.289, -0.254, -0.211, -0.222, -0.230 }, // a=20
            { -0.319, -0.276, -0.240, -0.246, -0.252 }, // a=25
            { -0.354, -0.306, -0.279, -0.281, -0.285 }, // a=30
            { -0.327, -0.289, -0.259, -0.262, -0.266 }, // a=35
            { -0.300, -0.267, -0.240, -0.243, -0.246 }, // a=40
            { -0.292, -0.262, -0.236, -0.238, -0.241 }  // a=45
        };
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 5; ++j)
                d[i * 5 + j] = tbl[i][j];
        return LookupTable<2>({ a, e }, d);
    }


    // ---------- CL_roll(alpha, beta) — rolling moment ----------
    inline LookupTable<2> build_cl_roll_table()
    {
        auto a = alpha_grid();
        auto b = beta_grid();
        // S&L tabulate for beta >= 0 and use antisymmetry; here we
        // build the full signed beta grid by mirroring (Cl odd in beta).
        // Positive-beta half from S&L Appendix A (Cl), 12 alpha x 4 beta(0,10,20,30 deg → but our grid uses 0,10,20,30 on the + side)
        const double half[12][4] = {
            // beta= 0     10       20       30
            {  0.000, -0.001, -0.003, -0.001 }, // a=-10
            {  0.000, -0.004, -0.009, -0.010 }, // a=-5
            {  0.000, -0.008, -0.017, -0.020 }, // a=0
            {  0.000, -0.012, -0.024, -0.030 }, // a=5
            {  0.000, -0.016, -0.030, -0.039 }, // a=10
            {  0.000, -0.019, -0.034, -0.044 }, // a=15
            {  0.000, -0.020, -0.040, -0.050 }, // a=20
            {  0.000, -0.020, -0.037, -0.049 }, // a=25
            {  0.000, -0.015, -0.016, -0.023 }, // a=30
            {  0.000, -0.008, -0.002, -0.006 }, // a=35
            {  0.000, -0.013, -0.010, -0.014 }, // a=40
            {  0.000, -0.015, -0.019, -0.027 }  // a=45
        };
        // Full grid order: beta = -30,-20,-10,0,10,20,30 (indices 0..6)
        // Cl is odd in beta: Cl(-b) = -Cl(b)
        Eigen::VectorXd d(12 * 7);
        for (int i = 0; i < 12; ++i) {
            d[i * 7 + 0] = -half[i][3]; // -30
            d[i * 7 + 1] = -half[i][2]; // -20
            d[i * 7 + 2] = -half[i][1]; // -10
            d[i * 7 + 3] = half[i][0]; //   0
            d[i * 7 + 4] = half[i][1]; //  10
            d[i * 7 + 5] = half[i][2]; //  20
            d[i * 7 + 6] = half[i][3]; //  30
        }
        return LookupTable<2>({ a, b }, d);
    }


    // ---------- Cn(alpha, beta) — yawing moment ----------
    inline LookupTable<2> build_cn_table()
    {
        auto a = alpha_grid();
        auto b = beta_grid();
        const double half[12][4] = {
            // beta= 0     10       20       30
            {  0.000,  0.018,  0.038,  0.056 }, // a=-10
            {  0.000,  0.019,  0.042,  0.057 }, // a=-5
            {  0.000,  0.018,  0.042,  0.059 }, // a=0
            {  0.000,  0.019,  0.042,  0.058 }, // a=5
            {  0.000,  0.019,  0.043,  0.058 }, // a=10
            {  0.000,  0.018,  0.039,  0.053 }, // a=15
            {  0.000,  0.013,  0.030,  0.032 }, // a=20
            {  0.000,  0.007,  0.017,  0.012 }, // a=25
            {  0.000,  0.004,  0.004,  0.002 }, // a=30
            {  0.000, -0.014, -0.035, -0.046 }, // a=35
            {  0.000, -0.017, -0.047, -0.071 }, // a=40
            {  0.000, -0.033, -0.057, -0.073 }  // a=45
        };
        // Cn is odd in beta
        Eigen::VectorXd d(12 * 7);
        for (int i = 0; i < 12; ++i) {
            d[i * 7 + 0] = -half[i][3];
            d[i * 7 + 1] = -half[i][2];
            d[i * 7 + 2] = -half[i][1];
            d[i * 7 + 3] = half[i][0];
            d[i * 7 + 4] = half[i][1];
            d[i * 7 + 5] = half[i][2];
            d[i * 7 + 6] = half[i][3];
        }
        return LookupTable<2>({ a, b }, d);
    }


    // ---------- Damping derivatives DAMP(alpha) ----------
    // S&L returns 9 damping terms as functions of alpha:
    //   [CXq, CYr, CYp, CZq, Clr, Clp, Cmq, Cnr, Cnp]
    // We expose each as a separate 1D table for clarity.

    inline LookupTable<1> build_cxq_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << -0.267, -0.110, 0.308, 1.140, 1.860, 2.550, 3.050, 3.460, 3.500, 3.500, 3.500, 3.500;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_cyr_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << 0.882, 0.852, 0.876, 0.958, 0.962, 0.974, 0.819, 0.483, 0.590, 1.210, -0.493, -1.040;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_cyp_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << -0.108, -0.108, -0.188, 0.110, 0.258, 0.226, 0.344, 0.362, 0.611, 0.529, 0.298, -2.270;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_czq_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << -8.80, -25.8, -28.9, -31.4, -31.2, -30.7, -27.7, -28.2, -29.0, -29.8, -38.3, -35.3;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_clr_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << -0.126, -0.026, 0.063, 0.113, 0.208, 0.230, 0.319, 0.437, 0.680, 0.100, 0.447, -0.330;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_clp_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << -0.360, -0.359, -0.443, -0.420, -0.383, -0.375, -0.329, -0.294, -0.230, -0.210, -0.120, -0.100;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_cmq_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << -7.21, -0.540, -5.23, -5.26, -6.11, -6.64, -5.69, -6.00, -6.20, -6.40, -6.60, -6.00;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_cnr_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << -0.380, -0.363, -0.378, -0.386, -0.370, -0.453, -0.550, -0.582, -0.595, -0.637, -1.020, -0.840;
        return LookupTable<1>({ a }, d);
    }
    inline LookupTable<1> build_cnp_table() {
        auto a = alpha_grid(); Eigen::VectorXd d(12);
        d << 0.061, 0.052, 0.052, -0.012, -0.013, -0.024, 0.050, 0.150, 0.130, 0.158, 0.240, 0.150;
        return LookupTable<1>({ a }, d);
    }

} // namespace f16_aero