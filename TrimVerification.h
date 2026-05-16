
#pragma once

#include "AircraftModel.h"
#include "AttitudeKinematics.h"
#include "GravityModel.h"
#include "SimulationEngine.h"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


// ============================================================
//  Trim verification — evaluate the EOM at the trim point and
//  confirm every derivative is below tolerance. Prints a table
//  with one row per condition and a final pass/fail count.
// ============================================================
inline void verifyStraightLevelTrimConditions(AircraftModel& ac,
    const AttitudeKinematics& kin,
    const Eigen::VectorXd& states,
    const Eigen::VectorXd& controls,
    const GravityModel& gravity_model)
{
    Eigen::VectorXd XDOT = Aircraft_Sim(ac, kin, states, controls, gravity_model);

    double u = states[0], v = states[1], w = states[2];
    double z_e = states[states.size() - 1];

    Eigen::VectorXd att = states.segment(kin.attitude_start_idx(), kin.n_attitude());
    Eigen::Vector3d eul = kin.to_euler(att);
    double phi = eul[0], theta = eul[1], psi = eul[2];

    double Va = std::sqrt(u * u + v * v + w * w);
    double alpha = std::atan2(w, u);
    double gamma = theta - alpha;
    double h = -z_e;

    constexpr double tol_d = 1e-5, tol_s = 1e-5, tol_a = 1e-3, tol_v = 1e-5;

    struct C { std::string name; double value, target; bool ok; };

    std::vector<C> conds = {
        {"udot",         XDOT[0],                 0.0,    std::abs(XDOT[0]) < tol_d},
        {"vdot",         XDOT[1],                 0.0,    std::abs(XDOT[1]) < tol_d},
        {"wdot",         XDOT[2],                 0.0,    std::abs(XDOT[2]) < tol_d},
        {"pdot",         XDOT[3],                 0.0,    std::abs(XDOT[3]) < tol_d},
        {"qdot",         XDOT[4],                 0.0,    std::abs(XDOT[4]) < tol_d},
        {"rdot",         XDOT[5],                 0.0,    std::abs(XDOT[5]) < tol_d},
        {"z_e_dot",      XDOT[XDOT.size() - 1],     0.0,    std::abs(XDOT[XDOT.size() - 1]) < tol_d},
        {"Va (m/s)",     Va,                      85.0,   std::abs(Va - 85.0) < tol_v},
        {"FPA (rad)",    gamma,                   0.0,    std::abs(gamma) < tol_s},
        {"v (m/s)",      v,                       0.0,    std::abs(v) < tol_s},
        {"phi (rad)",    phi,                     0.0,    std::abs(phi) < tol_s},
        {"psi (rad)",    psi,                     0.0,    std::abs(psi) < tol_s},
        {"Altitude (m)", h,                       3048.0, std::abs(h - 3048.0) < tol_a}
    };

    // Attitude-state derivatives (one per attitude component)
    for (int i = 0; i < kin.n_attitude(); ++i) {
        int idx = kin.attitude_start_idx() + i;
        std::string nm = "attdot[" + std::to_string(i) + "]";
        conds.push_back({ nm, XDOT[idx], 0.0, std::abs(XDOT[idx]) < tol_d });
    }

    std::cout << "\n=== Trim Verification (" << ac.name() << " / " << kin.name() << ") ===\n";
    std::cout << std::left << std::setw(20) << "Condition"
        << std::setw(15) << "Value" << std::setw(15) << "Target"
        << "Status\n" << std::string(60, '-') << "\n";

    int ok_count = 0;
    for (const auto& c : conds) {
        std::cout << std::setw(20) << c.name
            << std::fixed << std::setprecision(6) << std::setw(15) << c.value
            << std::setw(15) << c.target
            << (c.ok ? "OK" : "FAIL") << "\n";
        if (c.ok) ok_count++;
    }
    std::cout << ok_count << "/" << conds.size() << " satisfied.\n\n";
}
