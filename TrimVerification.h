#pragma once

#include "AircraftModel.h"
#include "AtmosphereModel.h"
#include "AttitudeKinematics.h"
#include "Constants.h"
#include "SimulationEngine.h"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


// ============================================================
//  Trim verification — evaluate the EOM at the trim point and
//  confirm every derivative is below tolerance. Now supports
//  non-zero flight path angle via gamma_target argument.
// ============================================================

inline void verifyStraightLevelTrimConditions(AircraftModel& ac,
    const AttitudeKinematics& kin,
    const Eigen::VectorXd& states,
    const Eigen::VectorXd& controls,
    const AtmosphereData& atm,
    double Va_target = 85.0,
    double h_target = 3048.0,
    double gam_target = 0.0)
{
    Eigen::VectorXd XDOT = Aircraft_Sim(ac, kin, states, controls, atm);

    double u = states[0], v = states[1], w = states[2];
    double z_e = states[states.size() - 1];

    Eigen::VectorXd att = states.segment(kin.attitude_start_idx(), kin.n_attitude());
    Eigen::Vector3d eul = kin.to_euler(att);
    double phi = eul[0], theta = eul[1], psi = eul[2];

    double Va = std::sqrt(u * u + v * v + w * w);
    double alpha = std::atan2(w, u);
    double gamma = theta - alpha;
    double h = -z_e;

    // Expected vertical velocity at steady climb/descent
    double ze_dot_expected = -Va_target * std::sin(gam_target);

    constexpr double tol_d = 1e-5, tol_s = 1e-5, tol_a = 1e-3, tol_v = 1e-5;

    struct C { std::string name; double value, target; bool ok; };

    std::vector<C> conds = {
        {"udot",         XDOT[0],                  0.0,               std::abs(XDOT[0]) < tol_d},
        {"vdot",         XDOT[1],                  0.0,               std::abs(XDOT[1]) < tol_d},
        {"wdot",         XDOT[2],                  0.0,               std::abs(XDOT[2]) < tol_d},
        {"pdot",         XDOT[3],                  0.0,               std::abs(XDOT[3]) < tol_d},
        {"qdot",         XDOT[4],                  0.0,               std::abs(XDOT[4]) < tol_d},
        {"rdot",         XDOT[5],                  0.0,               std::abs(XDOT[5]) < tol_d},
        {"z_e_dot",      XDOT[XDOT.size() - 1],      ze_dot_expected,   std::abs(XDOT[XDOT.size() - 1]
                                                                                - ze_dot_expected) < tol_d},
        {"Va (m/s)",     Va,                       Va_target,         std::abs(Va - Va_target) < tol_v},
        {"FPA (rad)",    gamma,                    gam_target,        std::abs(gamma - gam_target) < tol_s},
        {"v (m/s)",      v,                        0.0,               std::abs(v) < tol_s},
        {"phi (rad)",    phi,                      0.0,               std::abs(phi) < tol_s},
        {"psi (rad)",    psi,                      0.0,               std::abs(psi) < tol_s},
        {"Altitude (m)", h,                        h_target,          std::abs(h - h_target) < tol_a}
    };

    for (int i = 0; i < kin.n_attitude(); ++i) {
        int idx = kin.attitude_start_idx() + i;
        std::string nm = "attdot[" + std::to_string(i) + "]";
        conds.push_back({ nm, XDOT[idx], 0.0, std::abs(XDOT[idx]) < tol_d });
    }

    std::cout << "\n=== Trim Verification (" << ac.name() << " / " << kin.name() << ")"
        << "  gamma=" << gam_target * constants::RAD_TO_DEG << " deg ===\n";
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