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
#include <memory>


// ============================================================
//  TrimSolver
//
//  Weighted least-squares trim solver. Mirrors the MATLAB
//  formulation but generalized to non-zero flight path angle:
//
//      Z = [X(1:12); U(1:5)]                     (17 unknowns)
//      Q = [xdot(1:9);                           (9 dynamics)
//           Va - Va_target;                      (airspeed)
//           gamma - gamma_target;                (flight path angle)
//           v;                                   (sideslip vel)
//           phi;                                 (bank)
//           psi;                                 (heading)
//           z_e_dot + Va*sin(gamma_target);      (climb-rate consistency)
//           U(4) - U(5);                         (throttle symmetry)
//           y_e_dot;                             (lateral drift)
//           h - h_target]                        (altitude)     => 18 terms
//
//  At cruise (gamma_target = 0) this reduces to the original
//  level-flight cost. At gamma_target > 0 the aircraft climbs;
//  z_e_dot becomes -Va*sin(gamma) at the steady-state operating
//  point, so we enforce that exact value rather than zero.
//
//  Optimization: damped Gauss-Newton. Trim is solved internally
//  with the Euler attitude representation; the result is lifted
//  to the simulation's attitude representation via buildFullState.
// ============================================================


struct TrimResult {
    Eigen::VectorXd Z;             // 17-element [X_euler; U]
    Eigen::VectorXd states_euler;  // 12-element Euler-style state
    Eigen::VectorXd controls;      // NUM_CONTROLS-element controls
    Eigen::VectorXd Q_final;       // 18-element residual at convergence
    double F0_final;               // Q' * W * Q at convergence
    int    iterations;
    bool   converged;
};


inline Eigen::VectorXd defaultTrimInitialGuess(double Va_init, double h_init,
    double gam_init = 0.0)
{
    Eigen::VectorXd Z = Eigen::VectorXd::Zero(17);
    // For non-zero flight path angle, theta starts higher (theta ~ alpha + gamma)
    const double alpha_guess = 0.05;
    const double theta_guess = alpha_guess + gam_init;
    Z[0] = Va_init * std::cos(alpha_guess);  // u
    Z[2] = Va_init * std::sin(alpha_guess);  // w
    Z[7] = theta_guess;                      // theta
    Z[11] = -h_init;                          // z_e (NED)
    Z[13] = -0.10;                            // stabilizer
    // More throttle needed for climb
    const double throttle_guess = 0.20 + 0.5 * std::sin(gam_init);
    Z[15] = throttle_guess;
    Z[16] = throttle_guess;
    return Z;
}


inline Eigen::VectorXd defaultTrimWeights()
{
    Eigen::VectorXd H = Eigen::VectorXd::Ones(18);
    H[17] = 1000.0;                    // heavy weight on altitude
    return H;
}


namespace trim_detail {

    inline Eigen::VectorXd computeQ(AircraftModel& ac,
        const AttitudeKinematics& euler_kin,
        const AtmosphereData& atm,
        const Eigen::VectorXd& Z,
        double Va_target,
        double h_target,
        double gam_target)
    {
        Eigen::VectorXd X = Z.head(12);
        Eigen::VectorXd U = Z.tail(constants::NUM_CONTROLS);

        Eigen::VectorXd xdot = Aircraft_Sim(ac, euler_kin, X, U, atm);

        double u = X[0], v = X[1], w = X[2];
        double phi = X[6], theta = X[7], psi = X[8];
        double z_e = X[11];

        double Va = std::sqrt(u * u + v * v + w * w);
        double alpha = (Va > 1e-10) ? std::atan2(w, u) : 0.0;
        double gam = theta - alpha;           // FPA for symmetric flight
        double h = -z_e;

        double xe_dot_z = xdot[11];             // z_e_dot
        double xe_dot_y = xdot[10];             // y_e_dot

        // At steady climb, vertical velocity is -V·sin(gamma) in NED
        // (negative z_e_dot means altitude increasing).
        double ze_dot_target = -Va_target * std::sin(gam_target);

        Eigen::VectorXd Q(18);
        Q.head(9) = xdot.head(9);               // u,v,w,p,q,r,phi,theta,psi dot
        Q[9] = Va - Va_target;
        Q[10] = gam - gam_target;
        Q[11] = v;
        Q[12] = phi;
        Q[13] = psi;
        Q[14] = xe_dot_z - ze_dot_target;
        Q[15] = U[3] - U[4];
        Q[16] = xe_dot_y;
        Q[17] = h - h_target;
        return Q;
    }

    inline double weightedCost(const Eigen::VectorXd& Q, const Eigen::VectorXd& H_diag)
    {
        return (H_diag.array() * Q.array() * Q.array()).sum();
    }

}  // namespace trim_detail


inline TrimResult solveTrim(AircraftModel& ac,
    const AtmosphereData& atm,
    double Va_target,
    double h_target,
    double gam_target = 0.0,
    const Eigen::VectorXd& Z0 = Eigen::VectorXd{},
    const Eigen::VectorXd& H_diag = Eigen::VectorXd{},
    double tol = 1e-12,
    int    max_iter = 100,
    bool   verbose = true)
{
    using namespace constants;

    auto euler_kin = makeAttitudeKinematics(AttitudeMode::Euler);

    Eigen::VectorXd Z = Z0.size() == 17
        ? Z0
        : defaultTrimInitialGuess(Va_target, h_target, gam_target);
    Eigen::VectorXd H = H_diag.size() == 18 ? H_diag : defaultTrimWeights();

    Eigen::VectorXd Q;
    Eigen::MatrixXd J(18, 17);
    double F0 = 0.0;

    constexpr double fd_eps = 1e-7;

    if (verbose) {
        std::cout << "\n=== Trim Solver (" << ac.name() << ") ===\n";
        std::cout << "Target: Va = " << Va_target << " m/s, altitude = "
            << h_target << " m, gamma = "
            << gam_target * RAD_TO_DEG << " deg\n";
        std::cout << "17 unknowns / 18 cost terms\n";
        std::cout << std::string(72, '-') << "\n";
    }

    Q = trim_detail::computeQ(ac, *euler_kin, atm, Z, Va_target, h_target, gam_target);
    F0 = trim_detail::weightedCost(Q, H);

    int iter = 0;
    for (iter = 0; iter < max_iter; ++iter)
    {
        if (verbose) {
            std::cout << "  Iter " << std::setw(3) << iter
                << "  F0 = " << std::scientific << std::setprecision(4) << F0
                << "  ||Q|| = " << Q.norm() << "\n";
        }
        if (F0 < tol) break;

        for (int j = 0; j < 17; ++j) {
            Eigen::VectorXd Zp = Z, Zm = Z;
            Zp[j] += fd_eps;
            Zm[j] -= fd_eps;
            Eigen::VectorXd Qp = trim_detail::computeQ(ac, *euler_kin, atm,
                Zp, Va_target, h_target, gam_target);
            Eigen::VectorXd Qm = trim_detail::computeQ(ac, *euler_kin, atm,
                Zm, Va_target, h_target, gam_target);
            J.col(j) = (Qp - Qm) / (2.0 * fd_eps);
        }

        Eigen::MatrixXd W = H.asDiagonal();
        Eigen::MatrixXd JtWJ = J.transpose() * W * J;
        Eigen::VectorXd JtWQ = J.transpose() * W * Q;
        Eigen::VectorXd dZ = JtWJ.colPivHouseholderQr().solve(-JtWQ);

        // Damped step
        double step = 1.0;
        Eigen::VectorXd Z_trial, Q_trial;
        double F0_trial = F0;
        bool   accepted = false;
        for (int trial = 0; trial < 20; ++trial) {
            Z_trial = Z + step * dZ;
            Q_trial = trim_detail::computeQ(ac, *euler_kin, atm, Z_trial,
                Va_target, h_target, gam_target);
            F0_trial = trim_detail::weightedCost(Q_trial, H);
            if (F0_trial < F0) { accepted = true; break; }
            step *= 0.5;
        }
        if (!accepted) {
            if (verbose) {
                std::cout << "  Line search failed - no improvement found.\n";
            }
            break;
        }
        Z = Z_trial;
        Q = Q_trial;
        F0 = F0_trial;
    }

    TrimResult result;
    result.Z = Z;
    result.states_euler = Z.head(12);
    result.controls = Z.tail(NUM_CONTROLS);
    result.Q_final = Q;
    result.F0_final = F0;
    result.iterations = iter;
    result.converged = (F0 < tol);

    if (verbose) {
        std::cout << std::string(72, '-') << "\n";
        if (result.converged) {
            std::cout << "Converged in " << iter << " iterations.\n";
        }
        else {
            std::cout << "WARNING: did not converge (F0 = "
                << std::scientific << result.F0_final << ").\n";
        }
        double u = Z[0], w = Z[2];
        double Va_final = std::sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
        double alpha = (Va_final > 1e-10) ? std::atan2(w, u) : 0.0;
        double gam_final = Z[7] - alpha;
        std::cout << std::fixed << std::setprecision(5)
            << "  Va         = " << Va_final << " m/s\n"
            << "  alpha      = " << alpha << " rad ("
            << alpha * RAD_TO_DEG << " deg)\n"
            << "  theta      = " << Z[7] << " rad ("
            << Z[7] * RAD_TO_DEG << " deg)\n"
            << "  gamma      = " << gam_final << " rad ("
            << gam_final * RAD_TO_DEG << " deg)\n"
            << "  stabilizer = " << Z[13] << " rad ("
            << Z[13] * RAD_TO_DEG << " deg)\n"
            << "  throttle1  = " << Z[15] << "\n";
        if (ac.engines_.size() >= 2) {
            std::cout << "  throttle2  = " << Z[16] << "\n";
        }

        return result;
    }
}