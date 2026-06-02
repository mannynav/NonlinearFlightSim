#pragma once

#include "LQRController.h"
#include "AircraftModel.h"
#include "AttitudeKinematics.h"
#include "AtmosphereModel.h"
#include "SimulationEngine.h"
#include "ModeAnalysis.h"
#include "Constants.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>


// ============================================================
//  runLQRDemo
//
//  End-to-end LQR demonstration:
//    1. Extract longitudinal control columns from B_long
//       (stabilizer at index 1, throttle1 at index 3)
//    2. Build state and control weighting matrices Q, R using
//       Bryson's rule (1/max_allowed_deviation^2 per channel)
//    3. Solve the continuous-time Riccati equation
//    4. Compute the LQR gain K = R^{-1} B^T P
//    5. Print closed-loop eigenvalues (modes the controller
//       has shaped)
//    6. Run two nonlinear simulations from the same initial
//       condition (trim + pitch perturbation):
//         a. OPEN  loop — controls held at trim, the natural
//            dynamics oscillate at the short-period and phugoid
//         b. CLOSED loop — controls = u_trim - K*(x_long - x_long_trim)
//            at every step, the controller drives the state back
//            to trim
//
//  Two CSVs are written: lqr_open_loop.csv and lqr_closed_loop.csv.
//  Same number of rows, same time grid, ready for side-by-side plotting.
// ============================================================

inline void runLQRDemo(AircraftModel& ac,
    AttitudeKinematics& kin,
    const Eigen::VectorXd& trim_states,
    const Eigen::VectorXd& trim_states_euler,
    const Eigen::VectorXd& trim_controls,
    const Eigen::MatrixXd& A_long,
    const Eigen::MatrixXd& B_long,
    const Eigen::MatrixXd& S_perm,
    const AtmosphereData& atm,
    double sim_length = 10.0,
    int    sim_steps = 1000,
    double pitch_perturbation = 0.05)
{
    using namespace constants;

    std::cout << "\n############################################################\n";
    std::cout << "  LQR longitudinal autopilot\n";
    std::cout << "############################################################\n";

    // ---- Extract longitudinal control columns from B_long ----
    // Long subsystem inputs we'll use: stabilizer (idx 1), throttle1 (idx 3)
    Eigen::MatrixXd B_lon(4, 2);
    B_lon.col(0) = B_long.col(1);
    B_lon.col(1) = B_long.col(3);

    std::cout << "A_long (4x4):\n" << A_long << "\n\n";
    std::cout << "B_lon  (4x2)  [stabilizer, throttle1]:\n" << B_lon << "\n\n";

    // ---- LQR weights via Bryson's rule ----
    //   Q_ii = 1 / (max acceptable state deviation)^2
    //   R_jj = 1 / (max acceptable control deviation)^2
    Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
    Q(0, 0) = 1.0 / (10.0 * 10.0);    // u  : 10 m/s
    Q(1, 1) = 1.0 / (5.0 * 5.0);    // w  :  5 m/s
    Q(2, 2) = 1.0 / (0.1 * 0.1);    // q  : 0.1 rad/s
    Q(3, 3) = 1.0 / (0.1 * 0.1);    // theta-state: 0.1

    Eigen::Matrix2d R = Eigen::Matrix2d::Zero();
    const double de_max = 10.0 * DEG_TO_RAD;
    const double th_max = 0.3;
    R(0, 0) = 1.0 / (de_max * de_max);
    R(1, 1) = 1.0 / (th_max * th_max);

    std::cout << "Q (state weights, Bryson's rule):\n" << Q << "\n\n";
    std::cout << "R (control weights):\n" << R << "\n\n";

    // ---- Solve the Riccati equation ----
    Eigen::MatrixXd P = LQRController::solveCARE(A_long, B_lon, Q, R);
    Eigen::MatrixXd K = R.inverse() * B_lon.transpose() * P;

    std::cout << "Riccati solution P (4x4):\n" << P << "\n\n";
    std::cout << "LQR gain K (2x4)  [rows: stabilizer, throttle1]:\n"
        << K << "\n\n";

    // ---- Closed-loop eigenvalues (linear prediction) ----
    Eigen::MatrixXd A_cl = A_long - B_lon * K;
    printModeAnalysis("Open-loop  longitudinal", A_long);
    printModeAnalysis("Closed-loop longitudinal", A_cl);

    // ---- Reference longitudinal state at trim ----
    Eigen::Vector4d x_long_trim = (S_perm * trim_states).head(4);

    // ---- Initial condition: trim + pitch perturbation (rebuild from Euler) ----
    Eigen::VectorXd ic_euler = trim_states_euler;
    ic_euler[7] += pitch_perturbation;
    Eigen::VectorXd X0 = buildFullState(ic_euler, kin);

    const double dt = sim_length / sim_steps;
    const int    N = static_cast<int>(X0.size());
    const int    NU = NUM_CONTROLS;

    // ---- Helper to write a single trajectory to CSV ----
    auto writeCSV = [&](const std::string& filename,
        const Eigen::MatrixXd& T_states,
        const Eigen::MatrixXd& T_controls,
        const Eigen::VectorXd& T_time)
        {
            std::ofstream f(filename);
            f << "t";
            for (int i = 0; i < N; ++i) f << ",x" << i;
            for (int i = 0; i < NU; ++i) f << ",u" << i;
            f << "\n";
            f << std::setprecision(10);
            for (int k = 0; k <= sim_steps; ++k) {
                f << T_time[k];
                for (int i = 0; i < N; ++i) f << "," << T_states(i, k);
                for (int i = 0; i < NU; ++i) f << "," << T_controls(i, k);
                f << "\n";
            }
            std::cout << "  wrote " << filename << "\n";
        };

    // ============================================================
    //  OPEN LOOP — controls clamped to trim, watch the modes oscillate
    // ============================================================
    {
        Eigen::MatrixXd X_hist(N, sim_steps + 1);
        Eigen::MatrixXd U_hist(NU, sim_steps + 1);
        Eigen::VectorXd t_hist(sim_steps + 1);

        Eigen::VectorXd X = X0;
        Eigen::VectorXd U = trim_controls;

        X_hist.col(0) = X;
        U_hist.col(0) = U;
        t_hist[0] = 0.0;

        for (int k = 0; k < sim_steps; ++k) {
            Eigen::VectorXd k1 = Aircraft_Sim(ac, kin, X, U, atm);
            Eigen::VectorXd k2 = Aircraft_Sim(ac, kin, X + 0.5 * dt * k1, U, atm);
            Eigen::VectorXd k3 = Aircraft_Sim(ac, kin, X + 0.5 * dt * k2, U, atm);
            Eigen::VectorXd k4 = Aircraft_Sim(ac, kin, X + dt * k3, U, atm);
            X = X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

            X_hist.col(k + 1) = X;
            U_hist.col(k + 1) = U;
            t_hist[k + 1] = (k + 1) * dt;
        }
        std::cout << "Open-loop simulation done.\n";
        writeCSV("lqr_open_loop.csv", X_hist, U_hist, t_hist);
    }

    // ============================================================
    //  CLOSED LOOP — controls = trim + (-K * dx_long)
    // ============================================================
    {
        Eigen::MatrixXd X_hist(N, sim_steps + 1);
        Eigen::MatrixXd U_hist(NU, sim_steps + 1);
        Eigen::VectorXd t_hist(sim_steps + 1);

        Eigen::VectorXd X = X0;

        auto controlFromState = [&](const Eigen::VectorXd& Xs) {
            Eigen::Vector4d x_long = (S_perm * Xs).head(4);
            Eigen::Vector4d dx = x_long - x_long_trim;
            Eigen::Vector2d du = -K * dx;
            Eigen::VectorXd Uout = trim_controls;
            Uout[1] += du[0];                 // stabilizer
            Uout[3] += du[1];                 // throttle1
            Uout[4] = Uout[3];               // mirror for twin-engine aircraft;
            //   harmless for single-engine F-16
            return Uout;
            };

        X_hist.col(0) = X;
        U_hist.col(0) = controlFromState(X);
        t_hist[0] = 0.0;

        for (int k = 0; k < sim_steps; ++k) {
            Eigen::VectorXd U = controlFromState(X);
            Eigen::VectorXd k1 = Aircraft_Sim(ac, kin, X, U, atm);
            Eigen::VectorXd k2 = Aircraft_Sim(ac, kin, X + 0.5 * dt * k1, U, atm);
            Eigen::VectorXd k3 = Aircraft_Sim(ac, kin, X + 0.5 * dt * k2, U, atm);
            Eigen::VectorXd k4 = Aircraft_Sim(ac, kin, X + dt * k3, U, atm);
            X = X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

            X_hist.col(k + 1) = X;
            U_hist.col(k + 1) = controlFromState(X);
            t_hist[k + 1] = (k + 1) * dt;
        }
        std::cout << "Closed-loop simulation done.\n";
        writeCSV("lqr_closed_loop.csv", X_hist, U_hist, t_hist);
    }

    std::cout << "LQR demo complete.\n";
}