
#pragma once

#include "AircraftModel.h"
#include "AttitudeKinematics.h"
#include "CheckCases.h"
#include "GravityModel.h"
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <string>
#include <cmath>


// ============================================================
//  NESC Check-Case Validator
//
//  Drives the simulation with the rigid-body fixtures defined in
//  CheckCases.h and reports conservation diagnostics.
//
//    Case 1 (sphere)  : conserve total mechanical energy E = mgh + ½mv²
//                        and keep horizontal velocity at exactly zero
//    Case 2 (brick)   : conserve angular momentum |H_inertial| and
//                        rotational kinetic energy 0.5·ωᵀIω
//    Case 3 (damped)  : rotational KE strictly decreasing while the
//                        body continues to tumble
//
//  All diagnostics are computed against the RK4 simulation output;
//  drift over a 30 s sim is the signal of correctness.
// ============================================================


// Forward declaration — Aircraft_Sim lives in main.cpp; the validator
// below uses a local copy to keep the dependency one-way.

namespace nesc_internal {

    constexpr double PI = std::numbers::pi;

    struct CheckCaseConfig {
        double sim_time;
        int    n_steps;
        Eigen::VectorXd X12_euler_init;   // 12-element Euler IC: [u,v,w,p,q,r,phi,theta,psi,x,y,z]
    };

    // Build full state for the active attitude representation
    inline Eigen::VectorXd buildState(const Eigen::VectorXd& X12_euler,
        const AttitudeKinematics& kin)
    {
        Eigen::VectorXd X(kin.n_states());
        X.head(6) = X12_euler.head(6);
        Eigen::VectorXd att = kin.from_euler(X12_euler[6], X12_euler[7], X12_euler[8]);
        X.segment(kin.attitude_start_idx(), kin.n_attitude()) = att;
        X.tail(3) = X12_euler.tail(3);
        return X;
    }

    // Local Aircraft_Sim used only by the validator. Mirrors main.cpp's
    // Aircraft_Sim but with a guard against airspeed=0 since check-case
    // initial conditions can have V=0.
    inline Eigen::VectorXd checkcase_sim(AircraftModel& aircraft,
        const AttitudeKinematics& kin,
        const Eigen::VectorXd& X,
        const Eigen::VectorXd& U)
    {
        aircraft.initializeStates(X);
        aircraft.intializeControlInputs(U);

        double u = X[0], v = X[1], w = X[2];
        double p = X[3], q = X[4], r = X[5];
        Eigen::VectorXd att = X.segment(kin.attitude_start_idx(), kin.n_attitude());

        constexpr double rho = 1.225;
        constexpr double g = 9.81;

        double airSpeed = std::sqrt(u * u + v * v + w * w);
        double alpha = (airSpeed > 1e-10) ? std::atan2(w, u) : 0.0;
        double beta = (airSpeed > 1e-10) ? std::asin(v / airSpeed) : 0.0;
        double dynamicPressure = 0.5 * rho * airSpeed * airSpeed;

        Eigen::Vector3d omega_b(p, q, r);
        Eigen::Vector3d vel_b(u, v, w);

        Eigen::Matrix3d R_stab2body;
        R_stab2body << std::cos(alpha), 0, -std::sin(alpha),
            0, 1, 0,
            std::sin(alpha), 0, std::cos(alpha);

        Eigen::Vector3d F_grav = aircraft.mass * kin.gravity_body(g, att);
        Eigen::Vector3d F_aero = R_stab2body * aircraft.aerodynamic_forces_stability_axis(
            alpha, beta, airSpeed, dynamicPressure);
        Eigen::Vector3d F_eng = aircraft.engine_forces_body_frame(g);

        Eigen::Vector3d Vdot = (1.0 / aircraft.mass) * (F_grav + F_eng + F_aero)
            - omega_b.cross(vel_b);

        Eigen::Vector3d pqrdot = aircraft.cg_moments_body_frame(R_stab2body, omega_b,
            airSpeed, dynamicPressure, alpha, beta, g);

        Eigen::VectorXd att_dot = kin.attitude_derivative(omega_b, att);
        Eigen::Matrix3d R_n_b = kin.body_to_ned(att);
        Eigen::Vector3d nav = R_n_b * vel_b;

        Eigen::VectorXd XDOT(kin.n_states());
        XDOT.segment<3>(0) = Vdot;
        XDOT.segment<3>(3) = pqrdot;
        XDOT.segment(kin.attitude_start_idx(), kin.n_attitude()) = att_dot;
        XDOT.segment<3>(kin.position_start_idx()) = nav;
        return XDOT;
    }

    // RK4 sim with quaternion renormalization, mirrors runNonlinearSimulation
    inline Eigen::MatrixXd checkcase_rk4(AircraftModel& ac,
        const AttitudeKinematics& kin,
        const Eigen::VectorXd& X0,
        double sim_length, int steps)
    {
        const double dt = sim_length / steps;
        Eigen::MatrixXd sol(kin.n_states(), steps + 1);
        sol.setZero();
        sol.col(0) = X0;
        Eigen::VectorXd U = Eigen::VectorXd::Zero(5);

        for (int i = 0; i < steps; i++) {
            Eigen::VectorXd Xk = sol.col(i);
            Eigen::VectorXd k1 = checkcase_sim(ac, kin, Xk, U);
            Eigen::VectorXd k2 = checkcase_sim(ac, kin, (Xk + 0.5 * dt * k1).eval(), U);
            Eigen::VectorXd k3 = checkcase_sim(ac, kin, (Xk + 0.5 * dt * k2).eval(), U);
            Eigen::VectorXd k4 = checkcase_sim(ac, kin, (Xk + dt * k3).eval(), U);
            Eigen::VectorXd Xn = Xk + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

            Eigen::VectorXd att = Xn.segment(kin.attitude_start_idx(), kin.n_attitude());
            Xn.segment(kin.attitude_start_idx(), kin.n_attitude()) = kin.normalize(att);

            sol.col(i + 1) = Xn;
        }
        return sol;
    }
} // namespace nesc_internal


// ============================================================
//  Case 1 — Dropped sphere
// ============================================================
inline void runNESCCheckCase1(const AttitudeKinematics& kin)
{
    constexpr double sim_time = 30.0;
    constexpr int    n_steps = 1000;
    constexpr double g = 9.81;
    constexpr double h0 = 3048.0;   // 10,000 ft, matching CivilAircraft IC

    // IC: at altitude, zero velocity, level attitude, zero rates
    Eigen::VectorXd X12 = Eigen::VectorXd::Zero(12);
    X12[11] = -h0;   // z_e = -altitude (NED)

    Eigen::VectorXd U = Eigen::VectorXd::Zero(5);
    Eigen::VectorXd X0 = nesc_internal::buildState(X12, kin);
    DroppedSphere sphere(X0, U);

    Eigen::MatrixXd sol = nesc_internal::checkcase_rk4(sphere, kin, X0, sim_time, n_steps);

    // Conservation diagnostics
    const double dt = sim_time / n_steps;
    double E0 = sphere.mass * g * h0;   // initial total mechanical energy

    double E_min = E0, E_max = E0;
    double v_horiz_max = 0.0;
    int posZ = kin.position_start_idx() + 2;

    std::ofstream f("nesc_case1_sphere.csv");
    f << "time_s,h_m,v_vert_mps,v_horiz_mps,E_total_J,E_drift_rel\n";
    f << std::fixed << std::setprecision(8);

    for (int i = 0; i <= n_steps; ++i) {
        double t = i * dt;
        double u = sol(0, i), v = sol(1, i), w = sol(2, i);
        double z = sol(posZ, i);
        double h = -z;
        double Vsq = u * u + v * v + w * w;
        double E = sphere.mass * g * h + 0.5 * sphere.mass * Vsq;
        double E_drift = (E - E0) / E0;

        // In NED, positive w_NED is downward velocity, but body u,v,w are body frame.
        // For a falling sphere with no rotation, body and NED frames remain aligned,
        // so vertical velocity = w_body.
        double v_horiz = std::sqrt(u * u + v * v);
        v_horiz_max = std::max(v_horiz_max, v_horiz);
        E_min = std::min(E_min, E); E_max = std::max(E_max, E);

        f << t << "," << h << "," << w << "," << v_horiz << "," << E << "," << E_drift << "\n";
    }
    f.close();

    std::cout << "\n--- NESC Case 1: Dropped Sphere (no drag) ---\n"
        << "  sim time:                " << sim_time << " s\n"
        << "  initial altitude:        " << h0 << " m\n"
        << "  initial energy E0:       " << E0 << " J\n"
        << "  energy drift (rel):      " << ((E_max - E_min) / E0) << "\n"
        << "  max horizontal velocity: " << v_horiz_max << " m/s "
        << (v_horiz_max < 1e-9 ? "(OK)\n" : "(FAIL — should be 0)\n")
        << "  output: nesc_case1_sphere.csv\n";
}


// ============================================================
//  Case 2 — Tumbling brick (no damping)
// ============================================================
inline void runNESCCheckCase2(const AttitudeKinematics& kin)
{
    constexpr double sim_time = 30.0;
    constexpr int    n_steps = 6000;     // small dt for tight conservation

    // IC: at altitude, zero translational velocity, initial angular rates
    // designed to excite Euler's equations (rotation about all three axes)
    Eigen::VectorXd X12 = Eigen::VectorXd::Zero(12);
    X12[3] = 10.0 * nesc_internal::PI / 180.0;   // p = 10 deg/s
    X12[4] = 20.0 * nesc_internal::PI / 180.0;   // q = 20 deg/s
    X12[5] = 30.0 * nesc_internal::PI / 180.0;   // r = 30 deg/s
    X12[11] = -3048.0;              // 10,000 ft

    Eigen::VectorXd U = Eigen::VectorXd::Zero(5);
    Eigen::VectorXd X0 = nesc_internal::buildState(X12, kin);
    TumblingBrick brick(X0, U);

    Eigen::MatrixXd sol = nesc_internal::checkcase_rk4(brick, kin, X0, sim_time, n_steps);

    // Initial conservation quantities
    Eigen::Vector3d omega0(X12[3], X12[4], X12[5]);
    Eigen::Vector3d H0_body = brick.InertiaMatrix * omega0;
    double H0_mag = H0_body.norm();
    double T0_rot = 0.5 * omega0.dot(H0_body);

    const double dt = sim_time / n_steps;
    double H_min = H0_mag, H_max = H0_mag;
    double T_min = T0_rot, T_max = T0_rot;

    std::ofstream f("nesc_case2_brick.csv");
    f << "time_s,p_rads,q_rads,r_rads,H_inertial_mag,H_drift_rel,T_rot_J,T_drift_rel\n";
    f << std::fixed << std::setprecision(10);

    for (int i = 0; i <= n_steps; ++i) {
        double t = i * dt;
        Eigen::Vector3d omega(sol(3, i), sol(4, i), sol(5, i));
        Eigen::VectorXd att = sol.col(i).segment(kin.attitude_start_idx(), kin.n_attitude());

        Eigen::Matrix3d R_n_b = kin.body_to_ned(att);
        Eigen::Vector3d H_body = brick.InertiaMatrix * omega;
        Eigen::Vector3d H_inertial = R_n_b * H_body;
        double H_mag = H_inertial.norm();
        double T_rot = 0.5 * omega.dot(H_body);

        H_min = std::min(H_min, H_mag); H_max = std::max(H_max, H_mag);
        T_min = std::min(T_min, T_rot); T_max = std::max(T_max, T_rot);

        f << t << "," << omega[0] << "," << omega[1] << "," << omega[2]
            << "," << H_mag << "," << ((H_mag - H0_mag) / H0_mag)
            << "," << T_rot << "," << ((T_rot - T0_rot) / T0_rot) << "\n";
    }
    f.close();

    double H_drift = (H_max - H_min) / H0_mag;
    double T_drift = (T_max - T_min) / T0_rot;

    std::cout << "\n--- NESC Case 2: Tumbling Brick (no damping) ---\n"
        << "  sim time:                  " << sim_time << " s, dt = " << dt << " s\n"
        << "  initial omega (deg/s):     [" << X12[3] * 180 / nesc_internal::PI << ", "
        << X12[4] * 180 / nesc_internal::PI << ", "
        << X12[5] * 180 / nesc_internal::PI << "]\n"
        << "  initial |H| (kg·m²/s):     " << H0_mag << "\n"
        << "  initial T_rot (J):         " << T0_rot << "\n"
        << "  |H| drift (rel):           " << H_drift
        << (std::abs(H_drift) < 1e-4 ? " (OK)\n" : " (LARGE)\n")
        << "  T_rot drift (rel):         " << T_drift
        << (std::abs(T_drift) < 1e-4 ? " (OK)\n" : " (LARGE)\n")
        << "  output: nesc_case2_brick.csv\n";
}


// ============================================================
//  Case 3 — Tumbling brick with dynamic damping
// ============================================================
inline void runNESCCheckCase3(const AttitudeKinematics& kin)
{
    constexpr double sim_time = 30.0;
    constexpr int    n_steps = 6000;

    Eigen::VectorXd X12 = Eigen::VectorXd::Zero(12);
    X12[3] = 10.0 * nesc_internal::PI / 180.0;
    X12[4] = 20.0 * nesc_internal::PI / 180.0;
    X12[5] = 30.0 * nesc_internal::PI / 180.0;
    X12[11] = -3048.0;

    Eigen::VectorXd U = Eigen::VectorXd::Zero(5);
    Eigen::VectorXd X0 = nesc_internal::buildState(X12, kin);
    DampedTumblingBrick brick(X0, U);

    Eigen::MatrixXd sol = nesc_internal::checkcase_rk4(brick, kin, X0, sim_time, n_steps);

    Eigen::Vector3d omega0(X12[3], X12[4], X12[5]);
    Eigen::Vector3d H0_body = brick.InertiaMatrix * omega0;
    double H0_mag = H0_body.norm();
    double T0_rot = 0.5 * omega0.dot(H0_body);

    const double dt = sim_time / n_steps;

    std::ofstream f("nesc_case3_damped.csv");
    f << "time_s,p_rads,q_rads,r_rads,H_inertial_mag,T_rot_J,T_drift_rel\n";
    f << std::fixed << std::setprecision(10);

    bool monotonic = true;
    double T_prev = T0_rot;

    for (int i = 0; i <= n_steps; ++i) {
        double t = i * dt;
        Eigen::Vector3d omega(sol(3, i), sol(4, i), sol(5, i));
        Eigen::VectorXd att = sol.col(i).segment(kin.attitude_start_idx(), kin.n_attitude());

        Eigen::Matrix3d R_n_b = kin.body_to_ned(att);
        Eigen::Vector3d H_body = brick.InertiaMatrix * omega;
        Eigen::Vector3d H_inertial = R_n_b * H_body;
        double H_mag = H_inertial.norm();
        double T_rot = 0.5 * omega.dot(H_body);

        if (i > 0 && T_rot > T_prev + 1e-9) monotonic = false;
        T_prev = T_rot;

        f << t << "," << omega[0] << "," << omega[1] << "," << omega[2]
            << "," << H_mag << "," << T_rot << "," << ((T_rot - T0_rot) / T0_rot) << "\n";
    }
    f.close();

    double T_final = 0.5 * Eigen::Vector3d(sol(3, n_steps), sol(4, n_steps), sol(5, n_steps))
        .dot(brick.InertiaMatrix * Eigen::Vector3d(sol(3, n_steps), sol(4, n_steps), sol(5, n_steps)));

    std::cout << "\n--- NESC Case 3: Tumbling Brick (damped) ---\n"
        << "  sim time:                  " << sim_time << " s, dt = " << dt << " s\n"
        << "  initial T_rot (J):         " << T0_rot << "\n"
        << "  final   T_rot (J):         " << T_final << "\n"
        << "  fraction of energy lost:   " << (1.0 - T_final / T0_rot) << "\n"
        << "  T_rot monotonic decrease:  " << (monotonic ? "OK" : "FAIL") << "\n"
        << "  output: nesc_case3_damped.csv\n";
}


// ============================================================
//  Run all NESC check cases as a validation suite
// ============================================================
inline void runAllNESCCheckCases(const AttitudeKinematics& kin)
{
    std::cout << "\n========================================\n";
    std::cout << "  NESC 6-DOF Verification Check Cases\n";
    std::cout << "  Attitude representation: " << kin.name() << "\n";
    std::cout << "========================================\n";

    runNESCCheckCase1(kin);
    runNESCCheckCase2(kin);
    runNESCCheckCase3(kin);

    std::cout << "\nAll NESC check cases complete.\n";
    std::cout << "========================================\n\n";
}