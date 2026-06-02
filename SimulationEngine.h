#pragma once

#include "AircraftModel.h"
#include "AtmosphereModel.h"
#include "AttitudeKinematics.h"
#include "Constants.h"
#include <Eigen/Dense>
#include <cmath>


// ============================================================
//  SimulationEngine
//
//  Core nonlinear EOM, finite-difference linearization, and the
//  RK4 driver. None of these depend on which aircraft model or
//  attitude representation is in use — they take both as
//  parameters and dispatch through the AircraftModel and
//  AttitudeKinematics interfaces.
//
//  Environment (air density and gravity) is supplied via the
//  AtmosphereData table and looked up at the current altitude
//  on every evaluation.
// ============================================================


// Lift a 12-element Euler-style IC into the N-state vector
// required by the chosen attitude representation.
inline Eigen::VectorXd buildFullState(const Eigen::VectorXd& X12_euler,
    const AttitudeKinematics& kin)
{
    Eigen::VectorXd X(kin.n_states());
    X.head(6) = X12_euler.head(6);
    Eigen::VectorXd att = kin.from_euler(X12_euler[6], X12_euler[7], X12_euler[8]);
    X.segment(kin.attitude_start_idx(), kin.n_attitude()) = att;
    X.tail(3) = X12_euler.tail(3);
    return X;
}


// Full nonlinear EOM — returns Xdot for a single state evaluation
inline Eigen::VectorXd Aircraft_Sim(AircraftModel& aircraft,
    const AttitudeKinematics& kin,
    const Eigen::VectorXd& X,
    const Eigen::VectorXd& U,
    const AtmosphereData& atm)
{
    aircraft.initializeStates(X);
    aircraft.intializeControlInputs(U);

    double u = X[0], v = X[1], w = X[2];
    double p = X[3], q = X[4], r = X[5];
    Eigen::VectorXd att = X.segment(kin.attitude_start_idx(), kin.n_attitude());

    // ---- Atmosphere lookup at current altitude (NED: altitude = -z_pos) ----
    const double z_pos = X[kin.position_start_idx() + 2];
    const double altitude = -z_pos;
    const double rho = atm.airDensity(altitude);
    const double g = atm.gravity(altitude);

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


// Implicit residual: F(Xdot, X, U) = f(X, U) - Xdot
inline Eigen::VectorXd Implicit_Model(AircraftModel& ac,
    const AttitudeKinematics& kin,
    const Eigen::VectorXd& XDOT,
    const Eigen::VectorXd& X,
    const Eigen::VectorXd& U,
    const AtmosphereData& atm)
{
    return Aircraft_Sim(ac, kin, X, U, atm) - XDOT;
}


// A matrix via central differences in state space
inline Eigen::MatrixXd LinearizeSystem_A(AircraftModel& ac,
    const AttitudeKinematics& kin,
    const Eigen::VectorXd& XDOTo,
    const Eigen::VectorXd& Xo,
    const Eigen::VectorXd& Uo,
    const Eigen::MatrixXd& DX,
    const AtmosphereData& atm)
{
    int n = static_cast<int>(XDOTo.size());
    Eigen::MatrixXd A(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double dx = DX(i, j);
            Eigen::VectorXd x_plus = Xo;  x_plus(j) += dx;
            Eigen::VectorXd x_minus = Xo;  x_minus(j) -= dx;
            double Fp = Implicit_Model(ac, kin, XDOTo, x_plus, Uo, atm)[i];
            double Fm = Implicit_Model(ac, kin, XDOTo, x_minus, Uo, atm)[i];
            A(i, j) = (Fp - Fm) / (2.0 * dx);
        }
    }
    return A;
}


// B matrix via central differences in control space
inline Eigen::MatrixXd LinearizeSystem_B(AircraftModel& ac,
    const AttitudeKinematics& kin,
    const Eigen::VectorXd& XDOTo,
    const Eigen::VectorXd& Xo,
    const Eigen::VectorXd& Uo,
    const Eigen::MatrixXd& DU,
    const AtmosphereData& atm)
{
    int n = static_cast<int>(XDOTo.size());
    int m = static_cast<int>(Uo.size());
    Eigen::MatrixXd B(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double du = DU(i, j);
            Eigen::VectorXd u_plus = Uo;  u_plus(j) += du;
            Eigen::VectorXd u_minus = Uo;  u_minus(j) -= du;
            double Fp = Implicit_Model(ac, kin, XDOTo, Xo, u_plus, atm)[i];
            double Fm = Implicit_Model(ac, kin, XDOTo, Xo, u_minus, atm)[i];
            B(i, j) = (Fp - Fm) / (2.0 * du);
        }
    }
    return B;
}


// Linear state-space RHS: xdot = A·x + B·u
inline Eigen::VectorXd LinearStateSpace(const Eigen::VectorXd& x,
    const Eigen::VectorXd& u,
    const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& B)
{
    return A * x + B * u;
}


// RK4 nonlinear simulation with post-step attitude renormalization
inline Eigen::MatrixXd runNonlinearSimulation(AircraftModel& ac,
    const AttitudeKinematics& kin,
    const Eigen::VectorXd& X0,
    const Eigen::VectorXd& U0,
    const AtmosphereData& atm,
    double sim_length, int steps)
{
    const double dt = sim_length / steps;
    Eigen::MatrixXd sol(kin.n_states(), steps + 1);
    sol.setZero();
    sol.col(0) = X0;

    Eigen::VectorXd U = U0;

    for (int i = 0; i < steps; i++) {
        Eigen::VectorXd Xk = sol.col(i);
        Eigen::VectorXd k1 = Aircraft_Sim(ac, kin, Xk, U, atm);
        Eigen::VectorXd k2 = Aircraft_Sim(ac, kin, (Xk + 0.5 * dt * k1).eval(), U, atm);
        Eigen::VectorXd k3 = Aircraft_Sim(ac, kin, (Xk + 0.5 * dt * k2).eval(), U, atm);
        Eigen::VectorXd k4 = Aircraft_Sim(ac, kin, (Xk + dt * k3).eval(), U, atm);
        Eigen::VectorXd Xn = Xk + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

        Eigen::VectorXd att = Xn.segment(kin.attitude_start_idx(), kin.n_attitude());
        Xn.segment(kin.attitude_start_idx(), kin.n_attitude()) = kin.normalize(att);

        sol.col(i + 1) = Xn;
    }
    return sol;
}