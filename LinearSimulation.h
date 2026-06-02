
#pragma once

#include "AircraftModel.h"
#include "Constants.h"
#include "NumericalIntegration.h"
#include "SimulationEngine.h"   // for LinearStateSpace
#include <Eigen/Dense>
#include <functional>


// ============================================================
//  Linear subsystem simulator
//
//  Drives a precomputed (A_sub, B_sub) state-space block with a
//  control schedule built from doublet / sinusoidal commands on
//  the appropriate control surfaces. Returns the state-history
//  matrix and writes the control schedule to the output argument
//  so it can be saved alongside the response.
//
//  Used for both longitudinal (stabilizer doublet) and lateral
//  (aileron doublet) excitation campaigns.
// ============================================================

inline Eigen::MatrixXd runLinearSubsystemSim(AircraftModel& ac,
    const Eigen::MatrixXd& A_sub,
    const Eigen::MatrixXd& B_sub,
    const Eigen::VectorXd& X0,
    double sim_time, int steps,
    Eigen::MatrixXd& U_schedule_out,
    int stab_pos, int stab_neg,
    int ail_pos, int ail_neg,
    double thrust_amp, double thrust_period)
{
    const double dt = sim_time / steps;
    U_schedule_out.resize(constants::NUM_CONTROLS, steps + 1);
    U_schedule_out.setZero();

    ac.initialize_deflections(0, U_schedule_out, ail_pos, ail_neg, ail_pos, ail_neg, steps, dt);
    ac.initialize_deflections(1, U_schedule_out, stab_pos, stab_neg, stab_pos, stab_neg, steps, dt);
    ac.initialize_deflections(2, U_schedule_out, 0, 0, 0, 0, steps, dt);
    ac.initialize_thrusters(3, U_schedule_out, thrust_amp, thrust_period, steps, dt);
    ac.initialize_thrusters(4, U_schedule_out, thrust_amp, thrust_period, steps, dt);

    std::function<Eigen::VectorXd(const Eigen::VectorXd&, const Eigen::VectorXd&,
        const Eigen::MatrixXd&, const Eigen::MatrixXd&)>
        ss_func = LinearStateSpace;

    return rk4_simulate_ss(ss_func, X0, U_schedule_out, 0.0, sim_time, steps, A_sub, B_sub);
}