#pragma once

#include "AtmosphereModel.h"
#include "AttitudeKinematics.h"
#include "Constants.h"
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


// ============================================================
//  Post-processing of the nonlinear solution
//
//  Converts the raw state-history matrix into a 29-row table of
//  derived quantities (NED velocities, Euler angles in degrees,
//  altitude in feet, angular rates in deg/s, atmospheric state,
//  Mach, etc.) suitable for plotting or downstream analysis.
//
//  Atmospheric quantities (SoS, density, Mach) are looked up
//  from the AtmosphereData table at the current altitude.
// ============================================================

inline Eigen::MatrixXd postProcessNonlinearSolution(const Eigen::MatrixXd& X,
    const AttitudeKinematics& kin,
    const AtmosphereData& atm,
    double sim_length, int steps)
{
    using namespace constants;
    const int N = steps + 1;
    const double dt = sim_length / steps;

    Eigen::MatrixXd P(29, N);
    P.setZero();

    for (int i = 0; i < N; ++i)
    {
        double u = X(0, i), v = X(1, i), w = X(2, i);
        double z = X(X.rows() - 1, i);

        Eigen::VectorXd att = X.col(i).segment(kin.attitude_start_idx(), kin.n_attitude());
        Eigen::Matrix3d R_n_b = kin.body_to_ned(att);
        Eigen::Vector3d v_NED = R_n_b * Eigen::Vector3d(u, v, w);
        Eigen::Vector3d eul = kin.to_euler(att);

        double altitude_m = -z;
        double airSpeed = std::sqrt(u * u + v * v + w * w);
        double alpha = (std::abs(u) < 1e-12) ? 0.0 : std::atan(w / u);
        double beta = (airSpeed > 1e-10) ? std::asin(v / airSpeed) : 0.0;

        // Atmosphere lookups at current altitude
        double c_mps = atm.speedOfSound(altitude_m);
        double rho = atm.airDensity(altitude_m);
        double mach = (c_mps > 1e-10) ? airSpeed / c_mps : 0.0;

        P(0, i) = i * dt;
        P(1, i) = altitude_m;
        P(2, i) = altitude_m * METERS_TO_FEET;
        P(3, i) = c_mps * METERS_TO_FEET;        // SoS in ft/s
        P(4, i) = rho * KGM3_TO_SLUGFT3;       // density in slug/ft^3
        P(5, i) = airSpeed * METERS_TO_FEET;        // true airspeed in ft/s
        P(6, i) = mach;
        P(7, i) = alpha;
        P(8, i) = alpha * RAD_TO_DEG;
        P(9, i) = beta;
        P(10, i) = beta * RAD_TO_DEG;
        P(11, i) = v_NED[0];
        P(12, i) = v_NED[1];
        P(13, i) = v_NED[2];
        P(14, i) = v_NED[0] * METERS_TO_FEET;
        P(15, i) = v_NED[1] * METERS_TO_FEET;
        P(16, i) = v_NED[2] * METERS_TO_FEET;
        P(17, i) = eul[0];
        P(18, i) = eul[1];
        P(19, i) = eul[2];
        P(20, i) = eul[0] * RAD_TO_DEG;
        P(21, i) = eul[1] * RAD_TO_DEG;
        P(22, i) = eul[2] * RAD_TO_DEG;
        P(23, i) = X(3, i);
        P(24, i) = X(4, i);
        P(25, i) = X(5, i);
        P(26, i) = X(3, i) * RAD_TO_DEG;
        P(27, i) = X(4, i) * RAD_TO_DEG;
        P(28, i) = X(5, i) * RAD_TO_DEG;
    }
    return P;
}


inline void writePostProcessCsv(const Eigen::MatrixXd& P, const std::string& filename)
{
    std::ofstream f(filename);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return;
    }
    f << std::fixed << std::setprecision(10);
    f << "Time,Altitude_m,Altitude_ft,SoS_ft_s,AirDensity_slugft3,TransVel_ft_s,Mach,"
        "AoA_rad,AoA_deg,AoS_rad,AoS_deg,"
        "u_NED_ms,v_NED_ms,w_NED_ms,u_NED_fts,v_NED_fts,w_NED_fts,"
        "phi_rad,theta_rad,psi_rad,phi_deg,theta_deg,psi_deg,"
        "p_rads,q_rads,r_rads,p_degs,q_degs,r_degs\n";
    for (int j = 0; j < P.cols(); ++j) {
        for (int i = 0; i < P.rows(); ++i) {
            f << P(i, j);
            if (i < P.rows() - 1) f << ",";
        }
        f << "\n";
    }
    std::cout << "Wrote " << filename << "\n";
}