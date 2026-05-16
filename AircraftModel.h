#pragma once

#include <Eigen/Dense>
#include <numbers>
#include <map>
#include <stdexcept>
#include <string>


// ============================================================
//  AircraftModel — abstract base class for any aircraft.
//
//  Holds state and properties that are common to every airframe
//  (body-frame state, position, controls, mass, inertia tensor,
//  control-schedule helpers).
//
//  Derived classes implement the three pure virtual methods that
//  describe the airframe-specific aerodynamics, propulsion, and
//  moment buildup.
// ============================================================

class AircraftModel
{
public:
    virtual ~AircraftModel() = default;

    // Identifier for diagnostics
    virtual const char* name() const = 0;

    // ---- Common body / position / control state ----
    // Public for convenience; written via initializeStates / intializeControlInputs.
    double u{}, v{}, w{};
    double p{}, q{}, r{};
    double x_b{}, y_b{}, z_b{};

    double aileron_inp{}, stabilizer_inp{}, rudder_inp{};
    double throttle1_inp{}, throttle2_inp{};

    // ---- Common physical properties ----
    double mass{};
    Eigen::Matrix3d InertiaMatrix{};
    Eigen::Matrix3d InverseInertiaMatrix{};

    // pi shorthand for derived constructors
    double pi = std::numbers::pi;


    // ---- State / control update (same logic for every aircraft) ----
    // Body state at indices 0-5 and position at the last 3 indices,
    // regardless of attitude representation in between.
    void initializeStates(const Eigen::VectorXd& X)
    {
        const int n = static_cast<int>(X.size());
        u = X[0]; v = X[1]; w = X[2];
        p = X[3]; q = X[4]; r = X[5];
        x_b = X[n - 3]; y_b = X[n - 2]; z_b = X[n - 1];
    }

    void intializeControlInputs(const Eigen::VectorXd& U)
    {
        aileron_inp = U[0];
        stabilizer_inp = U[1];
        rudder_inp = U[2];
        throttle1_inp = U[3];
        throttle2_inp = U[4];
    }


    // ---- Force / moment interface (airframe-specific) ----
    virtual Eigen::Vector3d aerodynamic_forces_stability_axis(
        double alpha, double beta, double airSpeed, double dynamicPressure) = 0;

    virtual Eigen::Vector3d engine_forces_body_frame(double g) = 0;

    virtual Eigen::Vector3d cg_moments_body_frame(
        Eigen::Matrix3d& rot_stab_to_body,
        Eigen::Vector3d& omega_b,
        double airSpeed, double dynamicPressure,
        double alpha, double beta, double g) = 0;


    // Convenience: stability-to-body rotation applied to the aero force vector
    Eigen::Vector3d aerodynamic_forces_body_frame(
        double alpha, double beta, double airSpeed, double dynamicPressure,
        Eigen::Matrix3d& rot_stab_to_body)
    {
        return rot_stab_to_body * aerodynamic_forces_stability_axis(
            alpha, beta, airSpeed, dynamicPressure);
    }


    // ---- Control schedule helpers (common to all aircraft) ----
    void initialize_deflections(int U_index, Eigen::MatrixXd& U_inputs,
        double pos_def_deg, double neg_def_deg,
        double phase_time1, double phase_time2,
        double steps, double increment)
    {
        constexpr double DEG_TO_RAD = std::numbers::pi / 180.0;
        double pos_rad = pos_def_deg * DEG_TO_RAD;
        double neg_rad = neg_def_deg * DEG_TO_RAD;

        for (int i = 0; i <= steps; ++i) {
            double t = i * increment;
            if (t >= 0 && t < phase_time1)        U_inputs(U_index, i) = pos_rad;
            else if (t >= phase_time1 && t < 2 * phase_time2)    U_inputs(U_index, i) = neg_rad;
        }
    }

    void initialize_thrusters(int U_index, Eigen::MatrixXd& U_inputs,
        double thrust_value, double phase_time,
        double steps, double increment)
    {
        if (U_index != 3 && U_index != 4)
            throw std::out_of_range("U_index in initialize_thrusters must be 3 or 4");

        for (int i = 0; i <= steps; ++i) {
            double t = i * increment;
            if (t >= 0 && t < phase_time) U_inputs(U_index, i) = thrust_value;
        }
    }


    // ---- Optional metadata accessors (common implementations) ----
    virtual std::map<std::string, double> getAircraftSpecs() const = 0;

    std::map<std::string, Eigen::Matrix3d> getInertiaMatrices() const
    {
        std::map<std::string, Eigen::Matrix3d> m;
        m["inertia_matrix"] = InertiaMatrix;
        m["inverse_inertia_matrix"] = InverseInertiaMatrix;
        return m;
    }
};