
#pragma once

#include "EngineModel.h"
#include <Eigen/Dense>
#include <numbers>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>


// ============================================================
//  AircraftModel — abstract base class for any aircraft.
//
//  Holds state and properties common to every airframe (body
//  state, position, controls, mass, inertia tensor, control
//  schedule helpers) plus a vector of EngineModels.
//
//  Derived classes implement the airframe-specific aerodynamics
//  and moment buildup, and populate engines_ in their constructor.
//  Engine force and engine moment are then computed by the base
//  class by looping over engines_, so single- and twin-engine
//  aircraft share the same propulsion code.
// ============================================================

class AircraftModel
{
public:
    virtual ~AircraftModel() = default;

    virtual const char* name() const = 0;

    // ---- Common body / position / control state ----
    double u{}, v{}, w{};
    double p{}, q{}, r{};
    double x_b{}, y_b{}, z_b{};

    double aileron_inp{}, stabilizer_inp{}, rudder_inp{};
    double throttle1_inp{}, throttle2_inp{};

    // ---- Common physical properties ----
    double mass{};
    Eigen::Matrix3d InertiaMatrix{};
    Eigen::Matrix3d InverseInertiaMatrix{};

    double pi = std::numbers::pi;

    // ---- Propulsion ----
    // Populated by derived constructors. Each engine carries its own
    // position, thrust axis, and throttle index.
    std::vector<std::unique_ptr<EngineModel>> engines_;


    // ---- State / control update ----
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
        control_vector_ = U;   // keep the full vector for engine throttle lookup
    }


    // ---- Aerodynamics (airframe-specific) ----
    virtual Eigen::Vector3d aerodynamic_forces_stability_axis(
        double alpha, double beta, double airSpeed, double dynamicPressure) = 0;

    virtual Eigen::Vector3d cg_moments_body_frame(
        Eigen::Matrix3d& rot_stab_to_body,
        Eigen::Vector3d& omega_b,
        double airSpeed, double dynamicPressure,
        double alpha, double beta, double g) = 0;


    // ---- Propulsion (shared, computed from engines_) ----
    // Total body-frame thrust force, summed over all engines.
    // Throttle for each engine is read from its own control index.
    virtual Eigen::Vector3d engine_forces_body_frame(double g)
    {
        Eigen::Vector3d F = Eigen::Vector3d::Zero();
        const double V = std::sqrt(u * u + v * v + w * w);
        const double rho = last_rho_;   // set by the caller via setAirDensity (optional)
        for (const auto& eng : engines_) {
            double throttle = control_vector_.size() > eng->throttle_index()
                ? control_vector_[eng->throttle_index()]
                : 0.0;
            F += eng->force_body(throttle, V, rho, mass, g);
        }
        return F;
    }

    // Total engine moment about the CG, summed over all engines:
    //   M = sum_i  r_i x F_i
    Eigen::Vector3d engine_moments_body_frame(double g)
    {
        Eigen::Vector3d M = Eigen::Vector3d::Zero();
        const double V = std::sqrt(u * u + v * v + w * w);
        const double rho = last_rho_;
        for (const auto& eng : engines_) {
            double throttle = control_vector_.size() > eng->throttle_index()
                ? control_vector_[eng->throttle_index()]
                : 0.0;
            Eigen::Vector3d Fi = eng->force_body(throttle, V, rho, mass, g);
            M += eng->position_cg().cross(Fi);
        }
        return M;
    }

    // Let the simulation pass current air density to engines.
    void setAirDensity(double rho) { last_rho_ = rho; }


    // Convenience: stability-to-body rotation applied to the aero force vector
    Eigen::Vector3d aerodynamic_forces_body_frame(
        double alpha, double beta, double airSpeed, double dynamicPressure,
        Eigen::Matrix3d& rot_stab_to_body)
    {
        return rot_stab_to_body * aerodynamic_forces_stability_axis(
            alpha, beta, airSpeed, dynamicPressure);
    }


    // ---- Control schedule helpers ----
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
            if (t >= 0 && t < phase_time1)                    U_inputs(U_index, i) = pos_rad;
            else if (t >= phase_time1 && t < 2 * phase_time2) U_inputs(U_index, i) = neg_rad;
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


    // ---- Metadata ----
    virtual std::map<std::string, double> getAircraftSpecs() const = 0;

    std::map<std::string, Eigen::Matrix3d> getInertiaMatrices() const
    {
        std::map<std::string, Eigen::Matrix3d> m;
        m["inertia_matrix"] = InertiaMatrix;
        m["inverse_inertia_matrix"] = InverseInertiaMatrix;
        return m;
    }

protected:
    Eigen::VectorXd control_vector_{};   // full control vector from last update
    double          last_rho_ = 1.225;   // air density for thrust models that use it
};