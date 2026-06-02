
#pragma once

#include "AircraftModel.h"
#include <Eigen/Dense>
#include <cmath>


// ============================================================
//  NESC verification check cases
//
//  Reference: NESC-RP-12-00770 / NASA-TM-2015-218675
//             https://nescacademy.nasa.gov/flightsim
//             https://github.com/nasa/simupy-flight
//
//  These are rigid-body fixtures used
//  to verify the equations-of-motion implementation. They share
//  the AircraftModel interface so the same simulation pipeline
//  can run them with no special-casing in Aircraft_Sim.
//
//  The official NESC trajectories use J2 gravity and
//  the WGS-84 ellipsoid; this sim uses constant gravity and
//  flat-Earth navigation, so bit-for-bit match against the NESC
//  reference data is not possible. However the cases still test
//  the relevant physics:
//    Case 1 — translational EOM, gravity, integration
//    Case 2 — rotational EOM, conservation of H and rotational KE
//    Case 3 — rotational EOM with rate damping (energy dissipation)
// ============================================================


// ============================================================
//  Case 1 — Dropped sphere with no drag
//
//  No aerodynamics, no torques. Trajectory should be pure
//  free-fall with constant horizontal velocity.
// ============================================================
class DroppedSphere : public AircraftModel
{
public:
    DroppedSphere(Eigen::VectorXd& states, Eigen::VectorXd& controls)
    {
        initializeStates(states);
        intializeControlInputs(controls);

        // Mass cancels out for a no-drag sphere under gravity, so the
        // exact value doesn't affect the trajectory.
        mass = 1.0;

        // Inertia is irrelevant (no torques). Use unit-diagonal.
        InertiaMatrix = Eigen::Matrix3d::Identity();
        InverseInertiaMatrix = Eigen::Matrix3d::Identity();
    }

    const char* name() const override { return "NESC-Case1-Sphere"; }

    Eigen::Vector3d aerodynamic_forces_stability_axis(double, double, double, double) override
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d engine_forces_body_frame(double) override
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d cg_moments_body_frame(Eigen::Matrix3d&,
        Eigen::Vector3d& omega_b,
        double, double, double, double, double) override
    {
        // No external torque; only inertial coupling (which is zero
        // for a body with isotropic inertia and any ω).
        return InverseInertiaMatrix * (-omega_b.cross(InertiaMatrix * omega_b));
    }

    std::map<std::string, double> getAircraftSpecs() const override
    {
        return { {"mass", mass} };
    }
};


// ============================================================
//  Case 2 — Tumbling brick with no damping, no drag
//
//  Properties from arXiv:1905.09794 / SimuPy-Flight:
//    m   = 0.155404754 slug      = 2.268 kg
//    Ixx = 0.001894220 slug·ft²  = 0.002568 kg·m²
//    Iyy = 0.006211019 slug·ft²  = 0.008420 kg·m²
//    Izz = 0.007194665 slug·ft²  = 0.009754 kg·m²
//
//  Tests Euler's equations: angular momentum and rotational KE
//  must both be conserved exactly (to RK4 precision) with zero
//  external torque.
// ============================================================
class TumblingBrick : public AircraftModel
{
public:
    TumblingBrick(Eigen::VectorXd& states, Eigen::VectorXd& controls)
    {
        initializeStates(states);
        intializeControlInputs(controls);

        mass = 2.268;
        InertiaMatrix << 0.002568, 0, 0,
            0, 0.008420, 0,
            0, 0, 0.009754;
        InverseInertiaMatrix = InertiaMatrix.inverse();
    }

    const char* name() const override { return "NESC-Case2-TumblingBrick"; }

    Eigen::Vector3d aerodynamic_forces_stability_axis(double, double, double, double) override
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d engine_forces_body_frame(double) override
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d cg_moments_body_frame(Eigen::Matrix3d&,
        Eigen::Vector3d& omega_b,
        double, double, double, double, double) override
    {
        // Pure Euler dynamics: no external torque, only the
        //  -ω × (Iω)  inertial-coupling term inside Iω̇.
        return InverseInertiaMatrix * (-omega_b.cross(InertiaMatrix * omega_b));
    }

    std::map<std::string, double> getAircraftSpecs() const override
    {
        return { {"mass", mass}, {"Ixx", InertiaMatrix(0,0)},
                 {"Iyy", InertiaMatrix(1,1)}, {"Izz", InertiaMatrix(2,2)} };
    }
};


// ============================================================
//  Case 3 — Tumbling brick with dynamic damping, no drag
//
//  Same mass and inertia as Case 2, plus a linear rate-damping
//  torque  M_damp = -K_d · ω  in the body frame.
//
//  Tests that rate-damping torques are integrated correctly:
//  rotational KE should decrease monotonically while the body
//  continues to tumble.
// ============================================================
class DampedTumblingBrick : public AircraftModel
{
public:
    DampedTumblingBrick(Eigen::VectorXd& states, Eigen::VectorXd& controls)
    {
        initializeStates(states);
        intializeControlInputs(controls);

        mass = 2.268;
        InertiaMatrix << 0.002568, 0, 0,
            0, 0.008420, 0,
            0, 0, 0.009754;
        InverseInertiaMatrix = InertiaMatrix.inverse();

        // Linear rate-damping (N·m / (rad/s))
        // Roughly 10% of inertia magnitude → time constant ~10 s
        damping_matrix << 0.0003, 0, 0,
            0, 0.0008, 0,
            0, 0, 0.0010;
    }

    const char* name() const override { return "NESC-Case3-DampedBrick"; }

    Eigen::Vector3d aerodynamic_forces_stability_axis(double, double, double, double) override
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d engine_forces_body_frame(double) override
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d cg_moments_body_frame(Eigen::Matrix3d&,
        Eigen::Vector3d& omega_b,
        double, double, double, double, double) override
    {
        Eigen::Vector3d M_damp = -damping_matrix * omega_b;
        return InverseInertiaMatrix * (M_damp - omega_b.cross(InertiaMatrix * omega_b));
    }

    std::map<std::string, double> getAircraftSpecs() const override
    {
        return { {"mass", mass}, {"Ixx", InertiaMatrix(0,0)},
                 {"Iyy", InertiaMatrix(1,1)}, {"Izz", InertiaMatrix(2,2)} };
    }

    Eigen::Matrix3d damping_matrix = Eigen::Matrix3d::Zero();
};