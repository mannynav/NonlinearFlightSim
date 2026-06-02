
#pragma once

#include <Eigen/Dense>
#include <memory>
#include <vector>


// ============================================================
//  EngineModel
//
//  One physical engine. Knows:
//    - its mounting position relative to the CG (for moment arm)
//    - the body-frame unit vector its thrust acts along
//    - which throttle control index drives it
//    - how to map a throttle setting to a thrust magnitude
//
//  An aircraft owns a std::vector<std::unique_ptr<EngineModel>>.
//  Single-engine airframes push one; twin-engine push two. The
//  base AircraftModel can then compute total engine force and
//  moment by looping over the vector — no hard-coded engine count.
// ============================================================

class EngineModel
{
public:
    virtual ~EngineModel() = default;

    EngineModel(const Eigen::Vector3d& position_cg,
        int throttle_index,
        const Eigen::Vector3d& thrust_axis_body = Eigen::Vector3d(1, 0, 0))
        : position_cg_(position_cg)
        , throttle_index_(throttle_index)
        , thrust_axis_(thrust_axis_body.normalized())
    {
    }

    // Thrust magnitude (N) for a given throttle setting and flight condition.
    // airspeed, rho, mass, g are supplied so richer models can use them;
    // simple models ignore what they don't need.
    virtual double thrust_magnitude(double throttle, double airspeed,
        double rho, double mass, double g) const = 0;

    // Body-frame thrust force vector.
    Eigen::Vector3d force_body(double throttle, double airspeed,
        double rho, double mass, double g) const
    {
        return thrust_magnitude(throttle, airspeed, rho, mass, g) * thrust_axis_;
    }

    int                   throttle_index() const { return throttle_index_; }
    const Eigen::Vector3d& position_cg()   const { return position_cg_; }
    const Eigen::Vector3d& thrust_axis()   const { return thrust_axis_; }

protected:
    Eigen::Vector3d position_cg_;     // engine location relative to CG (body frame, m)
    int             throttle_index_;  // index into the control vector (3 or 4, typically)
    Eigen::Vector3d thrust_axis_;     // unit thrust direction in body frame
};


// ------------------------------------------------------------
//  GravityScaledEngine
//
//  The placeholder model used by CivilAircraft and NasaGTM:
//      T = throttle * mass * g
//  Thrust scales with weight, independent of airspeed/altitude.
//  Preserves the exact behavior those aircraft had before the
//  engine abstraction was introduced.
// ------------------------------------------------------------
class GravityScaledEngine : public EngineModel
{
public:
    using EngineModel::EngineModel;

    double thrust_magnitude(double throttle, double /*airspeed*/,
        double /*rho*/, double mass, double g) const override
    {
        return throttle * mass * g;
    }
};


// ------------------------------------------------------------
//  ConstantMaxEngine
//
//  Simple fixed-rating engine:
//      T = throttle * T_max
//  T_max is a constant sea-level-equivalent maximum thrust (N).
//  Used for the F-16 in the first pass (before any Mach/altitude
//  thrust table is added).
// ------------------------------------------------------------

class ConstantMaxEngine : public EngineModel
{
public:
    ConstantMaxEngine(const Eigen::Vector3d& position_cg,
        int throttle_index,
        double T_max_N,
        const Eigen::Vector3d& thrust_axis_body = Eigen::Vector3d(1, 0, 0))
        : EngineModel(position_cg, throttle_index, thrust_axis_body)
        , T_max_(T_max_N)
    {
    }

    double thrust_magnitude(double throttle, double /*airspeed*/,
        double /*rho*/, double /*mass*/, double /*g*/) const override
    {
        return throttle * T_max_;
    }

private:
    double T_max_;
};