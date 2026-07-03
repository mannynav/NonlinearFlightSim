#pragma once
#include <cmath>
#include <memory>


// ============================================================
//  GravityModel — pluggable gravity abstraction
//
//  Mirrors the EngineModel pattern: an abstract interface with
//  concrete implementations of increasing fidelity. The model
//  returns gravity MAGNITUDE (m/s^2) as a function of altitude;
//  direction remains local NED "down" (flat-Earth kinematics).
//
//  Latitude is a per-run configuration parameter, not a state:
//  for simulations covering < a few hundred km of ground track,
//  treating latitude as fixed is exact to well below the other
//  model errors. (A round-Earth navigation state would be the
//  next fidelity level and touches the EOM, not just gravity.)
//
//  Available models:
//
//  1. Tabulated (NESC CSV table)  — the existing behavior.
//     Selected by NOT attaching an override (nullptr). This
//     guarantees bit-perfect regression by construction.
//
//  2. InverseSquareGravity — analytic point-mass model
//         g(h) = g0 * ( r0 / (r0 + h) )^2
//     Defaults g0 = 9.80665 m/s^2, r0 = 6,356,766 m are the
//     US Standard Atmosphere 1976 constants; with these values
//     the model reproduces the NESC gravity table to ~6 digits
//     (the table was generated from this same law), which makes
//     it a useful cross-validation of the table data.
//
//  3. WGS84Gravity — ellipsoidal normal gravity (Somigliana):
//         gamma(phi) = gamma_e * (1 + k sin^2 phi)
//                      / sqrt(1 - e^2 sin^2 phi)
//     with the standard WGS-84 free-air altitude correction
//         g(phi,h) = gamma(phi) * [ 1
//                    - (2/a)(1 + f + m - 2 f sin^2 phi) h
//                    + (3/a^2) h^2 ]
//     Constants from NIMA TR8350.2 (WGS-84 definition).
//     Normal gravity already includes the centrifugal term for
//     an Earth-fixed observer, which is the correct "effective
//     gravity" for a flat-Earth NED simulation.
//     Surface gravity: 9.7803 m/s^2 (equator) -> 9.8322 (pole),
//     a 0.53% latitude variation the other models ignore.
// ============================================================

class GravityModel
{
public:
    virtual ~GravityModel() = default;

    // Gravity magnitude (m/s^2) at geometric altitude h (m).
    virtual double g(double h) const = 0;
    virtual const char* name() const = 0;
};


// ------------------------------------------------------------
//  InverseSquareGravity — spherical point-mass Earth
// ------------------------------------------------------------
class InverseSquareGravity : public GravityModel
{
public:
    explicit InverseSquareGravity(double g0 = 9.80665,      // m/s^2 (USSA76)
        double r0 = 6356766.0)     // m     (USSA76)
        : g0_(g0), r0_(r0)
    {
    }

    double g(double h) const override
    {
        const double ratio = r0_ / (r0_ + h);
        return g0_ * ratio * ratio;
    }

    const char* name() const override { return "InverseSquare (USSA76)"; }

private:
    double g0_, r0_;
};


// ------------------------------------------------------------
//  WGS84Gravity — Somigliana normal gravity + free-air correction
// ------------------------------------------------------------
class WGS84Gravity : public GravityModel
{
public:
    explicit WGS84Gravity(double latitude_rad)
        : lat_(latitude_rad)
    {
        // ---- WGS-84 defining/derived constants (NIMA TR8350.2) ----
        constexpr double A = 6378137.0;            // semi-major axis, m
        constexpr double GAMMA_E = 9.7803253359;         // equatorial normal gravity, m/s^2
        constexpr double K_SOM = 1.931852652458e-3;    // Somigliana's constant
        constexpr double E2 = 6.69437999014e-3;     // first eccentricity squared
        constexpr double F = 1.0 / 298.257223563;  // flattening
        constexpr double M_RATIO = 3.449786506841e-3;    // omega^2 a^2 b / GM

        const double s2 = std::sin(latitude_rad) * std::sin(latitude_rad);

        // Normal gravity on the ellipsoid surface at this latitude
        gamma_phi_ = GAMMA_E * (1.0 + K_SOM * s2) / std::sqrt(1.0 - E2 * s2);

        // Free-air altitude correction coefficients (precomputed)
        fa1_ = (2.0 / A) * (1.0 + F + M_RATIO - 2.0 * F * s2);
        fa2_ = 3.0 / (A * A);
    }

    double g(double h) const override
    {
        return gamma_phi_ * (1.0 - fa1_ * h + fa2_ * h * h);
    }

    const char* name() const override { return "WGS-84 (Somigliana)"; }

    double latitude_rad() const { return lat_; }

private:
    double lat_;
    double gamma_phi_;   // surface normal gravity at this latitude
    double fa1_, fa2_;   // free-air correction coefficients
};


// ------------------------------------------------------------
//  Selection enum + factory
//  Tabulated returns nullptr: "no override attached" means the
//  AtmosphereData table lookup is used, i.e. current behavior.
// ------------------------------------------------------------
enum class GravityModelType { Tabulated, InverseSquare, WGS84 };

inline std::shared_ptr<const GravityModel>
makeGravityModel(GravityModelType type, double latitude_deg)
{
    constexpr double D2R = 3.14159265358979323846 / 180.0;
    switch (type) {
    case GravityModelType::Tabulated:     return nullptr;
    case GravityModelType::InverseSquare: return std::make_shared<InverseSquareGravity>();
    case GravityModelType::WGS84:         return std::make_shared<WGS84Gravity>(latitude_deg * D2R);
    }
    return nullptr;
}