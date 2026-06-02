
#pragma once

#include "Input.h"     // readNumericCSVColumn
#include <algorithm>   // std::upper_bound
#include <cmath>       // std::isfinite
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


// ============================================================
//  Atmosphere Data
//
//  Tabulated US Standard Atmosphere 1976: altitude, speed of
//  sound, gravity, and air density, plus 1D-interpolated
//  lookups by altitude in meters. All quantities in SI.
//
//  Tables must be monotonically increasing in altitude. Queries
//  outside the table range are clamped to the nearest endpoint.
//
//  Queries with NaN/Inf altitude throw a runtime_error
// ============================================================

namespace atm_detail {

    inline double linear_interp_1d(const std::vector<double>& xs,
        const std::vector<double>& ys,
        double x)
    {
        if (!std::isfinite(x)) {
            throw std::runtime_error(
                "AtmosphereData lookup at non-finite altitude (NaN/Inf). "
                "Sim state has likely diverged - check timestep, trim, or aircraft scale.");
        }

        const size_t n = xs.size();
        if (n == 0) return 0.0;
        if (n == 1) return ys[0];
        if (x <= xs.front()) return ys.front();
        if (x >= xs.back())  return ys.back();

        auto it = std::upper_bound(xs.begin(), xs.end(), x);   // first xs[i] > x
        const size_t hi = static_cast<size_t>(it - xs.begin());
        const size_t lo = hi - 1;

        const double t = (x - xs[lo]) / (xs[hi] - xs[lo]);
        return ys[lo] + t * (ys[hi] - ys[lo]);
    }

}  // namespace atm_detail


struct AtmosphereData
{
    std::vector<double> altitude_m;          // m
    std::vector<double> speed_of_sound_mps;  // m/s
    std::vector<double> gravity_mps2;        // m/s^2
    std::vector<double> air_density_kgpm3;   // kg/m^3

    double speedOfSound(double alt_m) const {
        return atm_detail::linear_interp_1d(altitude_m, speed_of_sound_mps, alt_m);
    }
    double gravity(double alt_m) const {
        return atm_detail::linear_interp_1d(altitude_m, gravity_mps2, alt_m);
    }
    double airDensity(double alt_m) const {
        return atm_detail::linear_interp_1d(altitude_m, air_density_kgpm3, alt_m);
    }

    size_t size() const { return altitude_m.size(); }
};


inline AtmosphereData loadUSStandardAtmosphere(
    const std::string& alt_file = "gravity_atm_data_alt_m.csv",
    const std::string& c_file = "gravity_atm_data_c_mps.csv",
    const std::string& g_file = "gravity_atm_data_g_mps2.csv",
    const std::string& rho_file = "gravity_atm_data_rho_kgpm3.csv")
{
    AtmosphereData data;
    data.altitude_m = readNumericCSVColumn(alt_file);
    data.speed_of_sound_mps = readNumericCSVColumn(c_file);
    data.gravity_mps2 = readNumericCSVColumn(g_file);
    data.air_density_kgpm3 = readNumericCSVColumn(rho_file);

    const size_t n = data.altitude_m.size();
    if (n == 0) {
        throw std::runtime_error("Atmosphere altitude table is empty: " + alt_file);
    }
    if (data.speed_of_sound_mps.size() != n ||
        data.gravity_mps2.size() != n ||
        data.air_density_kgpm3.size() != n) {
        throw std::runtime_error("Atmosphere tables have inconsistent sizes.");
    }

    std::cout << "Loaded US Standard Atmosphere: " << n
        << " altitude points (" << data.altitude_m.front()
        << " - " << data.altitude_m.back() << " m)\n";
    return data;
}