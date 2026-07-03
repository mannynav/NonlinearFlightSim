
#pragma once
#include "AircraftFactory.h"        // AircraftType
#include "Attitudekinematics.h"   // AttitudeMode
#include "Constants.h"
#include "GravityModel.h"


// ============================================================
//  SimConfig
//
//  All run-time configuration for one simulation pass, in one
//  place. main() builds (or selects) a SimConfig and hands it
//  to the pipeline; nothing downstream needs hard-coded tuning
//  constants anymore.
//
//  Convenience factory functions cover the common scenarios.
//  For a one-off:
//
//      SimConfig cfg = SimConfig::gtm_cruise();
//      cfg.trim_gamma = 6.0 * constants::DEG_TO_RAD;   // make it a climb
//
//  Timing differs by aircraft because GTM modes are ~10x faster
//  than the civil transport and need a finer timestep for stable
//  RK4 integration.
// ============================================================

struct SimConfig
{
    // ---- What to simulate ----
    AircraftType ac_type = AircraftType::CivilAircraft;
    AttitudeMode att_mode = AttitudeMode::Quaternion;

    // ---- Trim target ----
    double trim_Va = 85.0;     // m/s (true airspeed)
    double trim_h = 3048.0;   // m   (altitude)
    double trim_gamma = 0.0;      // rad (flight path angle; 0 = level)

    // ---- Nonlinear sim ----
    double nl_length = 150.0;    // s
    int    nl_steps = 200;
    double pitch_perturbation = 0.01;  // rad, added to trim theta as the IC

    // ---- Linear longitudinal sim ----
    double long_length = 150.0;   // s
    int    long_steps = 150;

    // ---- Linear lateral sim ----
    double lat_length = 200.0;   // s
    int    lat_steps = 250;

    // ---- Output verbosity ----
    bool verbose_trim = true;
    bool run_lqr = false;

    // ---- Gravity ----
    GravityModelType gravity_model = GravityModelType::Tabulated;
    double latitude_deg = 45.0;

    // Convenience: derived timestep
    double nl_dt()   const { return nl_length / nl_steps; }
    double long_dt() const { return long_length / long_steps; }
    double lat_dt()  const { return lat_length / lat_steps; }


    // ========================================================
    //  Scenario factories
    // ========================================================

    static SimConfig civil_cruise()
    {
        SimConfig c;
        c.ac_type = AircraftType::CivilAircraft;
        c.trim_Va = 85.0;
        c.trim_h = 3048.0;
        c.trim_gamma = 0.0;
        // Civil transport: slow modes, coarse timestep
        c.nl_length = 150.0;   c.nl_steps = 200;
        c.long_length = 150.0; c.long_steps = 150;
        c.lat_length = 200.0;  c.lat_steps = 250;
        return c;
    }

    static SimConfig gtm_cruise()
    {
        SimConfig c;
        c.ac_type = AircraftType::NasaGTM;
        c.trim_Va = 85.0;
        c.trim_h = 3048.0;
        c.trim_gamma = 0.0;
        // GTM: ~10x faster modes, need fine timestep (dt = 0.005 s)
        c.nl_length = 5.0;    c.nl_steps = 1000;
        c.long_length = 10.0; c.long_steps = 2000;
        c.lat_length = 15.0;  c.lat_steps = 3000;
        return c;
    }

    // Steady climb for either aircraft. gamma_deg > 0 climbs,
    // < 0 descends. Inherits timing from the matching cruise config.
    static SimConfig climb(AircraftType type, double gamma_deg)
    {
        SimConfig c = (type == AircraftType::NasaGTM) ? gtm_cruise() : civil_cruise();
        c.trim_gamma = gamma_deg * constants::DEG_TO_RAD;
        return c;
    }


    static SimConfig f16_cruise()
    {
        SimConfig c;
        c.ac_type = AircraftType::F16;
        c.trim_Va = 180.0;     // m/s (~M 0.55 at altitude)
        c.trim_h = 9000.0;    // m  (about 30,000 ft)
        c.trim_gamma = 0.0;
        // F-16 timing: short period ~3-5 rad/s, period ~1-2 s
        c.nl_length = 20.0;   c.nl_steps = 1000;   // dt = 0.02 s
        c.long_length = 30.0; c.long_steps = 1500;
        c.lat_length = 40.0;  c.lat_steps = 2000;
        return c;
    }


    static SimConfig f16_lqr_demo() {
        SimConfig c = f16_cruise();
        c.run_lqr = true;
        return c;
    }

};

