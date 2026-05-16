
#pragma once
#include <numbers>


// ============================================================
//  Shared numeric constants used across the simulation modules.
// ============================================================
namespace constants {
    constexpr double PI = std::numbers::pi;
    constexpr double DEG_TO_RAD = PI / 180.0;
    constexpr double RAD_TO_DEG = 180.0 / PI;
    constexpr double METERS_TO_FEET = 3.28084;
    constexpr int    NUM_CONTROLS = 5;
}