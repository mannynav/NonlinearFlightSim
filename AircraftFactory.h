
#pragma once

#include "AircraftModel.h"
#include "Aircraft.h"
#include "NasaGTMAircraft.h"
#include "F16Aircraft.h"
#include <Eigen/Dense>
#include <memory>
#include <stdexcept>


// ============================================================
//  Aircraft factory
//
//  Single place that knows about every concrete aircraft type.
//  SimConfig and main() depend only on this header for the
//  enum and the factory function.
// ============================================================
enum class AircraftType { CivilAircraft, NasaGTM, F16 };

inline std::unique_ptr<AircraftModel> makeAircraft(AircraftType type,
    Eigen::VectorXd& X,
    Eigen::VectorXd& U)
{
    switch (type) {
    case AircraftType::CivilAircraft: return std::make_unique<CivilAircraft>(X, U);
    case AircraftType::NasaGTM:       return std::make_unique<NasaGTMAircraft>(X, U);
    case AircraftType::F16:           return std::make_unique<F16Aircraft>(X, U);
    }
    throw std::invalid_argument("Unknown aircraft type");
}