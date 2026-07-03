# Nonlinear Flight Dynamics Simulator

Six-degree-of-freedom nonlinear flight dynamics simulator in modern C++, supporting interchangeable airframes through a polymorphic aircraft model, with a Gauss-Newton trim solver, central-difference linearization, eigenvalue-based modal analysis, and an LQR autopilot demonstrated on the F-16 closed-loop nonlinear plant.

## Headline Results

Three aircraft trimmed at cruise, spanning four orders of magnitude in mass. All trims converge to machine precision in three iterations of damped Gauss-Newton:

| Aircraft | Mass | Short-period ωₙ | Dutch-roll ωₙ | Roll τ | Spiral |
|---|---|---|---|---|---|
| NASA GTM-T2 (subscale transport) | 26 kg | 14.6 rad/s | 11.4 rad/s | 0.085 s | stable, 24 s |
| F-16 (fighter) | 9,295 kg | 2.25 rad/s | 2.59 rad/s | 0.614 s | unstable, 177 s |
| Civil transport (heavy) | 120,000 kg | 1.59 rad/s | 0.696 rad/s | 1.02 s | stable, 7.5 s |

Each aircraft's modes match its class: well-damped slow modes for the heavy transport, fast lightly-damped modes for the subscale, and the famously unstable spiral mode of the F-16.

## What It Does

- 6-DOF nonlinear rigid-body dynamics with both quaternion and Euler-angle attitude representations, integrated via 4th-order Runge–Kutta
- Real NASA Generic Transport Model T2 polynomial aerodatabase driving the subscale aircraft, accessed through a templated N-dimensional lookup-table interpolator and an 8-table aerodynamic moment buildup
- Polynomial stability-derivative aerodynamic models for the civil transport and F-16, using published values from the standard flight-dynamics literature
- US Standard Atmosphere 1976 with interpolated density, speed of sound, and gravity over 1020 altitude points
- Steady-state trim via damped Gauss-Newton least squares with backtracking line search (17 unknowns, 18 residuals), supporting cruise and arbitrary flight-path-angle climb conditions
- Trim verification against 17 independent steady-flight conditions
- Central-difference Jacobian linearization at the trim point, with longitudinal/lateral mode decoupling via a permutation similarity transform
- Eigenvalue-based modal analysis (natural frequency, damping ratio, time constant) for the full system and decoupled longitudinal/lateral subsystems
- LQR autopilot design: continuous algebraic Riccati equation solved via Hamiltonian eigendecomposition, closed-loop response verified on the nonlinear plant

## Validation

The simulator passes all three NASA Engineering and Safety Center (NESC) 6-DOF benchmark check cases:

1. **Dropped sphere (no drag)** — mechanical energy conserved to machine precision over 30 s, zero horizontal velocity drift
2. **Tumbling brick (no damping)** — angular momentum and rotational kinetic energy conserved to machine precision
3. **Tumbling brick (damped)** — rotational energy decreases monotonically by the correct fraction

These are the standard verification suite for any 6-DOF simulator. Passing all three at machine precision confirms the integrator and rigid-body equations are correctly implemented.

## Architecture

The code is built around three layered abstractions:

**`AircraftModel`** (abstract base): holds body state, controls, mass, inertia tensor, and a `std::vector<std::unique_ptr<EngineModel>>`. Engine force and engine moment are computed by looping over the engine vector, so single- and twin-engine aircraft share the same propulsion code.

**`EngineModel`** (abstract): each engine knows its position relative to the CG, its thrust axis, its throttle index in the control vector, and how to map a throttle setting to a thrust force. Concrete: `GravityScaledEngine` (used by the transports) and `ConstantMaxEngine` (used by the F-16).

**`AttitudeKinematics`** (strategy pattern): swaps between quaternion (13-state) and Euler-angle (12-state) representations at runtime, with a longitudinal/lateral permutation matrix appropriate to each.

Run-time configuration is centralized in a `SimConfig` struct with factory functions (`civil_cruise()`, `gtm_cruise()`, `climb(type, gamma_deg)`, `f16_cruise()`, `f16_lqr_demo()`), so `main()` reduces to: load atmosphere, select a config, run.

## Gravity Models

Gravity is abstracted behind a `GravityModel` interface — the same pluggable pattern used for engines — with three implementations of increasing fidelity:

| Model | Law | Source |
|---|---|---|
| Tabulated (default) | NESC reference table, interpolated over 1020 altitude points | NESC check-case data |
| Inverse-square | g(h) = g₀·(r₀/(r₀+h))², g₀ = 9.80665 m/s², r₀ = 6,356,766 m | US Standard Atmosphere 1976 |
| WGS-84 | Somigliana normal gravity γ(φ) with second-order free-air altitude correction | NIMA TR8350.2 |

The tabulated model is selected by default and reproduces the pre-abstraction behavior exactly (regression-guaranteed by construction: no override attached means the original table lookup runs). The analytic inverse-square law reproduces the NESC table to ~6 significant digits, cross-validating the reference data against the physical law it was generated from.

The WGS-84 model adds the latitude dependence of gravity on the oblate Earth — a 0.53% variation from equator (9.7803 m/s²) to pole (9.8322 m/s²) that spherical models ignore. Normal gravity includes the centrifugal term for an Earth-fixed observer, which is the correct effective gravity for a flat-Earth NED simulation. Latitude is a per-run configuration parameter (`SimConfig::latitude_deg`); for trajectories under a few hundred kilometers of ground track, treating latitude as fixed introduces error well below the other model tolerances.

## Latitude sweep validation

F-16 trim at Va = 180 m/s, h = 9,000 m under WGS-84 gravity across three latitudes:

| Latitude | g (m/s²) | α (deg) | Stabilizer (deg) | Throttle |
|---|---|---|---|---|
| 0° (equator) | 9.75259 | 5.5083 | −1.8361 | 0.10596 |
| 45° | 9.77849 | 5.5228 | −1.8410 | 0.10623 |
| 90° (pole) | 9.80449 | 5.5374 | −1.8458 | 0.10650 |

The trim shifts are monotone and physically consistent: the aircraft is effectively 0.53% heavier at the pole, requiring ~0.5% more lift (higher α) and paying ~0.5% more induced drag (higher throttle). Trim converges quadratically in 3 iterations at every latitude.

The modal analysis provides a sharper theoretical check. Across the sweep, the phugoid natural frequency scales linearly with gravity to five significant digits (ω-ratio 1.00530 vs. g-ratio 1.00532, matching the classical approximation ωₙ ≈ √2·g/V), while the short-period frequency is invariant to the fifth digit — correct, since the short period is set by aerodynamic stiffness and pitch inertia, not weight. The gravity model propagates through trim, linearization, and modal structure exactly as classical flight dynamics predicts.


## LQR Autopilot

The F-16 longitudinal autopilot is designed by solving the continuous-time algebraic Riccati equation
```
Aᵀ P + P A − P B R⁻¹ Bᵀ P + Q = 0
```
via Hamiltonian eigendecomposition (Potter's method): the 2n × 2n Hamiltonian matrix has eigenvalues in stable/unstable pairs, and stacking the eigenvectors corresponding to the stable half-plane gives `P = V₂₁ V₁₁⁻¹`. The optimal feedback gain is `K = R⁻¹ Bᵀ P`. State and control weights follow Bryson's rule.

Result on the F-16 longitudinal subsystem at Va = 180 m/s, h = 9 km:

| Mode | Open loop | Closed loop |
|---|---|---|
| Short period | ωₙ = 2.25 rad/s, ζ = 0.23 | real poles at −24.7 and −3.75 |
| Phugoid | ωₙ = 0.075 rad/s, ζ = 0.07 | real poles at −0.22 and −0.11 |

All four closed-loop poles are on the real axis and stable. The controller is then applied to the nonlinear plant (not the linearized model) under a pitch perturbation. Open- and closed-loop responses are written to CSV for side-by-side comparison.

## Building and Running

Built with MSVC / Visual Studio 2022. Dependencies: Eigen (header-only). Open the solution, build, run.

In `main.cpp`, select one of the available scenarios:

```cpp
SimConfig cfg = SimConfig::gtm_cruise();
// SimConfig cfg = SimConfig::civil_cruise();
// SimConfig cfg = SimConfig::f16_cruise();
// SimConfig cfg = SimConfig::climb(AircraftType::NasaGTM, 6.0);
// SimConfig cfg = SimConfig::f16_lqr_demo();
runSimulation(cfg, atm);
```

Operating-point sweeps are a single loop:

```cpp
for (double gam : {0.0, 3.0, 6.0, 10.0, 15.0}) {
    runSimulation(SimConfig::climb(AircraftType::NasaGTM, gam), atm);
}
```

Throttle scales linearly with sin(γ) across this sweep, matching the first-principles relation T = D + W·sin(γ).

## File Layout

| File | Role |
|---|---|
| `AircraftModel.h` | Abstract base class |
| `Aircraft.h`, `NasaGTMAircraft.h`, `F16Aircraft.h` | Three concrete aircraft |
| `AircraftFactory.h` | Type enum and factory function |
| `EngineModel.h` | Engine abstraction with two concrete implementations |
| `AttitudeKinematics.h` | Quaternion / Euler strategy pattern |
| `LookupTable.h` | Templated N-dimensional interpolator |
| `GTMAeroTables.h` | NASA T2 polynomial aerodatabase coefficients |
| `AtmosphereModel.h` | US Standard Atmosphere 1976 |
| `SimulationEngine.h` | RK4 integration, linearization |
| `TrimSolver.h` | Damped Gauss-Newton trim |
| `TrimVerification.h` | Independent trim-condition checks |
| `ModeAnalysis.h` | Eigenvalue analysis |
| `LQRController.h` | CARE solver, optimal-gain computation |
| `LQRSimulation.h` | Open-/closed-loop nonlinear simulation |
| `SimConfig.h` | Configuration struct and scenario factories |
| `main.cpp` | Pipeline (`runSimulation`) and entry point |

## References

- *NASA Generic Transport Model T2 Polynomial Aerodatabase* (Cunis et al., v0.2, unlimited distribution).
- NASA Engineering and Safety Center 6-DOF Verification Check Cases.
- US Standard Atmosphere 1976, NASA-TM-X-74335.
- Stevens, Lewis, Johnson, *Aircraft Control and Simulation*, 3rd ed. — F-16 mass, geometry, and inertia.
- Nelson, *Flight Stability and Automatic Control* — published F-16 stability derivatives.
- Bryson, Ho, *Applied Optimal Control* — LQR formulation and weighting rule.
- Potter, J. E. (1966), "Matrix Quadratic Solutions," *SIAM Journal on Applied Mathematics* — Hamiltonian eigendecomposition for the Riccati equation.
- Department of Defense World Geodetic System 1984 (NIMA TR8350.2) — WGS-84 ellipsoid and normal-gravity constants.

## Future Work

- Faithful F-16 thrust model with Mach/altitude lookup and first-order spool-up dynamics
- Wind and turbulence disturbance models (Dryden, von Kármán) for control-robustness testing
- LQR lateral-directional autopilot (currently only longitudinal is implemented)
- Model Predictive Control extension with input/state constraint handling
- Adaptive control via recursive parameter estimation

---
