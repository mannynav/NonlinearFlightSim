#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <stdexcept>

// Aircraft models + attitude representation
#include "AircraftModel.h"
#include "Aircraft.h"
#include "NasaGTMAircraft.h"
#include "AttitudeKinematics.h"

// Simulation engine (split into modules)
#include "Constants.h"
#include "SimConfig.h"
#include "SimulationEngine.h"
#include "TrimRoutine.h"
#include "TrimVerification.h"
#include "PostProcessing.h"
#include "ModeAnalysis.h"
#include "LinearSimulation.h"

// Environment models
#include "AtmosphereModel.h"

// Validation / lookup-table tests
#include "CheckCases.h"
#include "CheckCaseValidator.h"
#include "LookupTableTest.h"
#include "GTMAeroTables.h"

// Project utilities
#include "Output.h"
#include "FlightModes.h"
#include "Input.h"
#include "Utility.h"

// LQR files
#include "LQRSimulation.h"
//#include "LQRController.h"


inline void runGTMAeroSanityCheck()
{
    using namespace constants;
    auto cm = gtm_aero::build_cm_table();
    auto cl_roll = gtm_aero::build_cl_roll_table();
    auto cn = gtm_aero::build_cn_table();
    std::cout << "GTM aero sanity check:\n";
    std::cout << "  Cm(0 deg, 0 deg)    = " << cm(0.0, 0.0) << "\n";
    std::cout << "  Cm(0 deg, -10 deg)  = " << cm(0.0, -10.0 * DEG_TO_RAD) << "\n";
    std::cout << "  Cl(10 deg, 0 deg)   = " << cl_roll(10.0 * DEG_TO_RAD, 0.0) << "\n";
    std::cout << "  Cn(10 deg, 0 deg)   = " << cn(10.0 * DEG_TO_RAD, 0.0) << "\n\n";
}


// ============================================================
//  runSimulation
//
//  The full pipeline for one configuration:
//    1. Build atmosphere, kinematics, aircraft
//    2. Solve trim
//    3. Verify trim
//    4. Nonlinear sim from trim + perturbation
//    5. Linearize about trim, decouple long/lat
//    6. Modal analysis
//    7. Linear longitudinal & lateral sims
// ============================================================

void runSimulation(const SimConfig& cfg, const AtmosphereData& atm)
{
    std::cout << "\n############################################################\n";
    std::cout
        << "  |  Va = " << cfg.trim_Va << " m/s"
        << "  |  h = " << cfg.trim_h << " m"
        << "  |  gamma = " << cfg.trim_gamma * constants::RAD_TO_DEG << " deg\n";
    std::cout << "############################################################\n";

    std::unique_ptr<AttitudeKinematics> kin = makeAttitudeKinematics(cfg.att_mode);
    std::cout << "Attitude mode: " << kin->name()
        << "  (n_states = " << kin->n_states() << ")\n";

    // ----- Build aircraft with a starting guess (trim will replace it) -----
    Eigen::VectorXd Z0 = defaultTrimInitialGuess(cfg.trim_Va, cfg.trim_h, cfg.trim_gamma);
    Eigen::VectorXd guessX_euler = Z0.head(12);
    Eigen::VectorXd guessU = Z0.tail(constants::NUM_CONTROLS);
    Eigen::VectorXd guessX = buildFullState(guessX_euler, *kin);

    std::unique_ptr<AircraftModel> ac = makeAircraft(cfg.ac_type, guessX, guessU);
    std::cout << "Aircraft model: " << ac->name() << "\n";
    std::cout << "Atmosphere @ " << cfg.trim_h << " m: "
        << "rho = " << atm.airDensity(cfg.trim_h) << " kg/m^3, "
        << "g = " << atm.gravity(cfg.trim_h) << " m/s^2\n";

    // ----- Trim -----
    TrimResult trim = solveTrim(*ac, atm, cfg.trim_Va, cfg.trim_h, cfg.trim_gamma,
        Eigen::VectorXd{}, Eigen::VectorXd{},
        1e-12, 100, cfg.verbose_trim);
    if (!trim.converged) {
        std::cerr << "ERROR: trim did not converge. Skipping this run.\n";
        return;
    }

    Eigen::VectorXd trim_states_euler = trim.states_euler;
    Eigen::VectorXd trim_controls = trim.controls;
    Eigen::VectorXd trim_states = buildFullState(trim_states_euler, *kin);

    verifyStraightLevelTrimConditions(*ac, *kin, trim_states, trim_controls, atm,
        cfg.trim_Va, cfg.trim_h, cfg.trim_gamma);

    // ----- Nonlinear sim from trim + pitch perturbation -----
    Eigen::VectorXd ic_euler = trim_states_euler;
    ic_euler[7] += cfg.pitch_perturbation;
    Eigen::VectorXd ic = buildFullState(ic_euler, *kin);

    try {
        Eigen::MatrixXd NL_sol = runNonlinearSimulation(
            *ac, *kin, ic, trim_controls, atm, cfg.nl_length, cfg.nl_steps);
        outputToFileWithHeaders(NL_sol, "Non Linear Solution.csv");
        std::cout << "Nonlinear sim complete (" << cfg.nl_length << " s, "
            << cfg.nl_steps << " steps, dt = " << cfg.nl_dt() << " s).\n";

        Eigen::MatrixXd PP = postProcessNonlinearSolution(
            NL_sol, *kin, atm, cfg.nl_length, cfg.nl_steps);
        writePostProcessCsv(PP, "postProcessSolution.csv");
    }
    catch (const std::exception& e) {
        std::cerr << "WARNING: nonlinear sim failed: " << e.what() << "\n";
    }

    // ----- Linearization about trim -----
    const int N = kin->n_states();
    Eigen::VectorXd Xdoto = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd Xo = trim_states;
    Eigen::VectorXd Uo = trim_controls;

    Eigen::MatrixXd dx(N, N);                       dx.setConstant(1e-5);
    Eigen::MatrixXd du(N, constants::NUM_CONTROLS); du.setConstant(1e-5);

    std::unique_ptr<AircraftModel> acl = makeAircraft(cfg.ac_type, Xo, Uo);

    std::cout << "Computing A and B (central differences)...\n";
    Eigen::MatrixXd A_full = LinearizeSystem_A(*acl, *kin, Xdoto, Xo, Uo, dx, atm);
    Eigen::MatrixXd B_full = LinearizeSystem_B(*acl, *kin, Xdoto, Xo, Uo, du, atm);

    outputToFileWithHeaders(A_full, "A_full.csv");
    outputToFileWithHeaders(B_full, "B_full.csv");

    Eigen::MatrixXd S = kin->longLatPermutation();
    Eigen::MatrixXd A_tf = S * A_full * S.transpose();
    Eigen::MatrixXd B_tf = S * B_full;

    Eigen::MatrixXd A_long = A_tf.block(0, 0, 4, 4);
    Eigen::MatrixXd B_long = B_tf.block(0, 0, 4, constants::NUM_CONTROLS);
    Eigen::MatrixXd A_lat = A_tf.block(4, 4, 4, 4);
    Eigen::MatrixXd B_lat = B_tf.block(4, 0, 4, constants::NUM_CONTROLS);

    outputToFileWithHeaders(A_long, "A_long.csv");
    outputToFileWithHeaders(B_long, "B_long.csv");
    outputToFileWithHeaders(A_lat, "A_lat.csv");
    outputToFileWithHeaders(B_lat, "B_lat.csv");

    printModeAnalysis("Full system", A_tf);
    printModeAnalysis("Longitudinal", A_long);
    printModeAnalysis("Lateral", A_lat);

    // ----- Linear longitudinal sim -----
    try {
        Eigen::VectorXd Xlong_euler = short_period_2();
        Eigen::VectorXd Xlong_o(4);
        Xlong_o << Xlong_euler[0], Xlong_euler[1], Xlong_euler[2],
            kin->euler_dtheta_to_attitude(Xlong_euler[3]);

        Eigen::MatrixXd long_ctrl;
        Eigen::MatrixXd sol_long = runLinearSubsystemSim(
            *acl, A_long, B_long, Xlong_o,
            cfg.long_length, cfg.long_steps, long_ctrl,
            /*stab*/ 10, -10, /*ail*/ 0, 0, /*thrust*/ 0.5, 60.0);

        outputToFileWithHeaders(sol_long, "Linear Longitudinal Solution.csv");
        outputToFileWithHeaders(long_ctrl, "Linear Longitudinal Control Inputs.csv");
        std::cout << "Longitudinal linear sim complete.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "WARNING: longitudinal linear sim failed: " << e.what() << "\n";
    }

    // ----- Linear lateral sim -----
    try {
        Eigen::VectorXd Xlat_o(4);
        Xlat_o << 0.5, 0.01, 0.01, kin->euler_dphi_to_attitude(0.01);

        Eigen::MatrixXd lat_ctrl;
        Eigen::MatrixXd sol_lat = runLinearSubsystemSim(
            *acl, A_lat, B_lat, Xlat_o,
            cfg.lat_length, cfg.lat_steps, lat_ctrl,
            /*stab*/ 10, -10, /*ail*/ 0, 0, /*thrust*/ 0.0, 60.0);

        outputToFileWithHeaders(sol_lat, "Linear Lateral Solution.csv");
        outputToFileWithHeaders(lat_ctrl, "Linear Lateral Control Inputs.csv");
        std::cout << "Lateral linear sim complete.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "WARNING: lateral linear sim failed: " << e.what() << "\n";
    }


    if (cfg.run_lqr) {
        runLQRDemo(*acl, *kin,
            trim_states, trim_states_euler, trim_controls,
            A_long, B_long, S,
            atm,
            /*sim_length=*/10.0,
            /*sim_steps=*/1000,
            /*pitch_perturbation=*/0.05);
    }

}


// ============================================================
//  main — pick configuration(s) and run
// ============================================================
int main()
{
    runGTMAeroSanityCheck();
    AtmosphereData atm = loadUSStandardAtmosphere();

    // ---- Single run: pick one ----
    //SimConfig cfg = SimConfig::gtm_cruise();
    //SimConfig cfg = SimConfig::civil_cruise();
    //SimConfig cfg = SimConfig::climb(AircraftType::NasaGTM, 6.0);
    SimConfig cfg = SimConfig::f16_lqr_demo();
    //SimConfig cfg = SimConfig::f16_cruise();
    runSimulation(cfg, atm);

    // ---- Or sweep several operating points (uncomment to use) ----
    // for (double gam : {0.0, 3.0, 6.0, 10.0, 15.0}) {
    //     runSimulation(SimConfig::climb(AircraftType::NasaGTM, gam), atm);
    // }

    // ----- NESC verification (representation-only, run once) -----
    auto kin = makeAttitudeKinematics(AttitudeMode::Quaternion);
    runAllNESCCheckCases(*kin);

    std::cout << "\nAll outputs written. Done.\n";
    return 0;
}