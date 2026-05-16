
#include <Eigen/Dense>
#include <iostream>
#include <memory>

// Aircraft models + attitude representation
#include "AircraftModel.h"
#include "Aircraft.h"
#include "NasaGTMAircraft.h"
#include "AttitudeKinematics.h"

// Simulation engine (split into modules)
#include "Constants.h"
#include "SimulationEngine.h"
#include "TrimVerification.h"
#include "PostProcessing.h"
#include "ModeAnalysis.h"
#include "LinearSimulation.h"

// Validation / lookup-table tests
#include "CheckCases.h"
#include "CheckCaseValidator.h"
#include "LookupTableTest.h"
#include "GTMAeroTables.h"

// Project utilities
#include "Output.h"
#include "FlightModes.h"
#include "GravityModel.h"
#include "Input.h"
#include "Utility.h"


// ============================================================
//  Initial conditions for the nonlinear sim
//  (10,000 ft, 85 m/s, straight & level, civil-aircraft trim)
// ============================================================
inline void initialStatesControls(Eigen::VectorXd& states, Eigen::VectorXd& controls)
{
    states[0] = 84.9993;
    states[1] = 0.0;
    states[2] = 1.2713;
    states[3] = 0.0;
    states[4] = 0.0;
    states[5] = 0.0;
    states[6] = 0.0;
    states[7] = 0.0150;
    states[8] = 0.0;
    states[9] = 0.0;
    states[10] = 0.0;
    states[11] = -3048.0;

    controls[0] = 0.0;
    controls[1] = -0.1780;
    controls[2] = 0.0;
    controls[3] = 0.0821;
    controls[4] = 0.0821;
}


// ============================================================
//  GTM aero lookup-table sanity check
//  Prints a handful of known points so you can spot-check the
//  tables before the full sim runs. Comment out the call in
//  main() once you trust the tables.
// ============================================================
inline void runGTMAeroSanityCheck()
{
    using namespace constants;
    auto cm = gtm_aero::build_cm_table();
    auto cl_roll = gtm_aero::build_cl_roll_table();
    auto cn = gtm_aero::build_cn_table();

    std::cout << "GTM aero sanity check:\n";
    std::cout << "  Cm(0 deg, 0 deg)    = " << cm(0.0, 0.0) << "   (expect -0.467)\n";
    std::cout << "  Cm(0 deg, -10 deg)  = " << cm(0.0, -10.0 * DEG_TO_RAD) << "   (expect ~0.041)\n";
    std::cout << "  Cl(10 deg, 0 deg)   = " << cl_roll(10.0 * DEG_TO_RAD, 0.0) << "   (expect ~-0.262)\n";
    std::cout << "  Cn(10 deg, 0 deg)   = " << cn(10.0 * DEG_TO_RAD, 0.0) << "   (expect ~0.175)\n\n";
}


// ============================================================
//  main
//
//  Pipeline:
//    1. (optional) GTM aero sanity check
//    2. Configure aircraft type and attitude representation
//    3. Build IC and run nonlinear RK4 simulation
//    4. Post-process to derived quantities
//    5. Load MATLAB trim point and verify
//    6. Linearize, decouple long/lat, run modal analysis
//    7. Drive linear longitudinal and lateral subsystems
//    8. Run NESC 6-DOF verification check cases
// ============================================================
int main()
{
    // ----- Optional pre-flight sanity checks -----
    // lookup_test::run_all();
    runGTMAeroSanityCheck();

    // ----- CONFIG -----
    AircraftType ac_type = AircraftType::CivilAircraft;        // or AircraftType::CivilAircraft
    AttitudeMode att_mode = AttitudeMode::Quaternion;     // or AttitudeMode::Euler

    std::unique_ptr<AttitudeKinematics> kin = makeAttitudeKinematics(att_mode);
    std::cout << "Attitude mode: " << kin->name()
        << "  (n_states = " << kin->n_states() << ")\n";

    // ----- Initial conditions -----
    Eigen::VectorXd initialX_euler(12);
    Eigen::VectorXd initialU(constants::NUM_CONTROLS);
    initialStatesControls(initialX_euler, initialU);
    Eigen::VectorXd initialX = buildFullState(initialX_euler, *kin);

    std::unique_ptr<AircraftModel> ac = makeAircraft(ac_type, initialX, initialU);
    std::cout << "Aircraft model: " << ac->name() << "\n";

    GravityModel gravity_model = standard_gravity_model;

    // ----- Nonlinear simulation -----
    constexpr double NL_SIM_LENGTH = 150.0;
    constexpr int    NL_STEPS = 200;

    Eigen::MatrixXd NL_sol = runNonlinearSimulation(*ac, *kin, initialX, initialU,
        gravity_model, NL_SIM_LENGTH, NL_STEPS);
    outputToFileWithHeaders(NL_sol, "Non Linear Solution.csv");
    std::cout << "Nonlinear simulation complete.\n";

    Eigen::MatrixXd PP = postProcessNonlinearSolution(NL_sol, *kin, NL_SIM_LENGTH, NL_STEPS);
    writePostProcessCsv(PP, "postProcessSolution.csv");

    // ----- Load trim point from MATLAB CSV -----
    std::cout << "\nLoading trim point from CSV...\n";
    Eigen::VectorXd trim_states_euler = readVectorFromCsv("trim_states.csv");
    Eigen::VectorXd trim_controls = readVectorFromCsv("trim_controls.csv");

    if (trim_states_euler.size() != 12 || trim_controls.size() != constants::NUM_CONTROLS) {
        std::cerr << "ERROR: Unexpected trim CSV sizes ("
            << trim_states_euler.size() << ", " << trim_controls.size() << ").\n";
        return 1;
    }

    Eigen::VectorXd trim_states = buildFullState(trim_states_euler, *kin);
    verifyStraightLevelTrimConditions(*ac, *kin, trim_states, trim_controls, gravity_model);

    // ----- Linearization about trim -----
    const int N = kin->n_states();
    Eigen::VectorXd Xdoto = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd Xo = trim_states;
    Eigen::VectorXd Uo = trim_controls;

    Eigen::MatrixXd dx(N, N);                       dx.setConstant(1e-5);
    Eigen::MatrixXd du(N, constants::NUM_CONTROLS); du.setConstant(1e-5);

    std::unique_ptr<AircraftModel> acl = makeAircraft(ac_type, Xo, Uo);

    std::cout << "Computing A and B (central differences)...\n";
    Eigen::MatrixXd A_full = LinearizeSystem_A(*acl, *kin, Xdoto, Xo, Uo, dx, gravity_model);
    Eigen::MatrixXd B_full = LinearizeSystem_B(*acl, *kin, Xdoto, Xo, Uo, du, gravity_model);

    outputToFileWithHeaders(A_full, "A_full.csv");
    outputToFileWithHeaders(B_full, "B_full.csv");

    // ----- Long/lat decoupling via similarity transform -----
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

    // ----- Longitudinal linear simulation (stabilizer doublet) -----
    constexpr double LONG_SIM_TIME = 150.0;
    constexpr int    LONG_STEPS = 150;

    Eigen::VectorXd Xlong_euler = short_period_2();
    Eigen::VectorXd Xlong_o(4);
    Xlong_o << Xlong_euler[0], Xlong_euler[1], Xlong_euler[2],
        kin->euler_dtheta_to_attitude(Xlong_euler[3]);

    Eigen::MatrixXd long_ctrl;
    Eigen::MatrixXd sol_long = runLinearSubsystemSim(
        *acl, A_long, B_long, Xlong_o, LONG_SIM_TIME, LONG_STEPS, long_ctrl,
        /*stab*/ 10, -10, /*ail*/ 0, 0, /*thrust*/ 0.5, 60.0);

    outputToFileWithHeaders(sol_long, "Linear Longitudinal Solution.csv");
    outputToFileWithHeaders(long_ctrl, "Linear Longitudinal Control Inputs.csv");
    std::cout << "Longitudinal linear simulation complete.\n";

    // ----- Lateral linear simulation -----
    constexpr double LAT_SIM_TIME = 200.0;
    constexpr int    LAT_STEPS = 250;

    Eigen::VectorXd Xlat_o(4);
    Xlat_o << 0.5, 0.01, 0.01, kin->euler_dphi_to_attitude(0.01);

    Eigen::MatrixXd lat_ctrl;
    Eigen::MatrixXd sol_lat = runLinearSubsystemSim(
        *acl, A_lat, B_lat, Xlat_o, LAT_SIM_TIME, LAT_STEPS, lat_ctrl,
        /*stab*/ 10, -10, /*ail*/ 0, 0, /*thrust*/ 0.0, 60.0);

    outputToFileWithHeaders(sol_lat, "Linear Lateral Solution.csv");
    outputToFileWithHeaders(lat_ctrl, "Linear Lateral Control Inputs.csv");
    std::cout << "Lateral linear simulation complete.\n";

    // ----- NESC 6-DOF Verification Check Cases -----
    runAllNESCCheckCases(*kin);

    std::cout << "\nAll outputs written. Done.\n";
    return 0;
}