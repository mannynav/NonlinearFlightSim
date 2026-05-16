
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>


// ============================================================
//  Mode analysis
//
//  Eigendecomposition of a system matrix A, with each eigenvalue
//  reported alongside natural frequency wn, damping ratio zeta,
//  time constant tau, and a stable/unstable/neutral classifier.
//
//  Works on any A: full 12x12 (Euler) / 13x13 (Quaternion) or the
//  4x4 longitudinal / lateral blocks after the long/lat similarity
//  transform.
// ============================================================
inline void printModeAnalysis(const std::string& label, const Eigen::MatrixXd& A_block)
{
    Eigen::EigenSolver<Eigen::MatrixXd> solver(A_block);
    auto evals = solver.eigenvalues();
    constexpr double NEUTRAL_TOL = 1e-6;

    std::cout << "\n--- " << label << " modes ---\n";
    std::cout << std::left
        << std::setw(6) << "#"
        << std::setw(28) << "Eigenvalue"
        << std::setw(14) << "wn (rad/s)"
        << std::setw(14) << "zeta"
        << std::setw(14) << "tau (s)"
        << "Status\n";
    std::cout << std::string(80, '-') << "\n";

    for (int i = 0; i < evals.size(); ++i) {
        double re = evals[i].real();
        double im = evals[i].imag();
        double wn = std::sqrt(re * re + im * im);
        double zeta = (wn > 1e-10) ? -re / wn : 0.0;
        std::string tau_str = (std::abs(re) > 1e-10) ? std::to_string(-1.0 / re) : "inf";
        std::string status = (re < -NEUTRAL_TOL) ? "stable"
            : (re > NEUTRAL_TOL) ? "UNSTABLE" : "neutral";

        std::cout << std::fixed << std::setprecision(6)
            << std::setw(6) << i
            << std::setw(28) << (std::to_string(re) + " + " + std::to_string(im) + "j")
            << std::setw(14) << wn
            << std::setw(14) << zeta
            << std::setw(14) << tau_str
            << status << "\n";
    }
}