
#pragma once

#include "LookupTable.h"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>


// ============================================================
//  Validation driver for LookupTable<N>
//
//  Builds synthetic tables from known analytical functions, then
//  exercises the interpolator at:
//    - grid points (should be EXACT, error < 1e-12)
//    - midpoints  (linear interp error ~ O(h^2 · max|f''|))
//    - out-of-range queries (should clamp to nearest edge)
//
//  Pass/fail thresholds are calibrated to the analytical second-
//  derivative of each test function and the chosen grid spacing.
// ============================================================
namespace lookup_test {

    inline std::string pass(bool ok) { return ok ? " [PASS]" : " [FAIL]"; }


    // ---- 1D test: f(x) = sin(x) on [0, 2π] ----
    inline bool test_1d()
    {
        const int N = 21;
        const double xmax = 2.0 * std::numbers::pi;
        Eigen::VectorXd x(N), y(N);
        for (int i = 0; i < N; ++i) {
            x[i] = i * xmax / (N - 1);
            y[i] = std::sin(x[i]);
        }
        LookupTable<1> table({ x }, y);

        // Test 1: queries exactly at grid points
        double err_grid = 0.0;
        for (int i = 0; i < N; ++i) {
            err_grid = std::max(err_grid, std::abs(table(x[i]) - y[i]));
        }
        const bool ok_grid = err_grid < 1e-12;

        // Test 2: queries at midpoints — linear interp error should be < h²·max|f''|/8
        // h = 2π/20 ≈ 0.314, max|f''| = 1, expected error < 0.012
        double err_mid = 0.0;
        for (int i = 0; i < N - 1; ++i) {
            double xq = 0.5 * (x[i] + x[i + 1]);
            err_mid = std::max(err_mid, std::abs(table(xq) - std::sin(xq)));
        }
        const bool ok_mid = err_mid < 0.02;

        // Test 3: clamping
        const double left = table(-100.0);
        const double right = table(+100.0);
        const bool ok_clamp = (std::abs(left - y[0]) < 1e-12)
            && (std::abs(right - y[N - 1]) < 1e-12);

        std::cout << "  1D f(x)=sin(x):\n";
        std::cout << "    at grid points: max err = " << err_grid << pass(ok_grid) << "\n";
        std::cout << "    midpoints:      max err = " << err_mid << pass(ok_mid) << "\n";
        std::cout << "    clamping:       L=" << left << " R=" << right << pass(ok_clamp) << "\n";

        return ok_grid && ok_mid && ok_clamp;
    }


    // ---- 2D test: f(x, y) = sin(x) · cos(y) on [0,3] x [0,3] ----
    inline bool test_2d()
    {
        const int Nx = 11, Ny = 11;
        Eigen::VectorXd ax(Nx), ay(Ny);
        for (int i = 0; i < Nx; ++i) ax[i] = i * 3.0 / (Nx - 1);
        for (int j = 0; j < Ny; ++j) ay[j] = j * 3.0 / (Ny - 1);

        // Row-major: data[i*Ny + j]
        Eigen::VectorXd data(Nx * Ny);
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                data[i * Ny + j] = std::sin(ax[i]) * std::cos(ay[j]);
            }
        }
        LookupTable<2> table({ ax, ay }, data);

        double err_grid = 0.0, err_mid = 0.0;

        // grid points
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                double got = table(ax[i], ay[j]);
                double want = std::sin(ax[i]) * std::cos(ay[j]);
                err_grid = std::max(err_grid, std::abs(got - want));
            }
        }

        // midpoints
        for (int i = 0; i < Nx - 1; ++i) {
            for (int j = 0; j < Ny - 1; ++j) {
                double xq = 0.5 * (ax[i] + ax[i + 1]);
                double yq = 0.5 * (ay[j] + ay[j + 1]);
                double got = table(xq, yq);
                double want = std::sin(xq) * std::cos(yq);
                err_mid = std::max(err_mid, std::abs(got - want));
            }
        }

        // clamping (corner)
        double lo = table(-10.0, -10.0);   // should clamp to f(0,0) = 0
        double hi = table(+10.0, +10.0);   // should clamp to f(3,3)
        bool ok_clamp = (std::abs(lo - 0.0) < 1e-12)
            && (std::abs(hi - std::sin(3.0) * std::cos(3.0)) < 1e-12);

        const bool ok_grid = err_grid < 1e-12;
        const bool ok_mid = err_mid < 0.05;

        std::cout << "  2D f(x,y)=sin(x)cos(y):\n";
        std::cout << "    at grid points: max err = " << err_grid << pass(ok_grid) << "\n";
        std::cout << "    midpoints:      max err = " << err_mid << pass(ok_mid) << "\n";
        std::cout << "    clamping:       lo=" << lo << " hi=" << hi << pass(ok_clamp) << "\n";

        return ok_grid && ok_mid && ok_clamp;
    }


    // ---- 3D test: f(x,y,z) = sin(x)·cos(y) + 0.5·z ----
    inline bool test_3d()
    {
        const int N = 7;
        Eigen::VectorXd a1(N), a2(N), a3(N);
        for (int i = 0; i < N; ++i) {
            a1[i] = i * 0.5;
            a2[i] = i * 0.5;
            a3[i] = i * 0.5;
        }

        // Row-major: data[i*N*N + j*N + k]
        Eigen::VectorXd data(N * N * N);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    data[i * N * N + j * N + k] =
                        std::sin(a1[i]) * std::cos(a2[j]) + 0.5 * a3[k];
                }
            }
        }
        LookupTable<3> table({ a1, a2, a3 }, data);

        double err_grid = 0.0, err_mid = 0.0;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    double got = table(a1[i], a2[j], a3[k]);
                    double want = std::sin(a1[i]) * std::cos(a2[j]) + 0.5 * a3[k];
                    err_grid = std::max(err_grid, std::abs(got - want));
                }
            }
        }

        for (int i = 0; i < N - 1; ++i) {
            for (int j = 0; j < N - 1; ++j) {
                for (int k = 0; k < N - 1; ++k) {
                    double xq = 0.5 * (a1[i] + a1[i + 1]);
                    double yq = 0.5 * (a2[j] + a2[j + 1]);
                    double zq = 0.5 * (a3[k] + a3[k + 1]);
                    double got = table(xq, yq, zq);
                    double want = std::sin(xq) * std::cos(yq) + 0.5 * zq;
                    err_mid = std::max(err_mid, std::abs(got - want));
                }
            }
        }

        const bool ok_grid = err_grid < 1e-12;
        const bool ok_mid = err_mid < 0.10;

        std::cout << "  3D f(x,y,z)=sin(x)cos(y)+0.5z:\n";
        std::cout << "    at grid points: max err = " << err_grid << pass(ok_grid) << "\n";
        std::cout << "    midpoints:      max err = " << err_mid << pass(ok_mid) << "\n";

        return ok_grid && ok_mid;
    }


    // ---- Edge case: constant function ----
    inline bool test_constant()
    {
        Eigen::VectorXd ax(5);  ax << -2, -1, 0, 1, 2;
        Eigen::VectorXd data(5); data.setConstant(3.14159);

        LookupTable<1> table({ ax }, data);

        // Any query should return the constant exactly
        bool ok = (std::abs(table(-100.0) - 3.14159) < 1e-12)
            && (std::abs(table(0.5) - 3.14159) < 1e-12)
            && (std::abs(table(+100.0) - 3.14159) < 1e-12);

        std::cout << "  Constant function:                          " << pass(ok) << "\n";
        return ok;
    }


    // ---- Top-level driver ----
    inline void run_all()
    {
        std::cout << "\n========== LookupTable validation ==========\n";
        std::cout << std::scientific << std::setprecision(3);

        int total = 0, passed = 0;
        auto record = [&](bool ok) { ++total; if (ok) ++passed; };

        record(test_1d());
        record(test_2d());
        record(test_3d());
        record(test_constant());

        std::cout << "\n  " << passed << "/" << total << " test groups passed.\n";
        std::cout << "============================================\n\n";
    }

} // namespace lookup_test