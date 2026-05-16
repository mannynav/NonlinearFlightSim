#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>


// ============================================================
//  LookupTable<N>
//
//  N-dimensional multilinear interpolation on a regular grid.
//  Each axis is an arbitrary monotonically-increasing breakpoint
//  vector; data is stored row-major (axis 0 outermost, axis N-1
//  innermost — i.e. C-style).
//
//  Out-of-range queries clamp to the nearest edge (standard
//  aerospace convention). Single-point axes are NOT supported;
//  every axis must have at least 2 breakpoints.
//
//  Lookup cost: O(N · log m + 2^N) where m is the average axis
//  size. With N=2 or N=3, that's effectively constant time and
//  the inner loop is small enough for the compiler to unroll.
//
//  Usage:
//      LookupTable<2> table({alpha_axis, beta_axis}, flat_data);
//      double CL = table(alpha, beta);     // variadic call
//      double CL = table.lookup({a, b});   // array call
// ============================================================
template <int N>
class LookupTable
{
    static_assert(N >= 1, "LookupTable requires at least 1 dimension");

public:
    using Axes = std::array<Eigen::VectorXd, N>;

    LookupTable() = default;

    LookupTable(const Axes& axes, const Eigen::VectorXd& data)
        : axes_(axes)
    {
        // Validate each axis
        int expected = 1;
        for (int a = 0; a < N; ++a) {
            if (axes_[a].size() < 2) {
                throw std::invalid_argument(
                    "LookupTable axis " + std::to_string(a) +
                    " must have at least 2 breakpoints");
            }
            for (int i = 1; i < axes_[a].size(); ++i) {
                if (axes_[a][i] <= axes_[a][i - 1]) {
                    throw std::invalid_argument(
                        "LookupTable axis " + std::to_string(a) +
                        " must be strictly increasing (at i=" +
                        std::to_string(i) + ")");
                }
            }
            sizes_[a] = static_cast<int>(axes_[a].size());
            expected *= sizes_[a];
        }

        if (data.size() != expected) {
            throw std::invalid_argument(
                "LookupTable data size mismatch: expected " +
                std::to_string(expected) + " got " +
                std::to_string(data.size()));
        }
        data_ = data;

        // Row-major strides: stride[N-1] = 1, stride[a] = stride[a+1] * size[a+1]
        strides_[N - 1] = 1;
        for (int a = N - 2; a >= 0; --a) {
            strides_[a] = strides_[a + 1] * sizes_[a + 1];
        }
    }


    // Variadic convenience overload
    template <typename... Args>
    double operator()(Args... xs) const
    {
        static_assert(sizeof...(Args) == N,
            "LookupTable: number of args must match dimension N");
        return lookup({ static_cast<double>(xs)... });
    }


    // Core lookup
    double lookup(const std::array<double, N>& x) const
    {
        std::array<int, N> base{};
        std::array<double, N> frac{};

        // For each axis: find bracketing pair, compute fractional position
        for (int a = 0; a < N; ++a) {
            const auto& ax = axes_[a];
            const int   n = sizes_[a];

            // upper_bound returns iterator to first element > x
            auto it = std::upper_bound(ax.data(), ax.data() + n, x[a]);
            int i = static_cast<int>(std::distance(ax.data(), it)) - 1;
            i = std::clamp(i, 0, n - 2);

            double f = (x[a] - ax[i]) / (ax[i + 1] - ax[i]);
            f = std::clamp(f, 0.0, 1.0);   // clamps out-of-range queries

            base[a] = i;
            frac[a] = f;
        }

        // Sum the 2^N corner contributions
        constexpr int n_corners = 1 << N;
        double result = 0.0;

        for (int c = 0; c < n_corners; ++c) {
            int    idx = 0;
            double weight = 1.0;
            for (int a = 0; a < N; ++a) {
                const bool upper = ((c >> a) & 1) != 0;
                idx += (base[a] + (upper ? 1 : 0)) * strides_[a];
                weight *= upper ? frac[a] : (1.0 - frac[a]);
            }
            result += weight * data_[idx];
        }

        return result;
    }


    // Accessors
    int   n_points()           const { return static_cast<int>(data_.size()); }
    int   axis_size(int a)     const { return sizes_[a]; }
    const Axes& axes() const { return axes_; }
    const Eigen::VectorXd& data() const { return data_; }


private:
    Axes              axes_;
    Eigen::VectorXd   data_;
    std::array<int, N> sizes_{};
    std::array<int, N> strides_{};
};