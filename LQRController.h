
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <stdexcept>


// ============================================================
//  LQRController
//
//  Solves the continuous-time algebraic Riccati equation
//      A^T P + P A  -  P B R^{-1} B^T P  +  Q  =  0
//  and computes the optimal LQR feedback gain
//      K = R^{-1} B^T P
//  for the infinite-horizon problem
//      J = integral_0^inf ( x^T Q x + u^T R u ) dt
//      u* = -K x
//
//  Method: Hamiltonian eigendecomposition (Potter, 1966).
//  The 2n x 2n Hamiltonian matrix
//      H = [ A          -B R^{-1} B^T ]
//          [ -Q                -A^T  ]
//  has eigenvalues in stable/unstable pairs (lambda, -lambda).
//  Stacking the n eigenvectors corresponding to the stable
//  half-plane into V = [V11; V21], the Riccati solution is
//      P = V21 * V11^{-1}.
// 
// ============================================================

class LQRController
{
public:

    // Solve CARE; returns P (symmetric positive semidefinite).
    static Eigen::MatrixXd solveCARE(const Eigen::MatrixXd& A,
        const Eigen::MatrixXd& B,
        const Eigen::MatrixXd& Q,
        const Eigen::MatrixXd& R)
    {
        const int n = static_cast<int>(A.rows());
        const Eigen::MatrixXd Rinv = R.inverse();

        // ---- Build the 2n x 2n Hamiltonian matrix ----
        Eigen::MatrixXd H(2 * n, 2 * n);
        H.topLeftCorner(n, n) = A;
        H.topRightCorner(n, n) = -B * Rinv * B.transpose();
        H.bottomLeftCorner(n, n) = -Q;
        H.bottomRightCorner(n, n) = -A.transpose();

        // ---- Eigendecomposition ----
        Eigen::EigenSolver<Eigen::MatrixXd> es(H);
        if (es.info() != Eigen::Success)
            throw std::runtime_error("CARE: Hamiltonian eigendecomposition failed");

        const Eigen::VectorXcd eigvals = es.eigenvalues();
        const Eigen::MatrixXcd eigvecs = es.eigenvectors();

        // ---- Select the n stable eigenvectors (Re(lambda) < 0) ----
        Eigen::MatrixXcd V(2 * n, n);
        int col = 0;
        for (int i = 0; i < 2 * n && col < n; ++i) {
            if (eigvals(i).real() < 0.0) {
                V.col(col++) = eigvecs.col(i);
            }
        }
        if (col != n)
            throw std::runtime_error("CARE: failed to find n stable eigenvalues "
                "(system may not be stabilizable)");

        // ---- P = V21 * V11^{-1}, taking the real part ----
        const Eigen::MatrixXcd V11 = V.topRows(n);
        const Eigen::MatrixXcd V21 = V.bottomRows(n);
        const Eigen::MatrixXcd Pc = V21 * V11.inverse();

        // Symmetrize to clean up tiny imaginary/asymmetric residuals
        Eigen::MatrixXd P = Pc.real();
        return 0.5 * (P + P.transpose());
    }


    // Compute K = R^{-1} B^T P directly.
    static Eigen::MatrixXd lqrGain(const Eigen::MatrixXd& A,
        const Eigen::MatrixXd& B,
        const Eigen::MatrixXd& Q,
        const Eigen::MatrixXd& R)
    {
        const Eigen::MatrixXd P = solveCARE(A, B, Q, R);
        return R.inverse() * B.transpose() * P;
    }
};