#pragma once

#include <Eigen/Dense>
#include <numbers>
#include <cmath>
#include <memory>


// ============================================================
//  Quaternion utilities (Hamilton, scalar-first: q = [q0, q1, q2, q3])
//
//  Convention: q rotates from NED to body, so trim => q = [1, 0, 0, 0]
//  Euler convention: 3-2-1 (yaw, pitch, roll)
// ============================================================

namespace quat {

    inline Eigen::Vector4d euler_to_quaternion(double phi, double theta, double psi)
    {
        double cp = std::cos(phi * 0.5), sp = std::sin(phi * 0.5);
        double ct = std::cos(theta * 0.5), st = std::sin(theta * 0.5);
        double cs = std::cos(psi * 0.5), ss = std::sin(psi * 0.5);

        Eigen::Vector4d q;
        q[0] = cp * ct * cs + sp * st * ss;   // q0
        q[1] = sp * ct * cs - cp * st * ss;   // q1
        q[2] = cp * st * cs + sp * ct * ss;   // q2
        q[3] = cp * ct * ss - sp * st * cs;   // q3
        return q;
    }

    inline Eigen::Vector3d quaternion_to_euler(double q0, double q1, double q2, double q3)
    {
        Eigen::Vector3d eul;
        eul[0] = std::atan2(2.0 * (q0 * q1 + q2 * q3), 1.0 - 2.0 * (q1 * q1 + q2 * q2));
        double s = 2.0 * (q0 * q2 - q3 * q1);
        s = std::max(-1.0, std::min(1.0, s));
        eul[1] = std::asin(s);
        eul[2] = std::atan2(2.0 * (q0 * q3 + q1 * q2), 1.0 - 2.0 * (q2 * q2 + q3 * q3));
        return eul;
    }

    inline Eigen::Vector4d normalize(const Eigen::Vector4d& q)
    {
        double n = q.norm();
        return (n > 1e-12) ? (q / n).eval() : Eigen::Vector4d(1, 0, 0, 0);
    }

    inline Eigen::Matrix3d body_to_ned(double q0, double q1, double q2, double q3)
    {
        Eigen::Matrix3d R;
        R(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
        R(0, 1) = 2.0 * (q1 * q2 - q0 * q3);
        R(0, 2) = 2.0 * (q1 * q3 + q0 * q2);
        R(1, 0) = 2.0 * (q1 * q2 + q0 * q3);
        R(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
        R(1, 2) = 2.0 * (q2 * q3 - q0 * q1);
        R(2, 0) = 2.0 * (q1 * q3 - q0 * q2);
        R(2, 1) = 2.0 * (q2 * q3 + q0 * q1);
        R(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
        return R;
    }

} // namespace quat


// ============================================================
//  AttitudeKinematics — strategy interface
//
//  Encapsulates everything that differs between Euler-angle
//  and quaternion attitude representations:
//    - state vector size (3 vs 4 attitude elements)
//    - gravity expressed in body frame
//    - body-to-NED rotation matrix
//    - attitude derivative given body angular velocity
//    - normalization (no-op for Euler, |q|=1 for quaternion)
//    - long/lat decoupling permutation
// ============================================================

class AttitudeKinematics
{
public:
    virtual ~AttitudeKinematics() = default;

    // Sizing — body (6) + attitude (n_att) + position (3)
    virtual int n_attitude() const = 0;
    int n_states() const { return 6 + n_attitude() + 3; }
    int attitude_start_idx() const { return 6; }
    int position_start_idx() const { return 6 + n_attitude(); }

    // (phi, theta, psi) -> attitude state segment
    virtual Eigen::VectorXd from_euler(double phi, double theta, double psi) const = 0;

    // attitude state segment -> (phi, theta, psi)
    virtual Eigen::Vector3d to_euler(const Eigen::VectorXd& att) const = 0;

    // Gravity in body frame
    virtual Eigen::Vector3d gravity_body(double g, const Eigen::VectorXd& att) const = 0;

    // Body-to-NED rotation matrix
    virtual Eigen::Matrix3d body_to_ned(const Eigen::VectorXd& att) const = 0;

    // Attitude derivative given body angular velocity ω = [p, q, r]
    virtual Eigen::VectorXd attitude_derivative(const Eigen::Vector3d& omega,
        const Eigen::VectorXd& att) const = 0;

    // Optional normalization (no-op for Euler; renormalizes to |q|=1 for quat)
    virtual Eigen::VectorXd normalize(const Eigen::VectorXd& att) const = 0;

    // Permutation S so that  S * x  has order [u, w, q, att_pitch | v, p, r, att_roll | rest...]
    virtual Eigen::MatrixXd longLatPermutation() const = 0;

    // Convert an Euler-angle perturbation to the corresponding attitude-state perturbation
    // Used to scale long/lat ICs uniformly across kinematics modes.
    virtual double euler_dtheta_to_attitude(double dtheta) const = 0;
    virtual double euler_dphi_to_attitude(double dphi)   const = 0;

    virtual const char* name() const = 0;
};


// ============================================================
//  EulerKinematics — 3-state attitude [phi, theta, psi]
// ============================================================
class EulerKinematics : public AttitudeKinematics
{
public:
    int n_attitude() const override { return 3; }

    Eigen::VectorXd from_euler(double phi, double theta, double psi) const override
    {
        Eigen::VectorXd a(3);
        a << phi, theta, psi;
        return a;
    }

    Eigen::Vector3d to_euler(const Eigen::VectorXd& att) const override
    {
        return Eigen::Vector3d(att[0], att[1], att[2]);
    }

    Eigen::Vector3d gravity_body(double g, const Eigen::VectorXd& att) const override
    {
        double phi = att[0], theta = att[1];
        return Eigen::Vector3d(
            -g * std::sin(theta),
            g * std::cos(theta) * std::sin(phi),
            g * std::cos(theta) * std::cos(phi));
    }

    Eigen::Matrix3d body_to_ned(const Eigen::VectorXd& att) const override
    {
        double phi = att[0], theta = att[1], psi = att[2];
        double cphi = std::cos(phi), sphi = std::sin(phi);
        double cth = std::cos(theta), sth = std::sin(theta);
        double cps = std::cos(psi), sps = std::sin(psi);

        Eigen::Matrix3d R;
        R(0, 0) = cth * cps;
        R(0, 1) = -cphi * sps + sphi * sth * cps;
        R(0, 2) = sphi * sps + cphi * sth * cps;
        R(1, 0) = cth * sps;
        R(1, 1) = cphi * cps + sphi * sth * sps;
        R(1, 2) = -sphi * cps + cphi * sth * sps;
        R(2, 0) = -sth;
        R(2, 1) = sphi * cth;
        R(2, 2) = cphi * cth;
        return R;
    }

    Eigen::VectorXd attitude_derivative(const Eigen::Vector3d& omega,
        const Eigen::VectorXd& att) const override
    {
        double phi = att[0], theta = att[1];
        double cphi = std::cos(phi), sphi = std::sin(phi);
        double cth = std::cos(theta), tth = std::tan(theta);

        // Standard 3-2-1 inverse kinematics/
        Eigen::Matrix3d H;
        H << 1.0, sphi* tth, cphi* tth,
            0.0, cphi, -sphi,
            0.0, sphi / cth, cphi / cth;
        return H * omega;
    }

    Eigen::VectorXd normalize(const Eigen::VectorXd& att) const override
    {
        return att;   // no-op
    }

    Eigen::MatrixXd longLatPermutation() const override
    {
        // 12x12 — original ordering [u,v,w,p,q,r,phi,theta,psi,x,y,z]
        // new ordering [u,w,q,theta | v,p,r,phi | psi,x,y,z]
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(12, 12);
        S(0, 0) = 1;  // u
        S(1, 2) = 1;  // w
        S(2, 4) = 1;  // q
        S(3, 7) = 1;  // theta
        S(4, 1) = 1;  // v
        S(5, 3) = 1;  // p
        S(6, 5) = 1;  // r
        S(7, 6) = 1;  // phi
        S(8, 8) = 1;  // psi
        S(9, 9) = 1;  // x_e
        S(10, 10) = 1;  // y_e
        S(11, 11) = 1;  // z_e
        return S;
    }

    double euler_dtheta_to_attitude(double dtheta) const override { return dtheta; }
    double euler_dphi_to_attitude(double dphi)   const override { return dphi; }

    const char* name() const override { return "Euler"; }
};


// ============================================================
//  QuaternionKinematics — 4-state attitude [q0, q1, q2, q3]
// ============================================================
class QuaternionKinematics : public AttitudeKinematics
{
public:
    int n_attitude() const override { return 4; }

    Eigen::VectorXd from_euler(double phi, double theta, double psi) const override
    {
        Eigen::Vector4d q = quat::euler_to_quaternion(phi, theta, psi);
        Eigen::VectorXd a(4);
        a << q[0], q[1], q[2], q[3];
        return a;
    }

    Eigen::Vector3d to_euler(const Eigen::VectorXd& att) const override
    {
        return quat::quaternion_to_euler(att[0], att[1], att[2], att[3]);
    }

    Eigen::Vector3d gravity_body(double g, const Eigen::VectorXd& att) const override
    {
        double q0 = att[0], q1 = att[1], q2 = att[2], q3 = att[3];
        return Eigen::Vector3d(
            g * 2.0 * (q1 * q3 - q0 * q2),
            g * 2.0 * (q2 * q3 + q0 * q1),
            g * (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3));
    }

    Eigen::Matrix3d body_to_ned(const Eigen::VectorXd& att) const override
    {
        return quat::body_to_ned(att[0], att[1], att[2], att[3]);
    }

    Eigen::VectorXd attitude_derivative(const Eigen::Vector3d& omega,
        const Eigen::VectorXd& att) const override
    {
        double q0 = att[0], q1 = att[1], q2 = att[2], q3 = att[3];
        double p = omega[0], q = omega[1], r = omega[2];
        Eigen::VectorXd qdot(4);
        qdot[0] = 0.5 * (-p * q1 - q * q2 - r * q3);
        qdot[1] = 0.5 * (p * q0 + r * q2 - q * q3);
        qdot[2] = 0.5 * (q * q0 - r * q1 + p * q3);
        qdot[3] = 0.5 * (r * q0 + q * q1 - p * q2);
        return qdot;
    }

    Eigen::VectorXd normalize(const Eigen::VectorXd& att) const override
    {
        Eigen::Vector4d q(att[0], att[1], att[2], att[3]);
        Eigen::Vector4d qn = quat::normalize(q);
        Eigen::VectorXd out(4);
        out << qn[0], qn[1], qn[2], qn[3];
        return out;
    }

    Eigen::MatrixXd longLatPermutation() const override
    {
        // 13x13 — original ordering [u,v,w,p,q,r,q0,q1,q2,q3,x,y,z]
        // new ordering [u,w,q,q2 | v,p,r,q1 | q0,q3,x,y,z]
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(13, 13);
        S(0, 0) = 1;  // u
        S(1, 2) = 1;  // w
        S(2, 4) = 1;  // q
        S(3, 8) = 1;  // q2  (≈ θ/2)
        S(4, 1) = 1;  // v
        S(5, 3) = 1;  // p
        S(6, 5) = 1;  // r
        S(7, 7) = 1;  // q1  (≈ φ/2)
        S(8, 6) = 1;  // q0
        S(9, 9) = 1;  // q3  (≈ ψ/2)
        S(10, 10) = 1;  // x_e
        S(11, 11) = 1;  // y_e
        S(12, 12) = 1;  // z_e
        return S;
    }

    // δq2 ≈ δθ/2  and  δq1 ≈ δφ/2 about q_o = [1,0,0,0]
    double euler_dtheta_to_attitude(double dtheta) const override { return 0.5 * dtheta; }
    double euler_dphi_to_attitude(double dphi)   const override { return 0.5 * dphi; }

    const char* name() const override { return "Quaternion"; }
};


// ============================================================
//  Factory helper
// ============================================================
enum class AttitudeMode { Euler, Quaternion };

inline std::unique_ptr<AttitudeKinematics> makeAttitudeKinematics(AttitudeMode mode)
{
    if (mode == AttitudeMode::Quaternion) return std::make_unique<QuaternionKinematics>();
    return std::make_unique<EulerKinematics>();
}
