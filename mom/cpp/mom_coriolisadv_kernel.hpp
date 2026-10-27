// SKILLS: 0.3.2
#pragma once
///
/// @file mom_coriolisadv_kernel.hpp
/// @brief Per-point device primitives for the MOM6 Coriolis and momentum advection kernels.
///
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_Math.H>
#include <AMReX_Algorithm.H>

namespace MOM {
using amrex::Real;

///
/// @brief Reconstruct the scalar (e.g., pv, vorticity) onto point i-1/2 using a third-order
///        upwind scheme
///
/// @param q4 Tracer values on points i-2, i-1, i, i+1 [A ~> a]
/// @param u  Velocity or thickness flux on point i-1/2
///           [l t-1 ~> m s-1] or [l2 t-1 ~> m2 s-1]
/// @param qr Reconstruction of tracer q at point i-1/2 [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void up3_reconstruction_point(Real const q4[4], Real const u, Real& qr) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_6 = 1.0_rt/6.0_rt;  // The ratio of 1/6 [nondim]

    if (u > 0.0_rt) {
        qr = ((2.0_rt*q4[2] + 5.0_rt*q4[1]) - q4[0]) * C1_6;
    } else {
        qr = ((2.0_rt*q4[1] + 5.0_rt*q4[2]) - q4[3]) * C1_6;
    }
}

///
/// @brief Reconstruct the scalar (e.g., PV, vorticity) onto point i-1/2
///        using a third-order upwind scheme with the Koren flux limiter
///
/// @param q4 Tracer values on points i-2, i-1, i, i+1 [A ~> a]
/// @param u  Velocity or thickness flux on point i-1/2
///           [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
/// @param qr Reconstruction of tracer q on point i-1/2 [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void up3_koren_limiter_reconstruction_point(Real const q4[4], Real const u, Real& qr) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_3 = 1.0_rt/3.0_rt;  // The ratio of 1/3 [nondim]
    constexpr Real C1_6 = 1.0_rt/6.0_rt;  // The ratio of 1/6 [nondim]
    Real theta;                     // Ratio of gradient [nondim]
    Real psi;                       // Limiter function [nondim]

    if (u > 0.0_rt) {
        if (q4[2] == q4[1]) {
            qr = q4[1];
        } else {
            theta = (q4[1] - q4[0])/(q4[2] - q4[1]);
            // limiter introduced by Koren (1993)
            psi = amrex::max(0.0_rt, amrex::min(1.0_rt, C1_3 + C1_6*theta, theta));
            qr = q4[1] + psi*(q4[2] - q4[1]);
        }
    } else {
        if (q4[2] == q4[1]) {
            qr = q4[2];
        } else {
            theta = (q4[3] - q4[2])/(q4[2] - q4[1]);
            psi = amrex::max(0.0_rt, amrex::min(1.0_rt, C1_3 + C1_6*theta, theta));
            qr = q4[2] + psi*(q4[1] - q4[2]);
        }
    }
}

///
/// @brief Compute the factor for the WENO weights
///
/// @param tau Difference of the smoothness indicator [A ~> a]
/// @param b   The smoothness indicator [A ~> a]
/// @return    The factor for the weight [nondim]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real fac_fn_point(Real const tau, Real const b) noexcept
{
    using namespace amrex::literals;

    Real fac = 1.0e40_rt;
    if (amrex::Math::abs(b) > 1.0e-20_rt*tau) {
        Real const r = 1.0_rt + tau / b;
        fac = r*r;
    }
    return fac;
}

///
/// @brief Compute the smoothness indicator for the two-point stencil of the third-order WENO scheme
///
/// @param q2 Tracer values on the two-point stencil [A ~> a]
/// @param w0 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_three_weight_point(Real const q2[2], Real& w0) noexcept
{
    Real const d = q2[0] - q2[1];
    w0 = d*d;
}

///
/// @brief Reconstruction in the second upwind stencil of the third-order WENO scheme
///
/// @param q2 Tracer values on the two-point stencil [A ~> a]
/// @param w0 Reconstruction of the quantity [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_three_reconstruction_0_point(Real const q2[2], Real& w0) noexcept
{
    using namespace amrex::literals;

    w0 = (q2[0] + q2[1]) * 0.5_rt;
}

///
/// @brief Reconstruction in the first upwind stencil for third-order WENO scheme
///
/// @param q2 Tracer values on the two-point stencil [A ~> a]
/// @param w0 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_three_reconstruction_1_point(Real const q2[2], Real& w0) noexcept
{
    using namespace amrex::literals;

    w0 = (- q2[0] + 3.0_rt * q2[1]) * 0.5_rt;
}

///
/// @brief Reconstruct the tracer (e.g., PV, vorticity) onto the point i-1/2 using a third-order
///        WENO scheme. This reconstruction is thickness-weighted
///
/// @param q4     Tracer value times thickness on points i-2, i-1, i, i+1 [A ~> a]
/// @param h4     Thickness values on points i-2, i-1, i, i+1 [L ~> m]
/// @param u4     Values on points i-2, i-1, i, i+1 used to compute the smoothness indicator:
///               velocity [L T-1 ~> m s-1] with velocity smoothing, otherwise q4 [A ~> a]
/// @param h_tiny A tiny thickness to prevent division by zero [L ~> m]
/// @param u      Velocity or thickness flux on point i-1/2
///               [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
/// @param qr     Reconstruction of tracer q on point i-1/2 [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_three_h_weight_reconstruction_point(Real const q4[4], Real const h4[4], Real const u4[4],
                                              Real const h_tiny, Real const u, Real& qr) noexcept
{
    using namespace amrex::literals;

    constexpr Real C2_3 = 2.0_rt/3.0_rt;  // The ratio of 2/3 [nondim]
    constexpr Real C1_3 = 1.0_rt/3.0_rt;  // The ratio of 1/3 [nondim]
    Real c0, c1;                    // Intermediate reconstruction of q [A ~> a]
    Real d0, d1;                    // Intermediate reconstruction of h [L ~> m]
    Real b0, b1;                    // Smoothness indicator [A ~> a]

    if (u > 0.0_rt) {
        weno_three_reconstruction_0_point(&q4[1], c0); // Reconstruction in the second upwind stencil
        weno_three_reconstruction_1_point(&q4[0], c1); // Reconstruction in the first upwind stencil

        weno_three_reconstruction_0_point(&h4[1], d0);
        weno_three_reconstruction_1_point(&h4[0], d1);
        weno_three_weight_point(&u4[1], b0);  // Smoothness indicator the second upwind stencil
        weno_three_weight_point(&u4[0], b1);  // Smoothness indicator the first upwind stencil
    } else {
        Real const q32[2] = {q4[2], q4[1]};
        Real const q43[2] = {q4[3], q4[2]};
        Real const h32[2] = {h4[2], h4[1]};
        Real const h43[2] = {h4[3], h4[2]};
        Real const u32[2] = {u4[2], u4[1]};
        Real const u43[2] = {u4[3], u4[2]};
        weno_three_reconstruction_0_point(q32, c0); // Reconstruction in the second upwind stencil
        weno_three_reconstruction_1_point(q43, c1); // Reconstruction in the first upwind stencil

        weno_three_reconstruction_0_point(h32, d0);
        weno_three_reconstruction_1_point(h43, d1);
        weno_three_weight_point(u32, b0);  // Smoothness indicator the second upwind stencil
        weno_three_weight_point(u43, b1);  // Smoothness indicator the first upwind stencil
    }

    Real const tau = amrex::Math::abs(b0-b1);  // Difference of smoothness indicator [A ~> a]
    Real w0 = C2_3 * fac_fn_point(tau, b0);    // Weights [nondim]
    Real w1 = C1_3 * fac_fn_point(tau, b1);

    Real const s = 1._rt / (w0 + w1);
    w0 = w0 * s;   // Weights of stencils
    w1 = w1 * s;

    Real const vr = (w0 * c0) + (w1 * c1);  // Reconstruction of hq [A ~> a]
    Real hr = (w0 * d0) + (w1 * d1);        // Reconstruction of h [L ~> m]
    hr = amrex::min(amrex::max(h4[2], h4[1]), hr); hr = amrex::max(amrex::min(h4[2], h4[1]), hr); // A monotonicity limiter

    qr = vr / amrex::max(hr, h_tiny);
}

///
/// @brief Compute the smoothness indicator for the third upwind stencil of the fifth-order WENO scheme
///
/// @param q3 Tracer values on the three-point stencil [A ~> a]
/// @param w0 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_five_weight_0_point(Real const q3[3], Real& w0) noexcept
{
    using namespace amrex::literals;

    w0 = (q3[0] * ((10.0_rt * q3[0] - 31.0_rt * q3[1]) + 11.0_rt * q3[2])) +
         ((q3[1] * (25.0_rt * q3[1] - 19.0_rt * q3[2])) + 4.0_rt * (q3[2] * q3[2]));
}

///
/// @brief Compute the smoothness indicator for the second upwind stencil of the fifth-order WENO scheme
///
/// @param q3 Tracer values on the three-point stencil [A ~> a]
/// @param w1 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_five_weight_1_point(Real const q3[3], Real& w1) noexcept
{
    using namespace amrex::literals;

    w1 = (q3[0] * ((4.0_rt * q3[0] - 13.0_rt * q3[1]) + 5.0_rt * q3[2])) +
         ((q3[1] * (13.0_rt * q3[1] - 13.0_rt * q3[2])) + 4.0_rt * (q3[2] * q3[2]));
}

///
/// @brief Compute the smoothness indicator for the first upwind stencil of the fifth-order WENO scheme
///
/// @param q3 Tracer values on the three-point stencil [A ~> a]
/// @param w2 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_five_weight_2_point(Real const q3[3], Real& w2) noexcept
{
    using namespace amrex::literals;

    w2 = (q3[0] * ((4.0_rt * q3[0] - 19.0_rt * q3[1]) + 11.0_rt * q3[2])) +
         ((q3[1] * (25.0_rt * q3[1] - 31.0_rt * q3[2])) + 10.0_rt * (q3[2] * q3[2]));
}

///
/// @brief Reconstruction in the third upwind stencil of the fifth-order WENO scheme
///
/// @param q3 Tracer values on three points [A ~> a]
/// @param p0 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_five_reconstruction_0_point(Real const q3[3], Real& p0) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_6 = 1.0_rt/6.0_rt;  // One sixth [nondim]
    p0 = ((2.0_rt*q3[0] + 5.0_rt*q3[1]) - q3[2]) * C1_6;
}

///
/// @brief Reconstruction in the second upwind stencil of the fifth-order WENO scheme
///
/// @param q3 Tracer values on the three-point stencil [A ~> a]
/// @param p1 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_five_reconstruction_1_point(Real const q3[3], Real& p1) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_6 = 1.0_rt/6.0_rt;  // One sixth [nondim]
    p1 = ((-q3[0] + 5.0_rt*q3[1]) + 2.0_rt*q3[2]) * C1_6;
}

///
/// @brief Reconstruction in the first upwind stencil of the fifth-order WENO scheme
///
/// @param q3 Tracer values on the three-point stencil [A ~> a]
/// @param p2 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_five_reconstruction_2_point(Real const q3[3], Real& p2) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_6 = 1.0_rt/6.0_rt;  // One sixth [nondim]
    p2 = ((2.0_rt*q3[0] - 7.0_rt*q3[1]) + 11.0_rt*q3[2]) * C1_6;
}

///
/// @brief Reconstruct the tracer (e.g., PV, vorticity) onto point i-1/2 using a fifth-order WENO
///        scheme. The reconstruction is weighted by the thickness
///
/// @param q6     Tracer values on points i-3, i-2, i-1, i, i+1, i+2 [A ~> a]
/// @param h6     Thickness values on points i-3, i-2, i-1, i, i+1, i+2 [L ~> m]
/// @param u6     Values on points i-3, i-2, i-1, i, i+1, i+2 used to compute the smoothness
///               indicator: velocity [L T-1 ~> m s-1] with velocity smoothing, otherwise q6 [A ~> a]
/// @param h_tiny A tiny thickness to prevent division by zero [L ~> m]
/// @param u      Velocity or thickness flux on point i-1/2
///               [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
/// @param qr     Reconstruction of tracer q on point i-1/2 [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_five_h_weight_reconstruction_point(Real const q6[6], Real const h6[6], Real const u6[6],
                                             Real const h_tiny, Real const u, Real& qr) noexcept
{
    using namespace amrex::literals;

    constexpr Real C3_10 = 3.0_rt/10.0_rt;  // The ratio of 3/10 [nondim]
    constexpr Real C3_5 = 3.0_rt/5.0_rt;    // The ratio of 3/5 [nondim]
    constexpr Real C1_10 = 1.0_rt/10.0_rt;  // The ratio of 1/10 [nondim]
    Real c0, c1, c2;                  // Intermediate reconstruction of hq[A ~> a]
    Real d0, d1, d2;                  // Intermediate reconstruction of h [L ~> m]
    Real b0, b1, b2;                  // Smoothness indicator [A ~> a]

    if (u > 0.0_rt) {
        weno_five_reconstruction_0_point(&q6[2], c0); // Reconstruction in the third upwind stencil
        weno_five_reconstruction_1_point(&q6[1], c1); // Reconstruction in the second upwind stencil
        weno_five_reconstruction_2_point(&q6[0], c2); // Reconstruction in the first upwind stencil

        weno_five_reconstruction_0_point(&h6[2], d0);
        weno_five_reconstruction_1_point(&h6[1], d1);
        weno_five_reconstruction_2_point(&h6[0], d2);
        weno_five_weight_0_point(&u6[2], b0);   // Smoothness indicator of the third upwind stencil
        weno_five_weight_1_point(&u6[1], b1);   // Smoothness indicator of the second upwind stencil
        weno_five_weight_2_point(&u6[0], b2);   // Smoothness indicator of the first upwind stencil
    } else {
        Real const q42[3] = {q6[3], q6[2], q6[1]};
        Real const q53[3] = {q6[4], q6[3], q6[2]};
        Real const q64[3] = {q6[5], q6[4], q6[3]};
        Real const h42[3] = {h6[3], h6[2], h6[1]};
        Real const h53[3] = {h6[4], h6[3], h6[2]};
        Real const h64[3] = {h6[5], h6[4], h6[3]};
        Real const u42[3] = {u6[3], u6[2], u6[1]};
        Real const u53[3] = {u6[4], u6[3], u6[2]};
        Real const u64[3] = {u6[5], u6[4], u6[3]};
        weno_five_reconstruction_0_point(q42, c0); // Reconstruction in the third upwind stencil
        weno_five_reconstruction_1_point(q53, c1); // Reconstruction in the second upwind stencil
        weno_five_reconstruction_2_point(q64, c2); // Reconstruction in the first upwind stencil

        weno_five_reconstruction_0_point(h42, d0);
        weno_five_reconstruction_1_point(h53, d1);
        weno_five_reconstruction_2_point(h64, d2);
        weno_five_weight_0_point(u42, b0); // Smoothness indicator of the third upwind stencil
        weno_five_weight_1_point(u53, b1); // Smoothness indicator of the second upwind stencil
        weno_five_weight_2_point(u64, b2); // Smoothness indicator of the first upwind stencil
    }

    Real const tau = amrex::Math::abs(b0 - b2);  // Difference of smoothness indicators [A ~> a]
    Real w0 = C3_10 * fac_fn_point(tau, b0);     // Weights [nondim]
    Real w1 = C3_5  * fac_fn_point(tau, b1);
    Real w2 = C1_10 * fac_fn_point(tau, b2);

    Real const s = 1._rt / ((w0 + w1) + w2);
    w0 = w0 * s;   // Weights of stencils
    w1 = w1 * s;
    w2 = w2 * s;

    Real const vr = ((w0 * c0) + (w1 * c1)) + (w2 * c2);  // Reconstruction of hq [A ~> a]
    Real hr = ((w0 * d0) + (w1 * d1)) + (w2 * d2);        // Reconstruction of h [L ~> m]
    hr = amrex::min(amrex::max(h6[2], h6[3]), hr); hr = amrex::max(amrex::min(h6[2], h6[3]), hr); //Impose a monotonicity limiter

    qr = vr / amrex::max(hr, h_tiny);
}

///
/// @brief Compute the smoothness indicator for the fourth upwind stencil of the seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param w0 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_weight_0_point(Real const q4[4], Real& w0) noexcept
{
    using namespace amrex::literals;

    // Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
    w0 = ((q4[0] * ((2.107_rt * q4[0] - 9.402_rt * q4[1]) + (7.042_rt * q4[2] - 1.854_rt * q4[3]))) +
       (q4[1] * ((11.003_rt * q4[1] - 17.246_rt * q4[2]) + 4.642_rt * q4[3]))) +
       ((q4[2] * (7.043_rt * q4[2] - 3.882_rt * q4[3])) + 0.547_rt * (q4[3] * q4[3]));
}

///
/// @brief Compute the smoothness indicator for the third upwind stencil of the seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param w1 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_weight_1_point(Real const q4[4], Real& w1) noexcept
{
    using namespace amrex::literals;

    // Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
    w1 = ((q4[0] * ((0.547_rt * q4[0] - 2.522_rt * q4[1]) + (1.922_rt * q4[2] - 0.494_rt * q4[3]))) +
       (q4[1] * ((3.443_rt * q4[1] - 5.966_rt * q4[2]) + 1.602_rt * q4[3]))) +
       ((q4[2] * (2.843_rt * q4[2] - 1.642_rt * q4[3])) + 0.267_rt * (q4[3] * q4[3]));
}

///
/// @brief Compute the smoothness indicator for the second upwind stencil of the seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param w2 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_weight_2_point(Real const q4[4], Real& w2) noexcept
{
    using namespace amrex::literals;

    // Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
    w2 = ((q4[0] * ((0.267_rt * q4[0] - 1.642_rt * q4[1]) + (1.602_rt * q4[2] - 0.494_rt * q4[3]))) +
       (q4[1] * ((2.843_rt * q4[1] - 5.966_rt * q4[2]) + 1.922_rt * q4[3]))) +
       ((q4[2] * (3.443_rt * q4[2] - 2.522_rt * q4[3])) + 0.547_rt * (q4[3] * q4[3]));
}

///
/// @brief Compute smoothness indicator for the first upwind stencil of the seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param w3 Smoothness indicator for this stencil [A2 ~> a2]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_weight_3_point(Real const q4[4], Real& w3) noexcept
{
    using namespace amrex::literals;

    // Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
    w3 = ((q4[0] * ((0.547_rt  * q4[0] - 3.882_rt * q4[1]) + (4.642_rt * q4[2] - 1.854_rt * q4[3]))) +
       (q4[1] * ((7.043_rt * q4[1] - 17.246_rt * q4[2]) + 7.042_rt * q4[3]))) +
       ((q4[2] * (11.003_rt * q4[2] - 9.402_rt * q4[3])) + 2.107_rt * (q4[3] * q4[3]));
}

///
/// @brief Reconstruction in the fourth upwind stencil for seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param p0 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_reconstruction_0_point(Real const q4[4], Real& p0) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_24 = 1.0_rt/24.0_rt;  // One twenty fourth [nondim]
    p0 = (((6.0_rt * q4[0] + 26.0_rt * q4[1]) - 10.0_rt * q4[2]) + 2.0_rt * q4[3]) * C1_24;
}

///
/// @brief Reconstruction in the third upwind stencil for seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param p1 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_reconstruction_1_point(Real const q4[4], Real& p1) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_24 = 1.0_rt/24.0_rt;  // One twenty fourth [nondim]
    p1 = (14.0_rt * (q4[1] + q4[2]) - 2.0_rt * (q4[0] + q4[3])) * C1_24;
}

///
/// @brief Reconstruction in the second upwind stencil for seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param p2 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_reconstruction_2_point(Real const q4[4], Real& p2) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_24 = 1.0_rt/24.0_rt;  // One twenty fourth [nondim]
    p2 = (((2.0_rt * q4[0] - 10.0_rt * q4[1]) + 26.0_rt * q4[2]) + 6.0_rt * q4[3]) * C1_24;
}

///
/// @brief Reconstruction in the first upwind stencil for seventh-order WENO scheme
///
/// @param q4 Tracer values on the four-point stencil [A ~> a]
/// @param p3 Reconstruction of the quantity [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_reconstruction_3_point(Real const q4[4], Real& p3) noexcept
{
    using namespace amrex::literals;

    constexpr Real C1_24 = 1.0_rt/24.0_rt;  // One twenty fourth [nondim]
    p3 = (((-6.0_rt * q4[0] + 26.0_rt * q4[1]) - 46.0_rt * q4[2]) + 50.0_rt * q4[3]) * C1_24;
}

///
/// @brief Reconstruct the tracer (e.g., PV, vorticity) onto point i-1/2 using a seventh-order WENO
///        scheme. This reconstruction computes a thickness weighted average of PV
///
/// @param q8     Tracer values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3
/// @param h8     Thickness on the same tracer points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 [L ~> m]
/// @param u8     Values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 used to compute the
///               smoothness indicator: velocity [L T-1 ~> m s-1] with velocity smoothing,
///               otherwise q8
/// @param h_tiny A tiny thickness to prevent division by zero [L ~> m]
/// @param u      Velocity or thickness flux on point i-1/2
///               [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
/// @param qr     Reconstruction of tracer q on point i-1/2 [A ~> a]
///
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_seven_h_weight_reconstruction_point(Real const q8[8], Real const h8[8], Real const u8[8],
                                              Real const h_tiny, Real const u, Real& qr) noexcept
{
    using namespace amrex::literals;

    constexpr Real C4_35 = 4.0_rt/35.0_rt;    // The ratio of 4/35 [nondim]
    constexpr Real C18_35 = 18.0_rt/35.0_rt;  // The ratio of 18/35 [nondim]
    constexpr Real C12_35 = 12.0_rt/35.0_rt;  // The ratio of 12/35 [nondim]
    constexpr Real C1_35 = 1.0_rt/35.0_rt;    // The ratio of 1/35 [nondim]
    Real c0, c1, c2, c3;                // Intermediate reconstruction of hq [A ~> a]
    Real d0, d1, d2, d3;                // Intermediate reconstruction of h [L ~> m]
    Real b0, b1, b2, b3;                // Smoothness indicator [A ~> a]

    if (u > 0.0_rt) {
        weno_seven_reconstruction_0_point(&q8[3], c0); // Reconstruction in the fourth upwind stencil
        weno_seven_reconstruction_1_point(&q8[2], c1); // Reconstruction in the third upwind stencil
        weno_seven_reconstruction_2_point(&q8[1], c2); // Reconstruction in the second upwind stencil
        weno_seven_reconstruction_3_point(&q8[0], c3); // Reconstruction in the first upwind stencil

        weno_seven_reconstruction_0_point(&h8[3], d0); // Reconstruction in the fourth upwind stencil
        weno_seven_reconstruction_1_point(&h8[2], d1); // Reconstruction in the third upwind stencil
        weno_seven_reconstruction_2_point(&h8[1], d2); // Reconstruction in the second upwind stencil
        weno_seven_reconstruction_3_point(&h8[0], d3); // Reconstruction in the first upwind stencil
        weno_seven_weight_0_point(&u8[3], b0);       // Smoothness indicator of the fourth upwind stencil
        weno_seven_weight_1_point(&u8[2], b1);       // Smoothness indicator of the third upwind stencil
        weno_seven_weight_2_point(&u8[1], b2);       // Smoothness indicator of the second upwind stencil
        weno_seven_weight_3_point(&u8[0], b3);       // Smoothness indicator of the first upwind stencil
    } else {
        Real const q52[4] = {q8[4], q8[3], q8[2], q8[1]};
        Real const q63[4] = {q8[5], q8[4], q8[3], q8[2]};
        Real const q74[4] = {q8[6], q8[5], q8[4], q8[3]};
        Real const q85[4] = {q8[7], q8[6], q8[5], q8[4]};
        Real const h52[4] = {h8[4], h8[3], h8[2], h8[1]};
        Real const h63[4] = {h8[5], h8[4], h8[3], h8[2]};
        Real const h74[4] = {h8[6], h8[5], h8[4], h8[3]};
        Real const h85[4] = {h8[7], h8[6], h8[5], h8[4]};
        Real const u52[4] = {u8[4], u8[3], u8[2], u8[1]};
        Real const u63[4] = {u8[5], u8[4], u8[3], u8[2]};
        Real const u74[4] = {u8[6], u8[5], u8[4], u8[3]};
        Real const u85[4] = {u8[7], u8[6], u8[5], u8[4]};
        weno_seven_reconstruction_0_point(q52, c0); // Reconstruction in the fourth upwind stencil
        weno_seven_reconstruction_1_point(q63, c1); // Reconstruction in the third upwind stencil
        weno_seven_reconstruction_2_point(q74, c2); // Reconstruction in the second upwind stencil
        weno_seven_reconstruction_3_point(q85, c3); // Reconstruction in the first upwind stencil

        weno_seven_reconstruction_0_point(h52, d0);
        weno_seven_reconstruction_1_point(h63, d1);
        weno_seven_reconstruction_2_point(h74, d2);
        weno_seven_reconstruction_3_point(h85, d3);
        weno_seven_weight_0_point(u52, b0);    // Smoothness indicator of the fourth upwind stencil
        weno_seven_weight_1_point(u63, b1);    // Smoothness indicator of the third upwind stencil
        weno_seven_weight_2_point(u74, b2);    // Smoothness indicator of the second upwind stencil
        weno_seven_weight_3_point(u85, b3);    // Smoothness indicator of the first upwind stencil
    }

    Real const tau = amrex::Math::abs((b0 - b3) + 3.0_rt * (b1 - b2));  // Difference of smoothness indicators [A ~> a]
    Real w0 = C4_35  * fac_fn_point(tau, b0);                        // Weights [nondim]
    Real w1 = C18_35 * fac_fn_point(tau, b1);
    Real w2 = C12_35 * fac_fn_point(tau, b2);
    Real w3 = C1_35  * fac_fn_point(tau, b3);

    Real const s = 1._rt / ((w0 + w1) + (w2 + w3));
    w0 = w0 * s;   // Weights of the stencils
    w1 = w1 * s;
    w2 = w2 * s;
    w3 = w3 * s;

    Real const vr = ((w0 * c0) + (w1 * c1)) + ((w2 * c2) + (w3 * c3));  // Reconstruction of hq [A ~> a]
    Real hr = ((w0 * d0) + (w1 * d1)) + ((w2 * d2) + (w3 * d3));        // Reconstruction of h [L ~> m]

    hr = amrex::min(amrex::max(h8[3], h8[4]), hr); hr = amrex::max(amrex::min(h8[3], h8[4]), hr); // Impose a monotonicity limiter

    qr = vr / amrex::max(hr, h_tiny);
}

}
