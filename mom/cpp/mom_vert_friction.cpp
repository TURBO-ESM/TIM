// SKILLS: 0.3.3
///
/// @file mom_vert_friction.cpp
/// @brief Box-level AMReX kernels for the MOM6 vertical friction (viscosity) terms.
///
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_Reduce.H>

#include "mom_vert_friction.hpp"

#include <string>

/// @brief Abort with @p msg annotated by the source file and line number.
/// @param msg Message to report before aborting.
#define AMREX_ABORT_LOC(msg) \
	amrex::Abort(std::string(msg) + " [" + __FILE__ + ":" + std::to_string(__LINE__) + "]")

using namespace amrex;

namespace MOM {

void vertvisc_limit_vel(const Box& bxU0,
                        const Box& bxV0,
                        Array4<Real>       const& u,
                        Array4<Real>       const& v,
                        const Real dt,
                        Array4<const Real> const& areaT,
                        Array4<const Real> const& IareaT,
                        Array4<const Real> const& dx_Cv,
                        Array4<const Real> const& dy_Cu,
                        const int ke,
                        const Real Angstrom_H,
                        const Real m_s_to_L_T,
                        const vertvisc_opts_CS_C& CS,
                        Array4<const Real> const& h_u,
                        Array4<const Real> const& h_v,
                        Array4<Real>       const& u_old,
                        Array4<Real>       const& v_old,
                        Array4<Real>       const& vel_report_u,
                        Array4<Real>       const& vel_report_v,
                        Array4<int>        const& dowrite_u,
                        Array4<int>        const& dowrite_v,
                        int& ntrunc_sum)
{
    BL_PROFILE("vertvisc_limit_vel");

    const Real vel_underflow = CS.vel_underflow;
    const Real CFL_trunc     = CS.CFL_trunc;
    const Real CFL_report    = CS.CFL_report;

    const int is = bxV0.smallEnd(0);
    const int ie = bxV0.bigEnd(0);
    const int js = bxU0.smallEnd(1);
    const int je = bxU0.bigEnd(1);
    const int nz = ke;

    /// A thickness below which not to report truncations [H ~> m or kg m-2]
    const Real H_report = 3.0 * Angstrom_H;

    Box colU = bxU0;
    colU.setSmall(2, 0);
    colU.setBig(2, 0);
    Box colV = bxV0;
    colV.setSmall(2, 0);
    colV.setBig(2, 0);

    if (CS.do_report_u) {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            dowrite_u(i,j,0) = 0;
            vel_report_u(i,j,0) = 3.0e8 * m_s_to_L_T;
        });

        ReduceOps<ReduceOpLogicalOr> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(colU, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept -> ReduceTuple
        {
            int trunc_any = 0;
            for (int k = 0; k < nz; ++k) {
                Real CFL;  ///< The local CFL number [nondim]
                if (Math::abs(u(i,j,k)) < vel_underflow) u(i,j,k) = 0.0;
                if (u(i,j,k) < 0.0) {
                    CFL = (-u(i,j,k) * dt) * (dy_Cu(i,j,0) * IareaT(i+1,j,0));
                } else {
                    CFL = (u(i,j,k) * dt) * (dy_Cu(i,j,0) * IareaT(i,j,0));
                }
                if (CFL > CFL_trunc) trunc_any = 1;
                if (CFL > CFL_report) {
                    dowrite_u(i,j,0) = 1;
                    vel_report_u(i,j,0) = amrex::min(vel_report_u(i,j,0), Math::abs(u(i,j,k)));
                }
            }
            return {trunc_any};
        });
        const bool trunc_any = (amrex::get<0>(reduce_data.value(reduce_op)) != 0);

        const int klo = Box(u_old).smallEnd(2);
        const int khi = Box(u_old).bigEnd(2);
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (dowrite_u(i,j,0) != 0) {
                for (int k = klo; k <= khi; ++k) u_old(i,j,k) = u(i,j,k);
            }
        });

        if (trunc_any) {
            ReduceOps<ReduceOpSum> sum_op;
            ReduceData<int> sum_data(sum_op);
            using SumTuple = typename decltype(sum_data)::Type;

            sum_op.eval(bxU0, sum_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> SumTuple
            {
                int ntrunc = 0;
                if ((u(i,j,k) * (dt * dy_Cu(i,j,0))) * IareaT(i+1,j,0) < -CFL_trunc) {
                    u(i,j,k) = (-0.9*CFL_trunc) * (areaT(i+1,j,0) / (dt * dy_Cu(i,j,0)));
                    if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                        (h_u(i,j,k) > H_report)) ntrunc = ntrunc + 1;
                } else if ((u(i,j,k) * (dt * dy_Cu(i,j,0))) * IareaT(i,j,0) > CFL_trunc) {
                    u(i,j,k) = (0.9*CFL_trunc) * (areaT(i,j,0) / (dt * dy_Cu(i,j,0)));
                    if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                        (h_u(i,j,k) > H_report)) ntrunc = ntrunc + 1;
                }
                return {ntrunc};
            });
            ntrunc_sum = ntrunc_sum + amrex::get<0>(sum_data.value(sum_op));
        }

    } else {  /// Do not report accelerations leading to large velocities.
        ReduceOps<ReduceOpSum> sum_op;
        ReduceData<int> sum_data(sum_op);
        using SumTuple = typename decltype(sum_data)::Type;

        sum_op.eval(bxU0, sum_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> SumTuple
        {
            int ntrunc = 0;
            if (Math::abs(u(i,j,k)) < vel_underflow) {
                u(i,j,k) = 0.0;
            } else if ((u(i,j,k) * (dt * dy_Cu(i,j,0))) * IareaT(i+1,j,0) < -CFL_trunc) {
                u(i,j,k) = (-0.9*CFL_trunc) * (areaT(i+1,j,0) / (dt * dy_Cu(i,j,0)));
                if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                    (h_u(i,j,k) > H_report)) ntrunc = ntrunc + 1;
            } else if ((u(i,j,k) * (dt * dy_Cu(i,j,0))) * IareaT(i,j,0) > CFL_trunc) {
                u(i,j,k) = (0.9*CFL_trunc) * (areaT(i,j,0) / (dt * dy_Cu(i,j,0)));
                if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                    (h_u(i,j,k) > H_report)) ntrunc = ntrunc + 1;
            }
            return {ntrunc};
        });
        ntrunc_sum = ntrunc_sum + amrex::get<0>(sum_data.value(sum_op));
    }

    if (CS.do_report_v) {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            dowrite_v(i,j,0) = 0;
            vel_report_v(i,j,0) = 3.0e8 * m_s_to_L_T;
        });

        ReduceOps<ReduceOpLogicalOr> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(colV, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept -> ReduceTuple
        {
            int trunc_any = 0;
            for (int k = 0; k < nz; ++k) {
                Real CFL;  ///< The local CFL number [nondim]
                if (Math::abs(v(i,j,k)) < vel_underflow) v(i,j,k) = 0.0;
                if (v(i,j,k) < 0.0) {
                    CFL = (-v(i,j,k) * dt) * (dx_Cv(i,j,0) * IareaT(i,j+1,0));
                } else {
                    CFL = (v(i,j,k) * dt) * (dx_Cv(i,j,0) * IareaT(i,j,0));
                }
                if (CFL > CFL_trunc) trunc_any = 1;
                if (CFL > CFL_report) {
                    dowrite_v(i,j,0) = 1;
                    vel_report_v(i,j,0) = amrex::min(vel_report_v(i,j,0), Math::abs(v(i,j,k)));
                }
            }
            return {trunc_any};
        });
        const bool trunc_any = (amrex::get<0>(reduce_data.value(reduce_op)) != 0);

        const int klo = Box(v_old).smallEnd(2);
        const int khi = Box(v_old).bigEnd(2);
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (dowrite_v(i,j,0) != 0) {
                for (int k = klo; k <= khi; ++k) v_old(i,j,k) = v(i,j,k);
            }
        });

        if (trunc_any) {
            ReduceOps<ReduceOpSum> sum_op;
            ReduceData<int> sum_data(sum_op);
            using SumTuple = typename decltype(sum_data)::Type;

            sum_op.eval(bxV0, sum_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> SumTuple
            {
                int ntrunc = 0;
                if ((v(i,j,k) * (dt * dx_Cv(i,j,0))) * IareaT(i,j+1,0) < -CFL_trunc) {
                    v(i,j,k) = (-0.9*CFL_trunc) * (areaT(i,j+1,0) / (dt * dx_Cv(i,j,0)));
                    if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                        (h_v(i,j,k) > H_report)) ntrunc = ntrunc + 1;
                } else if ((v(i,j,k) * (dt * dx_Cv(i,j,0))) * IareaT(i,j,0) > CFL_trunc) {
                    v(i,j,k) = (0.9*CFL_trunc) * (areaT(i,j,0) / (dt * dx_Cv(i,j,0)));
                    if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                        (h_v(i,j,k) > H_report)) ntrunc = ntrunc + 1;
                }
                return {ntrunc};
            });
            ntrunc_sum = ntrunc_sum + amrex::get<0>(sum_data.value(sum_op));
        }

    } else {  /// Do not report accelerations leading to large velocities.
        ReduceOps<ReduceOpSum> sum_op;
        ReduceData<int> sum_data(sum_op);
        using SumTuple = typename decltype(sum_data)::Type;

        sum_op.eval(bxV0, sum_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> SumTuple
        {
            int ntrunc = 0;
            if (Math::abs(v(i,j,k)) < vel_underflow) {
                v(i,j,k) = 0.0;
            } else if ((v(i,j,k) * (dt * dx_Cv(i,j,0))) * IareaT(i,j+1,0) < -CFL_trunc) {
                v(i,j,k) = (-0.9*CFL_trunc) * (areaT(i,j+1,0) / (dt * dx_Cv(i,j,0)));
                if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                    (h_v(i,j,k) > H_report)) ntrunc = ntrunc + 1;
            } else if ((v(i,j,k) * (dt * dx_Cv(i,j,0))) * IareaT(i,j,0) > CFL_trunc) {
                v(i,j,k) = (0.9*CFL_trunc) * (areaT(i,j,0) / (dt * dx_Cv(i,j,0)));
                if (((i >= is) && (i <= ie) && (j >= js) && (j <= je)) &&
                    (h_v(i,j,k) > H_report)) ntrunc = ntrunc + 1;
            }
            return {ntrunc};
        });
        ntrunc_sum = ntrunc_sum + amrex::get<0>(sum_data.value(sum_op));
    }
}

void vertvisc_remnant_tr(const Box& bxU0,
                         const Box& bxV0,
                         Array4<const Real> const& Ray_u,
                         Array4<const Real> const& Ray_v,
                         Array4<Real>       const& visc_rem_u,
                         Array4<Real>       const& visc_rem_v,
                         const Real dt,
                         Array4<const Real> const& mask2dCu,
                         Array4<const Real> const& mask2dCv,
                         const int ke,
                         Array4<const Real> const& a_u,
                         Array4<const Real> const& a_v,
                         Array4<const Real> const& h_u,
                         Array4<const Real> const& h_v)
{
    BL_PROFILE("vertvisc_remnant_tr");

    const int nz = ke;

    Box colU = bxU0;
    colU.setSmall(2, 0);
    colU.setBig(2, 0);
    Box colV = bxV0;
    colV.setSmall(2, 0);
    colV.setBig(2, 0);

    /// A variable used by the tridiagonal solver [nondim].
    Box c1U_bx = bxU0;
    c1U_bx.setSmall(2, 0);
    c1U_bx.setBig(2, nz-1);
    FArrayBox c1U_fab(c1U_bx, 1, The_Async_Arena());
    auto const& c1U = c1U_fab.array();
    Box c1V_bx = bxV0;
    c1V_bx.setSmall(2, 0);
    c1V_bx.setBig(2, nz-1);
    FArrayBox c1V_fab(c1V_bx, 1, The_Async_Arena());
    auto const& c1V = c1V_fab.array();

    /// Find the zonal viscous remnant using a modification of a standard tridagonal solver.
    if (Ray_u.p != nullptr) {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCu(i,j,0) > 0.0)) return;
            Real Ray = Ray_u(i,j,0);

            Real b_denom_1 = h_u(i,j,0) + dt * (Ray + a_u(i,j,0));
            Real b1 = 1.0 / (b_denom_1 + dt * a_u(i,j,1));
            Real d1 = b_denom_1 * b1;
            visc_rem_u(i,j,0) = b1 * h_u(i,j,0);

            for (int k = 1; k < nz; ++k) {
                Ray = Ray_u(i,j,k);

                c1U(i,j,k) = dt * a_u(i,j,k) * b1;
                b_denom_1 = h_u(i,j,k) + dt * (Ray + a_u(i,j,k) * d1);
                b1 = 1.0 / (b_denom_1 + dt * a_u(i,j,k+1));
                d1 = b_denom_1 * b1;
                visc_rem_u(i,j,k) = (h_u(i,j,k) + dt * a_u(i,j,k) * visc_rem_u(i,j,k-1)) * b1;

                /// Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
                /// XXX: Check GPU behavior
                if (dt < 0) break;
            }

            for (int k = nz-2; k >= 0; --k) {
                visc_rem_u(i,j,k) = visc_rem_u(i,j,k) + c1U(i,j,k+1) * visc_rem_u(i,j,k+1);
            }
        });
    } else {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCu(i,j,0) > 0.0)) return;
            const Real Ray = 0.0;

            Real b_denom_1 = h_u(i,j,0) + dt * (Ray + a_u(i,j,0));
            Real b1 = 1.0 / (b_denom_1 + dt * a_u(i,j,1));
            Real d1 = b_denom_1 * b1;
            visc_rem_u(i,j,0) = b1 * h_u(i,j,0);

            for (int k = 1; k < nz; ++k) {
                c1U(i,j,k) = dt * a_u(i,j,k) * b1;
                b_denom_1 = h_u(i,j,k) + dt * (Ray + a_u(i,j,k) * d1);
                b1 = 1.0 / (b_denom_1 + dt * a_u(i,j,k+1));
                d1 = b_denom_1 * b1;
                visc_rem_u(i,j,k) = (h_u(i,j,k) + dt * a_u(i,j,k) * visc_rem_u(i,j,k-1)) * b1;

                /// Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
                /// XXX: Check GPU behavior
                if (dt < 0) break;
            }

            for (int k = nz-2; k >= 0; --k) {
                visc_rem_u(i,j,k) = visc_rem_u(i,j,k) + c1U(i,j,k+1) * visc_rem_u(i,j,k+1);
            }
        });
    }

    /// Now find the meridional viscous remnant using the robust tridiagonal solver.
    if (Ray_v.p != nullptr) {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCv(i,j,0) > 0.0)) return;
            Real Ray = Ray_v(i,j,0);

            Real b_denom_1 = h_v(i,j,0) + dt * (Ray + a_v(i,j,0));
            Real b1 = 1.0 / (b_denom_1 + dt * a_v(i,j,1));
            Real d1 = b_denom_1 * b1;
            visc_rem_v(i,j,0) = b1 * h_v(i,j,0);

            for (int k = 1; k < nz; ++k) {
                Ray = Ray_v(i,j,k);

                c1V(i,j,k) = dt * a_v(i,j,k) * b1;
                b_denom_1 = h_v(i,j,k) + dt * (Ray + a_v(i,j,k) * d1);
                b1 = 1.0 / (b_denom_1 + dt * a_v(i,j,k+1));
                d1 = b_denom_1 * b1;
                visc_rem_v(i,j,k) = (h_v(i,j,k) + dt * a_v(i,j,k) * visc_rem_v(i,j,k-1)) * b1;

                /// Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
                /// XXX: Check GPU behavior
                if (dt < 0) break;
            }

            for (int k = nz-2; k >= 0; --k) {
                visc_rem_v(i,j,k) = visc_rem_v(i,j,k) + c1V(i,j,k+1) * visc_rem_v(i,j,k+1);
            }
        });
    } else {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCv(i,j,0) > 0.0)) return;
            const Real Ray = 0.0;

            Real b_denom_1 = h_v(i,j,0) + dt * (Ray + a_v(i,j,0));
            Real b1 = 1.0 / (b_denom_1 + dt * a_v(i,j,1));
            Real d1 = b_denom_1 * b1;
            visc_rem_v(i,j,0) = b1 * h_v(i,j,0);

            for (int k = 1; k < nz; ++k) {
                c1V(i,j,k) = dt * a_v(i,j,k) * b1;
                b_denom_1 = h_v(i,j,k) + dt * (Ray + a_v(i,j,k) * d1);
                b1 = 1.0 / (b_denom_1 + dt * a_v(i,j,k+1));
                d1 = b_denom_1 * b1;
                visc_rem_v(i,j,k) = (h_v(i,j,k) + dt * a_v(i,j,k) * visc_rem_v(i,j,k-1)) * b1;

                /// Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
                /// XXX: Check GPU behavior
                if (dt < 0) break;
            }

            for (int k = nz-2; k >= 0; --k) {
                visc_rem_v(i,j,k) = visc_rem_v(i,j,k) + c1V(i,j,k+1) * visc_rem_v(i,j,k+1);
            }
        });
    }
}

void vertfpmix_tr(const Box& bxU0,
                  const Box& bxV0,
                  Array4<Real>       const& ui,
                  Array4<Real>       const& vi,
                  Array4<const Real> const& hbl_h,
                  Array4<const Real> const& taux,
                  Array4<const Real> const& tauy,
                  Array4<const Real> const& Us_x,
                  Array4<const Real> const& Us_y,
                  const Real dt,
                  const Real Cemp_NL,
                  const Real Irho0,
                  Array4<const Real> const& mask2dCu,
                  Array4<const Real> const& mask2dCv,
                  Array4<const Real> const& mask2dT,
                  Array4<const Real> const& h_u,
                  Array4<const Real> const& h_v,
                  Array4<Real>       const& uE_u,
                  Array4<Real>       const& vE_v,
                  Array4<Real>       const& uInc_u,
                  Array4<Real>       const& vInc_v)
{
    BL_PROFILE("vertfpmix_tr");

    const int klo = 0;
    const int khi = bxU0.bigEnd(2);

    Box colU = bxU0;
    colU.setSmall(2, 0);
    colU.setBig(2, 0);
    Box colV = bxV0;
    colV.setSmall(2, 0);
    colV.setBig(2, 0);

    /// u-points
    ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        Real hbl_u;          ///< boundary layer depth (u-pts) [H ~> m]
        Real taux_u;         ///< kinematic zonal wind stress (u-pts) [L2 T-2 ~> m2 s-2]
        Real tmp_u;          ///< temporary ocean mask weights on u points [nondim]
        Real sigma;          ///< temporary normalize boundary layer coordinate [nondim]
        Real Gat1, Gsig, dGdsig; ///< Shape parameters [nondim]
        Real depth;          ///< Cumulative of thicknesses [H ~> m]

        taux_u  = taux(i,j,0) * Irho0;
        if ( (mask2dCu(i,j,0) > 0.5) ) {
            /// h to u-pts
            tmp_u  = amrex::max(1.0 ,(mask2dT(i,j,0) + mask2dT(i+1,j,0) ) );
            hbl_u = ((mask2dT(i,j,0) * hbl_h(i,j,0)) + (mask2dT(i+1,j,0) * hbl_h(i+1,j,0))) / tmp_u;
            depth   = 0.0;
            Gat1  = 0.0;
            for (int k = klo; k <= khi; ++k) {
                /// cell center
                depth = depth + 0.5*h_u(i,j,k);
                uE_u(i,j,k) = ui(i,j,k) - Us_x(i,j,k);
                if ( depth < hbl_u ) {
                    sigma = depth / hbl_u;
                    /// cell bottom
                    depth = depth + 0.5*h_u(i,j,k);
                    fpmix_gshape_point(sigma, Gat1, Gsig, dGdsig);
                    /// nonlocal boundary-layer increment
                    uInc_u(i,j,k)  = dt * Cemp_NL * taux_u * dGdsig / hbl_u;
                    ui(i,j,k) = ui(i,j,k) + uInc_u(i,j,k);
                } else {
                    uInc_u(i,j,k) = 0.0;
                }
            }
        } else {
            for (int k = klo; k <= khi; ++k) {
                uInc_u(i,j,k) = 0.0;
            }
        }
    });

    /// v-points
    ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        Real hbl_v;          ///< boundary layer depth (v-pts) [H ~> m]
        Real tauy_v;         ///< kinematic merid wind stress (v-pts) [L2 T-2 ~> m2 s-2]
        Real tmp_v;          ///< temporary ocean mask weights on v points [nondim]
        Real sigma;          ///< temporary normalize boundary layer coordinate [nondim]
        Real Gat1, Gsig, dGdsig; ///< Shape parameters [nondim]
        Real depth;          ///< Cumulative of thicknesses [H ~> m]

        tauy_v  = tauy(i,j,0) * Irho0;
        if ( (mask2dCv(i,j,0) > 0.5) ) {
            /// h to v-pts
            tmp_v  = amrex::max( 1.0 ,(mask2dT(i,j,0) + mask2dT(i,j+1,0)));
            hbl_v = (mask2dT(i,j,0) * hbl_h(i,j,0) + mask2dT(i,j+1,0) * hbl_h(i,j+1,0)) / tmp_v;
            depth = 0.0;
            Gat1  = 0.0;
            for (int k = klo; k <= khi; ++k) {
                /// cell center
                depth = depth + 0.5* h_v(i,j,k);
                vE_v(i,j,k) = vi(i,j,k) - Us_y(i,j,k);
                if ( depth < hbl_v ) {
                    sigma = depth / hbl_v;
                    /// cell bottom
                    depth = depth + 0.5* h_v(i,j,k);
                    fpmix_gshape_point(sigma, Gat1, Gsig, dGdsig);
                    /// nonlocal boundary-layer increment
                    vInc_v(i,j,k) = dt * Cemp_NL * tauy_v * dGdsig / hbl_v;
                    vi(i,j,k) = vi(i,j,k) + vInc_v(i,j,k);
                } else {
                    vInc_v(i,j,k)  = 0.0;
                }
            }
        } else {
            for (int k = klo; k <= khi; ++k) {
                vInc_v(i,j,k)  = 0.0;
            }
        }
    });
}

void vertvisc_coef_tr(const Box& bxU0,
                      const Box& bxV0,
                      Array4<const Real> const& u,
                      Array4<const Real> const& v,
                      Array4<const Real> const& h,
                      Array4<const Real> const& dz,
                      Array4<const Real> const& frac_shelf_u,
                      Array4<const Real> const& frac_shelf_v,
                      Array4<const Real> const& Ustar_2d,
                      Array4<const Real> const& Kv_bbl_u,
                      Array4<const Real> const& Kv_bbl_v,
                      Array4<const Real> const& bbl_thick_u,
                      Array4<const Real> const& bbl_thick_v,
                      Array4<const Real> const& tbl_thick_shelf_u,
                      Array4<const Real> const& tbl_thick_shelf_v,
                      Array4<const Real> const& Kv_tbl_shelf_u,
                      Array4<const Real> const& Kv_tbl_shelf_v,
                      Array4<const Real> const& nkml_visc_u,
                      Array4<const Real> const& nkml_visc_v,
                      Array4<const Real> const& Kv_shear,
                      Array4<const Real> const& Kv_shear_Bu,
                      Array4<const Real> const& SpV_avg,
                      const Real dt,
                      Array4<const Real> const& mask2dCu,
                      Array4<const Real> const& mask2dCv,
                      Array4<const Real> const& bathyT,
                      Array4<const Real> const& CoriolisBu,
                      const int ke,
                      const bool Boussinesq,
                      const int nkml,
                      const Real dZ_subroundoff,
                      const Real H_subroundoff,
                      const Real H_to_m,
                      const Real m_to_H,
                      const Real RZ_to_H,
                      const Real Z_to_H,
                      const Real T_to_s,
                      Array4<const Real> const& g_prime,
                      const coef_opts_CS_C& CS,
                      Array4<const Real> const& kappa_gl90_2d,
                      Array4<Real>       const& a_u,
                      Array4<Real>       const& a_u_gl90,
                      Array4<Real>       const& h_u,
                      Array4<Real>       const& a_v,
                      Array4<Real>       const& a_v_gl90,
                      Array4<Real>       const& h_v,
                      Array4<Real>       const& a1_shelf_u,
                      Array4<Real>       const& a1_shelf_v,
                      OceanOBC* OBC,
                      Array4<const Real> const& kdgl90_struct,
                      Array4<Real>       const& hML_u,
                      Array4<Real>       const& hML_v,
                      Array4<Real>       const& Kv_u,
                      Array4<Real>       const& Kv_v,
                      Array4<Real>       const& Kv_gl90_u,
                      Array4<Real>       const& Kv_gl90_v)
{
    BL_PROFILE("vertvisc_coef_tr");

    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    const int nz = ke;

    const coupling_CS_C cpl      = CS.coupling_CS;
    const Real Hbbl              = CS.Hbbl;
    const Real Hbbl_gl90         = CS.Hbbl_gl90;
    const bool use_GL90_in_SSW   = CS.use_GL90_in_SSW;
    const bool harmonic_visc     = CS.harmonic_visc;
    const Real harm_BL_val       = CS.harm_BL_val;
    const bool use_GL90_N2       = CS.GL90_CS.use_GL90_N2;
    const Real kappa_gl90        = CS.GL90_CS.kappa_gl90;
    const bool read_kappa_gl90   = CS.GL90_CS.read_kappa_gl90;
    const Real alpha_gl90        = CS.GL90_CS.alpha_gl90;

    /// A thickness that is so small it is usually lost in roundoff and can be neglected
    /// [H ~> m or kg m-2].
    const Real h_neglect = H_subroundoff;
    /// A vertical distance that is so small it is usually lost in roundoff and can be neglected
    /// [Z ~> m].
    const Real dz_neglect = dZ_subroundoff;
    /// The maximum drag coefficient across interfaces, set so that it will be representable as a
    /// 32-bit float in MKS units [H T-1 ~> m s-1 or Pa s m-1]
    const Real a_cpl_max = 1.0e37 * m_to_H * T_to_s;
    /// The inverse of a scaling factor determining when water is still within the boundary layer,
    /// as determined by the sum of the harmonic mean thicknesses [nondim].
    Real I_valBL = 0.0;
    if (harm_BL_val > 0.0) I_valBL = 1.0 / harm_BL_val;

    /// First do u-points
    {
        Box colU = bxU0;
        colU.setSmall(2, 0);
        colU.setBig(2, 0);
        Box sbx = bxU0;
        sbx.setSmall(2, 0);
        sbx.setBig(2, nz);
        FArrayBox hvel_fab        (sbx, 1, The_Async_Arena());
        FArrayBox dz_harm_fab     (sbx, 1, The_Async_Arena());
        FArrayBox dz_vel_fab      (sbx, 1, The_Async_Arena());
        FArrayBox hvel_shelf_fab  (sbx, 1, The_Async_Arena());
        FArrayBox dz_vel_shelf_fab(sbx, 1, The_Async_Arena());
        FArrayBox z_i_fab         (sbx, 1, The_Async_Arena());
        FArrayBox z_i_gl90_fab    (sbx, 1, The_Async_Arena());
        FArrayBox a_cpl_fab       (sbx, 1, The_Async_Arena());
        FArrayBox a_cpl_gl90_fab  (sbx, 1, The_Async_Arena());
        FArrayBox a_shelf_fab     (sbx, 1, The_Async_Arena());
        auto const& hvel         = hvel_fab.array();
        auto const& dz_harm      = dz_harm_fab.array();
        auto const& dz_vel       = dz_vel_fab.array();
        auto const& hvel_shelf   = hvel_shelf_fab.array();
        auto const& dz_vel_shelf = dz_vel_shelf_fab.array();
        auto const& z_i          = z_i_fab.array();
        auto const& z_i_gl90     = z_i_gl90_fab.array();
        auto const& a_cpl        = a_cpl_fab.array();
        auto const& a_cpl_gl90   = a_cpl_gl90_fab.array();
        auto const& a_shelf      = a_shelf_fab.array();

        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCu(i,j,0) > 0.0)) return;

            Real I_Hbbl, I_Hbbl_gl90 = 0.0, kv_bbl = 0.0, bbl_thick = 0.0, Dmin, zh = 0.0;
            Real zcol = 0.0, zcol_p1 = 0.0, h_harm, h_arith, h_delta, dz_arith, z2, botfn;
            Real z_clear, z2_wt, h_ml, I_HTbl = 0.0, Ztop_min = 0.0, topfn;
            int  zi_dir;
            bool do_any_shelf;

            I_Hbbl = 1.0 / (Hbbl + dz_neglect);
            if (use_GL90_in_SSW) {
                I_Hbbl_gl90 = 1.0 / (Hbbl_gl90 + dz_neglect);
            }

            if (cpl.bottomdraglaw) {
                kv_bbl = Kv_bbl_u(i,j,0);
                bbl_thick = bbl_thick_u(i,j,0) + dz_neglect;
                I_Hbbl = 1.0 / bbl_thick;
            }

            Dmin = amrex::min(bathyT(i,j,0), bathyT(i+1,j,0));
            zi_dir = 0;

            /// Project thickness outward across OBCs using a zero-gradient condition.
            /*
            if (OBC != nullptr) {
                if (OBC->u_E_OBCs_on_PE) {
                    if (OBC->segnum_u(i,j,0) > 0) {
                        Dmin = bathyT(i,j,0);
                        zi_dir = -1;
                    }
                }
                if (OBC->u_W_OBCs_on_PE) {
                    if (OBC->segnum_u(i,j,0) < 0) {
                        Dmin = bathyT(i+1,j,0);
                        zi_dir = 1;
                    }
                }
            }
            */

            /// The following block calculates the thicknesses at velocity grid points for
            /// the vertical viscosity (hvel and dz_vel).  Near the bottom an upwind biased
            /// thickness is used to control the effect of spurious Montgomery potential
            /// gradients at the bottom where nearly massless layers layers ride over the
            /// topography.

            z_i(i,j,nz) = 0.0;

            if (!harmonic_visc) {
                zh = 0.0;
                zcol = -bathyT(i,j,0);
                zcol_p1 = -bathyT(i+1,j,0);
            }

            if (use_GL90_in_SSW) {
                z_i_gl90(i,j,nz) = 0.0;
            }

            for (int k = nz-1; k >= 0; --k) {
                h_harm = 2.0 * h(i,j,k) * h(i+1,j,k) / (h(i,j,k) + h(i+1,j,k) + h_neglect);
                h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k));
                h_delta = h(i+1,j,k) - h(i,j,k);
                dz_harm(i,j,k) = 2.0 * dz(i,j,k) * dz(i+1,j,k) / (dz(i,j,k) + dz(i+1,j,k) + dz_neglect);
                dz_arith = 0.5 * (dz(i+1,j,k) + dz(i,j,k));

                /// Project thickness outward across OBCs using a zero-gradient condition.
                /*
                if (OBC != nullptr) {
                    if (OBC->u_E_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) > 0) {
                            h_harm = h(i,j,k);
                            h_arith = h(i,j,k);
                            h_delta = 0.0;
                            dz_harm(i,j,k) = dz(i,j,k);
                            dz_arith = dz(i,j,k);
                        }
                    }
                    if (OBC->u_W_OBCs_on_PE) {
                        if (OBC->segnum_u(i,j,0) < 0) {
                            h_harm = h(i+1,j,k);
                            h_arith = h(i+1,j,k);
                            h_delta = 0.0;
                            dz_harm(i,j,k) = dz(i+1,j,k);
                            dz_arith = dz(i+1,j,k);
                        }
                    }
                }
                */

                if (harmonic_visc) {
                    /// The following block calculates the thicknesses at velocity grid points
                    /// for the vertical viscosity (hvel and dz_vel).  Near the bottom an
                    /// upwind biased thickness is used to control the effect of spurious
                    /// Montgomery potential gradients at the bottom where nearly massless
                    /// layers ride over the topography.

                    hvel(i,j,k) = h_harm;
                    dz_vel(i,j,k) = dz_harm(i,j,k);

                    if (u(i,j,k) * h_delta < 0) {
                        z2 = z_i(i,j,k+1);
                        botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

                        hvel(i,j,k) = (1.0 - botfn) * h_harm + botfn * h_arith;
                        dz_vel(i,j,k) = (1.0 - botfn) * dz_harm(i,j,k) + botfn * dz_arith;
                    }

                    z_i(i,j,k) =  z_i(i,j,k+1) + dz_harm(i,j,k) * I_Hbbl;
                } else {
                    zcol = zcol + dz(i,j,k);
                    zcol_p1 = zcol_p1 + dz(i+1,j,k);

                    zh = zh + dz_harm(i,j,k);

                    z_clear = amrex::max(zcol, zcol_p1) + Dmin;
                    if (zi_dir < 0) z_clear = zcol + Dmin;
                    if (zi_dir > 0) z_clear = zcol_p1 + Dmin;

                    z_i(i,j,k) = amrex::max(zh, z_clear) * I_Hbbl;

                    hvel(i,j,k) = h_arith;
                    dz_vel(i,j,k) = dz_arith;

                    if (u(i,j,k) * h_delta > 0.0) {
                        if (zh * I_Hbbl < harm_BL_val) {
                            hvel(i,j,k) = h_harm;
                            dz_vel(i,j,k) = dz_harm(i,j,k);
                        } else {
                            z2_wt = 1.0;
                            if (zh * I_Hbbl < 2.0 * harm_BL_val)
                                z2_wt = amrex::max(0.0, amrex::min(1.0, zh * I_Hbbl * I_valBL - 1.0));

                            z2 = z2_wt * (amrex::max(zh, z_clear) * I_Hbbl);
                            botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

                            hvel(i,j,k) = (1.0 - botfn) * h_arith + botfn * h_harm;
                            dz_vel(i,j,k) = (1.0 - botfn) * dz_arith + botfn * dz_harm(i,j,k);
                        }
                    }
                }

                if (use_GL90_in_SSW) {
                    /// The following block calculates the normalized height above the GL90 BBL
                    /// (z_i_gl90), using a harmonic mean between layer thicknesses. For the
                    /// GL90 BBL we use simply a constant (Hbbl_gl90). The purpose isthat the
                    /// GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
                    /// that no momentum gets fluxed into vanished layers. The scheme is not
                    /// sensitive to the exact value of Hbbl_gl90, as long as it is in a
                    /// reasonable range (~1-20 m): large enough to capture vanished layers
                    /// over topography, small enough to not contaminate the interior.

                    z_i_gl90(i,j,k) = z_i_gl90(i,j,k+1) + dz_harm(i,j,k) * I_Hbbl_gl90;
                }
            }

            find_coupling_coef_point<true, false>(i, j, a_cpl, dz_vel, dz_harm, bbl_thick, kv_bbl, z_i,
                h_ml, dt, ke, Boussinesq, nkml, dZ_subroundoff, H_to_m, RZ_to_H, Z_to_H, CoriolisBu,
                cpl, Kv_shear, Kv_shear_Bu, Kv_tbl_shelf_u, Kv_tbl_shelf_v, tbl_thick_shelf_u,
                tbl_thick_shelf_v, nkml_visc_u, nkml_visc_v, Ustar_2d, SpV_avg);

            if (hML_u.p != nullptr) hML_u(i,j,0) = h_ml;

            if (use_GL90_in_SSW) {
                find_coupling_coef_gl90_point<true>(i, j, a_cpl_gl90, dz_vel, z_i_gl90, ke,
                    dZ_subroundoff, g_prime, CoriolisBu, use_GL90_N2, kappa_gl90, read_kappa_gl90,
                    alpha_gl90, kappa_gl90_2d, kdgl90_struct);
            }

            do_any_shelf = false;
            if (frac_shelf_u.p != nullptr) {
                a1_shelf_u(i,j,0) = 0.0;
                do_any_shelf = frac_shelf_u(i,j,0) > 0.0;

                if (do_any_shelf) {
                    if (!harmonic_visc) {
                        zh = 0.0;
                        Ztop_min = amrex::min(zcol, zcol_p1);
                        I_HTbl = 1.0 / (tbl_thick_shelf_u(i,j,0) + dz_neglect);
                    }

                    for (int k = 0; k < nz; ++k) {
                        if (harmonic_visc) {
                            hvel_shelf(i,j,k) = hvel(i,j,k);
                            dz_vel_shelf(i,j,k) = dz_vel(i,j,k);
                        } else {
                            /// Find upwind-biased thickness near the surface.
                            /// (Perhaps this needs to be done more carefully, via find_eta.)

                            h_harm = 2.0 * h(i,j,k) * h(i+1,j,k)
                                / (h(i,j,k) + h(i+1,j,k) + h_neglect);
                            h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k));
                            h_delta = h(i+1,j,k) - h(i,j,k);
                            dz_arith = 0.5 * (dz(i+1,j,k) + dz(i,j,k));

                            /*
                            if (OBC != nullptr) {
                                if (OBC->u_E_OBCs_on_PE) {
                                    if (OBC->segnum_u(i,j,0) > 0) {
                                        h_harm = h(i,j,k);
                                        h_arith = h(i,j,k);
                                        h_delta = 0.0;
                                        dz_arith = dz(i,j,k);
                                    }
                                }
                                if (OBC->u_W_OBCs_on_PE) {
                                    if (OBC->segnum_u(i,j,0) < 0) {
                                        h_harm = h(i+1,j,k);
                                        h_arith = h(i+1,j,k);
                                        h_delta = 0.0;
                                        dz_arith = dz(i+1,j,k);
                                    }
                                }
                            }
                            */

                            zcol = zcol - dz(i,j,k);
                            zcol_p1 = zcol_p1 - dz(i+1,j,k);

                            zh = zh + dz_harm(i,j,k);

                            hvel_shelf(i,j,k) = hvel(i,j,k);
                            dz_vel_shelf(i,j,k) = dz_vel(i,j,k);

                            if (u(i,j,k) * h_delta > 0.0) {
                                if (zh * I_HTbl < harm_BL_val) {
                                    hvel_shelf(i,j,k) = amrex::min(hvel(i,j,k), h_harm);
                                    dz_vel_shelf(i,j,k) = amrex::min(dz_vel(i,j,k), dz_harm(i,j,k));
                                } else {
                                    z2_wt = 1.0;
                                    if (zh * I_HTbl < 2.0 * harm_BL_val) {
                                        z2_wt = amrex::max(0.0, amrex::min(1.0, zh * I_HTbl * I_valBL - 1.0));
                                    }

                                    z2 = z2_wt * (amrex::max(zh, Ztop_min - amrex::min(zcol, zcol_p1)) * I_HTbl);
                                    /// TODO: replace **6 with multiply
                                    Real const z2_3 = (z2 * z2) * z2;
                                    topfn = 1.0 / (1.0 + 0.09 * (z2_3 * z2_3));

                                    hvel_shelf(i,j,k) = amrex::min(hvel(i,j,k), (1.0 - topfn) * h_arith + topfn * h_harm);
                                    dz_vel_shelf(i,j,k) = amrex::min(dz_vel(i,j,k), (1.0 - topfn) * dz_arith + topfn * dz_harm(i,j,k));
                                }
                            }
                        }
                    }

                    find_coupling_coef_point<true, true>(i, j, a_shelf, dz_vel_shelf, dz_harm, bbl_thick,
                        kv_bbl, z_i, h_ml, dt, ke, Boussinesq, nkml, dZ_subroundoff, H_to_m, RZ_to_H,
                        Z_to_H, CoriolisBu, cpl, Kv_shear, Kv_shear_Bu, Kv_tbl_shelf_u, Kv_tbl_shelf_v,
                        tbl_thick_shelf_u, tbl_thick_shelf_v, nkml_visc_u, nkml_visc_v, Ustar_2d, SpV_avg);

                    a1_shelf_u(i,j,0) = a_shelf(i,j,0);
                }
            }

            if (do_any_shelf) {
                if (use_GL90_in_SSW) {
                    for (int k = 0; k <= nz; ++k) {
                        a_u(i,j,k) = amrex::min(a_cpl_max, (frac_shelf_u(i,j,0) * a_shelf(i,j,k) +
                                                            (1.0 - frac_shelf_u(i,j,0)) * a_cpl(i,j,k)) + a_cpl_gl90(i,j,k));

                        /// This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
                        /// CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j) * max(a_shelf(K), a_cpl(K)) + &
                        ///                                (1. - forces%frac_shelf_u(I,j)) * a_cpl(K))

                        a_u_gl90(i,j,k) = amrex::min(a_cpl_max, a_cpl_gl90(i,j,k));
                    }
                } else {
                    for (int k = 0; k <= nz; ++k) {
                        a_u(i,j,k) = amrex::min(a_cpl_max, (frac_shelf_u(i,j,0) * a_shelf(i,j,k) +
                                                            (1.0 - frac_shelf_u(i,j,0)) * a_cpl(i,j,k)));
                        /// This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
                        /// CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j) * max(a_shelf(K), a_cpl(K)) + &
                        ///                                (1. - forces%frac_shelf_u(I,j)) * a_cpl(K))
                    }
                }

                for (int k = 0; k < nz; ++k) {
                    /// Should we instead take the inverse of the average of the inverses?
                    h_u(i,j,k) = frac_shelf_u(i,j,0) * hvel_shelf(i,j,k)
                        + (1.0 - frac_shelf_u(i,j,0)) * hvel(i,j,k) + h_neglect;
                }
            } else {
                if (use_GL90_in_SSW) {
                    for (int k = 0; k <= nz; ++k) {
                        a_cpl(i,j,k) = a_cpl(i,j,k) + a_cpl_gl90(i,j,k);
                    }

                    for (int k = 0; k <= nz; ++k) {
                        a_u_gl90(i,j,k) = amrex::min(a_cpl_max, a_cpl_gl90(i,j,k));
                    }
                }

                for (int k = 0; k <= nz; ++k) {
                    a_u(i,j,k) = amrex::min(a_cpl_max, a_cpl(i,j,k));
                }

                for (int k = 0; k < nz; ++k) {
                    h_u(i,j,k) = hvel(i,j,k) + h_neglect;
                }
            }

            /// Diagnose total Kv at u-points
            if (Kv_u.p != nullptr) {
                for (int k = 0; k < nz; ++k) {
                    Kv_u(i,j,k) = 0.5 * (a_u(i,j,k) + a_u(i,j,k+1)) * h_u(i,j,k);
                }
            }

            /// Diagnose GL90 Kv at u-points
            if (Kv_gl90_u.p != nullptr) {
                for (int k = 0; k < nz; ++k) {
                    Kv_gl90_u(i,j,k) = 0.5 * (a_u_gl90(i,j,k) + a_u_gl90(i,j,k+1)) * h_u(i,j,k);
                }
            }
        });
    }

    /// Now work on v-points.
    {
        Box colV = bxV0;
        colV.setSmall(2, 0);
        colV.setBig(2, 0);
        Box sbx = bxV0;
        sbx.setSmall(2, 0);
        sbx.setBig(2, nz);
        FArrayBox hvel_fab        (sbx, 1, The_Async_Arena());
        FArrayBox dz_harm_fab     (sbx, 1, The_Async_Arena());
        FArrayBox dz_vel_fab      (sbx, 1, The_Async_Arena());
        FArrayBox hvel_shelf_fab  (sbx, 1, The_Async_Arena());
        FArrayBox dz_vel_shelf_fab(sbx, 1, The_Async_Arena());
        FArrayBox z_i_fab         (sbx, 1, The_Async_Arena());
        FArrayBox z_i_gl90_fab    (sbx, 1, The_Async_Arena());
        FArrayBox a_cpl_fab       (sbx, 1, The_Async_Arena());
        FArrayBox a_cpl_gl90_fab  (sbx, 1, The_Async_Arena());
        FArrayBox a_shelf_fab     (sbx, 1, The_Async_Arena());
        auto const& hvel         = hvel_fab.array();
        auto const& dz_harm      = dz_harm_fab.array();
        auto const& dz_vel       = dz_vel_fab.array();
        auto const& hvel_shelf   = hvel_shelf_fab.array();
        auto const& dz_vel_shelf = dz_vel_shelf_fab.array();
        auto const& z_i          = z_i_fab.array();
        auto const& z_i_gl90     = z_i_gl90_fab.array();
        auto const& a_cpl        = a_cpl_fab.array();
        auto const& a_cpl_gl90   = a_cpl_gl90_fab.array();
        auto const& a_shelf      = a_shelf_fab.array();

        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCv(i,j,0) > 0.0)) return;

            Real I_Hbbl, I_Hbbl_gl90 = 0.0, kv_bbl = 0.0, bbl_thick = 0.0, Dmin, zh = 0.0;
            Real zcol = 0.0, zcol_p1 = 0.0, h_harm, h_arith, h_delta, dz_arith, z2, botfn;
            Real z_clear, z2_wt, h_ml, I_HTbl = 0.0, Ztop_min = 0.0, topfn;
            int  zi_dir;
            bool do_any_shelf;

            I_Hbbl = 1.0 / (Hbbl + dz_neglect);
            if (use_GL90_in_SSW) {
                I_Hbbl_gl90 = 1.0 / (Hbbl_gl90 + dz_neglect);
            }

            if (cpl.bottomdraglaw) {
                kv_bbl = Kv_bbl_v(i,j,0);
                bbl_thick = bbl_thick_v(i,j,0) + dz_neglect;
                I_Hbbl = 1.0 / bbl_thick;
            }

            Dmin = amrex::min(bathyT(i,j,0), bathyT(i,j+1,0));
            zi_dir = 0;

            /// Project thickness outward across OBCs using a zero-gradient condition.
            /*
            if (OBC != nullptr) {
                if (OBC->v_N_OBCs_on_PE) {
                    if (OBC->segnum_v(i,j,0) > 0) {
                        Dmin = bathyT(i,j,0);
                        zi_dir = -1;
                    }
                }
                if (OBC->v_S_OBCs_on_PE) {
                    if (OBC->segnum_v(i,j,0) < 0) {
                        Dmin = bathyT(i,j+1,0);
                        zi_dir = 1;
                    }
                }
            }
            */

            /// The following block calculates the thicknesses at velocity grid points for
            /// the vertical viscosity (hvel and dz_vel).  Near the bottom an upwind biased
            /// thickness is used to control the effect of spurious Montgomery potential
            /// gradients at the bottom where nearly massless layers layers ride over the
            /// topography.

            z_i(i,j,nz) = 0.0;

            if (!harmonic_visc) {
                zh = 0.0;
                zcol = -bathyT(i,j,0);
                zcol_p1 = -bathyT(i,j+1,0);
            }

            if (use_GL90_in_SSW) {
                z_i_gl90(i,j,nz) = 0.0;
            }

            for (int k = nz-1; k >= 0; --k) {
                h_harm = 2.0 * h(i,j,k) * h(i,j+1,k) / (h(i,j,k) + h(i,j+1,k) + h_neglect);
                h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k));
                h_delta = h(i,j+1,k) - h(i,j,k);
                dz_harm(i,j,k) = 2.0 * dz(i,j,k) * dz(i,j+1,k) / (dz(i,j,k) + dz(i,j+1,k) + dz_neglect);
                dz_arith = 0.5 * (dz(i,j+1,k) + dz(i,j,k));

                /// Project thickness outward across OBCs using a zero-gradient condition.
                /*
                if (OBC != nullptr) {
                    if (OBC->v_N_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) > 0) {
                            h_harm = h(i,j,k);
                            h_arith = h(i,j,k);
                            h_delta = 0.0;
                            dz_harm(i,j,k) = dz(i,j,k);
                            dz_arith = dz(i,j,k);
                        }
                    }
                    if (OBC->v_S_OBCs_on_PE) {
                        if (OBC->segnum_v(i,j,0) < 0) {
                            h_harm = h(i,j+1,k);
                            h_arith = h(i,j+1,k);
                            h_delta = 0.0;
                            dz_harm(i,j,k) = dz(i,j+1,k);
                            dz_arith = dz(i,j+1,k);
                        }
                    }
                }
                */

                if (harmonic_visc) {
                    /// The following block calculates the thicknesses at velocity grid points
                    /// for the vertical viscosity (hvel and dz_vel).  Near the bottom an
                    /// upwind biased thickness is used to control the effect of spurious
                    /// Montgomery potential gradients at the bottom where nearly massless
                    /// layers ride over the topography.

                    hvel(i,j,k) = h_harm;
                    dz_vel(i,j,k) = dz_harm(i,j,k);

                    if (v(i,j,k) * h_delta < 0) {
                        z2 = z_i(i,j,k+1);
                        botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

                        hvel(i,j,k) = (1.0 - botfn) * h_harm + botfn * h_arith;
                        dz_vel(i,j,k) = (1.0 - botfn) * dz_harm(i,j,k) + botfn * dz_arith;
                    }

                    z_i(i,j,k) = z_i(i,j,k+1) + dz_harm(i,j,k) * I_Hbbl;
                } else {
                    zcol = zcol + dz(i,j,k);
                    zcol_p1 = zcol_p1 + dz(i,j+1,k);

                    zh = zh + dz_harm(i,j,k);

                    z_clear = amrex::max(zcol, zcol_p1) + Dmin;
                    if (zi_dir < 0) z_clear = zcol + Dmin;
                    if (zi_dir > 0) z_clear = zcol_p1 + Dmin;

                    z_i(i,j,k) = amrex::max(zh, z_clear) * I_Hbbl;

                    hvel(i,j,k) = h_arith;
                    dz_vel(i,j,k) = dz_arith;

                    if (v(i,j,k) * h_delta > 0) {
                        if (zh * I_Hbbl < harm_BL_val) {
                            hvel(i,j,k) = h_harm;
                            dz_vel(i,j,k) = dz_harm(i,j,k);
                        } else {
                            z2_wt = 1.0;
                            if (zh * I_Hbbl < 2.0 * harm_BL_val)
                                z2_wt = amrex::max(0.0, amrex::min(1.0, zh * I_Hbbl * I_valBL - 1.0));

                            /// TODO: should z_clear be used here?
                            z2 = z2_wt * (amrex::max(zh, amrex::max(zcol, zcol_p1) + Dmin) * I_Hbbl);
                            botfn = 1.0 / (1.0 + 0.09 * z2 * z2 * z2 * z2 * z2 * z2);

                            hvel(i,j,k) = (1.0 - botfn) * h_arith + botfn * h_harm;
                            dz_vel(i,j,k) = (1.0 - botfn) * dz_arith + botfn * dz_harm(i,j,k);
                        }
                    }
                }

                if (use_GL90_in_SSW) {
                    /// The following block calculates the normalized height above the GL90 BBL
                    /// (z_i_gl90), using a harmonic mean between layer thicknesses. For the
                    /// GL90 BBL we use simply a constant (Hbbl_gl90). The purpose is that the
                    /// GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
                    /// that no momentum gets fluxed into vanished layers. The scheme is not
                    /// sensitive to the exact value of Hbbl_gl90, as long as it is in a
                    /// reasonable range (~1-20 m): large enough to capture vanished layers
                    /// over topography, small enough to not contaminate the interior.

                    z_i_gl90(i,j,k) = z_i_gl90(i,j,k+1) + dz_harm(i,j,k) * I_Hbbl_gl90;
                }
            }

            find_coupling_coef_point<false, false>(i, j, a_cpl, dz_vel, dz_harm, bbl_thick, kv_bbl, z_i,
                h_ml, dt, ke, Boussinesq, nkml, dZ_subroundoff, H_to_m, RZ_to_H, Z_to_H, CoriolisBu,
                cpl, Kv_shear, Kv_shear_Bu, Kv_tbl_shelf_u, Kv_tbl_shelf_v, tbl_thick_shelf_u,
                tbl_thick_shelf_v, nkml_visc_u, nkml_visc_v, Ustar_2d, SpV_avg);

            if (hML_v.p != nullptr) hML_v(i,j,0) = h_ml;

            if (use_GL90_in_SSW) {
                find_coupling_coef_gl90_point<false>(i, j, a_cpl_gl90, dz_vel, z_i_gl90, ke,
                    dZ_subroundoff, g_prime, CoriolisBu, use_GL90_N2, kappa_gl90, read_kappa_gl90,
                    alpha_gl90, kappa_gl90_2d, kdgl90_struct);
            }

            do_any_shelf = false;
            if (frac_shelf_v.p != nullptr) {
                a1_shelf_v(i,j,0) = 0.0;
                do_any_shelf = frac_shelf_v(i,j,0) > 0.0;

                if (do_any_shelf) {
                    /// Initialize non-harmonic depths
                    if (!harmonic_visc) {
                        zh = 0.0;
                        Ztop_min = amrex::min(zcol, zcol_p1);
                        I_HTbl = 1.0 / (tbl_thick_shelf_v(i,j,0) + dz_neglect);
                    }

                    for (int k = 0; k < nz; ++k) {
                        if (harmonic_visc) {
                            hvel_shelf(i,j,k) = hvel(i,j,k);
                            dz_vel_shelf(i,j,k) = dz_vel(i,j,k);
                        } else {
                            /// Find upwind-biased thickness near the surface.
                            /// Perhaps this needs to be done more carefully, via find_eta.
                            h_harm = 2.0 * h(i,j,k) * h(i,j+1,k)
                                / (h(i,j,k) + h(i,j+1,k) + h_neglect);
                            h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k));
                            h_delta = h(i,j+1,k) - h(i,j,k);
                            dz_arith = 0.5 * (dz(i,j+1,k) + dz(i,j,k));

                            /// Project thickness outward across OBCs using a zero-gradient condition.
                            /*
                            if (OBC != nullptr) {
                                if (OBC->v_N_OBCs_on_PE) {
                                    if (OBC->segnum_v(i,j,0) > 0) {
                                        h_harm = h(i,j,k);
                                        h_arith = h(i,j,k);
                                        h_delta = 0.0;
                                        dz_arith = dz(i,j,k);
                                    }
                                }
                                if (OBC->v_S_OBCs_on_PE) {
                                    if (OBC->segnum_v(i,j,0) < 0) {
                                        h_harm = h(i,j+1,k);
                                        h_arith = h(i,j+1,k);
                                        h_delta = 0.0;
                                        dz_arith = dz(i,j+1,k);
                                    }
                                }
                            }
                            */

                            zcol = zcol - dz(i,j,k);
                            zcol_p1 = zcol_p1 - dz(i,j+1,k);

                            zh = zh + dz_harm(i,j,k);

                            hvel_shelf(i,j,k) = hvel(i,j,k);
                            dz_vel_shelf(i,j,k) = dz_vel(i,j,k);

                            if (v(i,j,k) * h_delta > 0.0) {
                                if (zh * I_HTbl < harm_BL_val) {
                                    hvel_shelf(i,j,k) = amrex::min(hvel(i,j,k), h_harm);
                                    dz_vel_shelf(i,j,k) = amrex::min(dz_vel(i,j,k), dz_harm(i,j,k));
                                } else {
                                    z2_wt = 1.0;
                                    if (zh * I_HTbl < 2.0 * harm_BL_val) {
                                        z2_wt = amrex::max(0.0, amrex::min(1.0, zh * I_HTbl * I_valBL - 1.0));
                                    }

                                    z2 = z2_wt * (amrex::max(zh, Ztop_min - amrex::min(zcol, zcol_p1)) * I_HTbl);
                                    /// TODO: Replace **6
                                    Real const z2_3 = (z2 * z2) * z2;
                                    topfn = 1.0 / (1.0 + 0.09 * (z2_3 * z2_3));

                                    hvel_shelf(i,j,k) = amrex::min(hvel(i,j,k), (1.0 - topfn) * h_arith + topfn * h_harm);
                                    dz_vel_shelf(i,j,k) = amrex::min(dz_vel(i,j,k), (1.0 - topfn) * dz_arith + topfn * dz_harm(i,j,k));
                                }
                            }
                        }
                    }

                    find_coupling_coef_point<false, true>(i, j, a_shelf, dz_vel_shelf, dz_harm, bbl_thick,
                        kv_bbl, z_i, h_ml, dt, ke, Boussinesq, nkml, dZ_subroundoff, H_to_m, RZ_to_H,
                        Z_to_H, CoriolisBu, cpl, Kv_shear, Kv_shear_Bu, Kv_tbl_shelf_u, Kv_tbl_shelf_v,
                        tbl_thick_shelf_u, tbl_thick_shelf_v, nkml_visc_u, nkml_visc_v, Ustar_2d, SpV_avg);

                    a1_shelf_v(i,j,0) = a_shelf(i,j,0);
                }
            }

            if (do_any_shelf) {
                if (use_GL90_in_SSW) {
                    for (int k = 0; k <= nz; ++k) {
                        a_v(i,j,k) = amrex::min(a_cpl_max, (frac_shelf_v(i,j,0) * a_shelf(i,j,k) +
                                                            (1.0 - frac_shelf_v(i,j,0)) * a_cpl(i,j,k)) + a_cpl_gl90(i,j,k));

                        /// This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
                        /// CS%a_v(I,j,K) = min(a_cpl_max, forces%frac_shelf_v(I,j) * max(a_shelf(K), a_cpl(K)) + &
                        ///                                (1. - forces%frac_shelf_v(I,j)) * a_cpl(K))

                        a_v_gl90(i,j,k) = amrex::min(a_cpl_max, a_cpl_gl90(i,j,k));
                    }
                } else {
                    for (int k = 0; k <= nz; ++k) {
                        a_v(i,j,k) = amrex::min(a_cpl_max, (frac_shelf_v(i,j,0) * a_shelf(i,j,k) +
                                                            (1.0 - frac_shelf_v(i,j,0)) * a_cpl(i,j,k)));
                        /// This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
                        /// CS%a_v(I,j,K) = min(a_cpl_max, forces%frac_shelf_v(I,j) * max(a_shelf(K), a_cpl(K)) + &
                        ///                                (1. - forces%frac_shelf_v(I,j)) * a_cpl(K))
                    }
                }

                for (int k = 0; k < nz; ++k) {
                    /// Should we instead take the inverse of the average of the inverses?
                    h_v(i,j,k) = frac_shelf_v(i,j,0) * hvel_shelf(i,j,k)
                        + (1.0 - frac_shelf_v(i,j,0)) * hvel(i,j,k) + h_neglect;
                }
            } else {
                if (use_GL90_in_SSW) {
                    for (int k = 0; k <= nz; ++k) {
                        a_cpl(i,j,k) = a_cpl(i,j,k) + a_cpl_gl90(i,j,k);
                    }

                    for (int k = 0; k <= nz; ++k) {
                        a_v_gl90(i,j,k) = amrex::min(a_cpl_max, a_cpl_gl90(i,j,k));
                    }
                }

                for (int k = 0; k <= nz; ++k) {
                    a_v(i,j,k) = amrex::min(a_cpl_max, a_cpl(i,j,k));
                }

                for (int k = 0; k < nz; ++k) {
                    h_v(i,j,k) = hvel(i,j,k) + h_neglect;
                }
            }

            /// Diagnose total Kv at v-points
            if (Kv_v.p != nullptr) {
                for (int k = 0; k < nz; ++k) {
                    Kv_v(i,j,k) = 0.5 * (a_v(i,j,k) + a_v(i,j,k+1)) * h_v(i,j,k);
                }
            }

            /// Diagnose GL90 Kv at v-points
            if (Kv_gl90_v.p != nullptr) {
                for (int k = 0; k < nz; ++k) {
                    Kv_gl90_v(i,j,k) = 0.5 * (a_v_gl90(i,j,k) + a_v_gl90(i,j,k+1)) * h_v(i,j,k);
                }
            }
        });
    }
}

///
/// @brief Forward elimination and back substitution of the tridiagonal system for one velocity
///        component in each column of @p col.
///
/// @tparam has_Ray  If true, @p Ray is allocated.
/// @tparam has_str  If true, @p d_dt_str is allocated.
/// @param col             The velocity-point box of the whole column
/// @param vel             The velocity component [L T-1 ~> m s-1]
/// @param h_vel           The effective layer thickness at the velocity points [H ~> m or kg m-2]
/// @param a_vel           The drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
/// @param Ray             The Rayleigh drag velocity [H T-1 ~> m s-1 or Pa s m-1]
/// @param d_dt_str        The acceleration due to the surface stress [L T-2 ~> m s-2]
/// @param surface_stress  The surface stress times the time step, divided by the density
///                        [H L T-1 ~> m2 s-1 or kg m-1 s-1]
/// @param c1              A variable used by the tridiagonal solver [nondim]
/// @param mask            0 for boundary points and 1 for ocean points [nondim]
/// @param nz              The number of layers/levels in the vertical
/// @param dt              Time increment [T ~> s]
/// @param Idt             The inverse of the time step [T-1 ~> s-1]
/// @param accel_underflow An acceleration magnitude that is so small that values that are less
///                        than this are diagnosed as 0 [L T-2 ~> m s-2]
///
template <bool has_Ray, bool has_str>
void vertvisc_tridiag(const Box& col,
                      Array4<Real>       const& vel,
                      Array4<const Real> const& h_vel,
                      Array4<const Real> const& a_vel,
                      Array4<const Real> const& Ray,
                      Array4<Real>       const& d_dt_str,
                      Array4<const Real> const& surface_stress,
                      Array4<Real>       const& c1,
                      Array4<const Real> const& mask,
                      const int nz,
                      const Real dt,
                      const Real Idt,
                      const Real accel_underflow)
{
    ParallelFor(col, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        if (!(mask(i,j,0) > 0.0)) return;

        Real b1;         ///< A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
        Real d1;         ///< d1=1-c1 is used by the tridiagonal solver [nondim].
        Real Ray_k;      ///< Ray is the Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
        Real b_denom_1;  ///< The first term in the denominator of b1 [H ~> m or kg m-2].

        Ray_k = 0.0;
        if constexpr (has_Ray) Ray_k = Ray(i,j,0);

        b_denom_1 = h_vel(i,j,0) + dt * (Ray_k + a_vel(i,j,0));
        b1 = 1.0 / (b_denom_1 + dt * a_vel(i,j,1));
        d1 = b_denom_1 * b1;
        vel(i,j,0) = b1 * (h_vel(i,j,0) * vel(i,j,0) + surface_stress(i,j,0));

        if constexpr (has_str) {
            d_dt_str(i,j,0) = b1 * (h_vel(i,j,0) * d_dt_str(i,j,0) + surface_stress(i,j,0) * Idt);
        }

        for (int k = 1; k < nz; ++k) {
            if constexpr (has_Ray) Ray_k = Ray(i,j,k);

            c1(i,j,k) = dt * a_vel(i,j,k) * b1;
            b_denom_1 = h_vel(i,j,k) + dt * (Ray_k + a_vel(i,j,k) * d1);
            b1 = 1.0 / (b_denom_1 + dt * a_vel(i,j,k+1));
            d1 = b_denom_1 * b1;
            vel(i,j,k) = (h_vel(i,j,k) * vel(i,j,k) + dt * a_vel(i,j,k) * vel(i,j,k-1)) * b1;

            if constexpr (has_str) {
                d_dt_str(i,j,k) = (h_vel(i,j,k) * d_dt_str(i,j,k)
                    + dt * a_vel(i,j,k) * d_dt_str(i,j,k-1)) * b1;
            }

            /// Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
            /// XXX: Check GPU behavior
            if (dt < 0) break;
        }

        if constexpr (has_str) {
            if (Math::abs(d_dt_str(i,j,nz-1)) < accel_underflow)
                d_dt_str(i,j,nz-1) = 0.0;
        }

        for (int k = nz-2; k >= 0; --k) {
            vel(i,j,k) = vel(i,j,k) + c1(i,j,k+1) * vel(i,j,k+1);

            if constexpr (has_str) {
                d_dt_str(i,j,k) = d_dt_str(i,j,k) + c1(i,j,k+1) * d_dt_str(i,j,k+1);

                if (Math::abs(d_dt_str(i,j,k)) < accel_underflow)
                    d_dt_str(i,j,k) = 0.0;
            }
        }
    });
}

void vertvisc_tr(const Box& bxU0,
                 const Box& bxV0,
                 Array4<Real>       const& u,
                 Array4<Real>       const& v,
                 Array4<const Real> const& h,
                 Array4<const Real> const& taux,
                 Array4<const Real> const& tauy,
                 Array4<const Real> const& Ray_u,
                 Array4<const Real> const& Ray_v,
                 Array4<Real>       const& taux_shelf,
                 Array4<Real>       const& tauy_shelf,
                 const Real dt,
                 OceanOBC* OBC,
                 Array4<Real>       const& du_dt_visc,
                 Array4<Real>       const& dv_dt_visc,
                 Array4<Real>       const& du_dt_visc_gl90,
                 Array4<Real>       const& dv_dt_visc_gl90,
                 Array4<Real>       const& du_dt_str,
                 Array4<Real>       const& dv_dt_str,
                 Array4<const Real> const& mask2dCu,
                 Array4<const Real> const& mask2dCv,
                 Array4<const Real> const& areaCu,
                 Array4<const Real> const& areaCv,
                 Array4<const Real> const& areaT,
                 Array4<const Real> const& IareaT,
                 Array4<const Real> const& dx_Cv,
                 Array4<const Real> const& dy_Cu,
                 const int ke,
                 const Real H_subroundoff,
                 const Real Angstrom_H,
                 const Real H_to_RZ,
                 const Real m_s_to_L_T,
                 const vertvisc_opts_CS_C& CS,
                 Array4<const Real> const& a_u,
                 Array4<const Real> const& a_u_gl90,
                 Array4<const Real> const& h_u,
                 Array4<const Real> const& a_v,
                 Array4<const Real> const& a_v_gl90,
                 Array4<const Real> const& h_v,
                 Array4<const Real> const& a1_shelf_u,
                 Array4<const Real> const& a1_shelf_v,
                 Array4<Real>       const& taux_bot,
                 Array4<Real>       const& tauy_bot,
                 Array4<Real>       const& KE_u,
                 Array4<Real>       const& KE_v,
                 Array4<Real>       const& u_old,
                 Array4<Real>       const& v_old,
                 Array4<Real>       const& vel_report_u,
                 Array4<Real>       const& vel_report_v,
                 Array4<int>        const& dowrite_u,
                 Array4<int>        const& dowrite_v,
                 int& ntrunc_sum,
                 const bool fpmix,
                 Array4<const Real> const& Us_x,
                 Array4<const Real> const& Us_y)
{
    BL_PROFILE("vertvisc_tr");

    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    const int nz = ke;

    Real Hmix = 0.0;    ///< The mixed layer thickness over which stress
                        ///< is applied with direct_stress [H ~> m or kg m-2].
    Real I_Hmix = 0.0;  ///< The inverse of Hmix [H-1 ~> m-1 or m2 kg-1].

    if (CS.direct_stress) {
        Hmix = CS.Hmix_stress;
        I_Hmix = 1.0 / Hmix;
    }
    /// The time step divided by the mean density [T H Z-1 R-1 ~> s m3 kg-1 or s].
    const Real dt_Rho0 = dt / H_to_RZ;
    /// A thickness that is so small it is usually lost in roundoff and can be neglected
    /// [H ~> m or kg m-2].
    const Real h_neglect = H_subroundoff;
    /// The inverse of the time step [T-1 ~> s-1].
    const Real Idt = 1.0 / dt;

    /// An acceleration magnitude that is so small that values that are less than this are
    /// diagnosed as 0 [L T-2 ~> m s-2].
    const Real accel_underflow = CS.vel_underflow * Idt;

    const bool DoStokesMixing = (Us_x.p != nullptr);
    const bool lfpmix = fpmix;

    Box colU = bxU0;
    colU.setSmall(2, 0);
    colU.setBig(2, 0);
    Box colV = bxV0;
    colV.setSmall(2, 0);
    colV.setBig(2, 0);
    Box c1U_bx = bxU0;
    c1U_bx.setSmall(2, 0);
    c1U_bx.setBig(2, nz-1);
    Box c1V_bx = bxV0;
    c1V_bx.setSmall(2, 0);
    c1V_bx.setBig(2, nz-1);

    /// The same as stress, unless the wind stress is applied as a body force
    /// [H L T-1 ~> m2 s-1 or kg m-1 s-1].
    FArrayBox surface_stress_u_fab(colU, 1, The_Async_Arena());
    FArrayBox surface_stress_v_fab(colV, 1, The_Async_Arena());
    auto const& surface_stress_u = surface_stress_u_fab.array();
    auto const& surface_stress_v = surface_stress_v_fab.array();
    /// A variable used by the tridiagonal solver [nondim].
    FArrayBox c1U_fab(c1U_bx, 1, The_Async_Arena());
    FArrayBox c1V_fab(c1V_bx, 1, The_Async_Arena());
    auto const& c1U = c1U_fab.array();
    auto const& c1V = c1V_fab.array();

    ///   Update the zonal velocity component using a modification of a standard
    /// tridiagonal solver.

    /// WGL: Brandon Reichl says the following is obsolete. u(I,j,k) already
    /// includes Stokes.
    /// When mixing down Eulerian current + Stokes drift add before calling solver
    if (DoStokesMixing) {
        ParallelFor(bxU0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCu(i,j,0) > 0.0) u(i,j,k) = u(i,j,k) + Us_x(i,j,k);
        });
    }

    if (lfpmix) {
        ParallelFor(bxU0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCu(i,j,0) > 0.0) u(i,j,k) = u(i,j,k) - Us_x(i,j,k);
        });
    }

    if (du_dt_visc.p != nullptr) {
        ParallelFor(bxU0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            du_dt_visc(i,j,k) = u(i,j,k);
        });
    }

    if (du_dt_visc_gl90.p != nullptr) {
        ParallelFor(bxU0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            du_dt_visc_gl90(i,j,k) = u(i,j,k);
        });
    }

    if (du_dt_str.p != nullptr) {
        ParallelFor(bxU0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            du_dt_str(i,j,k) = 0.0;
        });
    }

    ///   One option is to have the wind stress applied as a body force
    /// over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
    /// the wind stress is applied as a stress boundary condition.
    if (CS.direct_stress) {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCu(i,j,0) > 0.0)) return;
            Real zDS, h_avg;  ///< Temporary thickness variables used with direct_stress
                              ///< [H ~> m or kg m-2]
            Real hfr;         ///< Temporary ratio of thicknesses used with direct_stress [nondim]
            Real stress;      ///<   The surface stress times the time step, divided
                              ///< by the density [H L T-1 ~> m2 s-1 or kg m-1 s-1].
            surface_stress_u(i,j,0) = 0.0;
            zDS = 0.0;
            stress = dt_Rho0 * taux(i,j,0);
            for (int k = 0; k < nz; ++k) {
                h_avg = 0.5 * (h(i,j,k) + h(i+1,j,k)) + h_neglect;
                hfr = 1.0; if ((zDS+h_avg) > Hmix) hfr = (Hmix - zDS) / h_avg;
                u(i,j,k) = u(i,j,k) + I_Hmix * hfr * stress;
                if (du_dt_str.p != nullptr) du_dt_str(i,j,k) = (I_Hmix * hfr * stress) * Idt;
                zDS = zDS + h_avg; if (zDS >= Hmix) break;
            }
        });
    } else {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            surface_stress_u(i,j,0) = dt_Rho0 * (mask2dCu(i,j,0)*taux(i,j,0));
        });
    }

    /// perform forward elimination on the tridiagonal system
    ///
    /// denote the diagonal of the system as b_k, the subdiagonal as a_k
    /// and the superdiagonal as c_k. The right-hand side terms are d_k.
    ///
    /// ignoring the Rayleigh drag contribution,
    /// we have a_k = -dt * a_u(k)
    ///         b_k = h_u(k) + dt * (a_u(k) + a_u(k+1))
    ///         c_k = -dt * a_u(k+1)
    ///
    /// for forward elimination, we want to:
    /// calculate c'_k = - c_k                / (b_k + a_k c'_(k-1))
    /// and       d'_k = (d_k - a_k d'_(k-1)) / (b_k + a_k c'_(k-1))
    /// where c'_1 = c_1/b_1 and d'_1 = d_1/b_1
    ///
    /// This form is mathematically equivalent to Thomas' tridiagonal matrix algorithm, but it
    /// does not suffer from the acute sensitivity to truncation errors of the Thomas algorithm
    /// because it involves no subtraction, as discussed by Schopf & Loughe, MWR, 1995.
    ///
    /// b1 is the denominator term 1 / (b_k + a_k c'_(k-1))
    /// b_denom_1 is (b_k + a_k + c_k) - a_k(1 - c'_(k-1))
    ///            = (b_k + c_k + c'_(k-1))
    /// this is done so that d1 = b1 * b_denom_1 = 1 - c'_(k-1)
    /// c1(k) is -c'_(k - 1)
    /// and the right-hand-side is destructively updated to be d'_k

    if (Ray_u.p != nullptr) {
        if (du_dt_str.p != nullptr) {
            vertvisc_tridiag<true, true>(colU, u, h_u, a_u, Ray_u, du_dt_str, surface_stress_u, c1U,
                                         mask2dCu, nz, dt, Idt, accel_underflow);
        } else {
            vertvisc_tridiag<true, false>(colU, u, h_u, a_u, Ray_u, du_dt_str, surface_stress_u, c1U,
                                          mask2dCu, nz, dt, Idt, accel_underflow);
        }
    } else {
        if (du_dt_str.p != nullptr) {
            vertvisc_tridiag<false, true>(colU, u, h_u, a_u, Ray_u, du_dt_str, surface_stress_u, c1U,
                                          mask2dCu, nz, dt, Idt, accel_underflow);
        } else {
            vertvisc_tridiag<false, false>(colU, u, h_u, a_u, Ray_u, du_dt_str, surface_stress_u, c1U,
                                           mask2dCu, nz, dt, Idt, accel_underflow);
        }
    }

    /// compute vertical velocity tendency that arises from GL90 viscosity;
    /// follow tridiagonal solve method as above; to avoid corrupting u,
    /// use ADp%du_dt_visc_gl90 as a placeholder for updated u (due to GL90) until last do loop
    if ((CS.do_du_dt_visc_gl90) || (CS.do_GLwork)) {
        if (du_dt_visc_gl90.p != nullptr) {
            const bool do_GLwork = CS.do_GLwork;
            ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                if (!(mask2dCu(i,j,0) > 0.0)) return;
                Real b1, d1, b_denom_1;

                b_denom_1 = h_u(i,j,0);  /// CS%a_u_gl90(I,j,1) is zero
                b1 = 1.0 / (b_denom_1 + dt * a_u_gl90(i,j,1));
                d1 = b_denom_1 * b1;

                du_dt_visc_gl90(i,j,0) = b1 * (h_u(i,j,0) * du_dt_visc_gl90(i,j,0));

                for (int k = 1; k < nz; ++k) {
                    c1U(i,j,k) = dt * a_u_gl90(i,j,k) * b1;
                    b_denom_1 = h_u(i,j,k) + dt * (a_u_gl90(i,j,k)*d1);
                    b1 = 1.0 / (b_denom_1 + dt * a_u_gl90(i,j,k+1));
                    d1 = b_denom_1 * b1;

                    du_dt_visc_gl90(i,j,k) = (h_u(i,j,k) * du_dt_visc_gl90(i,j,k)
                        + dt * a_u_gl90(i,j,k) * du_dt_visc_gl90(i,j,k-1)) * b1;
                }

                /// back substitute to solve for new velocities, held by ADp%du_dt_visc_gl90
                for (int k = nz-2; k >= 0; --k) {
                    du_dt_visc_gl90(i,j,k) =
                        du_dt_visc_gl90(i,j,k) + c1U(i,j,k+1) * du_dt_visc_gl90(i,j,k+1);
                }

                for (int k = 0; k < nz; ++k) {
                    /// now fill ADp%du_dt_visc_gl90(I,j,k) with actual velocity tendency due to GL90;
                    /// note that on RHS: ADp%du_dt_visc(I,j,k) holds the original velocity value u(I,j,k)
                    /// and ADp%du_dt_visc_gl90(I,j,k) the updated velocity due to GL90
                    du_dt_visc_gl90(i,j,k) =
                        (du_dt_visc_gl90(i,j,k) - du_dt_visc(i,j,k)) * Idt;

                    if (Math::abs(du_dt_visc_gl90(i,j,k)) < accel_underflow) {
                        du_dt_visc_gl90(i,j,k) = 0.0;
                    }
                }

                /// to compute energetics, we need to multiply by u*h, where u is original velocity before
                /// velocity update; note that ADp%du_dt_visc(I,j,k) holds the original velocity value u(I,j,k)
                if (do_GLwork) {
                    for (int k = 0; k < nz; ++k) {
                        KE_u(i,j,k) = du_dt_visc(i,j,k) * h_u(i,j,k) * areaCu(i,j,0) *
                                      du_dt_visc_gl90(i,j,k);
                    }
                }
            });
        }
    }

    if (du_dt_visc.p != nullptr) {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            for (int k = 0; k < nz; ++k) {
                du_dt_visc(i,j,k) = (u(i,j,k) - du_dt_visc(i,j,k)) * Idt;

                if (Math::abs(du_dt_visc(i,j,k)) < accel_underflow)
                    du_dt_visc(i,j,k) = 0.0;
            }
        });
    }

    if (taux_shelf.p != nullptr) {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            taux_shelf(i,j,0) = -H_to_RZ * a1_shelf_u(i,j,0) * u(i,j,0); /// - u_shelf?
        });
    }

    if (taux_bot.p != nullptr) {
        ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            taux_bot(i,j,0) = H_to_RZ * (u(i,j,nz-1) * a_u(i,j,nz));
        });

        if (Ray_u.p != nullptr) {
            ParallelFor(colU, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int k = 0; k < nz; ++k) {
                    taux_bot(i,j,0) = taux_bot(i,j,0) + H_to_RZ * (Ray_u(i,j,k) * u(i,j,k));
                }
            });
        }
    }

    /// When mixing down Eulerian current + Stokes drift subtract after calling solver
    if (DoStokesMixing) {
        ParallelFor(bxU0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCu(i,j,0) > 0.0) u(i,j,k) = u(i,j,k) - Us_x(i,j,k);
        });
    }

    if (lfpmix) {
        ParallelFor(bxU0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCu(i,j,0) > 0.0) u(i,j,k) = u(i,j,k) + Us_x(i,j,k);
        });
    }

    /// == Now work on the meridional velocity component.

    /// When mixing down Eulerian current + Stokes drift add before calling solver
    if (DoStokesMixing) {
        ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCv(i,j,0) > 0.0) v(i,j,k) = v(i,j,k) + Us_y(i,j,k);
        });
    }

    if (lfpmix) {
        ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCv(i,j,0) > 0.0) v(i,j,k) = v(i,j,k) - Us_y(i,j,k);
        });
    }

    if (dv_dt_visc.p != nullptr) {
        ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dv_dt_visc(i,j,k) = v(i,j,k);
        });
    }

    if (dv_dt_visc_gl90.p != nullptr) {
        ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dv_dt_visc_gl90(i,j,k) = v(i,j,k);
        });
    }

    if (dv_dt_str.p != nullptr) {
        ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dv_dt_str(i,j,k) = 0.0;
        });
    }

    ///   One option is to have the wind stress applied as a body force
    /// over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
    /// the wind stress is applied as a stress boundary condition.
    if (CS.direct_stress) {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (!(mask2dCv(i,j,0) > 0.0)) return;
            Real zDS, h_avg;  ///< Temporary thickness variables used with direct_stress
                              ///< [H ~> m or kg m-2]
            Real hfr;         ///< Temporary ratio of thicknesses used with direct_stress [nondim]
            Real stress;      ///<   The surface stress times the time step, divided
                              ///< by the density [H L T-1 ~> m2 s-1 or kg m-1 s-1].
            surface_stress_v(i,j,0) = 0.0;
            zDS = 0.0;
            stress = dt_Rho0 * tauy(i,j,0);
            for (int k = 0; k < nz; ++k) {
                h_avg = 0.5 * (h(i,j,k) + h(i,j+1,k)) + h_neglect;
                hfr = 1.0; if ((zDS+h_avg) > Hmix) hfr = (Hmix - zDS) / h_avg;
                v(i,j,k) = v(i,j,k) + I_Hmix * hfr * stress;
                if (dv_dt_str.p != nullptr) dv_dt_str(i,j,k) = (I_Hmix * hfr * stress) * Idt;
                zDS = zDS + h_avg; if (zDS >= Hmix) break;
            }
        });
    } else {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            surface_stress_v(i,j,0) = dt_Rho0 * (mask2dCv(i,j,0) * tauy(i,j,0));
        });
    }

    if (Ray_v.p != nullptr) {
        if (dv_dt_str.p != nullptr) {
            vertvisc_tridiag<true, true>(colV, v, h_v, a_v, Ray_v, dv_dt_str, surface_stress_v, c1V,
                                         mask2dCv, nz, dt, Idt, accel_underflow);
        } else {
            vertvisc_tridiag<true, false>(colV, v, h_v, a_v, Ray_v, dv_dt_str, surface_stress_v, c1V,
                                          mask2dCv, nz, dt, Idt, accel_underflow);
        }
    } else {
        if (dv_dt_str.p != nullptr) {
            vertvisc_tridiag<false, true>(colV, v, h_v, a_v, Ray_v, dv_dt_str, surface_stress_v, c1V,
                                          mask2dCv, nz, dt, Idt, accel_underflow);
        } else {
            vertvisc_tridiag<false, false>(colV, v, h_v, a_v, Ray_v, dv_dt_str, surface_stress_v, c1V,
                                           mask2dCv, nz, dt, Idt, accel_underflow);
        }
    }

    /// compute vertical velocity tendency that arises from GL90 viscosity;
    /// follow tridiagonal solve method as above; to avoid corrupting v,
    /// use ADp%dv_dt_visc_gl90 as a placeholder for updated v (due to GL90) until last do loop
    if ((CS.do_dv_dt_visc_gl90) || (CS.do_GLwork)) {
        if (dv_dt_visc_gl90.p != nullptr) {
            ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                if (!(mask2dCv(i,j,0) > 0.0)) return;
                Real b1, d1, b_denom_1;

                b_denom_1 = h_v(i,j,0);  /// CS%a_v_gl90(i,J,1) is zero
                b1 = 1.0 / (b_denom_1 + dt*a_v_gl90(i,j,1));
                d1 = b_denom_1 * b1;
                dv_dt_visc_gl90(i,j,0) = b1 * (h_v(i,j,0) * dv_dt_visc_gl90(i,j,0));

                for (int k = 1; k < nz; ++k) {
                    c1V(i,j,k) = dt * a_v_gl90(i,j,k) * b1;
                    b_denom_1 = h_v(i,j,k) + dt * (a_v_gl90(i,j,k) * d1);
                    b1 = 1.0 / (b_denom_1 + dt * a_v_gl90(i,j,k+1));
                    d1 = b_denom_1 * b1;
                    dv_dt_visc_gl90(i,j,k) = (h_v(i,j,k) * dv_dt_visc_gl90(i,j,k)
                        + dt * a_v_gl90(i,j,k) * dv_dt_visc_gl90(i,j,k-1)) * b1;
                }

                /// back substitute to solve for new velocities, held by ADp%dv_dt_visc_gl90
                for (int k = nz-2; k >= 0; --k) {
                    dv_dt_visc_gl90(i,j,k) = dv_dt_visc_gl90(i,j,k) + c1V(i,j,k+1) * dv_dt_visc_gl90(i,j,k+1);
                }
            });

            ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (!(mask2dCv(i,j,0) > 0.0)) return;
                /// now fill ADp%dv_dt_visc_gl90(i,J,k) with actual velocity tendency due to GL90;
                /// note that on RHS: ADp%dv_dt_visc(i,J,k) holds the original velocity value v(i,J,k)
                /// and ADp%dv_dt_visc_gl90(i,J,k) the updated velocity due to GL90
                dv_dt_visc_gl90(i,j,k) = (dv_dt_visc_gl90(i,j,k) - dv_dt_visc(i,j,k)) * Idt;

                if (Math::abs(dv_dt_visc_gl90(i,j,k)) < accel_underflow)
                    dv_dt_visc_gl90(i,j,k) = 0.0;
            });

            /// to compute energetics, we need to multiply by v*h, where u is original velocity before
            /// velocity update; note that ADp%dv_dt_visc(I,j,k) holds the original velocity value v(i,J,k)
            if (CS.do_GLwork) {
                ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (!(mask2dCv(i,j,0) > 0.0)) return;
                    /// note that on RHS: ADp%dv_dt_visc(I,j,k) holds the original velocity value v(I,j,k)
                    KE_v(i,j,k) = dv_dt_visc(i,j,k) * h_v(i,j,k) * areaCv(i,j,0) *
                                  dv_dt_visc_gl90(i,j,k);
                });
            }
        }
    }

    if (dv_dt_visc.p != nullptr) {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            for (int k = 0; k < nz; ++k) {
                dv_dt_visc(i,j,k) = (v(i,j,k) - dv_dt_visc(i,j,k))*Idt;
                if (Math::abs(dv_dt_visc(i,j,k)) < accel_underflow) dv_dt_visc(i,j,k) = 0.0;
            }
        });
    }

    if (tauy_shelf.p != nullptr) {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            tauy_shelf(i,j,0) = -H_to_RZ * a1_shelf_v(i,j,0) * v(i,j,0); /// - v_shelf?
        });
    }

    /// JORGE TODO: this has to be malloced
    if (tauy_bot.p != nullptr) {
        ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            tauy_bot(i,j,0) = H_to_RZ * (v(i,j,nz-1) * a_v(i,j,nz));
        });

        if (Ray_v.p != nullptr) {
            ParallelFor(colV, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int k = 0; k < nz; ++k) {
                    tauy_bot(i,j,0) = tauy_bot(i,j,0) + H_to_RZ * (Ray_v(i,j,k)*v(i,j,k));
                }
            });
        }
    }

    /// When mixing down Eulerian current + Stokes drift subtract after calling solver
    if (DoStokesMixing) {
        ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCv(i,j,0) > 0.0) v(i,j,k) = v(i,j,k) - Us_y(i,j,k);
        });
    }

    if (lfpmix) {
        ParallelFor(bxV0, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask2dCv(i,j,0) > 0.0) v(i,j,k) = v(i,j,k) + Us_y(i,j,k);
        });
    }

    vertvisc_limit_vel(bxU0, bxV0, u, v, dt, areaT, IareaT, dx_Cv, dy_Cu, ke, Angstrom_H,
                       m_s_to_L_T, CS, h_u, h_v, u_old, v_old, vel_report_u, vel_report_v,
                       dowrite_u, dowrite_v, ntrunc_sum);

    /// Here the velocities associated with open boundary conditions are applied.
    /*
    if (OBC != nullptr) {
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto& segment = OBC->segment[n];
            if (segment.specified) {
                if (segment.is_N_or_S) {
                    int J = segment.HI.JsdB;
                    ParallelFor(Box(IntVect(segment.HI.isd, J, 0),
                                    IntVect(segment.HI.ied, J, nz-1)),
                    [=] AMREX_GPU_DEVICE (int i, int jj, int k)
                    {
                        v(i,J,k) = segment.normal_vel(i,J,k);
                    });
                } else if (segment.is_E_or_W) {
                    int I = segment.HI.IsdB;
                    ParallelFor(Box(IntVect(I, segment.HI.jsd, 0),
                                    IntVect(I, segment.HI.jed, nz-1)),
                    [=] AMREX_GPU_DEVICE (int ii, int j, int k)
                    {
                        u(I,j,k) = segment.normal_vel(I,j,k);
                    });
                }
            }
        }
    }
    */
}

} // namespace MOM
