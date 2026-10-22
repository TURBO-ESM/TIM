// SKILLS: 0.3.2
///
/// @file mom_coriolisadv.cpp
/// @brief Box-level AMReX kernels for the MOM6 Coriolis and momentum advection terms.
///
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#include "mom_coriolisadv.hpp"

/// @brief Abort with @p msg annotated by the source file and line number.
/// @param msg Message to report before aborting.
#define AMREX_ABORT_LOC(msg) \
	amrex::Abort(std::string(msg) + " [" + __FILE__ + ":" + std::to_string(__LINE__) + "]")

using namespace amrex;

namespace MOM {

void gradke(const Box& bxQ,
            const Box& bxU,
            const Box& bxV,
            Array4<const Real> const& u,
            Array4<const Real> const& v,
            Array4<Real>       const& KE,
            Array4<Real>       const& KEx,
            Array4<Real>       const& KEy,
            Array4<const Real> const& areaCu,
            Array4<const Real> const& areaCv,
            Array4<const Real> const& IareaT,
            Array4<const Real> const& mask2dCu,
            Array4<const Real> const& mask2dCv,
            Array4<const Real> const& IdxCu_OBCmask,
            Array4<const Real> const& IdyCv_OBCmask,
            const int KE_Scheme,
            const bool KE_use_limiter)
{
    constexpr Real C1_12 = 1.0/12.0;  // The ratio of 1/12 [nondim]

    // bxQ grown by one point at the end in i and j, [IscB:IecB+1, JscB:JecB+1, ksc:kec]
    Box bxQ_pij = bxQ;
    bxQ_pij.growHi(0, 1).growHi(1, 1);

    // Calculate KE (Kinetic energy for use in the -grad(KE) acceleration term).
    if (KE_Scheme == KE_ARAKAWA) {
        // The following calculation of Kinetic energy includes the metric terms
        // identified in Arakawa & Lamb 1982 as important for KE conservation.  It
        // also includes the possibility of partially-blocked tracer cell faces.
        ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            KE(i,j,k) = ( ( (areaCu(i  ,j,0)*(u(i  ,j,k)*u(i  ,j,k))) +
                            (areaCu(i-1,j,0)*(u(i-1,j,k)*u(i-1,j,k))) ) +
                          ( (areaCv(i,j  ,0)*(v(i,j  ,k)*v(i,j  ,k))) +
                            (areaCv(i,j-1,0)*(v(i,j-1,k)*v(i,j-1,k))) ) )*0.25*IareaT(i,j,0);
        });
    } else if (KE_Scheme == KE_SIMPLE_GUDONOV) {
        // The following discretization of KE is based on the one-dimensional Gudonov
        // scheme which does not take into account any geometric factors
        ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const up = 0.5*( u(i-1,j,k) + Math::abs( u(i-1,j,k) ) ); Real const up2 = up*up;
            Real const um = 0.5*( u(i  ,j,k) - Math::abs( u(i  ,j,k) ) ); Real const um2 = um*um;
            Real const vp = 0.5*( v(i,j-1,k) + Math::abs( v(i,j-1,k) ) ); Real const vp2 = vp*vp;
            Real const vm = 0.5*( v(i,j  ,k) - Math::abs( v(i,j  ,k) ) ); Real const vm2 = vm*vm;
            KE(i,j,k) = ( amrex::max(up2,um2) + amrex::max(vp2,vm2) ) *0.5;
        });
    } else if (KE_Scheme == KE_GUDONOV) {
        // The following discretization of KE is based on the one-dimensional Gudonov
        // scheme but has been adapted to take horizontal grid factors into account
        ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const up = 0.5*( u(i-1,j,k) + Math::abs( u(i-1,j,k) ) ); Real const up2a = up*up*areaCu(i-1,j,0);
            Real const um = 0.5*( u(i  ,j,k) - Math::abs( u(i  ,j,k) ) ); Real const um2a = um*um*areaCu(i  ,j,0);
            Real const vp = 0.5*( v(i,j-1,k) + Math::abs( v(i,j-1,k) ) ); Real const vp2a = vp*vp*areaCv(i,j-1,0);
            Real const vm = 0.5*( v(i,j  ,k) - Math::abs( v(i,j  ,k) ) ); Real const vm2a = vm*vm*areaCv(i,j  ,0);
            KE(i,j,k) = ( amrex::max(um2a,up2a) + amrex::max(vm2a,vp2a) )*0.5*IareaT(i,j,0);
        });
    } else if (KE_Scheme == KE_UP3) {
        // The following discretization of KE is based on the one-dimensional third-order
        // upwind scheme which does not take horizontal grid factors into account

        if (KE_use_limiter) {
            ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real up, um, vp, vm;
                // compute the masking to make sure that inland values are not used
                Real const third_order_u = (mask2dCu(i-2,j,0) * mask2dCu(i-1,j,0)*
                                            mask2dCu(i,j,0) * mask2dCu(i+1,j,0));

                if (third_order_u == 1.0) {
                    up = (7.0 * (u(i-1,j,k) + u(i,j,k)) - (u(i-2,j,k) + u(i+1,j,k))) * C1_12;
                    Real const u4[4] = {u(i-2,j,k), u(i-1,j,k), u(i,j,k), u(i+1,j,k)};
                    up3_koren_limiter_reconstruction_point(u4, up, um);
                } else {
                    up = (u(i-1,j,k) + u(i,j,k))*0.5;
                    if (up > 0.0) {
                        um = u(i-1,j,k);
                    } else if (up < 0.0) {
                        um = u(i,j,k);
                    } else {
                        um = up;
                    }
                }

                Real const third_order_v = (mask2dCv(i,j-2,0) * mask2dCv(i,j-1,0)*
                                            mask2dCv(i,j,0) * mask2dCv(i,j+1,0));
                if (third_order_v == 1.0) {
                    vp = (7.0 * (v(i,j-1,k) + v(i,j,k)) - (v(i,j-2,k) + v(i,j+1,k))) * C1_12;
                    Real const v4[4] = {v(i,j-2,k), v(i,j-1,k), v(i,j,k), v(i,j+1,k)};
                    up3_koren_limiter_reconstruction_point(v4, vp, vm);
                } else {
                    vp = (v(i,j-1,k) + v(i,j,k))*0.5;
                    if (vp > 0.0) {
                        vm = v(i,j-1,k);
                    } else if (vp < 0.0) {
                        vm = v(i,j,k);
                    } else {
                        vm = vp;
                    }
                }

                KE(i,j,k) = ( (um*um) + (vm*vm) )*0.5;
            });
        } else {
            ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real up, um, vp, vm;
                // compute the masking to make sure that inland values are not used
                Real const third_order_u = (mask2dCu(i-2,j,0) * mask2dCu(i-1,j,0)*
                                            mask2dCu(i,j,0) * mask2dCu(i+1,j,0));

                if (third_order_u == 1.0) {
                    up = (7.0 * (u(i-1,j,k) + u(i,j,k)) - (u(i-2,j,k) + u(i+1,j,k))) * C1_12;
                    Real const u4[4] = {u(i-2,j,k), u(i-1,j,k), u(i,j,k), u(i+1,j,k)};
                    up3_reconstruction_point(u4, up, um);
                } else {
                    up = (u(i-1,j,k) + u(i,j,k))*0.5;
                    if (up > 0.0) {
                        um = u(i-1,j,k);
                    } else if (up < 0.0) {
                        um = u(i,j,k);
                    } else {
                        um = up;
                    }
                }

                Real const third_order_v = (mask2dCv(i,j-2,0) * mask2dCv(i,j-1,0)*
                                            mask2dCv(i,j,0) * mask2dCv(i,j+1,0));
                if (third_order_v == 1.0) {
                    vp = (7.0 * (v(i,j-1,k) + v(i,j,k)) - (v(i,j-2,k) + v(i,j+1,k))) * C1_12;
                    Real const v4[4] = {v(i,j-2,k), v(i,j-1,k), v(i,j,k), v(i,j+1,k)};
                    up3_reconstruction_point(v4, vp, vm);
                } else {
                    vp = (v(i,j-1,k) + v(i,j,k))*0.5;
                    if (vp > 0.0) {
                        vm = v(i,j-1,k);
                    } else if (vp < 0.0) {
                        vm = v(i,j,k);
                    } else {
                        vm = vp;
                    }
                }

                KE(i,j,k) = ( (um*um) + (vm*vm) )*0.5;
            });
        }
    }

    // Term - d(KE)/dx.
    ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        KEx(i,j,k) = (KE(i+1,j,k) - KE(i,j,k)) * IdxCu_OBCmask(i,j,0);
    });

    // Term - d(KE)/dy.
    ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        KEy(i,j,k) = (KE(i,j+1,k) - KE(i,j,k)) * IdyCv_OBCmask(i,j,0);
    });
}

void coradv_setup(const Box& bxH,
                  const Box& bxQ,
                  const Box& bxU,
                  const Box& bxV,
                  const Box& bxQs,
                  Array4<const Real> const& u,
                  Array4<const Real> const& v,
                  Array4<const Real> const& h,
                  OceanOBC* OBC,
                  const bool rv_x_v_associated,
                  const bool rv_x_u_associated,
                  Array4<const Real> const& por_face_areaU,
                  Array4<const Real> const& por_face_areaV,
                  WaveParametersCS* Waves,
                  Array4<const Real> const& dyCv,
                  Array4<const Real> const& dxCu,
                  Array4<const Real> const& dy_Cu,
                  Array4<const Real> const& dx_Cv,
                  Array4<const Real> const& dyBu,
                  Array4<const Real> const& dxBu,
                  Array4<const Real> const& mask2dBu,
                  Array4<const Real> const& IareaBu,
                  Array4<const Real> const& CoriolisBu,
                  Array4<const Real> const& Area_h,
                  Array4<const Real> const& Area_q,
                  const bool Stokes_VF,
                  const bool use_weno,
                  const Real vol_neglect,
                  const Real area_neglect,
                  Array4<Real>       const& q,
                  Array4<Real>       const& Ih_q,
                  Array4<Real>       const& abs_vort,
                  Array4<Real>       const& h_q,
                  Array4<Real>       const& qS,
                  Array4<Real>       const& q2,
                  Array4<Real>       const& uh_center,
                  Array4<Real>       const& vh_center,
                  Array4<Real>       const& RV,
                  Array4<Real>       const& PV,
                  const CorAdv_opts_CS_C& CS)
{
    BL_PROFILE("coradv_setup");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    // NOTE: Stokes drift support temporarily disabled.
    // WaveParametersCS is forward-declared only.
    if (Waves != nullptr) {
        AMREX_ABORT_LOC("Waves pointer provided but not yet implemented");
    }

    const int ksc = bxH.smallEnd(2);
    const int kec = bxH.bigEnd(2);

    // bxU grown by one point at the end in i and the start in j,
    // [IscB:IecB+1, jsc-1:jec, ksc:kec]
    Box bxUx = bxU;
    bxUx.growHi(0, 1).growLo(1, 1);
    // bxV grown by one point at the start in i and the end in j,
    // [isc-1:iec, JscB:JecB+1, ksc:kec]
    Box bxVx = bxV;
    bxVx.growLo(0, 1).growHi(1, 1);
    // bxH grown by one point at the start in i and j,
    // [isc-1:iec, jsc-1:jec, ksc:kec]
    Box bxH_mij = bxH;
    bxH_mij.growLo(0, 1).growLo(1, 1);
    // bxQ grown by one point on both sides in i and j,
    // [IscB-1:IecB+1, JscB-1:JecB+1, ksc:kec]
    Box bxQ_gij = bxQ;
    bxQ_gij.grow(0, 1).grow(1, 1);
    // bxQs grown by one point at the end in i,
    // [Is_q:Ie_q+1, Js_q:Je_q, ksc:kec]
    Box bxQs_pi = bxQs;
    bxQs_pi.growHi(0, 1);
    // bxQs grown by one point at the end in j,
    // [Is_q:Ie_q, Js_q:Je_q+1, ksc:kec]
    Box bxQs_pj = bxQs;
    bxQs_pj.growHi(1, 1);

    const Box q_box(IntVect(Box(q).smallEnd(0), Box(q).smallEnd(1), ksc),
                    IntVect(Box(q).bigEnd(0), Box(q).bigEnd(1), kec));
    const Box u_box(IntVect(Box(uh_center).smallEnd(0), Box(uh_center).smallEnd(1), ksc),
                    IntVect(Box(uh_center).bigEnd(0), Box(uh_center).bigEnd(1), kec));
    const Box v_box(IntVect(Box(vh_center).smallEnd(0), Box(vh_center).smallEnd(1), ksc),
                    IntVect(Box(vh_center).bigEnd(0), Box(vh_center).bigEnd(1), kec));

    // Contributions to the circulation around q-points [L2 T-1 ~> m2 s-1]
    FArrayBox dvdx_fab(q_box, 1, The_Async_Arena());
    FArrayBox dudy_fab(q_box, 1, The_Async_Arena());
    // Relative vorticity at q-points [T-1 ~> s-1].
    FArrayBox rel_vort_fab(q_box, 1, The_Async_Arena());
    // The cell area weighted thickness interpolated to u points
    // times the effective areas [H L2 ~> m3 or kg].
    FArrayBox hArea_u_fab(u_box, 1, The_Async_Arena());
    // The cell area weighted thickness interpolated to v points
    // times the effective areas [H L2 ~> m3 or kg].
    FArrayBox hArea_v_fab(v_box, 1, The_Async_Arena());

    Array4<Real> const dvdx     = dvdx_fab.array();
    Array4<Real> const dudy     = dudy_fab.array();
    Array4<Real> const rel_vort = rel_vort_fab.array();
    Array4<Real> const hArea_u  = hArea_u_fab.array();
    Array4<Real> const hArea_v  = hArea_v_fab.array();

    /*
    // idem. for Stokes drift [L2 T-1 ~> m2 s-1]
    FArrayBox dvSdx_fab(q_box, 1, The_Async_Arena());
    FArrayBox duSdy_fab(q_box, 1, The_Async_Arena());
    // Stokes vorticity at q-points [T-1 ~> s-1].
    FArrayBox stk_vort_fab(q_box, 1, The_Async_Arena());

    Array4<Real> const dvSdx    = dvSdx_fab.array();
    Array4<Real> const duSdy    = duSdy_fab.array();
    Array4<Real> const stk_vort = stk_vort_fab.array();

    Array4<const Real> const us_x = Waves->us_x;
    Array4<const Real> const us_y = Waves->us_y;
    */

    // Here the second order accurate layer potential vorticities, q,
    // are calculated.  hq is  second order accurate in space.  Relative
    // vorticity is second order accurate everywhere with free slip b.c.s,
    // but only first order accurate at boundaries with no slip b.c.s.
    // First calculate the contributions to the circulation around the q-point.
    if (Stokes_VF) {
        /*
        if (CS.do_CAS_diag) {
            ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dvSdx(i,j,k) = (-us_y(i+1,j,k)*dyCv(i+1,j,0)) -
                               (-us_y(i,j,k)*dyCv(i,j,0));
                duSdy(i,j,k) = (-us_x(i,j+1,k)*dxCu(i,j+1,0)) -
                               (-us_x(i,j,k)*dxCu(i,j,0));
            });
        }
        if (!Waves->Passive_Stokes_VF) {
            ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dvdx(i,j,k) = ((v(i+1,j,k)-us_y(i+1,j,k))*dyCv(i+1,j,0)) -
                              ((v(i,j,k)-us_y(i,j,k))*dyCv(i,j,0));
                dudy(i,j,k) = ((u(i,j+1,k)-us_x(i,j+1,k))*dxCu(i,j+1,0)) -
                              ((u(i,j,k)-us_x(i,j,k))*dxCu(i,j,0));
            });
        } else {
            ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dvdx(i,j,k) = (v(i+1,j,k)*dyCv(i+1,j,0)) - (v(i,j,k)*dyCv(i,j,0));
                dudy(i,j,k) = (u(i,j+1,k)*dxCu(i,j+1,0)) - (u(i,j,k)*dxCu(i,j,0));
            });
        }
        */
    } else {
        ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dvdx(i,j,k) = (v(i+1,j,k)*dyCv(i+1,j,0)) - (v(i,j,k)*dyCv(i,j,0));
            dudy(i,j,k) = (u(i,j+1,k)*dxCu(i,j+1,0)) - (u(i,j,k)*dxCu(i,j,0));
        });
    }

    ParallelFor(bxQs_pi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        hArea_v(i,j,k) = 0.5*((Area_h(i,j,0) * h(i,j,k)) + (Area_h(i,j+1,0) * h(i,j+1,k)));
    });

    ParallelFor(bxQs_pj, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        hArea_u(i,j,k) = 0.5*((Area_h(i,j,0) * h(i,j,k)) + (Area_h(i+1,j,0) * h(i+1,j,k)));
    });

    if (CS.Coriolis_scheme_CS.Coriolis_En_Dis) {
        ParallelFor(bxVx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            uh_center(i,j,k) = 0.5 * ((dy_Cu(i,j,0)*por_face_areaU(i,j,k)) * u(i,j,k)) * (h(i,j,k) + h(i+1,j,k));
        });

        ParallelFor(bxUx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vh_center(i,j,k) = 0.5 * ((dx_Cv(i,j,0)*por_face_areaV(i,j,k)) * v(i,j,k)) * (h(i,j,k) + h(i,j+1,k));
        });
    }

    /*
    // Adjust circulation components to relative vorticity and thickness projected onto
    // velocity points on open boundaries.
    if (OBC != nullptr) {
        const int Isq  = bxQ.smallEnd(0);
        const int Ieq  = bxQ.bigEnd(0);
        const int Jsq  = bxQ.smallEnd(1);
        const int Jeq  = bxQ.bigEnd(1);
        const int Is_q = bxQs.smallEnd(0);
        const int Ie_q = bxQs.bigEnd(0);
        const int Js_q = bxQs.smallEnd(1);
        const int Je_q = bxQs.bigEnd(1);
        const bool Coriolis_En_Dis = CS.Coriolis_scheme_CS.Coriolis_En_Dis;

        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto& segment = OBC->segment[n];
            if (!segment.on_pe) continue;

            const int I = segment.HI.IsdB;
            const int J = segment.HI.JsdB;
            const int direction = segment.direction;
            Array4<const Real> const tangential_vel  = segment.tangential_vel;
            Array4<const Real> const tangential_grad = segment.tangential_grad;

            if (segment.is_N_or_S && (J >= Js_q) && (J <= Je_q)) {
                const Box bxI(IntVect(segment.HI.IsdB, J, ksc), IntVect(segment.HI.IedB, J, kec));
                switch (OBC->vorticity_config) {
                case OBC_VORTICITY_ZERO:
                    ParallelFor(bxI, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                    {
                        dvdx(i,J,k) = 0.0; dudy(i,J,k) = 0.0;
                    });
                    break;
                case OBC_VORTICITY_FREESLIP:
                    ParallelFor(bxI, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                    {
                        dudy(i,J,k) = 0.0;
                    });
                    break;
                case OBC_VORTICITY_COMPUTED:
                    ParallelFor(bxI, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_N) {
                            dudy(i,J,k) = 2.0*(tangential_vel(i,J,k) - u(i,J,k))*dxCu(i,J,0);
                        } else { // (direction == OBC_DIRECTION_S)
                            dudy(i,J,k) = 2.0*(u(i,J+1,k) - tangential_vel(i,J,k))*dxCu(i,J+1,0);
                        }
                    });
                    break;
                case OBC_VORTICITY_SPECIFIED:
                    ParallelFor(bxI, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_N) {
                            dudy(i,J,k) = tangential_grad(i,J,k)*dxCu(i,J,0)*dyBu(i,J,0);
                        } else { // (direction == OBC_DIRECTION_S)
                            dudy(i,J,k) = tangential_grad(i,J,k)*dxCu(i,J+1,0)*dyBu(i,J,0);
                        }
                    });
                    break;
                }

                // Project thicknesses across OBC points with a no-gradient condition.
                const Box bxh(IntVect(amrex::max(Is_q, segment.HI.isd), J, ksc),
                              IntVect(amrex::min(Ie_q+1, segment.HI.ied), J, kec));
                if (bxh.ok()) {
                    ParallelFor(bxh, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_N) {
                            hArea_v(i,J,k) = 0.5 * (Area_h(i,J,0) + Area_h(i,J+1,0)) * h(i,J,k);
                        } else { // (direction == OBC_DIRECTION_S)
                            hArea_v(i,J,k) = 0.5 * (Area_h(i,J,0) + Area_h(i,J+1,0)) * h(i,J+1,k);
                        }
                    });
                }

                if (Coriolis_En_Dis) {
                    const Box bxc(IntVect(amrex::max(Isq, segment.HI.isd), J, ksc),
                                  IntVect(amrex::min(Ieq+1, segment.HI.ied), J, kec));
                    if (bxc.ok()) {
                        ParallelFor(bxc, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                        {
                            if (direction == OBC_DIRECTION_N) {
                                vh_center(i,J,k) = (dx_Cv(i,J,0)*por_face_areaV(i,J,k)) * v(i,J,k) * h(i,J,k);
                            } else { // (direction == OBC_DIRECTION_S)
                                vh_center(i,J,k) = (dx_Cv(i,J,0)*por_face_areaV(i,J,k)) * v(i,J,k) * h(i,J+1,k);
                            }
                        });
                    }
                }
            } else if (segment.is_E_or_W && (I >= Is_q) && (I <= Ie_q)) {
                const Box bxJ(IntVect(I, segment.HI.JsdB, ksc), IntVect(I, segment.HI.JedB, kec));
                switch (OBC->vorticity_config) {
                case OBC_VORTICITY_ZERO:
                    ParallelFor(bxJ, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                    {
                        dvdx(I,j,k) = 0.0; dudy(I,j,k) = 0.0;
                    });
                    break;
                case OBC_VORTICITY_FREESLIP:
                    ParallelFor(bxJ, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                    {
                        dvdx(I,j,k) = 0.0;
                    });
                    break;
                case OBC_VORTICITY_COMPUTED:
                    ParallelFor(bxJ, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_E) {
                            dvdx(I,j,k) = 2.0*(tangential_vel(I,j,k) - v(I,j,k))*dyCv(I,j,0);
                        } else { // (direction == OBC_DIRECTION_W)
                            dvdx(I,j,k) = 2.0*(v(I+1,j,k) - tangential_vel(I,j,k))*dyCv(I+1,j,0);
                        }
                    });
                    break;
                case OBC_VORTICITY_SPECIFIED:
                    ParallelFor(bxJ, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_E) {
                            dvdx(I,j,k) = tangential_grad(I,j,k)*dyCv(I,j,0)*dxBu(I,j,0);
                        } else { // (direction == OBC_DIRECTION_W)
                            dvdx(I,j,k) = tangential_grad(I,j,k)*dyCv(I+1,j,0)*dxBu(I,j,0);
                        }
                    });
                    break;
                }

                // Project thicknesses across OBC points with a no-gradient condition.
                const Box bxh(IntVect(I, amrex::max(Js_q, segment.HI.jsd), ksc),
                              IntVect(I, amrex::min(Je_q+1, segment.HI.jed), kec));
                if (bxh.ok()) {
                    ParallelFor(bxh, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_E) {
                            hArea_u(I,j,k) = 0.5*(Area_h(I,j,0) + Area_h(I+1,j,0)) * h(I,j,k);
                        } else { // (direction == OBC_DIRECTION_W)
                            hArea_u(I,j,k) = 0.5*(Area_h(I,j,0) + Area_h(I+1,j,0)) * h(I+1,j,k);
                        }
                    });
                }

                if (Coriolis_En_Dis) {
                    const Box bxc(IntVect(I, amrex::max(Jsq, segment.HI.jsd), ksc),
                                  IntVect(I, amrex::min(Jeq+1, segment.HI.jed), kec));
                    if (bxc.ok()) {
                        ParallelFor(bxc, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                        {
                            if (direction == OBC_DIRECTION_E) {
                                uh_center(I,j,k) = (dy_Cu(I,j,0)*por_face_areaU(I,j,k)) * u(I,j,k) * h(I,j,k);
                            } else { // (direction == OBC_DIRECTION_W)
                                uh_center(I,j,k) = (dy_Cu(I,j,0)*por_face_areaU(I,j,k)) * u(I,j,k) * h(I+1,j,k);
                            }
                        });
                    }
                }
            }
        }

        // Now project thicknesses across cell-corner points in the OBCs.  The two
        // projections have to occur in sequence and can not be combined easily.
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto& segment = OBC->segment[n];
            if (!segment.on_pe) continue;

            const int I = segment.HI.IsdB;
            const int J = segment.HI.JsdB;
            const int direction = segment.direction;

            if (segment.is_N_or_S && (J >= Js_q) && (J <= Je_q)) {
                const Box bxI(IntVect(amrex::max(Is_q, segment.HI.IsdB), J, ksc),
                              IntVect(amrex::min(Ie_q, segment.HI.IedB), J, kec));
                if (bxI.ok()) {
                    ParallelFor(bxI, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_N) {
                            if (Area_h(i,J,0) + Area_h(i+1,J,0) > 0.0) {
                                hArea_u(i,J+1,k) = hArea_u(i,J,k) * ((Area_h(i,J+1,0) + Area_h(i+1,J+1,0)) /
                                                                     (Area_h(i,J,0) + Area_h(i+1,J,0)));
                            } else { hArea_u(i,J+1,k) = 0.0; }
                        } else { // (direction == OBC_DIRECTION_S)
                            if (Area_h(i,J+1,0) + Area_h(i+1,J+1,0) > 0.0) {
                                hArea_u(i,J,k) = hArea_u(i,J+1,k) * ((Area_h(i,J,0) + Area_h(i+1,J,0)) /
                                                                     (Area_h(i,J+1,0) + Area_h(i+1,J+1,0)));
                            } else { hArea_u(i,J,k) = 0.0; }
                        }
                    });
                }
            } else if (segment.is_E_or_W && (I >= Is_q) && (I <= Ie_q)) {
                const Box bxJ(IntVect(I, amrex::max(Js_q, segment.HI.JsdB), ksc),
                              IntVect(I, amrex::min(Je_q, segment.HI.JedB), kec));
                if (bxJ.ok()) {
                    ParallelFor(bxJ, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                    {
                        if (direction == OBC_DIRECTION_E) {
                            if (Area_h(I,j,0) + Area_h(I,j+1,0) > 0.0) {
                                hArea_v(I+1,j,k) = hArea_v(I,j,k) * ((Area_h(I+1,j,0) + Area_h(I+1,j+1,0)) /
                                                                     (Area_h(I,j,0) + Area_h(I,j+1,0)));
                            } else { hArea_v(I+1,j,k) = 0.0; }
                        } else { // (direction == OBC_DIRECTION_W)
                            hArea_v(I,j,k) = 0.5 * (Area_h(I,j,0) + Area_h(I,j+1,0)) * h(I,j+1,k);
                            if (Area_h(I+1,j,0) + Area_h(I+1,j+1,0) > 0.0) {
                                hArea_v(I,j,k) = hArea_v(I+1,j,k) * ((Area_h(I,j,0) + Area_h(I,j+1,0)) /
                                                                     (Area_h(I+1,j,0) + Area_h(I+1,j+1,0)));
                            } else { hArea_v(I,j,k) = 0.0; }
                        }
                    });
                }
            }
        }
    }
    */

    if (CS.no_slip) {
        ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rel_vort(i,j,k) = (2.0 - mask2dBu(i,j,0)) * (dvdx(i,j,k) - dudy(i,j,k)) * IareaBu(i,j,0);
        });

        if (Stokes_VF) {
            /*
            if (CS.do_CAS_diag) {
                ParallelFor(bxQ_gij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    stk_vort(i,j,k) = (2.0 - mask2dBu(i,j,0)) * (dvSdx(i,j,k) - duSdy(i,j,k)) * IareaBu(i,j,0);
                });
            }
            */
        }
    } else {
        ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rel_vort(i,j,k) = mask2dBu(i,j,0) * (dvdx(i,j,k) - dudy(i,j,k)) * IareaBu(i,j,0);
        });

        if (Stokes_VF) {
            /*
            if (CS.do_CAS_diag) {
                ParallelFor(bxQ_gij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    stk_vort(i,j,k) = (2.0 - mask2dBu(i,j,0)) * (dvSdx(i,j,k) - duSdy(i,j,k)) * IareaBu(i,j,0);
                });
            }
            */
        }
    }

    ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        abs_vort(i,j,k) = CoriolisBu(i,j,0) + rel_vort(i,j,k);
    });

    ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // The sum of area times thickness of the cells
        // surrounding a q point [H L2 ~> m3 or kg].
        Real const hArea_q = (hArea_u(i,j,k) + hArea_u(i,j+1,k)) + (hArea_v(i,j,k) + hArea_v(i+1,j,k));
        Ih_q(i,j,k) = Area_q(i,j,0) / (hArea_q + vol_neglect);
        q(i,j,k) = abs_vort(i,j,k) * Ih_q(i,j,k);
    });

    // NOTE: `h_q` is only used by WENO and was pulled out of the above loop to
    //   improve GPU performance, but it may need to be moved back.
    if (use_weno) {
        ParallelFor(bxQs, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const hArea_q = (hArea_u(i,j,k) + hArea_u(i,j+1,k)) + (hArea_v(i,j,k) + hArea_v(i+1,j,k));
            h_q(i,j,k) = hArea_q / amrex::max(Area_q(i,j,0), area_neglect);
        });
    }

    if (Stokes_VF) {
        /*
        if (CS.do_CAS_diag) {
            ParallelFor(bxH_mij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                qS(i,j,k) = stk_vort(i,j,k) * Ih_q(i,j,k);
            });
        }
        */
    }

    if (CS.do_RV_diag) {
        ParallelFor(bxQ_gij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            RV(i,j,k) = rel_vort(i,j,k);
        });
    }

    if (CS.do_PV_diag) {
        ParallelFor(bxQ_gij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            PV(i,j,k) = q(i,j,k);
        });
    }

    if (rv_x_v_associated || rv_x_u_associated) {
        ParallelFor(bxQ_gij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            q2(i,j,k) = rel_vort(i,j,k) * Ih_q(i,j,k);
        });
    }
}

void coradv_common_terms(const Box& bxU,
                         const Box& bxV,
                         Array4<const Real> const& u,
                         Array4<const Real> const& v,
                         Array4<const Real> const& uh,
                         Array4<const Real> const& vh,
                         Array4<const Real> const& abs_vort,
                         Array4<const Real> const& qS,
                         Array4<const Real> const& q2,
                         Array4<const Real> const& KEx,
                         Array4<const Real> const& KEy,
                         const bool Stokes_VF,
                         Array4<const Real> const& IdxCu,
                         Array4<const Real> const& IdyCv,
                         Array4<Real>       const& CAu,
                         Array4<Real>       const& CAv,
                         Array4<Real>       const& CAuS,
                         Array4<Real>       const& CAvS,
                         Array4<Real>       const& gradKEu,
                         Array4<Real>       const& gradKEv,
                         Array4<Real>       const& rv_x_v,
                         Array4<Real>       const& rv_x_u,
                         const CorAdv_opts_CS_C& CS)
{
    BL_PROFILE("coradv_common_terms");

    constexpr Real C1_12 = 1.0 / 12.0;  // C1_12 = 1/12 [nondim]

    const bool gradKEu_associated = (gradKEu.p != nullptr);
    const bool gradKEv_associated = (gradKEv.p != nullptr);
    const bool rv_x_u_associated  = (rv_x_u.p != nullptr);
    const bool rv_x_v_associated  = (rv_x_v.p != nullptr);

    if (Stokes_VF) {
        if (CS.do_CAS_diag) {
            // Computing the diagnostic Stokes contribution to CAu
            ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                CAuS(i,j,k) = 0.25 *
                      ((qS(i,j,k) * (vh(i+1,j,k) + vh(i,j,k))) +
                       (qS(i,j-1,k) * (vh(i,j-1,k) + vh(i+1,j-1,k)))) * IdxCu(i,j,0);
            });
        }
    }

    if (CS.bound_Coriolis) {
        ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // (f+rv)*v at the 4 points surrounding a u points[L T-2 ~> m s-2]
            Real const fv1 = abs_vort(i,j,k) * v(i+1,j,k);
            Real const fv2 = abs_vort(i,j,k) * v(i,j,k);
            Real const fv3 = abs_vort(i,j-1,k) * v(i+1,j-1,k);
            Real const fv4 = abs_vort(i,j-1,k) * v(i,j-1,k);

            // The maximum and minimum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
            Real const max_fv = amrex::max(fv1, fv2, fv3, fv4);
            Real const min_fv = amrex::min(fv1, fv2, fv3, fv4);

            CAu(i,j,k) = amrex::min(CAu(i,j,k), max_fv);
            CAu(i,j,k) = amrex::max(CAu(i,j,k), min_fv);
        });
    }

    // Term - d(KE)/dx.
    ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        CAu(i,j,k) = CAu(i,j,k) - KEx(i,j,k);
    });

    if (gradKEu_associated) {
        ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            gradKEu(i,j,k) = -KEx(i,j,k);
        });
    }

    if (Stokes_VF) {
        if (CS.do_CAS_diag) {
            // Computing the diagnostic Stokes contribution to CAv
            ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                CAvS(i,j,k) = 0.25 *
                      ((qS(i,j,k) * (uh(i,j+1,k) + uh(i,j,k))) +
                       (qS(i-1,j,k) * (uh(i-1,j,k) + uh(i-1,j+1,k)))) * IdyCv(i,j,0);
            });
        }
    }

    if (CS.bound_Coriolis) {
        ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // -(f+rv)*u at the 4 points surrounding a v point [L T-2 ~> m s-2]
            Real const fu1 = -abs_vort(i,j,k) * u(i,j+1,k);
            Real const fu2 = -abs_vort(i,j,k) * u(i,j,k);
            Real const fu3 = -abs_vort(i-1,j,k) * u(i-1,j+1,k);
            Real const fu4 = -abs_vort(i-1,j,k) * u(i-1,j,k);

            // The maximum and minimum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
            Real const max_fu = amrex::max(fu1, fu2, fu3, fu4);
            Real const min_fu = amrex::min(fu1, fu2, fu3, fu4);

            CAv(i,j,k) = amrex::min(CAv(i,j,k), max_fu);
            CAv(i,j,k) = amrex::max(CAv(i,j,k), min_fu);
        });
    }

    // Term - d(KE)/dy.
    ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        CAv(i,j,k) = CAv(i,j,k) - KEy(i,j,k);
    });
    if (gradKEv_associated) {
        ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            gradKEv(i,j,k) = -KEy(i,j,k);
        });
    }

    if (rv_x_u_associated || rv_x_v_associated) {
        // Calculate the Coriolis-like acceleration due to relative vorticity.
        if (CS.Coriolis_scheme_CS.Coriolis_Scheme == SADOURNY75_ENERGY) {
            if (rv_x_u_associated) {
                ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rv_x_u(i,j,k) = - 0.25*
                        ((q2(i-1,j,k)*(uh(i-1,j,k) + uh(i-1,j+1,k))) +
                         (q2(i,j,k)*(uh(i,j,k) + uh(i,j+1,k)))) * IdyCv(i,j,0);
                });
            }

            if (rv_x_v_associated) {
                ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rv_x_v(i,j,k) = 0.25 *
                        ((q2(i,j,k) * (vh(i+1,j,k) + vh(i,j,k))) +
                         (q2(i,j-1,k) * (vh(i,j-1,k) + vh(i+1,j-1,k)))) * IdxCu(i,j,0);
                });
            }
        } else {
            if (rv_x_u_associated) {
                ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rv_x_u(i,j,k) = -IdyCv(i,j,0) * C1_12 *
                        (((((q2(i,j,k) + q2(i-1,j-1,k)) + q2(i-1,j,k)) * uh(i-1,j,k)) +
                          (((q2(i-1,j,k) + q2(i,j+1,k)) + q2(i,j,k)) * uh(i,j+1,k))) +
                         ((((q2(i-1,j,k) + q2(i,j-1,k)) + q2(i,j,k)) * uh(i,j,k))+
                          (((q2(i,j,k) + q2(i-1,j+1,k)) + q2(i-1,j,k)) * uh(i-1,j+1,k))));
                });
            }

            if (rv_x_v_associated) {
                ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rv_x_v(i,j,k) = IdxCu(i,j,0) * C1_12 *
                        (((((q2(i+1,j,k) + q2(i,j-1,k)) + q2(i,j,k)) * vh(i+1,j,k)) +
                          (((q2(i-1,j-1,k) + q2(i,j,k)) + q2(i,j-1,k)) * vh(i,j-1,k))) +
                         ((((q2(i-1,j,k) + q2(i,j-1,k)) + q2(i,j,k)) * vh(i,j,k)) +
                          (((q2(i+1,j-1,k) + q2(i,j,k)) + q2(i,j-1,k)) * vh(i+1,j-1,k))));
                });
            }
        }
    }
}

void coradv_arakawa(const Box& bxH,
                    const Box& bxQ,
                    const Box& bxU,
                    const Box& bxV,
                    Array4<const Real> const& uh,
                    Array4<const Real> const& vh,
                    Array4<const Real> const& q,
                    Array4<const Real> const& Ih_q,
                    const Real Fe_m2,
                    const Real rat_lin,
                    Array4<const Real> const& IdxCu,
                    Array4<const Real> const& IdyCv,
                    Array4<Real>       const& CAu,
                    Array4<Real>       const& CAv,
                    const Real wt_lin_blend,
                    const Coriolis_scheme_CS_C& CS)
{
    BL_PROFILE("coradv_arakawa");

    constexpr Real C1_12 = 1.0 / 12.0;  // C1_12 = 1/12 [nondim]
    constexpr Real C1_24 = 1.0 / 24.0;  // C1_24 = 1/24 [nondim]

    const int ksc = bxH.smallEnd(2);
    const int kec = bxH.bigEnd(2);

    // bxV grown by one point at the start in i and the end in j,
    // [isc-1:iec, JscB:JecB+1, ksc:kec]
    Box bxVx = bxV;
    bxVx.growLo(0, 1).growHi(1, 1);
    // bxQ grown by one point at the end in i and j,
    // [IscB:IecB+1, JscB:JecB+1, ksc:kec]
    Box bxQ_pij = bxQ;
    bxQ_pij.growHi(0, 1).growHi(1, 1);
    // bxQ grown by one point at the end in j,
    // [IscB:IecB, JscB:JecB+1, ksc:kec]
    Box bxQ_pj = bxQ;
    bxQ_pj.growHi(1, 1);

    const Box abcd_box(IntVect(Box(uh).smallEnd(0), Box(uh).smallEnd(1), ksc),
                       IntVect(Box(uh).bigEnd(0), Box(uh).bigEnd(1), kec));
    const Box ep_box(IntVect(Box(vh).smallEnd(0), Box(uh).smallEnd(1), ksc),
                     IntVect(Box(vh).bigEnd(0), Box(uh).bigEnd(1), kec));

    // a, b, c, & d are combinations of the potential vorticities
    // surrounding an h grid point.  At small scales, a = q/4,
    // b = q/4, etc.  All are in [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1],
    // and use the indexing of the corresponding u point.
    FArrayBox a_fab(abcd_box, 1, The_Async_Arena());
    FArrayBox b_fab(abcd_box, 1, The_Async_Arena());
    FArrayBox c_fab(abcd_box, 1, The_Async_Arena());
    FArrayBox d_fab(abcd_box, 1, The_Async_Arena());
    // Additional pseudo-Coriolis terms in the Arakawa and Lamb
    // discretization [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    FArrayBox ep_u_fab(ep_box, 1, The_Async_Arena());
    FArrayBox ep_v_fab(ep_box, 1, The_Async_Arena());

    Array4<Real> const a    = a_fab.array();
    Array4<Real> const b    = b_fab.array();
    Array4<Real> const c    = c_fab.array();
    Array4<Real> const d    = d_fab.array();
    Array4<Real> const ep_u = ep_u_fab.array();
    Array4<Real> const ep_v = ep_v_fab.array();

    //   a, b, c, and d are combinations of neighboring potential
    // vorticities which form the Arakawa and Hsu vorticity advection
    // scheme.  All are defined at u grid points.

    if (CS.Coriolis_Scheme == ARAKAWA_HSU90) {
        ParallelFor(bxVx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) = (q(i,j,k) + (q(i+1,j,k) + q(i,j-1,k))) * C1_12;
            d(i,j,k) = ((q(i,j,k) + q(i+1,j-1,k)) + q(i,j-1,k)) * C1_12;
        });

        ParallelFor(bxQ_pj, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            b(i,j,k) = (q(i,j,k) + (q(i-1,j,k) + q(i,j-1,k))) * C1_12;
            c(i,j,k) = ((q(i,j,k) + q(i-1,j-1,k)) + q(i,j-1,k)) * C1_12;
        });
    } else if (CS.Coriolis_Scheme == ARAKAWA_LAMB81) {
        ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i-1,j,k) = (2.0*(q(i,j,k) + q(i-1,j-1,k)) + (q(i-1,j,k) + q(i,j-1,k))) * C1_24;
            d(i-1,j,k) = ((q(i,j,k) + q(i-1,j-1,k)) + 2.0*(q(i-1,j,k) + q(i,j-1,k))) * C1_24;
            b(i,j,k) =   ((q(i,j,k) + q(i-1,j-1,k)) + 2.0*(q(i-1,j,k) + q(i,j-1,k))) * C1_24;
            c(i,j,k) =   (2.0*(q(i,j,k) + q(i-1,j-1,k)) + (q(i-1,j,k) + q(i,j-1,k))) * C1_24;
            ep_u(i,j,k) = ((q(i,j,k) - q(i-1,j-1,k)) + (q(i-1,j,k) - q(i,j-1,k))) * C1_24;
            ep_v(i,j,k) = (-(q(i,j,k) - q(i-1,j-1,k)) + (q(i-1,j,k) - q(i,j-1,k))) * C1_24;
        });
    } else if (CS.Coriolis_Scheme == AL_BLEND) {
        // Fe_m2 and rat_lin are k-independent; computed in CorAdCalc_TR before the tile loop.
        ParallelFor(bxQ_pij, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // The maximum and minimum of the nearby Ihq [H-1 ~> m-1 or m2 kg-1].
            Real const min_Ihq = amrex::min(Ih_q(i-1,j-1,k), Ih_q(i,j-1,k), Ih_q(i-1,j,k), Ih_q(i,j,k));
            Real const max_Ihq = amrex::max(Ih_q(i-1,j-1,k), Ih_q(i,j-1,k), Ih_q(i-1,j,k), Ih_q(i,j,k));
            // The ratio of the maximum neighboring inverse thickness
            // to the minimum inverse thickness minus 1 [nondim]. rat_m1 >= 0.
            Real rat_m1 = 1.0e15;
            if (max_Ihq < 1.0e15*min_Ihq) rat_m1 = max_Ihq / min_Ihq - 1.0;
            // The weights used here are designed to keep the effective Coriolis
            // acceleration from any one point on its neighbors within a factor
            // of F_eff_max.  The minimum permitted value is 2 (the factor for
            // Sadourny's energy conserving scheme).

            // The relative weight of the Arakawa & Lamb scheme to the
            // Arakawa & Hsu scheme [nondim], between 0 and 1.
            Real AL_wt;
            // The relative weight of the Sadourny energy scheme to
            // the other two with the ARAKAWA_LAMB_BLEND scheme [nondim],
            // between 0 and 1.
            Real Sad_wt;

            // Determine the relative weights of Arakawa & Lamb vs. Arakawa and Hsu.
            if (rat_m1 <= Fe_m2) { AL_wt = 1.0; }
            else if (rat_m1 < 1.5*Fe_m2) { AL_wt = 3.0*Fe_m2 / rat_m1 - 2.0; }
            else { AL_wt = 0.0; }

            // Determine the relative weights of Sadourny Energy vs. the other two.
            if (rat_m1 <= 1.5*Fe_m2) { Sad_wt = 0.0; }
            else if (rat_m1 <= rat_lin) {
                Sad_wt = 1.0 - (1.5*Fe_m2) / rat_m1;
            } else if (rat_m1 < 2.0*rat_lin) {
                Sad_wt = 1.0 - (wt_lin_blend / rat_lin) * (rat_m1 - 2.0*rat_lin);
            } else { Sad_wt = 1.0; }

            a(i-1,j,k) = Sad_wt * 0.25 * q(i-1,j,k) + (1.0 - Sad_wt) *
                          ( ((2.0-AL_wt)* q(i-1,j,k) + AL_wt*q(i,j-1,k)) +
                             2.0 * (q(i,j,k) + q(i-1,j-1,k)) ) * C1_24;
            d(i-1,j,k) = Sad_wt * 0.25 * q(i-1,j-1,k) + (1.0 - Sad_wt) *
                          ( ((2.0-AL_wt)* q(i-1,j-1,k) + AL_wt*q(i,j,k)) +
                             2.0 * (q(i-1,j,k) + q(i,j-1,k)) ) * C1_24;
            b(i,j,k) =   Sad_wt * 0.25 * q(i,j,k) + (1.0 - Sad_wt) *
                          ( ((2.0-AL_wt)* q(i,j,k) + AL_wt*q(i-1,j-1,k)) +
                             2.0 * (q(i-1,j,k) + q(i,j-1,k)) ) * C1_24;
            c(i,j,k) =   Sad_wt * 0.25 * q(i,j-1,k) + (1.0 - Sad_wt) *
                          ( ((2.0-AL_wt)* q(i,j-1,k) + AL_wt*q(i-1,j,k)) +
                             2.0 * (q(i,j,k) + q(i-1,j-1,k)) ) * C1_24;
            ep_u(i,j,k) = AL_wt  * ((q(i,j,k) - q(i-1,j-1,k)) + (q(i-1,j,k) - q(i,j-1,k))) * C1_24;
            ep_v(i,j,k) = AL_wt * (-(q(i,j,k) - q(i-1,j-1,k)) + (q(i-1,j,k) - q(i,j-1,k))) * C1_24;
        });
    }

    // Calculate the tendencies of zonal velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAu =  q * vh - d(KE)/dx.
    // (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
    ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        CAu(i,j,k) = (((a(i,j,k) * vh(i+1,j,k)) +  (c(i,j,k) * vh(i,j-1,k)))  +
                      ((b(i,j,k) * vh(i,j,k)) +  (d(i,j,k) * vh(i+1,j-1,k)))) * IdxCu(i,j,0);
    });

    // Add in the additional terms with Arakawa & Lamb.
    if ((CS.Coriolis_Scheme == ARAKAWA_LAMB81) ||
        (CS.Coriolis_Scheme == AL_BLEND)) {
        ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            CAu(i,j,k) = CAu(i,j,k) +
                  ((ep_u(i,j,k)*uh(i-1,j,k)) - (ep_u(i+1,j,k)*uh(i+1,j,k))) * IdxCu(i,j,0);
        });
    }

    // Calculate the tendencies of meridional velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAv = - q * uh - d(KE)/dy.
    // (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
    ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        CAv(i,j,k) = - (((a(i-1,j,k)   * uh(i-1,j,k)) +
                         (c(i,j+1,k)   * uh(i,j+1,k)))
                      + ((b(i,j,k)     * uh(i,j,k)) +
                         (d(i-1,j+1,k) * uh(i-1,j+1,k)))) * IdyCv(i,j,0);
    });
    // Add in the additonal terms with Arakawa & Lamb.
    if ((CS.Coriolis_Scheme == ARAKAWA_LAMB81) ||
        (CS.Coriolis_Scheme == AL_BLEND)) {
        ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            CAv(i,j,k) = CAv(i,j,k) +
                  ((ep_v(i,j,k)*vh(i,j-1,k)) - (ep_v(i,j+1,k)*vh(i,j+1,k))) * IdyCv(i,j,0);
        });
    }
}

void coradv_robust_enstro(const Box& bxU,
                          const Box& bxV,
                          Array4<const Real> const& u,
                          Array4<const Real> const& v,
                          Array4<const Real> const& h,
                          Array4<const Real> const& uh,
                          Array4<const Real> const& vh,
                          Array4<const Real> const& abs_vort,
                          const Real eps_vel,
                          const Real h_tiny,
                          Array4<const Real> const& IdxCv,
                          Array4<const Real> const& IdxCu,
                          Array4<const Real> const& IdyCu,
                          Array4<const Real> const& IdyCv,
                          Array4<Real>       const& CAu,
                          Array4<Real>       const& CAv,
                          const int PV_Adv_Scheme)
{
    BL_PROFILE("coradv_robust_enstro");

    // Calculate the tendencies of zonal velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAu =  q * vh - d(KE)/dx.
    // An enstrophy conserving scheme robust to vanishing layers
    // Note: Heffs are in lieu of h_at_v that should be returned by the
    //       continuity solver. AJA
    if (PV_Adv_Scheme == PV_ADV_CENTERED) {
        ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Temporary effective H at U or V points [H ~> m or kg m-2].
            Real Heff1 = amrex::Math::abs(vh(i,j,k) * IdxCv(i,j,0)) / (eps_vel+amrex::Math::abs(v(i,j,k)));
            Heff1 = amrex::max(Heff1, amrex::min(h(i,j,k),h(i,j+1,k)));
            Heff1 = amrex::min(Heff1, amrex::max(h(i,j,k),h(i,j+1,k)));
            Real Heff2 = amrex::Math::abs(vh(i,j-1,k) * IdxCv(i,j-1,0)) / (eps_vel+amrex::Math::abs(v(i,j-1,k)));
            Heff2 = amrex::max(Heff2, amrex::min(h(i,j-1,k),h(i,j,k)));
            Heff2 = amrex::min(Heff2, amrex::max(h(i,j-1,k),h(i,j,k)));
            Real Heff3 = amrex::Math::abs(vh(i+1,j,k) * IdxCv(i+1,j,0)) / (eps_vel+amrex::Math::abs(v(i+1,j,k)));
            Heff3 = amrex::max(Heff3, amrex::min(h(i+1,j,k),h(i+1,j+1,k)));
            Heff3 = amrex::min(Heff3, amrex::max(h(i+1,j,k),h(i+1,j+1,k)));
            Real Heff4 = amrex::Math::abs(vh(i+1,j-1,k) * IdxCv(i+1,j-1,0)) / (eps_vel+amrex::Math::abs(v(i+1,j-1,k)));
            Heff4 = amrex::max(Heff4, amrex::min(h(i+1,j-1,k),h(i+1,j,k)));
            Heff4 = amrex::min(Heff4, amrex::max(h(i+1,j-1,k),h(i+1,j,k)));
            CAu(i,j,k) = 0.5*(abs_vort(i,j,k)+abs_vort(i,j-1,k)) *
                         ((vh(i,j,k) + vh(i+1,j-1,k)) + (vh(i,j-1,k) + vh(i+1,j,k)) ) /
                         (h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) * IdxCu(i,j,0);
        });
    } else if (PV_Adv_Scheme == PV_ADV_UPWIND1) {
        ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Temporary effective H at U or V points [H ~> m or kg m-2].
            Real Heff1 = amrex::Math::abs(vh(i,j,k) * IdxCv(i,j,0)) / (eps_vel+amrex::Math::abs(v(i,j,k)));
            Heff1 = amrex::max(Heff1, amrex::min(h(i,j,k),h(i,j+1,k)));
            Heff1 = amrex::min(Heff1, amrex::max(h(i,j,k),h(i,j+1,k)));
            Real Heff2 = amrex::Math::abs(vh(i,j-1,k) * IdxCv(i,j-1,0)) / (eps_vel+amrex::Math::abs(v(i,j-1,k)));
            Heff2 = amrex::max(Heff2, amrex::min(h(i,j-1,k),h(i,j,k)));
            Heff2 = amrex::min(Heff2, amrex::max(h(i,j-1,k),h(i,j,k)));
            Real Heff3 = amrex::Math::abs(vh(i+1,j,k) * IdxCv(i+1,j,0)) / (eps_vel+amrex::Math::abs(v(i+1,j,k)));
            Heff3 = amrex::max(Heff3, amrex::min(h(i+1,j,k),h(i+1,j+1,k)));
            Heff3 = amrex::min(Heff3, amrex::max(h(i+1,j,k),h(i+1,j+1,k)));
            Real Heff4 = amrex::Math::abs(vh(i+1,j-1,k) * IdxCv(i+1,j-1,0)) / (eps_vel+amrex::Math::abs(v(i+1,j-1,k)));
            Heff4 = amrex::max(Heff4, amrex::min(h(i+1,j-1,k),h(i+1,j,k)));
            Heff4 = amrex::min(Heff4, amrex::max(h(i+1,j-1,k),h(i+1,j,k)));
            // More temporary variables [H L2 T-1 ~> m3 s-1 or kg s-1].
            Real const VHeff = ((vh(i,j,k) + vh(i+1,j-1,k)) + (vh(i,j-1,k) + vh(i+1,j,k)) );
            // More temporary variables [H L2 T-2 ~> m3 s-2 or kg s-2].
            Real const QVHeff = 0.5*( ((abs_vort(i,j,k)+abs_vort(i,j-1,k))*VHeff)
                                    - ((abs_vort(i,j,k)-abs_vort(i,j-1,k))*amrex::Math::abs(VHeff)) );
            CAu(i,j,k) = (QVHeff / ( h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) ) * IdxCu(i,j,0);
        });
    }

    // Calculate the tendencies of meridional velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAv = - q * uh - d(KE)/dy.
    // An enstrophy conserving scheme robust to vanishing layers
    // Note: Heffs are in lieu of h_at_u that should be returned by the
    //       continuity solver. AJA
    if (PV_Adv_Scheme == PV_ADV_CENTERED) {
        ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Temporary effective H at U or V points [H ~> m or kg m-2].
            Real Heff1 = amrex::Math::abs(uh(i,j,k) * IdyCu(i,j,0)) / (eps_vel+amrex::Math::abs(u(i,j,k)));
            Heff1 = amrex::max(Heff1, amrex::min(h(i,j,k),h(i+1,j,k)));
            Heff1 = amrex::min(Heff1, amrex::max(h(i,j,k),h(i+1,j,k)));
            Real Heff2 = amrex::Math::abs(uh(i-1,j,k) * IdyCu(i-1,j,0)) / (eps_vel+amrex::Math::abs(u(i-1,j,k)));
            Heff2 = amrex::max(Heff2, amrex::min(h(i-1,j,k),h(i,j,k)));
            Heff2 = amrex::min(Heff2, amrex::max(h(i-1,j,k),h(i,j,k)));
            Real Heff3 = amrex::Math::abs(uh(i,j+1,k) * IdyCu(i,j+1,0)) / (eps_vel+amrex::Math::abs(u(i,j+1,k)));
            Heff3 = amrex::max(Heff3, amrex::min(h(i,j+1,k),h(i+1,j+1,k)));
            Heff3 = amrex::min(Heff3, amrex::max(h(i,j+1,k),h(i+1,j+1,k)));
            Real Heff4 = amrex::Math::abs(uh(i-1,j+1,k) * IdyCu(i-1,j+1,0)) / (eps_vel+amrex::Math::abs(u(i-1,j+1,k)));
            Heff4 = amrex::max(Heff4, amrex::min(h(i-1,j+1,k),h(i,j+1,k)));
            Heff4 = amrex::min(Heff4, amrex::max(h(i-1,j+1,k),h(i,j+1,k)));
            CAv(i,j,k) = - 0.5*(abs_vort(i,j,k)+abs_vort(i-1,j,k)) *
                           ((uh(i  ,j  ,k)+uh(i-1,j+1,k)) +
                            (uh(i-1,j  ,k)+uh(i  ,j+1,k)) ) /
                        (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * IdyCv(i,j,0);
        });
    } else if (PV_Adv_Scheme == PV_ADV_UPWIND1) {
        ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Temporary effective H at U or V points [H ~> m or kg m-2].
            Real Heff1 = amrex::Math::abs(uh(i,j,k) * IdyCu(i,j,0)) / (eps_vel+amrex::Math::abs(u(i,j,k)));
            Heff1 = amrex::max(Heff1, amrex::min(h(i,j,k),h(i+1,j,k)));
            Heff1 = amrex::min(Heff1, amrex::max(h(i,j,k),h(i+1,j,k)));
            Real Heff2 = amrex::Math::abs(uh(i-1,j,k) * IdyCu(i-1,j,0)) / (eps_vel+amrex::Math::abs(u(i-1,j,k)));
            Heff2 = amrex::max(Heff2, amrex::min(h(i-1,j,k),h(i,j,k)));
            Heff2 = amrex::min(Heff2, amrex::max(h(i-1,j,k),h(i,j,k)));
            Real Heff3 = amrex::Math::abs(uh(i,j+1,k) * IdyCu(i,j+1,0)) / (eps_vel+amrex::Math::abs(u(i,j+1,k)));
            Heff3 = amrex::max(Heff3, amrex::min(h(i,j+1,k),h(i+1,j+1,k)));
            Heff3 = amrex::min(Heff3, amrex::max(h(i,j+1,k),h(i+1,j+1,k)));
            Real Heff4 = amrex::Math::abs(uh(i-1,j+1,k) * IdyCu(i-1,j+1,0)) / (eps_vel+amrex::Math::abs(u(i-1,j+1,k)));
            Heff4 = amrex::max(Heff4, amrex::min(h(i-1,j+1,k),h(i,j+1,k)));
            Heff4 = amrex::min(Heff4, amrex::max(h(i-1,j+1,k),h(i,j+1,k)));
            // More temporary variables [H L2 T-1 ~> m3 s-1 or kg s-1].
            Real const UHeff = ((uh(i  ,j  ,k)+uh(i-1,j+1,k)) +
                                (uh(i-1,j  ,k)+uh(i  ,j+1,k)) );
            // More temporary variables [H L2 T-2 ~> m3 s-2 or kg s-2].
            Real const QUHeff = 0.5*( ((abs_vort(i,j,k)+abs_vort(i-1,j,k))*UHeff)
                                    - ((abs_vort(i,j,k)-abs_vort(i-1,j,k))*amrex::Math::abs(UHeff)) );
            CAv(i,j,k) = - QUHeff /
                         (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * IdyCv(i,j,0);
        });
    }
}

void coradv_sadourny(const Box& bxH,
                     const Box& bxU,
                     const Box& bxV,
                     Array4<const Real> const& u,
                     Array4<const Real> const& v,
                     Array4<const Real> const& uh,
                     Array4<const Real> const& vh,
                     Array4<const Real> const& q,
                     Array4<const Real> const& uh_center,
                     Array4<const Real> const& vh_center,
                     Array4<const Real> const& dy_Cu,
                     Array4<const Real> const& dx_Cv,
                     Array4<const Real> const& IdxCu,
                     Array4<const Real> const& IdyCv,
                     Array4<Real>       const& CAu,
                     Array4<Real>       const& CAv,
                     const Coriolis_scheme_CS_C& CS)
{
    BL_PROFILE("coradv_sadourny");

    const int ksc = bxH.smallEnd(2);
    const int kec = bxH.bigEnd(2);

    // bxU grown by one point at the end in i and the start in j,
    // [IscB:IecB+1, jsc-1:jec, ksc:kec]
    Box bxUx = bxU;
    bxUx.growHi(0, 1).growLo(1, 1);
    // bxV grown by one point at the start in i and the end in j,
    // [isc-1:iec, JscB:JecB+1, ksc:kec]
    Box bxVx = bxV;
    bxVx.growLo(0, 1).growHi(1, 1);

    const Box minmax_box(IntVect(Box(v).smallEnd(0), Box(u).smallEnd(1), ksc),
                         IntVect(Box(v).bigEnd(0), Box(u).bigEnd(1), kec));

    // The smallest and largest estimates of the zonal volume fluxes through
    // the faces (i.e. u*h*dy) [H L2 T-1 ~> m3 s-1 or kg s-1]
    FArrayBox uh_min_fab(minmax_box, 1, The_Async_Arena());
    FArrayBox uh_max_fab(minmax_box, 1, The_Async_Arena());
    // The smallest and largest estimates of the meridional volume fluxes through
    // the faces (i.e. v*h*dx) [H L2 T-1 ~> m3 s-1 or kg s-1]
    FArrayBox vh_min_fab(minmax_box, 1, The_Async_Arena());
    FArrayBox vh_max_fab(minmax_box, 1, The_Async_Arena());

    Array4<Real> const uh_min = uh_min_fab.array();
    Array4<Real> const uh_max = uh_max_fab.array();
    Array4<Real> const vh_min = vh_min_fab.array();
    Array4<Real> const vh_max = vh_max_fab.array();

    // The energy-dissipating biased limiter (Coriolis_En_Dis) is only used by these schemes.
    // .and. SADOURNEY75_ENERGY ??
    if (CS.Coriolis_En_Dis) {
        //  c1 = 1.0-1.5*RANGE ; c2 = 1.0-RANGE ; c3 = 2.0 ; slope = 0.5
        // Nondimensional parameters for the Coriolis limiter scheme [nondim]
        constexpr Real c1 = 1.0-1.5*0.5;
        constexpr Real c2 = 1.0-0.5;
        constexpr Real c3 = 2.0;
        constexpr Real slope = 0.5;

        ParallelFor(bxVx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Centered and input estimates of uh [H L2 T-1 ~> m3 s-1 or kg s-1].
            Real uhc = uh_center(i,j,k);
            Real uhm = uh(i,j,k);
            // This sometimes matters with some types of open boundary conditions.
            if (dy_Cu(i,j,0) == 0.0) uhc = uhm;

            if (amrex::Math::abs(uhc) < 0.1*amrex::Math::abs(uhm)) {
                uhm = 10.0*uhc;
            } else if (amrex::Math::abs(uhc) > c1*amrex::Math::abs(uhm)) {
                if (amrex::Math::abs(uhc) < c2*amrex::Math::abs(uhm)) { uhc = (3.0*uhc+(1.0-c2*3.0)*uhm); }
                else if (amrex::Math::abs(uhc) <= c3*amrex::Math::abs(uhm)) { uhc = uhm; }
                else { uhc = slope*uhc+(1.0-c3*slope)*uhm; }
            }

            if (uhc > uhm) {
                uh_min(i,j,k) = uhm; uh_max(i,j,k) = uhc;
            } else {
                uh_max(i,j,k) = uhm; uh_min(i,j,k) = uhc;
            }
        });

        ParallelFor(bxUx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Centered and input estimates of vh [H L2 T-1 ~> m3 s-1 or kg s-1].
            Real vhc = vh_center(i,j,k);
            Real vhm = vh(i,j,k);
            // This sometimes matters with some types of open boundary conditions.
            if (dx_Cv(i,j,0) == 0.0) vhc = vhm;

            if (amrex::Math::abs(vhc) < 0.1*amrex::Math::abs(vhm)) {
                vhm = 10.0*vhc;
            } else if (amrex::Math::abs(vhc) > c1*amrex::Math::abs(vhm)) {
                if (amrex::Math::abs(vhc) < c2*amrex::Math::abs(vhm)) { vhc = (3.0*vhc+(1.0-c2*3.0)*vhm); }
                else if (amrex::Math::abs(vhc) <= c3*amrex::Math::abs(vhm)) { vhc = vhm; }
                else { vhc = slope*vhc+(1.0-c3*slope)*vhm; }
            }

            if (vhc > vhm) {
                vh_min(i,j,k) = vhm; vh_max(i,j,k) = vhc;
            } else {
                vh_max(i,j,k) = vhm; vh_min(i,j,k) = vhc;
            }
        });
    }

    // Calculate the tendencies of zonal velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAu =  q * vh - d(KE)/dx.
    if (CS.Coriolis_Scheme == SADOURNY75_ENERGY) {
        if (CS.Coriolis_En_Dis) {
            // Energy dissipating biased scheme, Hallberg 200x
            ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Temporary variables [L2 T-2 ~> m2 s-2].
                Real temp1, temp2;
                if (q(i,j,k)*u(i,j,k) == 0.0) {
                    temp1 = q(i,j,k) * ( (vh_max(i,j,k)+vh_max(i+1,j,k))
                                       + (vh_min(i,j,k)+vh_min(i+1,j,k)) )*0.5;
                } else if (q(i,j,k)*u(i,j,k) < 0.0) {
                    temp1 = q(i,j,k) * (vh_max(i,j,k)+vh_max(i+1,j,k));
                } else {
                    temp1 = q(i,j,k) * (vh_min(i,j,k)+vh_min(i+1,j,k));
                }
                if (q(i,j-1,k)*u(i,j,k) == 0.0) {
                    temp2 = q(i,j-1,k) * ( (vh_max(i,j-1,k)+vh_max(i+1,j-1,k))
                                         + (vh_min(i,j-1,k)+vh_min(i+1,j-1,k)) )*0.5;
                } else if (q(i,j-1,k)*u(i,j,k) < 0.0) {
                    temp2 = q(i,j-1,k) * (vh_max(i,j-1,k)+vh_max(i+1,j-1,k));
                } else {
                    temp2 = q(i,j-1,k) * (vh_min(i,j-1,k)+vh_min(i+1,j-1,k));
                }
                CAu(i,j,k) = 0.25 * IdxCu(i,j,0) * (temp1 + temp2);
            });
        } else {
            // Energy conserving scheme, Sadourny 1975
            ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                CAu(i,j,k) = 0.25 *
                  ((q(i,j,k) * (vh(i+1,j,k) + vh(i,j,k))) +
                   (q(i,j-1,k) * (vh(i,j-1,k) + vh(i+1,j-1,k)))) * IdxCu(i,j,0);
            });
        }
    } else if (CS.Coriolis_Scheme == SADOURNY75_ENSTRO) {
        ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            CAu(i,j,k) = 0.125 * (IdxCu(i,j,0) * (q(i,j,k) + q(i,j-1,k))) *
                         ((vh(i+1,j,k) + vh(i,j,k)) + (vh(i,j-1,k) + vh(i+1,j-1,k)));
        });
    }

    // Calculate the tendencies of meridional velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAv = - q * uh - d(KE)/dy.
    if (CS.Coriolis_Scheme == SADOURNY75_ENERGY) {
        if (CS.Coriolis_En_Dis) {
            // Energy dissipating biased scheme, Hallberg 200x
            ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Temporary variables [L2 T-2 ~> m2 s-2].
                Real temp1, temp2;
                if (q(i-1,j,k)*v(i,j,k) == 0.0) {
                    temp1 = q(i-1,j,k) * ( (uh_max(i-1,j,k)+uh_max(i-1,j+1,k))
                                         + (uh_min(i-1,j,k)+uh_min(i-1,j+1,k)) )*0.5;
                } else if (q(i-1,j,k)*v(i,j,k) > 0.0) {
                    temp1 = q(i-1,j,k) * (uh_max(i-1,j,k)+uh_max(i-1,j+1,k));
                } else {
                    temp1 = q(i-1,j,k) * (uh_min(i-1,j,k)+uh_min(i-1,j+1,k));
                }
                if (q(i,j,k)*v(i,j,k) == 0.0) {
                    temp2 = q(i,j,k) * ( (uh_max(i,j,k)+uh_max(i,j+1,k))
                                       + (uh_min(i,j,k)+uh_min(i,j+1,k)) )*0.5;
                } else if (q(i,j,k)*v(i,j,k) > 0.0) {
                    temp2 = q(i,j,k) * (uh_max(i,j,k)+uh_max(i,j+1,k));
                } else {
                    temp2 = q(i,j,k) * (uh_min(i,j,k)+uh_min(i,j+1,k));
                }
                CAv(i,j,k) = -0.25 * IdyCv(i,j,0) * (temp1 + temp2);
            });
        } else {
            // Energy conserving scheme, Sadourny 1975
            ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                CAv(i,j,k) = - 0.25*
                    ((q(i-1,j,k)*(uh(i-1,j,k) + uh(i-1,j+1,k))) +
                     (q(i,j,k)*(uh(i,j,k) + uh(i,j+1,k)))) * IdyCv(i,j,0);
            });
        }
    } else if (CS.Coriolis_Scheme == SADOURNY75_ENSTRO) {
        ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            CAv(i,j,k) = -0.125 * (IdyCv(i,j,0) * (q(i-1,j,k) + q(i,j,k))) *
                         ((uh(i-1,j,k) + uh(i-1,j+1,k)) + (uh(i,j,k) + uh(i,j+1,k)));
        });
    }
}

namespace {

/// The zonal Coriolis and momentum advection acceleration at one u point with the WENO PV
/// reconstruction of order up to MaxOrder (7, 5 or 3).
template <int MaxOrder, bool VelocitySmoothing>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void coradv_weno_cau_cell(int i, int j, int k,
                          Array4<const Real> const& u,
                          Array4<const Real> const& vh,
                          Array4<const Real> const& q,
                          Array4<const Real> const& abs_vort,
                          Array4<const Real> const& h_q,
                          Array4<const Real> const& IdxCu,
                          Array4<const Real> const& mask2dCu,
                          Real const H_subroundoff,
                          Array4<Real> const& CAu) noexcept
{
    // v_u is the v velocity at u point [L T-1 ~> m s-1]
    Real const v_u = 0.25*IdxCu(i,j,0)*((vh(i+1,j,k) + vh(i,j,k)) + (vh(i,j-1,k) + vh(i+1,j-1,k)));
    // PV at the u point [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
    Real q_u;
    // Order of accuracy for the WENO calculations
    int fifth_order = 0, seventh_order = 0;

    // check whether there is masked land points in the stencil
    int const third_order = static_cast<int>(mask2dCu(i,j-2,0) * mask2dCu(i,j-1,0) * mask2dCu(i,j,0) *
                                             mask2dCu(i,j+1,0) * mask2dCu(i,j+2,0));
    if constexpr (MaxOrder >= 5) {
        fifth_order   = static_cast<int>(third_order * mask2dCu(i,j-3,0) * mask2dCu(i,j+3,0));
    }
    if constexpr (MaxOrder >= 7) {
        seventh_order = static_cast<int>(fifth_order * mask2dCu(i,j-4,0) * mask2dCu(i,j+4,0));
    }

    // compute the masking to make sure that inland values are not used
    if (seventh_order == 1) {
        // all values are valid, we use seventh order reconstruction
        Real q8[8], h8[8];
        for (int n = 0; n < 8; ++n) { q8[n] = abs_vort(i,j-4+n,k); h8[n] = h_q(i,j-4+n,k); }
        if constexpr (VelocitySmoothing) {
            // Eight-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
            Real u_q8[8];
            for (int n = 0; n < 8; ++n) { u_q8[n] = (u(i,j-4+n,k) + u(i,j-3+n,k)) * 0.5; }
            weno_seven_h_weight_reconstruction_point(q8, h8, u_q8, H_subroundoff, v_u, q_u);
        } else {
            weno_seven_h_weight_reconstruction_point(q8, h8, q8, H_subroundoff, v_u, q_u);
        }
        CAu(i,j,k) = (q_u * v_u);

    } else if (fifth_order == 1) {
        // all values are valid, we use fifth order reconstruction
        Real q6[6], h6[6];
        for (int n = 0; n < 6; ++n) { q6[n] = abs_vort(i,j-3+n,k); h6[n] = h_q(i,j-3+n,k); }
        if constexpr (VelocitySmoothing) {
            // Six-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
            Real u_q6[6];
            for (int n = 0; n < 6; ++n) { u_q6[n] = (u(i,j-3+n,k) + u(i,j-2+n,k)) * 0.5; }
            weno_five_h_weight_reconstruction_point(q6, h6, u_q6, H_subroundoff, v_u, q_u);
        } else {
            weno_five_h_weight_reconstruction_point(q6, h6, q6, H_subroundoff, v_u, q_u);
        }
        CAu(i,j,k) = (q_u * v_u);

    } else if (third_order == 1) {
        // only the middle values are valid, we use third order reconstruction
        Real q4[4], h4[4];
        for (int n = 0; n < 4; ++n) { q4[n] = abs_vort(i,j-2+n,k); h4[n] = h_q(i,j-2+n,k); }
        if constexpr (VelocitySmoothing) {
            // Four-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
            Real u_q4[4];
            for (int n = 0; n < 4; ++n) { u_q4[n] = (u(i,j-2+n,k) + u(i,j-1+n,k)) * 0.5; }
            weno_three_h_weight_reconstruction_point(q4, h4, u_q4, H_subroundoff, v_u, q_u);
        } else {
            weno_three_h_weight_reconstruction_point(q4, h4, q4, H_subroundoff, v_u, q_u);
        }
        CAu(i,j,k) = (q_u * v_u);

    } else { // Upwind first order
        if (v_u > 0.0) {
            q_u = q(i,j-1,k);
        } else {
            q_u = q(i,j,k);
        }
        CAu(i,j,k) = (q_u * v_u);
    }
}

/// The meridional Coriolis and momentum advection acceleration at one v point with the WENO PV
/// reconstruction of order up to MaxOrder (7, 5 or 3).
template <int MaxOrder, bool VelocitySmoothing>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void coradv_weno_cav_cell(int i, int j, int k,
                          Array4<const Real> const& v,
                          Array4<const Real> const& uh,
                          Array4<const Real> const& q,
                          Array4<const Real> const& abs_vort,
                          Array4<const Real> const& h_q,
                          Array4<const Real> const& IdyCv,
                          Array4<const Real> const& mask2dCv,
                          Real const H_subroundoff,
                          Array4<Real> const& CAv) noexcept
{
    // u_v is the u velocity at v point [L T-1 ~> m s-1]
    Real const u_v = 0.25*IdyCv(i,j,0)*((uh(i-1,j,k) + uh(i-1,j+1,k)) + (uh(i,j,k) + uh(i,j+1,k)));
    // PV at the v point [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
    Real q_v;
    // Order of accuracy for the WENO calculations
    int fifth_order = 0, seventh_order = 0;

    // check whether there is any masked land values within the stencils
    int const third_order = static_cast<int>(mask2dCv(i-2,j,0) * mask2dCv(i-1,j,0) * mask2dCv(i,j,0) *
                                             mask2dCv(i+1,j,0) * mask2dCv(i+2,j,0));
    if constexpr (MaxOrder >= 5) {
        fifth_order   = static_cast<int>(third_order * mask2dCv(i-3,j,0) * mask2dCv(i+3,j,0));
    }
    if constexpr (MaxOrder >= 7) {
        seventh_order = static_cast<int>(fifth_order * mask2dCv(i-4,j,0) * mask2dCv(i+4,j,0));
    }

    // compute the masking to make sure that inland values are not used
    if (seventh_order == 1) {
        // all values are valid, we use seventh order reconstruction
        Real q8[8], h8[8];
        for (int n = 0; n < 8; ++n) { q8[n] = abs_vort(i-4+n,j,k); h8[n] = h_q(i-4+n,j,k); }
        if constexpr (VelocitySmoothing) {
            // Eight-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
            Real v_q8[8];
            for (int n = 0; n < 8; ++n) { v_q8[n] = (v(i-4+n,j,k) + v(i-3+n,j,k)) * 0.5; }
            weno_seven_h_weight_reconstruction_point(q8, h8, v_q8, H_subroundoff, u_v, q_v);
        } else {
            weno_seven_h_weight_reconstruction_point(q8, h8, q8, H_subroundoff, u_v, q_v);
        }
        CAv(i,j,k) = - (q_v * u_v);

    } else if (fifth_order == 1) {
        // all values are valid, we use fifth order reconstruction
        Real q6[6], h6[6];
        for (int n = 0; n < 6; ++n) { q6[n] = abs_vort(i-3+n,j,k); h6[n] = h_q(i-3+n,j,k); }
        if constexpr (VelocitySmoothing) {
            // Six-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
            Real v_q6[6];
            for (int n = 0; n < 6; ++n) { v_q6[n] = (v(i-3+n,j,k) + v(i-2+n,j,k)) * 0.5; }
            weno_five_h_weight_reconstruction_point(q6, h6, v_q6, H_subroundoff, u_v, q_v);
        } else {
            weno_five_h_weight_reconstruction_point(q6, h6, q6, H_subroundoff, u_v, q_v);
        }
        CAv(i,j,k) = - (q_v * u_v);

    } else if (third_order == 1) {
        // only the middle values are valid, we use third order reconstruction
        Real q4[4], h4[4];
        for (int n = 0; n < 4; ++n) { q4[n] = abs_vort(i-2+n,j,k); h4[n] = h_q(i-2+n,j,k); }
        if constexpr (VelocitySmoothing) {
            // Four-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
            Real v_q4[4];
            for (int n = 0; n < 4; ++n) { v_q4[n] = (v(i-2+n,j,k) + v(i-1+n,j,k)) * 0.5; }
            weno_three_h_weight_reconstruction_point(q4, h4, v_q4, H_subroundoff, u_v, q_v);
        } else {
            weno_three_h_weight_reconstruction_point(q4, h4, q4, H_subroundoff, u_v, q_v);
        }
        CAv(i,j,k) = - (q_v * u_v);

    } else { // Upwind first order!
        if (u_v > 0.0) {
            q_v = q(i-1,j,k);
        } else {
            q_v = q(i,j,k);
        }
        CAv(i,j,k) = - (q_v * u_v);
    }
}

/// Launches the zonal and meridional WENO cell updates over bxU and bxV.
template <int MaxOrder, bool VelocitySmoothing>
void coradv_weno_launch(const Box& bxU,
                        const Box& bxV,
                        Array4<const Real> const& u,
                        Array4<const Real> const& v,
                        Array4<const Real> const& uh,
                        Array4<const Real> const& vh,
                        Array4<const Real> const& q,
                        Array4<const Real> const& abs_vort,
                        Array4<const Real> const& h_q,
                        Array4<const Real> const& IdxCu,
                        Array4<const Real> const& mask2dCu,
                        Array4<const Real> const& IdyCv,
                        Array4<const Real> const& mask2dCv,
                        Real const H_subroundoff,
                        Array4<Real> const& CAu,
                        Array4<Real> const& CAv)
{
    // Calculate the tendencies of zonal velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAu =  q * vh - d(KE)/dx.
    ParallelFor(bxU, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        coradv_weno_cau_cell<MaxOrder, VelocitySmoothing>(i, j, k, u, vh, q, abs_vort, h_q,
                                                          IdxCu, mask2dCu, H_subroundoff, CAu);
    });

    // Calculate the tendencies of meridional velocity due to the Coriolis
    // force and momentum advection.  On a Cartesian grid, this is
    //     CAv = - q * uh - d(KE)/dy.
    ParallelFor(bxV, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        coradv_weno_cav_cell<MaxOrder, VelocitySmoothing>(i, j, k, v, uh, q, abs_vort, h_q,
                                                          IdyCv, mask2dCv, H_subroundoff, CAv);
    });
}

} // namespace

void coradv_weno(const Box& bxU,
                 const Box& bxV,
                 Array4<const Real> const& u,
                 Array4<const Real> const& v,
                 Array4<const Real> const& uh,
                 Array4<const Real> const& vh,
                 Array4<const Real> const& q,
                 Array4<const Real> const& abs_vort,
                 Array4<const Real> const& h_q,
                 Array4<const Real> const& IdxCu,
                 Array4<const Real> const& mask2dCu,
                 Array4<const Real> const& IdyCv,
                 Array4<const Real> const& mask2dCv,
                 const Real H_subroundoff,
                 Array4<Real>       const& CAu,
                 Array4<Real>       const& CAv,
                 const bool weno_velocity_smooth,
                 const Coriolis_scheme_CS_C& CS)
{
    BL_PROFILE("coradv_weno");

    if (CS.Coriolis_Scheme == wenovi7th_PV_ENSTRO) {
        if (weno_velocity_smooth) {
            coradv_weno_launch<7, true>(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q,
                                        IdxCu, mask2dCu, IdyCv, mask2dCv, H_subroundoff, CAu, CAv);
        } else {
            coradv_weno_launch<7, false>(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q,
                                         IdxCu, mask2dCu, IdyCv, mask2dCv, H_subroundoff, CAu, CAv);
        }
    } else if (CS.Coriolis_Scheme == wenovi5th_PV_ENSTRO) {
        if (weno_velocity_smooth) {
            coradv_weno_launch<5, true>(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q,
                                        IdxCu, mask2dCu, IdyCv, mask2dCv, H_subroundoff, CAu, CAv);
        } else {
            coradv_weno_launch<5, false>(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q,
                                         IdxCu, mask2dCu, IdyCv, mask2dCv, H_subroundoff, CAu, CAv);
        }
    } else if (CS.Coriolis_Scheme == wenovi3rd_PV_ENSTRO) {
        if (weno_velocity_smooth) {
            coradv_weno_launch<3, true>(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q,
                                        IdxCu, mask2dCu, IdyCv, mask2dCv, H_subroundoff, CAu, CAv);
        } else {
            coradv_weno_launch<3, false>(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q,
                                         IdxCu, mask2dCu, IdyCv, mask2dCv, H_subroundoff, CAu, CAv);
        }
    }
}

namespace {

/// The halo stencil size for the Coriolis advection scheme
int coriolis_scheme_stencil(const Coriolis_scheme_CS_C& CS)
{
    int stencil = 2;
    if (CS.Coriolis_Scheme == wenovi7th_PV_ENSTRO) stencil = 4;
    if (CS.Coriolis_Scheme == wenovi5th_PV_ENSTRO) stencil = 3;
    return stencil;
}

/// A view of fab with its k range set to [k_start, k_end]
Array4<Real> rebound_k(FArrayBox& fab, int k_start, int k_end)
{
    const Box& b = fab.box();
    return Array4<Real>(fab.dataPtr(),
                        Dim3{b.smallEnd(0), b.smallEnd(1), k_start},
                        Dim3{b.bigEnd(0)+1, b.bigEnd(1)+1, k_end+1}, 1);
}

} // namespace

void coradcalc_tr(const Box& bxH0,
                  const Box& bxQ0,
                  Array4<const Real> const& u,
                  Array4<const Real> const& v,
                  Array4<const Real> const& h,
                  Array4<const Real> const& uh,
                  Array4<const Real> const& vh,
                  Array4<Real>       const& CAu,
                  Array4<Real>       const& CAv,
                  Array4<Real>       const& RV,
                  Array4<Real>       const& PV,
                  Array4<Real>       const& CAuS,
                  Array4<Real>       const& CAvS,
                  OceanOBC* OBC,
                  Array4<Real>       const& gradKEu,
                  Array4<Real>       const& gradKEv,
                  Array4<Real>       const& rv_x_v,
                  Array4<Real>       const& rv_x_u,
                  Array4<const Real> const& mask2dT,
                  Array4<const Real> const& areaT,
                  Array4<const Real> const& areaCu,
                  Array4<const Real> const& areaCv,
                  Array4<const Real> const& IareaT,
                  Array4<const Real> const& mask2dCu,
                  Array4<const Real> const& mask2dCv,
                  Array4<const Real> const& dyCv,
                  Array4<const Real> const& dxCu,
                  Array4<const Real> const& dy_Cu,
                  Array4<const Real> const& dx_Cv,
                  Array4<const Real> const& dyBu,
                  Array4<const Real> const& dxBu,
                  Array4<const Real> const& mask2dBu,
                  Array4<const Real> const& IareaBu,
                  Array4<const Real> const& CoriolisBu,
                  Array4<const Real> const& IdxCu,
                  Array4<const Real> const& IdyCv,
                  Array4<const Real> const& IdxCv,
                  Array4<const Real> const& IdyCu,
                  Array4<const Real> const& IdxCu_OBCmask,
                  Array4<const Real> const& IdyCv_OBCmask,
                  const int nz,
                  const Real H_subroundoff,
                  const Real Angstrom_H,
                  const Real m_to_L,
                  const Real m_s_to_L_T,
                  const CorAdv_opts_CS_C& CS,
                  Array4<const Real> const& por_face_areaU,
                  Array4<const Real> const& por_face_areaV,
                  WaveParametersCS* Waves,
                  const bool Stokes_VF)
{
    BL_PROFILE("coradcalc_tr");

    // NOTE: OBC support temporarily disabled.
    // OceanOBC is forward-declared only.
    if (OBC != nullptr) {
        AMREX_ABORT_LOC("OBC pointer provided but not yet implemented");
    }

    // NOTE: Stokes drift support temporarily disabled.
    // WaveParametersCS is forward-declared only.
    if (Waves != nullptr) {
        AMREX_ABORT_LOC("Waves pointer provided but not yet implemented");
    }

    // To work, the following fields must be set outside of the usual
    // is to ie range before this subroutine is called:
    //   v(is-1:ie+2,js-1:je+1), u(is-1:ie+1,js-1:je+2), h(is-1:ie+2,js-1:je+2),
    //   uh(is-1,ie,js:je+1) and vh(is:ie+1,js-1:je).

    const int is = bxH0.smallEnd(0);
    const int ie = bxH0.bigEnd(0);
    const int js = bxH0.smallEnd(1);
    const int je = bxH0.bigEnd(1);
    const int nkblock = (CS.nkblock == 0) ? nz : CS.nkblock;
    const Real neglect_len = 1e-4 * m_to_L;
    // A volume so small that is expected to be
    // lost in roundoff [H L2 ~> m3 or kg].
    const Real vol_neglect = H_subroundoff * (neglect_len*neglect_len);
    // An area so small that is expected to be
    // lost in roundoff [L2 ~> m2].
    const Real area_neglect = (neglect_len*neglect_len);
    // A tiny, positive velocity [L T-1 ~> m s-1].
    const Real eps_vel = 1.0e-10*m_s_to_L_T;
    // A very small thickness [H ~> m or kg m-2].
    const Real h_tiny = Angstrom_H;  // Perhaps this should be set to h_neglect instead.

    // Stencil size of WENO scheme
    const int stencil = coriolis_scheme_stencil(CS.Coriolis_scheme_CS);

    // True if using one of the WENO schemes
    const bool use_weno = CS.Coriolis_scheme_CS.Coriolis_Scheme == wenovi7th_PV_ENSTRO
        || CS.Coriolis_scheme_CS.Coriolis_Scheme == wenovi5th_PV_ENSTRO
        || CS.Coriolis_scheme_CS.Coriolis_Scheme == wenovi3rd_PV_ENSTRO;

    // The scheme-dependent range of values at which vorticity is set.
    int Is_q, Ie_q, Js_q, Je_q;
    if (use_weno) {
        Is_q = is - stencil; Ie_q = ie + stencil - 1; Js_q = js - stencil; Je_q = je + stencil - 1;
    } else {
        Is_q = bxQ0.smallEnd(0) - 1; Ie_q = bxQ0.bigEnd(0) + 1;
        Js_q = bxQ0.smallEnd(1) - 1; Je_q = bxQ0.bigEnd(1) + 1;
    }

    // The sum of the ocean areas at the 4 adjacent thickness points [L2 ~> m2].
    FArrayBox Area_q_fab(Box(IntVect(Box(u).smallEnd(0), Box(v).smallEnd(1), 0), IntVect(Box(u).bigEnd(0), Box(v).bigEnd(1), 0)),
                         1, The_Async_Arena());
    // The ocean area at h points [L2 ~> m2].  Area_h is used to find the
    // average thickness in the denominator of q.  0 for land points.
    FArrayBox Area_h_fab(Box(IntVect(Box(h).smallEnd(0), Box(h).smallEnd(1), 0), IntVect(Box(h).bigEnd(0), Box(h).bigEnd(1), 0)),
                         1, The_Async_Arena());
    Array4<Real> const Area_q = Area_q_fab.array();
    Array4<Real> const Area_h = Area_h_fab.array();

    // The scheme-dependent vorticity-point box of one tile,
    // [Is_q:Ie_q, Js_q:Je_q, k_start:k_end]
    Box bxQs(IntVect(Is_q, Js_q, 0), IntVect(Ie_q, Je_q, nz - 1));
    // bxQs grown by one point at the end in i and j, before the tile loop,
    // [Is_q:Ie_q+1, Js_q:Je_q+1, 1:nz]
    Box bxQs_pij = bxQs;
    bxQs_pij.growHi(0, 1).growHi(1, 1);

    const Box bxQs_pij_2d(IntVect(bxQs_pij.smallEnd(0), bxQs_pij.smallEnd(1), 0),
                          IntVect(bxQs_pij.bigEnd(0), bxQs_pij.bigEnd(1), 0));
    ParallelFor(bxQs_pij_2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        Area_h(i,j,0) = mask2dT(i,j,0) * areaT(i,j,0);
    });

    /*
    if (OBC != nullptr) {
        for (int n = 0; n < OBC->number_of_segments; ++n) {
            auto& segment = OBC->segment[n];
            if (!segment.on_pe) continue;
            const int I = segment.HI.IsdB;
            const int J = segment.HI.JsdB;
            const int direction = segment.direction;
            if (segment.is_N_or_S && (J >= Js_q) && (J <= Je_q)) {
                const Box bxI(IntVect(amrex::max(Is_q, segment.HI.isd), J, 0),
                              IntVect(amrex::min(Ie_q+1, segment.HI.ied), J, 0));
                if (bxI.ok()) {
                    ParallelFor(bxI, [=] AMREX_GPU_DEVICE (int i, int, int) noexcept
                    {
                        if (direction == OBC_DIRECTION_N) {
                            Area_h(i,J+1,0) = Area_h(i,J,0);
                        } else { // (direction == OBC_DIRECTION_S)
                            Area_h(i,J,0) = Area_h(i,J+1,0);
                        }
                    });
                }
            } else if (segment.is_E_or_W && (I >= Is_q) && (I <= Ie_q)) {
                const Box bxJ(IntVect(I, amrex::max(Js_q, segment.HI.jsd), 0),
                              IntVect(I, amrex::min(Je_q+1, segment.HI.jed), 0));
                if (bxJ.ok()) {
                    ParallelFor(bxJ, [=] AMREX_GPU_DEVICE (int, int j, int) noexcept
                    {
                        if (direction == OBC_DIRECTION_E) {
                            Area_h(I+1,j,0) = Area_h(I,j,0);
                        } else { // (direction == OBC_DIRECTION_W)
                            Area_h(I,j,0) = Area_h(I+1,j,0);
                        }
                    });
                }
            }
        }
    }
    */

    const Box bxQs_2d(IntVect(Is_q, Js_q, 0), IntVect(Ie_q, Je_q, 0));
    ParallelFor(bxQs_2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        Area_q(i,j,0) = (Area_h(i,j,0) + Area_h(i+1,j+1,0)) +
                        (Area_h(i+1,j,0) + Area_h(i,j+1,0));
    });

    // Hoist AL_BLEND k-independent scalars out of the tile loop.
    // Temporary variables associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
    Real Fe_m2 = 0.0, rat_lin = 0.0;
    if (CS.Coriolis_scheme_CS.Coriolis_Scheme == AL_BLEND) {
        Fe_m2 = CS.F_eff_max_blend - 2.0;
        rat_lin = 1.5 * Fe_m2 / amrex::max(CS.wt_lin_blend, 1.0e-16);
        if (CS.F_eff_max_blend <= 2.0) { Fe_m2 = -1.; rat_lin = -1.0; }
    }

    // Loop over iteration tiles.  Each tile covers the whole horizontal compute domain and
    // nkblock layers; the tile boxes carry all of the index ranges used within the loop.
    const int nk_tile = amrex::min(nkblock, nz);
    auto tile_box = [nk_tile] (int ilo, int jlo, int ihi, int jhi) {
        return Box(IntVect(ilo, jlo, 0), IntVect(ihi, jhi, nk_tile - 1));
    };
    const Box& Area_q_box = Area_q_fab.box();
    const Box& Area_h_box = Area_h_fab.box();
    const Box q_tile  = tile_box(Area_q_box.smallEnd(0), Area_q_box.smallEnd(1),
                                 Area_q_box.bigEnd(0), Area_q_box.bigEnd(1));
    const Box u_tile  = tile_box(Box(u).smallEnd(0), Box(u).smallEnd(1), Box(u).bigEnd(0), Box(u).bigEnd(1));
    const Box v_tile  = tile_box(Box(v).smallEnd(0), Box(v).smallEnd(1), Box(v).bigEnd(0), Box(v).bigEnd(1));
    const Box h_tile  = tile_box(Area_h_box.smallEnd(0), Area_h_box.smallEnd(1),
                                 Area_h_box.bigEnd(0), Area_h_box.bigEnd(1));

    // Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    FArrayBox q_fab(q_tile, 1, The_Async_Arena());
    // Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    FArrayBox qS_fab(q_tile, 1, The_Async_Arena());
    // The inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    FArrayBox Ih_q_fab(q_tile, 1, The_Async_Arena());
    // The thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    FArrayBox h_q_fab(q_tile, 1, The_Async_Arena());
    // Absolute vorticity at q-points [T-1 ~> s-1].
    FArrayBox abs_vort_fab(q_tile, 1, The_Async_Arena());
    // Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    FArrayBox q2_fab(q_tile, 1, The_Async_Arena());
    // The zonal gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
    // KEx = d/dx KE.
    FArrayBox KEx_fab(u_tile, 1, The_Async_Arena());
    // Transport based on arithmetic mean h at u-points [H L2 T-1 ~> m3 s-1 or kg s-1]
    FArrayBox uh_center_fab(u_tile, 1, The_Async_Arena());
    // Kinetic energy per unit mass [L2 T-2 ~> m2 s-2], KE = (u^2 + v^2)/2.
    FArrayBox KE_fab(h_tile, 1, The_Async_Arena());
    // The meridional gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
    // KEy = d/dy KE.
    FArrayBox KEy_fab(v_tile, 1, The_Async_Arena());
    // Transport based on arithmetic mean h at v-points [H L2 T-1 ~> m3 s-1 or kg s-1]
    FArrayBox vh_center_fab(v_tile, 1, The_Async_Arena());

    const bool rv_x_v_associated = (rv_x_v.p != nullptr);
    const bool rv_x_u_associated = (rv_x_u.p != nullptr);

    for (int k_start = 0; k_start < nz; k_start += nkblock) {
        const int k_end = amrex::min(k_start+nkblock-1, nz-1);
        // The h-point iteration box of one tile, [isc:iec, jsc:jec, k_start:k_end]
        const Box bxH(IntVect(bxH0.smallEnd(0), bxH0.smallEnd(1), k_start),
                      IntVect(bxH0.bigEnd(0), bxH0.bigEnd(1), k_end));
        // The B-grid-index iteration box of one tile, [IscB:IecB, JscB:JecB, k_start:k_end]
        const Box bxQ(IntVect(bxQ0.smallEnd(0), bxQ0.smallEnd(1), k_start),
                      IntVect(bxQ0.bigEnd(0), bxQ0.bigEnd(1), k_end));
        bxQs = Box(IntVect(Is_q, Js_q, k_start), IntVect(Ie_q, Je_q, k_end));
        // The u-point iteration box of one tile, [IscB:IecB, jsc:jec, k_start:k_end]
        const Box bxU(IntVect(bxQ0.smallEnd(0), bxH0.smallEnd(1), k_start),
                      IntVect(bxQ0.bigEnd(0), bxH0.bigEnd(1), k_end));
        // The v-point iteration box of one tile, [isc:iec, JscB:JecB, k_start:k_end]
        const Box bxV(IntVect(bxH0.smallEnd(0), bxQ0.smallEnd(1), k_start),
                      IntVect(bxH0.bigEnd(0), bxQ0.bigEnd(1), k_end));

        Array4<Real> const q         = rebound_k(q_fab, k_start, k_end);
        Array4<Real> const qS        = rebound_k(qS_fab, k_start, k_end);
        Array4<Real> const Ih_q      = rebound_k(Ih_q_fab, k_start, k_end);
        Array4<Real> const h_q       = rebound_k(h_q_fab, k_start, k_end);
        Array4<Real> const abs_vort  = rebound_k(abs_vort_fab, k_start, k_end);
        Array4<Real> const q2        = rebound_k(q2_fab, k_start, k_end);
        Array4<Real> const KEx       = rebound_k(KEx_fab, k_start, k_end);
        Array4<Real> const uh_center = rebound_k(uh_center_fab, k_start, k_end);
        Array4<Real> const KE        = rebound_k(KE_fab, k_start, k_end);
        Array4<Real> const KEy       = rebound_k(KEy_fab, k_start, k_end);
        Array4<Real> const vh_center = rebound_k(vh_center_fab, k_start, k_end);

        // Potential vorticity and the related quantities at q points used by every scheme.
        coradv_setup(bxH, bxQ, bxU, bxV, bxQs, u, v, h, OBC, rv_x_v_associated, rv_x_u_associated,
                     por_face_areaU, por_face_areaV, Waves, dyCv, dxCu, dy_Cu, dx_Cv, dyBu, dxBu,
                     mask2dBu, IareaBu, CoriolisBu, Area_h, Area_q, Stokes_VF, use_weno,
                     vol_neglect, area_neglect, q, Ih_q, abs_vort, h_q, qS, q2,
                     uh_center, vh_center, RV, PV, CS);

        // Calculate KE and the gradient of KE
        gradke(bxQ, bxU, bxV, u, v, KE, KEx, KEy, areaCu, areaCv, IareaT, mask2dCu, mask2dCv,
               IdxCu_OBCmask, IdyCv_OBCmask, CS.KE_Scheme, CS.KE_use_limiter);

        // Calculate the Coriolis and momentum advection accelerations, CAu and CAv, with the selected
        // scheme.  On a Cartesian grid, CAu =  q * vh - d(KE)/dx and CAv = - q * uh - d(KE)/dy.
        switch (CS.Coriolis_scheme_CS.Coriolis_Scheme) {
        case SADOURNY75_ENERGY:
        case SADOURNY75_ENSTRO:
            coradv_sadourny(bxH, bxU, bxV, u, v, uh, vh, q, uh_center, vh_center,
                            dy_Cu, dx_Cv, IdxCu, IdyCv, CAu, CAv, CS.Coriolis_scheme_CS);
            break;
        case ARAKAWA_HSU90:
        case ARAKAWA_LAMB81:
        case AL_BLEND:
            coradv_arakawa(bxH, bxQ, bxU, bxV, uh, vh, q, Ih_q, Fe_m2, rat_lin,
                           IdxCu, IdyCv, CAu, CAv, CS.wt_lin_blend, CS.Coriolis_scheme_CS);
            break;
        case ROBUST_ENSTRO:
            coradv_robust_enstro(bxU, bxV, u, v, h, uh, vh, abs_vort, eps_vel, h_tiny,
                                 IdxCv, IdxCu, IdyCu, IdyCv, CAu, CAv, CS.PV_Adv_Scheme);
            break;
        case wenovi7th_PV_ENSTRO:
        case wenovi5th_PV_ENSTRO:
        case wenovi3rd_PV_ENSTRO:
            coradv_weno(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q,
                        IdxCu, mask2dCu, IdyCv, mask2dCv, H_subroundoff, CAu, CAv,
                        CS.weno_velocity_smooth, CS.Coriolis_scheme_CS);
            break;
        }

        // The Stokes-drift diagnostic, bounding, the kinetic energy gradient, and diagnostics.
        coradv_common_terms(bxU, bxV, u, v, uh, vh, abs_vort, qS, q2, KEx, KEy, Stokes_VF,
                            IdxCu, IdyCv, CAu, CAv, CAuS, CAvS, gradKEu, gradKEv, rv_x_v, rv_x_u, CS);
    } // end of tile loop.
}

}
