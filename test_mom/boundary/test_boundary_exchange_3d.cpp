#include <format>
#include <iostream>

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>

#include "amrex_assertions.hpp"
#include "captured_io.hpp"
#include "data_dir.hpp"

using test_mom::expect_arrays_equal;
using test_mom::to_host_fab;

TEST(BoundaryExchange, BoundaryMatchesAfterOneTimeStep) {
    int myRank = amrex::ParallelDescriptor::MyProc();
    test_mom::CapturedFile captured(test_mom::data_dir /
                                    ("MOM_dynamics_split_RK2_" + std::format("{:04}", myRank)));

    int  nihalo       = captured.integer("G%Domain%nihalo");
    int  njhalo       = captured.integer("G%Domain%njhalo");
    int  niglobal     = captured.integer("G%Domain%niglobal");
    int  njglobal     = captured.integer("G%Domain%njglobal");
    int  isd_global   = captured.integer("G%isd_global");
    int  jsd_global   = captured.integer("G%jsd_global");
    int  ke           = captured.integer("G%ke");
    auto h_tmp_before = captured.fab_host("h_tmp_TIM_before", isd_global, jsd_global);
    auto h_tmp_after  = captured.fab_host("h_tmp_TIM_after" , isd_global, jsd_global);

    amrex::BoxList     boxList;
    amrex::Vector<int> processBoxList;
    int iLocalGridSize = h_tmp_before.box().length(0) - (2*nihalo);
    int jLocalGridSize = h_tmp_before.box().length(1) - (2*njhalo);

    // Assumes equally sized grids
    for(int i = 0; i < amrex::ParallelDescriptor::NProcs(); i++)
    {
        test_mom::CapturedFile rankc(test_mom::data_dir /
                                    ("MOM_dynamics_split_RK2_" + std::format("{:04}", i)));
        int iGridStartLoc = rankc.integer("G%isd_global") + nihalo;
        int jGridStartLoc = rankc.integer("G%jsd_global") + njhalo;
        amrex::IntVect lo(iGridStartLoc, jGridStartLoc, 0);
        amrex::IntVect hi(lo[0] + iLocalGridSize - 1, lo[1] + jLocalGridSize - 1, ke - 1);
        amrex::Box b(lo, hi);
        boxList.push_back(b);
        processBoxList.push_back(i);
    }

    // Setup multifab
    amrex::BoxArray boxes(boxList);
    amrex::DistributionMapping dm(processBoxList);
    amrex::IntVect fourRowTwo2DGhostCells{nihalo, njhalo, 0};
    amrex::MultiFab mf(boxes, dm, /*ncomp=*/1, fourRowTwo2DGhostCells);

    // Setup box
    amrex::Box globalComputationalDomain(amrex::IntVect(           0,            0,  0),
                                         amrex::IntVect(niglobal - 1, njglobal - 1, ke-1));
    amrex::RealBox realBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    amrex::Geometry geom(globalComputationalDomain, &realBox);

    // Verify that the loops iterate only the expected number of times
    int mf_iterations = 0;
    for (amrex::MFIter mfi(mf, /*tiling=*/false); mfi.isValid(); ++mfi)
    {
        amrex::FArrayBox& localValidFab  = mf[mfi];
        ASSERT_EQ(localValidFab.box(), h_tmp_before.box());

        localValidFab.copy(h_tmp_before);
        mf_iterations++;
    }
    ASSERT_EQ(1, mf_iterations);

    mf.FillBoundary(geom.periodicity());

    mf_iterations = 0;
    for (amrex::MFIter mfi(mf, /*tiling=*/false); mfi.isValid(); ++mfi)
    {
        amrex::FArrayBox& localValidFab  = mf[mfi];
        ASSERT_EQ(localValidFab.box(), h_tmp_after.box());

        expect_arrays_equal(h_tmp_after, localValidFab, "h_tmp");
        mf_iterations++;
    }
    ASSERT_EQ(1, mf_iterations);
}
