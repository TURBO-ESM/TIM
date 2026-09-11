#include <format>

#include <gtest/gtest.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>

#include "amrex_assertions.hpp"
#include "captured_io.hpp"
#include "data_dir.hpp"

using test_mom::expect_arrays_equal;

TEST(BoundaryExchange, BoundaryMatchesAfterOneTimeStep) {
    int myRank = amrex::ParallelDescriptor::MyProc();
    test_mom::CapturedFile captured(test_mom::data_dir /
                                    ("MOM_fixed_initialization_" + std::format("{:04}", myRank)));

    // Number of halo cells in x,y directions
    int  nihalo       = captured.integer("G%Domain%nihalo");
    int  njhalo       = captured.integer("G%Domain%njhalo");
    // Total computational domain not including ghost cells
    int  niglobal     = captured.integer("G%Domain%niglobal");
    int  njglobal     = captured.integer("G%Domain%njglobal");
    // Global location of current grid with ghost cells counting negative at origin
    // using 0-based indexing
    int  isd_global   = captured.integer("G%isd_global");
    int  jsd_global   = captured.integer("G%jsd_global");

    auto bathyTBefore = captured.fab_host("G%bathyT_brefore", isd_global, jsd_global);
    auto bathyTAfter  = captured.fab_host("G%bathyT_after"  , isd_global, jsd_global);

    // Sets up data domain region
    int iLocalDataGridSize = bathyTBefore.box().length(0) - (2*nihalo);
    int jLocalDataGridSize = bathyTBefore.box().length(1) - (2*njhalo);

    amrex::BoxList     boxList;
    amrex::Vector<int> processBoxList;

    // Assumes equally sized grids
    for(int i = 0; i < amrex::ParallelDescriptor::NProcs(); i++)
    {
        test_mom::CapturedFile rankc(test_mom::data_dir /
                                    ("MOM_fixed_initialization_" + std::format("{:04}", i)));

        // Each rank needs to map the region every other rank convers.
        int iGridStartLoc = rankc.integer("G%isd_global") + nihalo;
        int jGridStartLoc = rankc.integer("G%jsd_global") + njhalo;
        amrex::IntVect lo(iGridStartLoc, jGridStartLoc, 0);
        amrex::IntVect hi(lo[0] + iLocalDataGridSize - 1, lo[1] + jLocalDataGridSize - 1, 0);

        boxList.push_back(amrex::Box(lo, hi));
        processBoxList.push_back(i);
    }

    // Setup multifab
    amrex::BoxArray            boxes(boxList);
    amrex::DistributionMapping dm(processBoxList);
    amrex::IntVect             fourRowTwo2DGhostCells{nihalo, njhalo, 0};
    amrex::MultiFab            mf(boxes, dm, /*ncomp=*/1, fourRowTwo2DGhostCells);

    // Setup geometry
    amrex::Box globalComputationalDomain(amrex::IntVect(           0,            0, 0),
                                         amrex::IntVect(niglobal - 1, njglobal - 1, 0));
    amrex::RealBox  realBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    amrex::Geometry geom(globalComputationalDomain, &realBox);

    // Verify that the loops iterate only the expected number of times
    int mf_iterations = 0;
    for (amrex::MFIter mfi(mf, /*tiling=*/false); mfi.isValid(); ++mfi)
    {
        amrex::FArrayBox& localValidFab  = mf[mfi];
        ASSERT_EQ(bathyTBefore.box(), localValidFab.box());

        localValidFab.copy(bathyTBefore);
        mf_iterations++;
    }
    ASSERT_EQ(1, mf_iterations);

    mf.FillBoundary(geom.periodicity());

    mf_iterations = 0;
    for (amrex::MFIter mfi(mf, /*tiling=*/false); mfi.isValid(); ++mfi)
    {
        amrex::FArrayBox& localValidFab  = mf[mfi];
        ASSERT_EQ(bathyTAfter.box(), localValidFab.box());

        expect_arrays_equal(bathyTAfter, localValidFab, "G%bathyT");
        mf_iterations++;
    }
    ASSERT_EQ(1, mf_iterations);
}
