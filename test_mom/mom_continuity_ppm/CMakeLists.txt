add_executable(test_mom_continuity_ppm test_mom_continuity_ppm.cpp)
target_sources(test_mom_continuity_ppm PRIVATE $<TARGET_OBJECTS:data_driven_shared_main>)
target_link_libraries(test_mom_continuity_ppm PRIVATE common mom_continuity)
gtest_discover_tests(test_mom_continuity_ppm
    EXTRA_ARGS --data-dir=${TEST_MOM_DATA_DIR}
    # Discover at test-run time rather than build time. With the CUDA backend
    # the binary calls amrex::Initialize() which aborts on hosts without a GPU
    # (e.g. Derecho login nodes), so running it at POST_BUILD would break the
    # build.
    DISCOVERY_MODE PRE_TEST
)
