set(COMMON_SRC
    amrex_assertions.cpp
    captured_io.cpp)

add_library(common OBJECT ${COMMON_SRC})

target_link_libraries(common PUBLIC GTest::gtest AMReX::amrex)
target_include_directories(common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
