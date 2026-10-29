# Precision-based Fortran compiler flags
set(r8_flags "-fdefault-real-8 -fdefault-double-8")
# flang's default REAL is already 4 bytes, and flang has no -fdefault-real-4.
set(r4_flags "")

# LLVM Flang
# flang has no -fcray-pointer, -ffree-line-length-none or -fno-range-check:
# its defaults already behave the way those flags make gfortran behave.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fconvert=big-endian")
# GNU's Release level; flang has no -funroll-all-loops or -finline-functions.
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
# Unlike GNU's Debug flags, these neither check array bounds nor trap
# floating-point exceptions: flang 21 has no -fcheck=bounds or -ffpe-trap.
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
set(CMAKE_Fortran_LINK_FLAGS "")

# UFS flags are not known for LLVM Flang.
