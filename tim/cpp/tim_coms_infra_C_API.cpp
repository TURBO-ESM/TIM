#include "tim_coms_infra_C_API.h"
#include "tim_coms_infra.hpp"

int64_t tim_chksum_c(double* field, size_t field_size, double* mask_val) {
    int field_size_int = static_cast<int>(field_size);
    
    amrex::Dim3 begin = {0, 0, 0};
    amrex::Dim3 end   = {field_size_int, 1, 1};
    amrex::Array4<amrex::Real> arr(field, begin, end, 1);
    amrex::Box bx(amrex::IntVect{0,0,0}, amrex::IntVect{field_size_int-1, 0, 0});

    if(mask_val)
        return TIM::checksum(bx, arr, *mask_val);
    else
        return TIM::checksum(bx, arr);
}
