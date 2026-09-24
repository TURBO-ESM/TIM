module tim_coms_infra_interface

use iso_fortran_env,         only : int64
use iso_c_binding,           only : c_int64_t, c_size_t, c_ptr, c_null_ptr, c_loc, c_int
use turbotmp_bridge_c_types, only : RealArray_c
implicit none
private

public :: tim_chksum

interface tim_chksum_c
  function tim_chksum_c(field, mask_ptr) bind(c, name="tim_chksum_c")
    import c_ptr, c_int64_t, c_size_t, RealArray_C
    integer(c_int64_t)                        :: tim_chksum_c
    type(RealArray_C),             intent(in) :: field
    type(c_ptr),            value, intent(in) :: mask_ptr
  end function tim_chksum_c
end interface tim_chksum_c

interface tim_chksum
  module procedure tim_chksum_real_0d
  module procedure tim_chksum_real_nd
end interface tim_chksum

contains

function tim_chksum_real_0d(field, mask_val) result(chksum)
  real,              target, intent(in) :: field               !< Input scalar
  real,    optional, target, intent(in) :: mask_val            !< FMS mask value
  type(RealArray_C)                     :: field_in
  type(c_ptr)                           :: mask_loc !< c pointers to field and mask
  integer(kind=int64)                   :: chksum              !< checksum of array
  integer(c_int), target :: shp(3), lb(3), ub(3)

  shp=1; lb=1; ub=1
  field_in%data  = c_loc(field)
  field_in%shape = c_loc(shp)
  field_in%lb    = c_loc(lb)
  field_in%ub    = c_loc(ub)
  field_in%rank  = 3

  mask_loc = c_null_ptr
  if(present(mask_val)) mask_loc = c_loc(mask_val)

  chksum = tim_chksum_c(field_in, mask_loc)
end function tim_chksum_real_0d

function tim_chksum_real_nd(field, mask_val) result(chksum)
  type(RealArray_C),          intent(in) :: field               !< Input array
  real,     optional, target, intent(in) :: mask_val            !< FMS mask value
  type(c_ptr)                            :: mask_loc !< c pointers to field and mask
  integer(kind=int64)                    :: chksum              !< checksum of array

  mask_loc = c_null_ptr
  if(present(mask_val)) mask_loc = c_loc(mask_val)

  chksum = tim_chksum_c(field, mask_loc)
end function tim_chksum_real_nd

end module tim_coms_infra_interface

