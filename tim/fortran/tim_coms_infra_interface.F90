module tim_coms_infra_interface

use iso_fortran_env, only : int32, int64
use iso_c_binding,   only : c_int64_t, c_double, c_size_t, c_ptr, c_null_ptr, c_loc
use array_mod, only : RealArray_c
use box_mod,   only : Box_c
implicit none
private

public :: tim_chksum

interface tim_chksum_c
  function tim_chksum_c(field, mask_ptr, global_chksum) bind(c, name="tim_chksum_c")
    import c_ptr, c_int64_t, c_size_t
    integer(c_int64_t)                        :: tim_chksum_c
    type(RealArray_C),      value, intent(in) :: field
    type(c_ptr),            value, intent(in) :: mask_ptr
  end function tim_chksum_c
end interface tim_chksum_c

interface tim_chksum
  module procedure tim_chksum_real_0d
  module procedure tim_chksum_real_1d
  module procedure tim_chksum_real_2d
  module procedure tim_chksum_real_3d
  module procedure tim_chksum_real_4d
end interface tim_chksum

contains

function tim_chksum_real_0d(field, pelist, mask_val) result(chksum)
  real,              target, intent(in) :: field               !< Input scalar
  integer, optional, target, intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,    optional, target, intent(in) :: mask_val            !< FMS mask value
  type(RealArray_C)                     :: field_in
  type(c_ptr)                           :: mask_loc !< c pointers to field and mask
  integer(kind=int64)                   :: chksum              !< checksum of array

  call field_in%alloc(lb=[1], ub=[1], source=[field])

  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  chksum = tim_chksum_c(field_in, mask_loc)

  field_in%free()
end function tim_chksum_real_0d

function tim_chksum_real_1d(field, pelist, mask_val) result(chksum)
  real, dimension(:), target, intent(in) :: field               !< Input array
  integer,  optional, target, intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,     optional, target, intent(in) :: mask_val            !< FMS mask value
  type(RealArray_C)                      :: field_in
  type(c_ptr)                            :: mask_loc !< c pointers to field and mask
  integer(kind=int64)                    :: chksum              !< checksum of array

  call field_in%alloc(lb=LBOUND(field), ub=UBOUND(field), source=field)

  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  chksum = tim_chksum_c(field_in, mask_loc)

  field_in%free()
end function tim_chksum_real_1d

function tim_chksum_real_2d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:), target, intent(in) :: field               !< Unrotated input field
  integer,    optional, target, intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,       optional, target, intent(in) :: mask_val            !< FMS mask value
  type(RealArray_C)                        :: field_in
  type(c_ptr)                              :: mask_loc !< c pointers to field and mask
  type(c_bool)                             :: global_chksum
  integer(kind=int64)                      :: chksum              !< checksum of array

  call field_in%alloc(lb=LBOUND(field), ub=UBOUND(field), source=field)

  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  global_chksum = present(pelist)

  chksum = tim_chksum_c(field_in, mask_loc, global_chksum)

  field_in%free()
end function tim_chksum_real_2d

function tim_chksum_real_3d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:), target, intent(in) :: field               !< Unrotated input field
  integer,      optional, target, intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,         optional, target, intent(in) :: mask_val            !< FMS mask value
  type(RealArray_C)                          :: field_in
  type(c_ptr)                                :: mask_loc !< c pointers to field and mask
  type(c_bool)                               :: global_chksum
  integer(kind=int64)                        :: chksum              !< checksum of array

  call field_in%alloc(lb=LBOUND(field), ub=UBOUND(field), source=field)

  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  global_chksum = present(pelist)

  chksum = tim_chksum_c(field_in, mask_loc, global_chksum)

  field_in%free()
end function tim_chksum_real_3d

function tim_chksum_real_4d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:,:), target, intent(in) :: field               !< Unrotated input field
  integer,        optional, target, intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,           optional, target, intent(in) :: mask_val            !< FMS mask value
  type(RealArray_C)                            :: field_in
  type(c_ptr)                                  :: mask_loc !< c pointers to field and mask
  type(c_bool)                                 :: global_chksum
  integer(kind=int64)                          :: chksum              !< checksum of array

  call field_in%alloc(lb=LBOUND(field), ub=UBOUND(field), source=field)
  
  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  global_chksum = present(pelist)

  chksum = tim_chksum_c(field_in, mask_loc, global_chksum)

  field_in%free()
end function tim_chksum_real_4d

end module tim_coms_infra_interface

