module turbotmp_bridge_c_types
  use iso_c_binding, only : c_ptr, c_int
  implicit none
  private
  public :: RealArray_C, IntArray_C, LogicalArray_C

  !< Type IntArray_C struct for C++ bridge layer
  type, bind(C) :: IntArray_C
     type(c_ptr) :: data               !< Storage pointer for array container
     type(c_ptr) :: shape              !< An array of dimension extents
     type(c_ptr) :: lb                 !< Lower bounds
     type(c_ptr) :: ub                 !< Upper bounds
     integer(c_int) :: rank            !< The number of dimensions
  end type IntArray_C

  !< RealArray struct for C bridge
  type, bind(C) :: RealArray_C
     type(c_ptr) :: data               !< Storage pointer for array container
     type(c_ptr) :: shape              !< An array of dimension extents
     type(c_ptr) :: lb                 !< Lower bounds
     type(c_ptr) :: ub                 !< Upper bounds
     integer(c_int) :: rank            !< The number of dimensions
  end type RealArray_C

  !< LogicalArray struct for C bridge. The data pointer is integer-encoded
  type, bind(C) :: LogicalArray_C
     type(c_ptr) :: data               !< Storage pointer for array container
     type(c_ptr) :: shape              !< An array of dimension extents
     type(c_ptr) :: lb                 !< Lower bounds
     type(c_ptr) :: ub                 !< Upper bounds
     integer(c_int) :: rank            !< The number of dimensions
  end type LogicalArray_C

end module turbotmp_bridge_c_types
