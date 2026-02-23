! -*-f90-*-
!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!> @brief Routines to retrieve data used in @ref mpp_io_mod

!> @addtogroup mpp_io_mod
!> @{

    !> @brief Get some general information about a file.
    !!
    !> <br>Example usage:
    !! @code{.F90}
    !! call mpp_get_info( unit, ndim, nvar, natt, ntime )
    !! @endcode
    subroutine mpp_get_info( unit, ndim, nvar, natt, ntime )

      integer, intent(in) :: unit
      integer, intent(out) :: ndim, nvar, natt, ntime


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_INFO: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )&
           call mpp_error(FATAL, 'MPP_GET_INFO: invalid unit number, file '//trim(mpp_file(unit)%name))

      ndim = mpp_file(unit)%ndim
      nvar = mpp_file(unit)%nvar
      natt = mpp_file(unit)%natt
      ntime = mpp_file(unit)%time_level

      return

    end subroutine mpp_get_info

    !> @brief Copy global file attributes for use by user
    subroutine mpp_get_global_atts( unit, global_atts )

      integer,       intent(in)    :: unit
      type(atttype), intent(inout) :: global_atts(:) !< an attribute type which is allocated from the
                                                     !!  calling routine

      integer :: natt,i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_INFO: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )&
           call mpp_error( FATAL, 'MPP_GET_INFO: invalid unit number,file '//trim(mpp_file(unit)%name))

      if (size(global_atts(:)).lt.mpp_file(unit)%natt) &
           call mpp_error(FATAL, 'MPP_GET_ATTS: atttype not dimensioned properly in calling routine, file '// &
           trim(mpp_file(unit)%name))

      natt = mpp_file(unit)%natt
      global_atts = default_att

      do i=1,natt
         global_atts(i) = mpp_file(unit)%Att(i)
      enddo

      return
   end subroutine mpp_get_global_atts

  !#####################################################################
   subroutine mpp_get_field_atts(field, name, units, longname, min, max, missing, ndim, siz, axes, atts, &
                                 valid, scale, add, checksum)

     type(fieldtype),    intent(in)                            :: field
     character(len=*),   intent(out),                 optional :: name, units
     character(len=*),   intent(out),                 optional :: longname
     real,               intent(out),                 optional :: min,max,missing
     integer,            intent(out),                 optional :: ndim
     integer,            intent(out),   dimension(:), optional :: siz
     type(validtype),    intent(out),                 optional :: valid
     real,               intent(out),                 optional :: scale
     real,               intent(out),                 optional :: add
     integer(i8_kind), intent(out),   dimension(:), optional :: checksum

     type(atttype),      intent(inout), dimension(:), optional :: atts
     type(axistype),     intent(inout), dimension(:), optional :: axes

     integer :: n,m, check_exist

     if (PRESENT(name)) name = field%name
     if (PRESENT(units)) units = field%units
     if (PRESENT(longname)) longname = field%longname
     if (PRESENT(min)) min = field%min
     if (PRESENT(max)) max = field%max
     if (PRESENT(missing)) missing = field%missing
     if (PRESENT(ndim)) ndim = field%ndim
     if (PRESENT(atts)) then
        atts = default_att
        n = size(atts(:));m=size(field%Att(:))
        if (n.LT.m)&
             call mpp_error(FATAL,'attribute array not large enough in mpp_get_field_atts, field '//&
             trim(field%name))
        do n=1,m
          atts(n) = field%Att(n)
        end do
     end if
     if (PRESENT(axes)) then
        axes = default_axis
        n = size(axes(:));m=field%ndim
        if (n.LT.m) &
             call mpp_error(FATAL,'axis array not large enough in mpp_get_field_atts, field '//&
             trim(field%name))
        do n=1,m
          axes(n) = field%axes(n)
        end do
     end if
     if (PRESENT(siz)) then
        siz = -1
        n = size(siz(:));m=field%ndim
        if (n.LT.m) &
             call mpp_error(FATAL,'size array not large enough in mpp_get_field_atts, field '//&
             trim(field%name))
        do n=1,m
          siz(n) = field%size(n)
        end do
     end if

     if(PRESENT(valid)) then
       call mpp_get_valid(field,valid)
     endif

     if(PRESENT(scale))    scale    = field%scale
     if(present(add))      add      = field%add
     if(present(checksum)) then
       checksum = 0
       check_exist = mpp_find_att(field%Att(:),"checksum")
       if ( check_exist >= 0 ) then
         if(size(checksum(:)) >size(field%checksum(:))) call mpp_error(FATAL, &
            & "size(checksum(:)) >size(field%checksum(:))")
         checksum = field%checksum(1:size(checksum(:)))
       endif
     endif

     return
   end subroutine mpp_get_field_atts

  !#####################################################################
   subroutine mpp_get_axis_atts( axis, name, units, longname, cartesian, &
                                 calendar, sense, len, natts, atts, compressed )

     type(axistype), intent(in) :: axis
     character(len=*), intent(out) , optional :: name, units
     character(len=*), intent(out), optional :: longname, cartesian
     character(len=*), intent(out), optional :: compressed, calendar
     integer,intent(out), optional :: sense, len , natts
     type(atttype), intent(inout), optional, dimension(:) :: atts

     integer :: n,m

     if (PRESENT(name)) name = axis%name
     if (PRESENT(units)) units = axis%units
     if (PRESENT(longname)) longname = axis%longname
     if (PRESENT(cartesian)) cartesian = axis%cartesian
     if (PRESENT(compressed)) compressed = axis%compressed
     if (PRESENT(calendar)) calendar = axis%calendar
     if (PRESENT(sense)) sense = axis%sense
     if (PRESENT(len)) len = axis%len
     if (PRESENT(atts)) then
        atts = default_att
        n = size(atts(:));m=size(axis%Att(:))
        if (n.LT.m) &
             call mpp_error(FATAL,'attribute array not large enough in mpp_get_field_atts, axis '//&
             trim(axis%name))
        do n=1,m
          atts(n) = axis%Att(n)
        end do
     end if
     if (PRESENT(natts)) natts = size(axis%Att(:))

     return
   end subroutine mpp_get_axis_atts


  !#####################################################################
    !> @brief Copy variable information from file (excluding data)
    subroutine mpp_get_fields( unit, variables )
      integer,         intent(in)    :: unit
      type(fieldtype), intent(inout) :: variables(:)

      integer :: nvar,i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_FIELDS: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_GET_FIELDS: invalid unit number.' )

      if (size(variables(:)).ne.mpp_file(unit)%nvar) &
          call mpp_error(FATAL,'MPP_GET_FIELDS: fieldtype not dimensioned properly in calling routine, file '//&
          trim(mpp_file(unit)%name))

      nvar = mpp_file(unit)%nvar

      do i=1,nvar
         variables(i) = mpp_file(unit)%Var(i)
      enddo

      return
   end subroutine mpp_get_fields



  !#####################################################################
    !> @brief Copy variable information from file (excluding data)
    subroutine mpp_get_axes( unit, axes, time_axis )
      integer, intent(in) :: unit
      type(axistype), intent(inout) :: axes(:)
      type(axistype), intent(inout), optional :: time_axis
      integer :: ndim,i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GET_AXES: must first call mpp_io_init.' )
      if( .NOT.mpp_file(unit)%opened )&
           call mpp_error( FATAL, 'MPP_GET_AXES: invalid unit number, file '//trim(mpp_file(unit)%name))

      if (size(axes(:)).ne.mpp_file(unit)%ndim) &
           call mpp_error(FATAL, 'MPP_GET_AXES: axistype not dimensioned properly in calling routine, file '//&
           trim(mpp_file(unit)%name))


      if (PRESENT(time_axis)) time_axis = default_axis
      ndim = mpp_file(unit)%ndim

      do i=1,ndim
         axes(i)=mpp_file(unit)%Axis(i)

         if (PRESENT(time_axis) &
             .AND. .NOT. ASSOCIATED(mpp_file(unit)%Axis(i)%data) &
             .AND. mpp_file(unit)%Axis(i)%type /= -1) then
            time_axis = mpp_file(unit)%Axis(i)
         endif
      enddo

      return
   end subroutine mpp_get_axes

  !#####################################################################
   function mpp_get_field_size(field)

     type(fieldtype) :: field
     integer :: mpp_get_field_size(4)

     mpp_get_field_size = -1

     mpp_get_field_size(1) = field%size(1)
     mpp_get_field_size(2) = field%size(2)
     mpp_get_field_size(3) = field%size(3)
     mpp_get_field_size(4) = field%size(4)

     return
   end function mpp_get_field_size

  !#####################################################################
    !> @brief Set the mpp_io_stack variable to be at least n LONG words long
    subroutine mpp_io_set_stack_size(n)
      integer, intent(in) :: n
      character(len=10) :: text

      if( n.GT.mpp_io_stack_size .AND. allocated(mpp_io_stack) )deallocate(mpp_io_stack)
      if( .NOT.allocated(mpp_io_stack) )then
          allocate( mpp_io_stack(n) )
          mpp_io_stack_size = n
          write( text,'(i10)' )n
          if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, 'MPP_IO_SET_STACK_SIZE: stack size set to '//text//'.' )
      end if

      return
    end subroutine mpp_io_set_stack_size

  !#####################################################################
  !> Based on presence/absence of attributes, defines valid range or missing
  ! value. For details, see section 8.1 of NetCDF User Guide
  subroutine mpp_get_valid(f,v)
     type(fieldtype),intent(in)  :: f ! field
     type(validtype),intent(out) :: v ! validator

     integer :: irange,imin,imax,ifill,imissing,iscale
     integer :: valid_T, scale_T ! types of attributes

     v%is_range = .true.
     v%min = -HUGE(v%min); v%max = HUGE(v%max)
     if (f%natt == 0) return
     ! find indices of relevant attributes
     irange   = mpp_find_att(f%att,'valid_range')
     imin     = mpp_find_att(f%att,'valid_min')
     imax     = mpp_find_att(f%att,'valid_max')
     ifill    = mpp_find_att(f%att,'_FillValue')
     imissing = mpp_find_att(f%att,'missing_value')

     ! find the widest type of scale and offset; note that the code
     ! uses assumption that NetCDF types are arranged in th order of rank,
     ! that is NF_BYTE < NF_SHORT < NF_INT < NF_FLOAT < NF_DOUBLE
     scale_T = 0
     iscale   = mpp_find_att(f%att,'scale_factor')
     if(iscale>0) scale_T = f%att(iscale)%type
     iscale = mpp_find_att(f%att,'add_offset')
     if(iscale>0) scale_T = max(scale_T,f%att(iscale)%type)


     ! examine possible range attributes
     valid_T = 0
     if (irange>0) then
        v%min = f%att(irange)%fatt(1)
        v%max = f%att(irange)%fatt(2)
        valid_T = f%att(irange)%type
     else if (imax>0.or.imin>0) then
        if(imax>0) then
           v%max = f%att(imax)%fatt(1)
           valid_T = max(valid_T,f%att(imax)%type)
        endif
        if(imin>0) then
           v%min = f%att(imin)%fatt(1)
           valid_T = max(valid_T,f%att(imin)%type)
        endif
     else if (imissing > 0) then
        v%is_range = .false.
        ! here we always scale, since missing_value is supposed to be in
        ! external representation
        v%min = f%att(imissing)%fatt(1)*f%scale + f%add
     else if (ifill>0) then
     !z1l ifdef is added in to be able to compile without using use_netCDF.
#ifdef use_netCDF
        ! define min and max according to _FillValue
        if(f%att(ifill)%fatt(1)>0) then
            ! if _FillValue is positive, then it defines valid maximum
            v%max = f%att(ifill)%fatt(1)
            select case(f%type)
            case (NF_BYTE,NF_SHORT,NF_INT)
               v%max = v%max-1
            case (NF_FLOAT)
               v%max = nearest(nearest(real(v%max,4),-1.0),-1.0)
            case (NF_DOUBLE)
               v%max = nearest(nearest(real(v%max,8),-1.0),-1.0)
            end select
            ! always do the scaling, as the _FillValue is in external
            ! representation
            v%max = v%max*f%scale + f%add
        else
            ! if _FillValue is negative or zero, then it defines valid minimum
            v%min = f%att(ifill)%fatt(1)
            select case(f%type)
            case (NF_BYTE,NF_SHORT,NF_INT)
               v%min = v%min+1
            case (NF_FLOAT)
               v%min = nearest(nearest(real(v%min,4),+1.0),+1.0)
            case (NF_DOUBLE)
               v%min = nearest(nearest(real(v%min,8),+1.0),+1.0)
            end select
            ! always do the scaling, as the _FillValue is in external
            ! representation
            v%min = v%min*f%scale + f%add
        endif
#endif
    endif
   ! If valid_range is the same type as scale_factor (actually the wider of
   ! scale_factor and add_offset) and this is wider than the external data, then it
   ! will be interpreted as being in the units of the internal (unpacked) data.
   ! Otherwise it is in the units of the external (packed) data.
    ! Note that it is not relevant if we went through the missing_data of _FillValue
    ! brances, because in this case all irange, imin, and imax are less then 0
    if(.not.((valid_T == scale_T).and.(scale_T>f%type))) then
       if(irange>0 .or. imin>0) then
          v%min = v%min*f%scale + f%add
       endif
       if(irange>0 .or. imax>0) then
          v%max = v%max*f%scale + f%add
       endif
    endif

   end subroutine mpp_get_valid

   !#####################################################################
   ! finds an attribute by name in the array; returns -1 if it is not
   ! found
   function mpp_find_att(atts, name)
     integer                   :: mpp_find_att
     type(atttype), intent(in) :: atts(:) ! array of attributes
     character(len=*)          :: name ! name of the attributes

     integer :: i

     mpp_find_att = -1
     do i = 1, size(atts)
        if (trim(name)==trim(atts(i)%name)) then
           mpp_find_att=i
           exit
        endif
     enddo
   end function mpp_find_att
   !#####################################################################

   !> @brief return the name of an attribute.
   function mpp_get_att_name(att)
      type(atttype),    intent(in) :: att
      character(len=len(att%name)) :: mpp_get_att_name

      mpp_get_att_name = att%name
      return

   end function mpp_get_att_name

   !#####################################################################

   !> @brief return the char value of an attribute.
   function mpp_get_att_char(att)
      type(atttype), intent(in) :: att
      character(len=att%len)    :: mpp_get_att_char

      mpp_get_att_char = att%catt
      return

   end function mpp_get_att_char

!#######################################################################
logical function mpp_is_dist_ioroot(ssize,ioroot,lsize)
  integer, intent(in)  :: ssize   ! Dist io set size
  integer, intent(out), optional :: ioroot, lsize
  integer :: pe, npes, mypos, maxpe, d_ioroot, d_lsize, last_ioroot
  integer :: rootpe

  if(ssize < 1) call mpp_error(FATAL,'mpp_is_dist_ioroot: I/O stripe size < 1')

  mpp_is_dist_ioroot = .false.
  rootpe = mpp_root_pe()
  d_lsize = ssize
  pe = mpp_pe()
  mypos = modulo(pe-rootpe,ssize)  ! Which PE am I in the io group?
  d_ioroot = pe - mypos            ! What is the io root for the group?
  npes = mpp_npes()
  maxpe = min(d_ioroot+ssize,npes+rootpe) - 1  ! Handle end case
  d_lsize = maxpe - d_ioroot + 1
  if(mod(npes,ssize) == 1)then  ! Ensure there are no sets with 1 member
    last_ioroot = (npes-1) - ssize
    if(pe >= last_ioroot) then
       d_ioroot = last_ioroot
       d_lsize = ssize + 1
    endif
  endif
  if(pe == d_ioroot) mpp_is_dist_ioroot = .true.
  if(PRESENT(ioroot)) ioroot = d_ioroot
  if(PRESENT(lsize)) lsize = d_lsize
end function mpp_is_dist_ioroot
!> @}
