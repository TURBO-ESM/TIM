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
!> @brief This series of routines is used to describe the contents of the file
!! being written on <unit>.

!> @addtogroup mpp_io_mod
!> @{

!> Writes a global metadata attribute to unit <unit>
!! attribute <name> can be an real, integer or character
!! one and only one of rval, ival, and cval should be present
!! the first found will be used
!! for a non-netCDF file, it is encoded into a string "GLOBAL <name> <val>"
    subroutine mpp_write_meta_global( unit, name, rval, ival, cval, pack)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer, intent(in), optional :: pack

!      call mpp_clock_begin(mpp_write_clock)
      if( .NOT.module_is_initialized    )call mpp_error( FATAL, 'MPP_WRITE_META: must first call mpp_io_init.' )
      if( .NOT. mpp_file(unit)%write_on_this_pe) then
!         call mpp_clock_end(mpp_write_clock)
         return
      endif
      if( .NOT.mpp_file(unit)%opened )call mpp_error( FATAL, 'MPP_WRITE_META: invalid unit number.' )
      if( mpp_file(unit)%initialized ) &
           call mpp_error( FATAL, 'MPP_WRITE_META: cannot write metadata to file after an mpp_write.' )

      if( mpp_file(unit)%format.EQ.MPP_NETCDF )then
#ifdef use_netCDF
          call write_attribute_netcdf( unit, NF_GLOBAL, name, rval, ival, cval, pack )
#endif
      else
          call write_attribute( unit, 'GLOBAL '//trim(name), rval, ival, cval, pack )
      end if
!      call mpp_clock_end(mpp_write_clock)

      return
    end subroutine mpp_write_meta_global

    subroutine mpp_write_meta_global_scalar_i( unit, name, ival, pack )
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      integer, intent(in) :: ival
      integer, intent(in), optional :: pack

      call mpp_write_meta_global( unit, name, ival=(/ival/), pack=pack )
      return
    end subroutine mpp_write_meta_global_scalar_i

    subroutine mpp_def_dim_nodata(unit,name,size)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      integer, intent(in) :: size
      integer :: error,did

    ! This routine assumes the file is in define mode
      if(.NOT. mpp_file(unit)%write_on_this_pe) return
#ifdef use_netCDF
      error = NF_DEF_DIM(mpp_file(unit)%ncid,name,size,did)
      call netcdf_err(error, mpp_file(unit),string='Axis='//trim(name))
#endif
    end subroutine mpp_def_dim_nodata

    subroutine mpp_def_dim_int(unit,name,dsize,longname,data)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      integer, intent(in) :: dsize
      character(len=*), intent(in) :: longname
      integer, intent(in) :: data(:)
      integer :: error,did,id

    ! This routine assumes the file is in define mode
#ifdef use_netCDF
      if(.NOT. mpp_file(unit)%write_on_this_pe) return
      error = NF_DEF_DIM(mpp_file(unit)%ncid,name,dsize,did)
      call netcdf_err(error, mpp_file(unit),string='Axis='//trim(name))

    ! Write dimension data.
      error = NF_DEF_VAR( mpp_file(unit)%ncid, name, NF_INT, 1, (/did/), id )
      call netcdf_err( error, mpp_file(unit), string=' axis varable '//trim(name))

      error = NF_PUT_ATT_TEXT( mpp_file(unit)%ncid, id, 'long_name', len_trim(longname), longname )
      call netcdf_err( error, mpp_file(unit), string=' Attribute=long_name' )

      if( mpp_file(unit)%action.EQ.MPP_WRONLY )then
         if(header_buffer_val>0) then
            error = NF__ENDDEF(mpp_file(unit)%ncid,header_buffer_val,4,0,4)
         else
            error = NF_ENDDEF(mpp_file(unit)%ncid)
         endif
      endif
      call netcdf_err( error, mpp_file(unit), string=' subroutine mpp_def_dim')
      error = NF_PUT_VARA_INT ( mpp_file(unit)%ncid, id, (/1/), (/size(data)/), data )
      call netcdf_err( error, mpp_file(unit), string=' axis varable '//trim(name))
      error = NF_REDEF(mpp_file(unit)%ncid)
      call netcdf_err( error, mpp_file(unit), string=' subroutine mpp_def_dim')
#endif
    return
  end subroutine mpp_def_dim_int

    subroutine mpp_def_dim_real(unit,name,dsize,longname,data)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      integer, intent(in) :: dsize
      character(len=*), intent(in) :: longname
      real, intent(in) :: data(:)
      integer :: error,did,id

    ! This routine assumes the file is in define mode
#ifdef use_netCDF
      if(.NOT. mpp_file(unit)%write_on_this_pe) return
      error = NF_DEF_DIM(mpp_file(unit)%ncid,name,dsize,did)
      call netcdf_err(error, mpp_file(unit),string='Axis='//trim(name))

    ! Write dimension data.
      error = NF_DEF_VAR( mpp_file(unit)%ncid, name, NF_INT, 1, (/did/), id )
      call netcdf_err( error, mpp_file(unit), string=' axis varable '//trim(name))

      error = NF_PUT_ATT_TEXT( mpp_file(unit)%ncid, id, 'long_name', len_trim(longname), longname )
      call netcdf_err( error, mpp_file(unit), string=' Attribute=long_name' )

      if( mpp_file(unit)%action.EQ.MPP_WRONLY )then
         if(header_buffer_val>0) then
            error = NF__ENDDEF(mpp_file(unit)%ncid,header_buffer_val,4,0,4)
         else
            error = NF_ENDDEF(mpp_file(unit)%ncid)
         endif
      endif
      call netcdf_err( error, mpp_file(unit), string=' subroutine mpp_def_dim')
      error = NF90_PUT_VAR ( mpp_file(unit)%ncid, id, data, start=(/1/), count=(/size(data)/) )
      call netcdf_err( error, mpp_file(unit), string=' axis varable '//trim(name))
      error = NF_REDEF(mpp_file(unit)%ncid)
      call netcdf_err( error, mpp_file(unit), string=' subroutine mpp_def_dim')
#endif
    return
  end subroutine mpp_def_dim_real


    subroutine write_attribute( unit, name, rval, ival, cval, pack )
!called to write metadata for non-netCDF I/O
      integer, intent(in) :: unit
      character(len=*), intent(in) :: name
      real, intent(in), optional :: rval(:)
      integer, intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
!pack is currently ignored in this routine: only used by netCDF I/O
      integer, intent(in), optional :: pack

      if( mpp_file(unit)%nohdrs )return
!encode text string
      if( PRESENT(rval) )then
          write( text,* )trim(name)//'=', rval
      else if( PRESENT(ival) )then
          write( text,* )trim(name)//'=', ival
      else if( PRESENT(cval) )then
          text = ' '//trim(name)//'='//trim(cval)
      else
          call mpp_error( FATAL, 'WRITE_ATTRIBUTE: one of rval, ival, cval must be present.' )
      end if
      if( mpp_file(unit)%format.EQ.MPP_ASCII )then
!implies sequential access
          write( unit,fmt='(a)' )trim(text)//char(10)
      else                      !MPP_IEEE32 or MPP_NATIVE
          if( mpp_file(unit)%access.EQ.MPP_SEQUENTIAL )then
              write(unit)trim(text)//char(10)
          else                  !MPP_DIRECT
              write( unit,rec=mpp_file(unit)%record )trim(text)//char(10)
              if( verbose )print '(a,i6,a,i3)', 'WRITE_ATTRIBUTE: PE=', pe, ' wrote record ', mpp_file(unit)%record
              mpp_file(unit)%record = mpp_file(unit)%record + 1
          end if
      end if
      return
    end subroutine write_attribute

    subroutine write_attribute_netcdf( unit, id, name, rval, ival, cval, pack )
!called to write metadata for netCDF I/O
      integer, intent(in) :: unit
      integer, intent(in) :: id
      character(len=*), intent(in) :: name
      real,             intent(in), optional :: rval(:)
      integer,          intent(in), optional :: ival(:)
      character(len=*), intent(in), optional :: cval
      integer, intent(in), optional :: pack
      integer, allocatable :: rval_i(:)
#ifdef use_netCDF
      if( PRESENT(rval) )then
!pack was only meaningful for FP numbers, but is now extended by the ival branch of this routine
          if( PRESENT(pack) )then
             if( pack== 0 ) then !! here be dragons, use ival branch!...
                if( KIND(rval).EQ.r8_kind )then
                   call mpp_error( FATAL, &
                        'WRITE_ATTRIBUTE_NETCDF: attempting to write internal NF_INT, currently int32, as double.' )
                else if( KIND(rval).EQ.r4_kind )then
                   call mpp_error( FATAL, &
                        'WRITE_ATTRIBUTE_NETCDF: attempting to write internal NF_INT, currently int32, as float.' )
                end if
             else if( pack.EQ.1 )then
                  if( KIND(rval).EQ.r8_kind )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_DOUBLE, &
                                                size(rval(:)), real(rval, kind=r8_kind))
                  else if( KIND(rval).EQ.r4_kind )then
                      call mpp_error( WARNING, &
                           'WRITE_ATTRIBUTE_NETCDF: attempting to write internal 32-bit real as external 64-bit.' )
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_DOUBLE, &
                                                size(rval(:)), real(rval, kind=r4_kind))
                  end if
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
              else if( pack.EQ.2 )then
                  if( KIND(rval).EQ.r8_kind )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_FLOAT, &
                                                size(rval(:)), real(rval, kind=r8_kind))
                  else if( KIND(rval).EQ.r4_kind )then
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_FLOAT, &
                                                size(rval(:)), real(rval, kind=r4_kind))
                  end if
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
              else if( pack.EQ.4 )then
                  allocate( rval_i(size(rval(:))) )
                  rval_i = rval
                  if( KIND(rval).EQ.r8_kind )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_SHORT, &
                                                size(rval_i(:)), real(rval, kind=r8_kind))
                  else if( KIND(rval).EQ.r4_kind )then
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_SHORT, &
                                                size(rval_i(:)), real(rval, kind=r4_kind))
                  end if
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
                  deallocate(rval_i)
              else if( pack.EQ.8 )then
                  allocate( rval_i(size(rval(:))) )
                  rval_i = rval
                  if( KIND(rval).EQ.r8_kind )then
                      error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_BYTE, &
                                                size(rval_i(:)), real(rval, kind=r8_kind))
                  else if( KIND(rval).EQ.r4_kind )then
                      error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_BYTE, &
                                                size(rval_i(:)), real(rval, kind=r4_kind))
                  end if
                  call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
                  deallocate(rval_i)
              else
                  call mpp_error( FATAL, 'WRITE_ATTRIBUTE_NETCDF: only legal packing values are 1,2,4,8.' )
              end if
          else
!default is to write FLOATs (32-bit)
              if( KIND(rval).EQ.r8_kind )then
                  error = NF_PUT_ATT_DOUBLE( mpp_file(unit)%ncid, id, name, NF_FLOAT, &
                                            size(rval(:)), real(rval, kind=r8_kind))
              else if( KIND(rval).EQ.r4_kind )then
                  error = NF_PUT_ATT_REAL  ( mpp_file(unit)%ncid, id, name, NF_FLOAT, &
                                            size(rval(:)), real(rval, kind=r4_kind))
              end if
              call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
          end if
      else if( PRESENT(ival) )then
          if( PRESENT(pack) ) then
             if (pack ==0) then
                if (KIND(ival).EQ.i8_kind) then
                   call mpp_error(FATAL,'only use NF_INTs with pack=0 for now')
                end if
                error = NF_PUT_ATT_INT( mpp_file(unit)%ncid, id, name, &
                                       NF_INT, size(ival(:)), ival ) !!XXX int32_t
                call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
             else
                call mpp_error( FATAL, 'WRITE_ATTRIBUTE_NETCDF: only implimented ints when pack=0, else use reals.' )
             endif
          else
          error = NF_PUT_ATT_INT ( mpp_file(unit)%ncid, id, name, NF_INT, size(ival(:)), ival )
          call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
          end if
      else if( present(cval) )then
          if (.NOT.cf_compliance .or. trim(name).NE.'calendar') then
            error = NF_PUT_ATT_TEXT( mpp_file(unit)%ncid, id, name, len_trim(cval), cval )
          else
            error = NF_PUT_ATT_TEXT( mpp_file(unit)%ncid, id, name, len_trim(cval), lowercase(cval) )
          endif
          call netcdf_err( error, mpp_file(unit), string=' Attribute='//name )
      else
          call mpp_error( FATAL, 'WRITE_ATTRIBUTE_NETCDF: one of rval, ival, cval must be present.' )
      end if
#endif /* use_netCDF */
      return
    end subroutine write_attribute_netcdf

    subroutine mpp_modify_axis_meta( axis, name, units, longname, cartesian, data )

      type(axistype), intent(inout) :: axis
      character(len=*), intent(in), optional :: name, units, longname, cartesian
      real, dimension(:), intent(in), optional :: data

      if (PRESENT(name)) axis%name = trim(name)
      if (PRESENT(units)) axis%units = trim(units)
      if (PRESENT(longname)) axis%longname = trim(longname)
      if (PRESENT(cartesian)) axis%cartesian = trim(cartesian)
      if (PRESENT(data)) then
         axis%len = size(data(:))
         if (ASSOCIATED(axis%data)) deallocate(axis%data)
         allocate(axis%data(axis%len))
         axis%data = data
      endif

      return
    end subroutine mpp_modify_axis_meta

    subroutine mpp_modify_field_meta( field, name, units, longname, min, max, missing, axes )

      type(fieldtype), intent(inout) :: field
      character(len=*), intent(in), optional :: name, units, longname
      real, intent(in), optional :: min, max, missing
      type(axistype), dimension(:), intent(inout), optional :: axes

      if (PRESENT(name)) field%name = trim(name)
      if (PRESENT(units)) field%units = trim(units)
      if (PRESENT(longname)) field%longname = trim(longname)
      if (PRESENT(min)) field%min = min
      if (PRESENT(max)) field%max = max
      if (PRESENT(missing)) field%missing = missing
!      if (PRESENT(axes)) then
!         axis%len = size(data(:))
!         deallocate(axis%data)
!         allocate(axis%data(axis%len))
!         axis%data = data
!      endif

      return
    end subroutine mpp_modify_field_meta


!> \brief Fills in a fieldtype variable, and is used with the diag_manager when using
!! fms2_io
    subroutine fillin_fieldtype(field, axes, name, units, longname,&
         min, max, missing, fill, scale, add, pack, time_method, standard_name, checksum)
!define field: must have already called mpp_write_meta(axis) for each axis
      type(fieldtype), intent(inout) :: field
      type(axistype), intent(in) :: axes(:)
      character(len=*), intent(in) :: name, units, longname
      real, intent(in), optional :: min, max, missing, fill, scale, add
      integer, intent(in), optional :: pack
      character(len=*), intent(in), optional :: time_method
      character(len=*), intent(in), optional :: standard_name
      integer(i8_kind), dimension(:), intent(in), optional :: checksum
!this array is required because of f77 binding on netCDF interface
      integer, allocatable :: axis_id(:)
      real :: a, b
      integer :: i, istat, ishift, jshift
      character(len=64) :: checksum_char

!      call mpp_clock_begin(mpp_write_clock)

      !--- figure out the location of data, this is needed in mpp_write.
      !--- for NON-symmetry domain, the position is not an issue.
      !--- we may need to rethink how to address the symmetric issue.
      ishift = 0; jshift = 0
      do i = 1, size(axes(:))
         select case ( lowercase( axes(i)%cartesian ) )
         case ( 'x' )
            ishift = axes(i)%shift
         case ( 'y' )
            jshift = axes(i)%shift
         end select
      end do

      field%position = CENTER
      if(ishift == 1 .AND. jshift == 1) then
         field%position = CORNER
      else if(ishift == 1) then
         field%position = EAST
      else if(jshift == 1) then
         field%position = NORTH
      endif

!pre-existing pointers need to be nullified
      if( ASSOCIATED(field%axes) ) DEALLOCATE(field%axes, stat=istat)
      if( ASSOCIATED(field%size) ) DEALLOCATE(field%size, stat=istat)
!fill in field metadata
      field%name = name
      field%units = units
      field%longname = longname
      allocate( field%axes(size(axes(:))) )
      field%axes = axes
      field%ndim = size(axes(:))
      field%time_axis_index = -1 !this value will never match any axis index
!size is buffer area for the corresponding axis info: it is required to buffer this info in the fieldtype
!because axis might be reused in different files
      allocate( field%size(size(axes(:))) )
      do i = 1,size(axes(:))
         if( ASSOCIATED(axes(i)%data) )then !space axis
             field%size(i) = size(axes(i)%data(:))
         else               !time
             field%size(i) = 1
             field%time_axis_index = i
         end if
      end do
!attributes
      if( PRESENT(min) )          field%min           = min
      if( PRESENT(max) )          field%max           = max
      if( PRESENT(scale) )        field%scale         = scale
      if( PRESENT(add) )          field%add           = add
      if( PRESENT(standard_name)) field%standard_name = standard_name
      if( PRESENT(missing) )      field%missing       = missing
      if( PRESENT(fill) )         field%fill          = fill
      field%checksum      = 0
      if( PRESENT(checksum) )     field%checksum(1:size(checksum)) = checksum(:)

      ! Issue warning if fill and missing are different
      if ( (present(fill).and.present(missing)) .and. (field%missing .ne. field%fill) ) then
         call mpp_error(WARNING, 'MPP_WRITE_META: NetCDF attributes _FillValue and missing_value should be equal.')
      end if
!pack is currently used only for netCDF
      field%pack = 2        !default write 32-bit floats
      if( PRESENT(pack) )field%pack = pack

      return
    end subroutine fillin_fieldtype
!> @}
