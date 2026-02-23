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

!> @}
