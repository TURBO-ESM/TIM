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
!> @defgroup fms_io_mod fms_io_mod
!> @ingroup fms
!> @brief Module for writing and reading restart data via NetCDF files
!> @author M.J. Harrison, Zhi Liang
!!
!! This module is for writing and reading restart data in NetCDF format.
!! fms_io_init must be called before the first write_data/read_data call
!! For writing, fms_io_exit must be called after ALL write calls have
!! been made. Typically, fms_io_init and fms_io_exit are placed in the
!! main (driver) program while read_data and write_data can be called where needed.
!! Presently, two combinations of threading and fileset are supported, users can choose
!! one line of the following by setting namelist:
!!
!! With the introduction of netCDF restart files, there is a need for a global
!! switch to turn on/off netCDF restart options in all of the modules that deal with
!! restart files. Here two more namelist variables (logical type) are introduced to fms_io
!!
!! - fms_netcdf_override
!! - fms_netcdf_restart
!!
!! because default values of both flags are .true., the default behavior of the entire model is
!! to use netCDF IO mode. To turn off netCDF restart, simply set fms_netcdf_restart to .false.

! <NAMELIST NAME="fms_io_nml">
! <DATA NAME="threading_read" TYPE="character">
! threading_read can be 'single' or 'multi'
! </DATA>
! <DATA NAME="fms_netcdf_override" TYPE="logical">
!   .true. : fms_netcdf_restart overrides individual do_netcdf_restart value (default behavior)
!   .false.: individual module settings has a precedence over the global setting, therefore
!   fms_netcdf_restart is ignored
! </DATA>
! <DATA NAME="fms_netcdf_restart" TYPE="logical">
!   .true. : all modules deal with restart files will operate under netCDF mode (default behavior)
!   .false.: all modules deal with restart files will operate under binary mode
!   This flag is effective only when fms_netcdf_override is .true. When fms_netcdf_override is .false., individual
!   module setting takes over.
! </DATA>
! <DATA NAME="time_stamped_restart" TYPE="logical">
!   .true. : time_stamp will be added to the restart file name as a prefix when
!            optional argument time_stamp is passed into routine save_restart.
!   .false.: time_stmp will not be added to the restart file name even though
!            time_stamp is passed into save_restart.
!    default is true.
! </DATA>
! <DATA NAME="print_chksum" TYPE="logical">
!    set print_chksum (default is false) to true to print out chksum of fields that are
!    read and written through save_restart/restore_state. The chksum is accross all the
!    processors, so there will be only one chksum even there are multiple-tiles in the
!    grid. For the multiple case, the filename appeared in the message will contain
!    tile1 because the message is print out from root pe and on root pe the tile id is tile1.
! </DATA>
! <DATA NAME="debug_mask_list" TYPE="logical">
!    set debug_mask_list (default is false) to true to print out mask_list reading from mask_table.
! </DATA>
! <DATA NAME="checksum_required" TYPE="logical">
!    Set checksum_required (default is true) to true to compare checksums stored in the attribute of a
!    field against the checksum after reading in the data. This check mitigates the possibility of data
!    that gets corrupted on write or read from being used in a n ongoing fashion. The checksum is across
!    all the  processors, so there will be only one checksum even if there are multiple-tiles in the
!    grid. For the decomposed file case, the filename appearing in the message will contain tile1
!    because the message is printed out from the root pe and on root pe the tile id is tile1.
!
!    Set checksum_required to false if you do not want to compare checksums.
! </DATA>
!</NAMELIST>

!> @addtogroup fms_io_mod
!> @{
module fms_io_mod

#include <fms_platform.h>

use mpp_io_mod,      only: mpp_open, mpp_close, mpp_io_init, mpp_io_exit, mpp_read, mpp_write
use mpp_io_mod,      only: mpp_write_meta, mpp_get_info, mpp_get_atts, mpp_get_fields
use mpp_io_mod,      only: mpp_read_compressed, mpp_write_compressed, mpp_def_dim
use mpp_io_mod,      only: mpp_write_unlimited_axis, mpp_read_distributed_ascii
use mpp_io_mod,      only: mpp_get_axes, mpp_get_axis_data, mpp_get_att_char, mpp_get_att_name
use mpp_io_mod,      only: mpp_get_att_real_scalar, mpp_attribute_exist, mpp_is_dist_ioroot
use mpp_io_mod,      only: fieldtype, axistype, atttype, default_field, default_axis, default_att
use mpp_io_mod,      only: MPP_NETCDF, MPP_ASCII, MPP_MULTI, MPP_SINGLE, MPP_OVERWR, MPP_RDONLY
use mpp_io_mod,      only: MPP_IEEE32, MPP_NATIVE, MPP_DELETE, MPP_APPEND, MPP_SEQUENTIAL, MPP_DIRECT
use mpp_io_mod,      only: MAX_FILE_SIZE, mpp_get_att_value
use mpp_io_mod,      only: mpp_get_dimension_length
use mpp_domains_mod, only: domain2d, domain1d, NULL_DOMAIN1D, NULL_DOMAIN2D, operator( .EQ. )
use mpp_domains_mod, only: CENTER, EAST, WEST, NORTH, SOUTH, CORNER
use mpp_domains_mod, only: mpp_get_domain_components, mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: mpp_get_domain_shift, mpp_get_global_domain, mpp_global_field, mpp_domain_is_tile_root_pe
use mpp_domains_mod, only: mpp_get_ntile_count, mpp_get_current_ntile, mpp_get_tile_id
use mpp_domains_mod, only: mpp_get_pelist, mpp_get_io_domain, mpp_get_domain_npes
use mpp_mod,         only: mpp_error, FATAL, NOTE, WARNING, mpp_pe, mpp_root_pe, mpp_npes, stdlog, stdout
use mpp_mod,         only: mpp_broadcast, ALL_PES, mpp_chksum, mpp_get_current_pelist, mpp_npes, lowercase
use mpp_mod,         only: input_nml_file, mpp_get_current_pelist_name, uppercase
use mpp_mod,         only: mpp_gather, mpp_scatter, mpp_send, mpp_recv, mpp_sync_self, COMM_TAG_1, EVENT_RECV
use mpp_mod,         only: MPP_FILL_DOUBLE,MPP_FILL_INT

use platform_mod, only: r8_kind

implicit none
private


integer, parameter, private :: max_split_file = 50
integer, parameter, private :: max_fields=400
integer, parameter, private :: max_axes=40
integer, parameter, private :: max_atts=20
integer, parameter, private :: max_domains = 100
integer, parameter, private :: MAX_TIME_LEVEL_REGISTER = 2
integer, parameter, private :: MAX_TIME_LEVEL_WRITE = 20
integer, parameter          :: max_axis_size=10000

! Index postions for axes in restart_file_type
! This is done so the user may define the axes
! in any order but a check can be performed
! to ensure no registration of duplicate axis

integer(INT_KIND),parameter,public :: XIDX = 1
integer(INT_KIND),parameter,public :: YIDX = 2
integer(INT_KIND),parameter,public :: CIDX = 3
integer(INT_KIND),parameter,public :: ZIDX = 4
integer(INT_KIND),parameter,public :: HIDX = 5
integer(INT_KIND),parameter,public :: TIDX = 6
integer(INT_KIND),parameter,public :: UIDX = 7
integer(INT_KIND),parameter,public :: CCIDX = 8

integer, parameter, private :: NIDX=8

logical, private :: warn_string_function = .true.

!> @}
!> @ingroup fms_io_mod
type, private :: meta_type
  type(meta_type), pointer :: prev=>null(), next=>null()
!!$ Gfortran on gaea does not yet support deferred length character strings
!!$  character(len=:),allocatable  :: name
  character(len=256)   :: name
  real,    allocatable :: rval(:)
  integer, allocatable :: ival(:)
!!$ Gfortran on gaea does not yet support deferred length character strings
!!$  character(len=:), allocatable :: cval
  character(len=256)   :: cval
end type meta_type

!> @ingroup fms_io_mod
type, private :: ax_type
   private
   character(len=128) :: name = ''
   character(len=128) :: units = ''
   character(len=128) :: longname = ''
   character(len=8)   :: cartesian = ''
   character(len=256) :: compressed = ''
   character(len=128) :: dimlen_name = ''
   character(len=128) :: dimlen_lname = ''
   character(len=128) :: calendar = ''
   integer            :: sense              !< Orientation of z axis definition
   integer            :: dimlen             !< max dim of elements across global domain
   real               :: min             !< valid min for real axis data
   integer            :: imin            !< valid min for integer axis data
   integer,allocatable :: idx(:)         !< compressed io-domain index vector
   integer,allocatable :: nelems(:)      !< num elements for each rank in io domain
   real, pointer      :: data(:) =>NULL()    !< real axis values (not used if time axis)
   type(domain2d),pointer :: domain =>NULL() !< domain associated with compressed axis

end type ax_type

!> @ingroup fms_io_mod
type, private :: var_type
   private
   character(len=128)                     :: name = ''
   character(len=128)                     :: longname = ''
   character(len=128)                     :: units = ''
   real, dimension(:,:,:,:), allocatable :: buffer
   logical                                :: domain_present = .FALSE.
   integer                                :: domain_idx = -1
   logical                                :: is_dimvar = .FALSE.
   logical                                :: read_only = .FALSE.
   logical                                :: owns_data = .FALSE. !< if true, restart owns the
                                                                 !! data and will deallocate them when freed
   type(fieldtype)                        :: field
   type(axistype)                         :: axis
   integer                                :: position
   integer                                :: ndim
   integer                                :: siz(5)      !< X/Y/Z/T/A extent of fields (data domain
                                                         !< size for distributed writes;global size for reads)
   integer                                :: gsiz(4)     !< global X/Y/Z/A extent of fields
   integer                                :: id_axes(4)  !< store index for x/y/z/a axistype.
   logical                                :: initialized !< indicate if the field is read or not in routine save_state.
   logical                                :: mandatory   !< indicate if the field is mandatory to be when restart.
   integer                                :: is, ie, js, je  !< index of the data in compute domain
   real                                   :: default_data
   character(len=8)                       :: compressed_axis !<< If on a compressed axis, which axis
   integer, dimension(:), allocatable     :: pelist
   integer                                :: ishift, jshift !< can be used to shift indices when no_domain=T
   integer                                :: x_halo, y_halo !< can be used to indicate halo size when no_domain=T

    integer(INT_KIND),dimension(5)    :: field_dimension_order !< Array telling the ordering
                                                               !! of the dimensions for the field.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !< Array of sizes of the dimensions for the field.

end type var_type

!> @ingroup fms_io_mod
type Ptr0Dr
   real,                   pointer :: p => NULL()
end type Ptr0Dr

!> @ingroup fms_io_mod
type Ptr1Dr
   real, dimension(:),     pointer :: p => NULL()
end type Ptr1Dr

!> @ingroup fms_io_mod
type Ptr2Dr
   real, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Dr

!> @ingroup fms_io_mod
type Ptr3Dr
   real, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Dr

!> @ingroup fms_io_mod
type Ptr2Dr8
   real(DOUBLE_KIND), dimension(:,:),   pointer :: p => NULL()
end type Ptr2Dr8

!> @ingroup fms_io_mod
type Ptr3Dr8
   real(DOUBLE_KIND), dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Dr8

!> @ingroup fms_io_mod
type Ptr4Dr
   real, dimension(:,:,:,:), pointer :: p => NULL()
end type Ptr4Dr

!> @ingroup fms_io_mod
type Ptr0Di
   integer,                   pointer :: p => NULL()
end type Ptr0Di

!> @ingroup fms_io_mod
type Ptr1Di
   integer, dimension(:),     pointer :: p => NULL()
end type Ptr1Di

!> @ingroup fms_io_mod
type Ptr2Di
   integer, dimension(:,:),   pointer :: p => NULL()
end type Ptr2Di

!> @ingroup fms_io_mod
type Ptr3Di
   integer, dimension(:,:,:), pointer :: p => NULL()
end type Ptr3Di

!> @ingroup fms_io_mod
type restart_file_type
   private
   integer                                  :: unit = -1 ! mpp_io unit for netcdf file
   character(len=256)                       :: name = ''
   integer                                  :: register_id = 0
   integer                                  :: nvar = 0
   integer                                  :: natt = 0
   integer                                  :: max_ntime = 0
   logical                                  :: is_root_pe = .FALSE.
   logical                                  :: is_compressed = .FALSE.
   logical                                  :: unlimited_axis = .FALSE.
   integer                                  :: tile_count = 1
   type(ax_type),  allocatable              :: axes(:)  ! Currently define X,Y,Compressed, unlimited and maybe Z
   type(meta_type),                pointer  :: first =>NULL() ! pointer to first additional global metadata element
   type(var_type), dimension(:),   pointer  :: var  => NULL()
   type(Ptr0Dr),   dimension(:,:), pointer  :: p0dr => NULL()
   type(Ptr1Dr),   dimension(:,:), pointer  :: p1dr => NULL()
   type(Ptr2Dr),   dimension(:,:), pointer  :: p2dr => NULL()
   type(Ptr3Dr),   dimension(:,:), pointer  :: p3dr => NULL()
   type(Ptr2Dr8),  dimension(:,:), pointer  :: p2dr8 => NULL()
   type(Ptr3Dr8),  dimension(:,:), pointer  :: p3dr8 => NULL()
   type(Ptr4Dr),   dimension(:,:), pointer  :: p4dr => NULL()
   type(Ptr0Di),   dimension(:,:), pointer  :: p0di => NULL()
   type(Ptr1Di),   dimension(:,:), pointer  :: p1di => NULL()
   type(Ptr2Di),   dimension(:,:), pointer  :: p2di => NULL()
   type(Ptr3Di),   dimension(:,:), pointer  :: p3di => NULL()
end type restart_file_type

!> Read data from a file
!> @ingroup fms_io_mod
interface read_data
   module procedure read_data_text
end interface

!> @ingroup fms_io_mod
interface query_initialized
   module procedure query_initialized_id
   module procedure query_initialized_name
   module procedure query_initialized_r2d
   module procedure query_initialized_r3d
   module procedure query_initialized_r4d
end interface

!> @ingroup fms_io_mod
interface get_global_att_value
  module procedure get_global_att_value_text
  module procedure get_global_att_value_real
end interface

!> @ingroup fms_io_mod
interface get_var_att_value
  module procedure get_var_att_value_text
end interface

!> @ingroup fms_io_mod
interface parse_mask_table
  module procedure parse_mask_table_2d
  module procedure parse_mask_table_3d
end interface

!> @ingroup fms_io_mod
interface get_mosaic_tile_file
  module procedure get_mosaic_tile_file_sg
end interface

!> @addtogroup fms_io_mod
!> @{

integer :: num_files_r = 0 !< number of currently opened files for reading
integer :: num_files_w = 0 !< number of currently opened files for writing
integer :: num_domains = 0 !< number of domains in array_domain
integer :: num_registered_files = 0 !< mumber of files registered by calling register_restart_file

integer :: thread_r, form
logical :: module_is_initialized = .FALSE.

character(len=128):: error_msg
logical           :: great_circle_algorithm=.FALSE.

!------ private data, pointer to current 2d domain ------
! entrained from fms_mod.  This will be deprecated in the future.
type(domain2D), pointer, private :: Current_domain =>NULL()

integer, private :: is,ie,js,je      !< compute domain
integer, private :: isd,ied,jsd,jed  !< data domain
integer, private :: isg,ieg,jsg,jeg  !< global domain
character(len=128),      dimension(:), allocatable         :: registered_file !< file names
                                                                            !! registered through register_restart_file
type(restart_file_type), dimension(:), allocatable         :: files_read  !< store files that are read
                                                                          !! through read_data
type(restart_file_type), dimension(:), allocatable, target :: files_write !< store files that
                                                                          !! are written through write_data
type(domain2d), dimension(max_domains), target, save  :: array_domain
type(domain1d), dimension(max_domains), save       :: domain_x, domain_y
public  :: read_data
public  :: fms_io_init, fms_io_exit
public  :: open_namelist_file, open_restart_file, open_ieee32_file, close_file
public  :: open_file, open_direct_file
public  :: string
public  :: get_mosaic_tile_file, get_file_name
public  :: get_global_att_value, get_var_att_value
public  :: file_exist, field_exist
public  :: restart_file_type, query_initialized
private :: lookup_field_r, lookup_axis, unique_axes
public  :: set_filename_appendix, get_instance_filename
public  :: get_filename_appendix, nullify_filename_appendix
public  :: parse_mask_table
public  :: get_great_circle_algorithm
character(len=32), save :: filename_appendix = ''

!--- public interface ---
!> @}
!> @ingroup fms_io_mod
interface string
   module procedure string_from_integer
   module procedure string_from_real
end interface
!> @addtogroup fms_io_mod
!> @{

!--- namelist interface
logical           :: fms_netcdf_override = .true.
logical           :: fms_netcdf_restart  = .true.
character(len=32) :: threading_read      = 'multi'
character(len=32) :: format              = 'netcdf'
logical           :: read_all_pe         = .TRUE.
integer           :: max_files_w         = 40
integer           :: max_files_r         = 40
integer           :: dr_set_size         = 10
logical           :: read_data_bug       = .false.
logical           :: time_stamp_restart  = .true.
logical           :: print_chksum        = .false.
logical           :: show_open_namelist_file_warning = .false.
logical           :: debug_mask_list     = .false.
logical           :: checksum_required   = .true.
  namelist /fms_io_nml/ fms_netcdf_override, fms_netcdf_restart, &
       threading_read, format, read_all_pe, max_files_w,max_files_r, &
       read_data_bug, time_stamp_restart, print_chksum, show_open_namelist_file_warning, &
       debug_mask_list, checksum_required, dr_set_size

integer            :: pack_size  ! = 1 for double = 2 for float

! Include variable "version" to be written to log file.
#include<file_version.h>

! make version public so it can be written in fms_init()
character(len=*), parameter, public :: fms_io_version = version


!> @addtogroup fms_io_mod
!> @{
contains

!.....................................................................
! <SUBROUTINE NAME="fms_io_init">
!   <DESCRIPTION>
! Initialize fms_io module
!   </DESCRIPTION>
!   <TEMPLATE>
! call fms_io_init()
!   </TEMPLATE>
subroutine fms_io_init()

  integer                            :: i, unit, io_status, logunit
  integer, allocatable, dimension(:) :: pelist
  real(DOUBLE_KIND)                  :: doubledata = 0
  real                               :: realarray(4)
  character(len=256)                 :: grd_file, filename
  logical                            :: is_mosaic_grid
  character(len=4096)                :: attvalue

  if (module_is_initialized) return
  call mpp_io_init()

  read (input_nml_file, fms_io_nml, iostat=io_status)
  if (io_status > 0) then
     call mpp_error(FATAL,'=>fms_io_init: Error reading input nml file')
  endif

! take namelist options if present
! read_data_bug is no longer supported.
  if (read_data_bug) then
    call mpp_error(FATAL, "fms_io_init: You have overridden the default value of " // &
       "read_data_bug and set it to .true. in fms_io_nml. This was a temporary workaround " // &
       "that is no longer supported. Please remove this namelist variable.")
  endif

! determine packsize
  pack_size = size(transfer(doubledata, realarray))
  if( pack_size .NE. 1 .AND. pack_size .NE. 2) call mpp_error(FATAL,'=>fms_io_init: pack_size should be 1 or 2')

  select case (threading_read)
  case ('multi')
     thread_r = MPP_MULTI
  case ('single')
     thread_r = MPP_SINGLE
  case default
     call mpp_error(FATAL,'fms_io_init: threading_read should be multi/single but you chose'//trim(threading_read))
  end select
! take namelist options if present

  select case(format)
  case ('netcdf')
     form=MPP_NETCDF
  case default
     call mpp_error(FATAL,'fms_io_init: only NetCDF format currently supported in fms_io')
  end select

! Initially allocate  files_write and files_read
  if (.not. allocated(files_write) ) allocate(files_write(max_files_w))
  if (.not. allocated(files_read) ) allocate(files_read(max_files_r))
  if (.not. allocated(registered_file)) allocate(registered_file(max_files_w))

  do i = 1, max_domains
     array_domain(i) = NULL_DOMAIN2D
  enddo

  !---- initialize module domain2d pointer ----
  nullify (Current_domain)

  !This is set here instead of at the end of the routine to prevent the read_data call below from stopping the model
  module_is_initialized = .TRUE.

  !--- read INPUT/grid_spec.nc to decide the value of great_circle_algorithm
  !--- great_circle_algorithm could be true only for mosaic grid.
  great_circle_algorithm = .false.
  grd_file = "INPUT/grid_spec.nc"

  is_mosaic_grid = .FALSE.
  if (file_exist(grd_file)) then
     if(field_exist(grd_file, 'atm_mosaic_file')) then  ! coupled grid
        is_mosaic_grid = .TRUE.
     else if(field_exist(grd_file, "gridfiles")) then
        call read_data(grd_file, "gridfiles", filename, level=1)
        grd_file = 'INPUT/'//trim(filename)
        is_mosaic_grid = .TRUE.
     endif
  endif

  if(is_mosaic_grid) then
     if( get_global_att_value(grd_file, "great_circle_algorithm", attvalue) ) then
        if(trim(attvalue) == "TRUE") then
           great_circle_algorithm = .true.
        else if(trim(attvalue) == "FALSE") then
           great_circle_algorithm = .false.
        else
           call mpp_error(FATAL, "fms_io(fms_io_init: value of global attribute great_circle_algorithm in file"// &
             trim(grd_file)//" should be TRUE of FALSE")
        endif
     endif
  endif

  if(great_circle_algorithm .AND. (mpp_pe() == mpp_root_pe()) ) then
     call mpp_error(NOTE,"fms_io_mod: great_circle algorithm will be used in the model run")
  endif

end subroutine fms_io_init

! </SUBROUTINE>
! <SUBROUTINE NAME="fms_io_exit">
!   <DESCRIPTION>
! This routine is called after ALL fields have been written to temporary files
! The result NETCDF files are created here.
!   </DESCRIPTION>
!   <TEMPLATE>
! call fms_io_exit
!   </TEMPLATE>

subroutine fms_io_exit()
    integer                             :: num_x_axes, num_y_axes, num_z_axes
    integer                             :: unit
    real, dimension(max_axis_size)      :: axisdata
    real                                :: tlev
    integer,        dimension(max_axes) :: id_x_axes, siz_x_axes
    integer,        dimension(max_axes) :: id_y_axes, siz_y_axes
    integer,        dimension(max_axes) :: id_z_axes, siz_z_axes
    type(axistype), dimension(max_axes) :: x_axes, y_axes, z_axes
    type(axistype)                      :: t_axes
    type(var_type), pointer, save       :: cur_var=>NULL()
    integer                             :: i, j, k, kk
    character(len=256)                  :: filename
    character(len=10)                   :: axisname
    logical                             :: domain_present
    logical                             :: write_on_this_pe
    type(domain2d), pointer :: io_domain =>NULL()

    if( .NOT.module_is_initialized )return !make sure it's only called once per PE

    do i=1,max_axis_size
       axisdata(i) = i
    enddo

    ! each field has an associated domain type (may be undefined).
    ! each file only needs to write unique axes (i.e. if 2 fields share an identical axis,
    ! then only write the axis once)
    ! unique axes are defined by the global size and domain decomposition (i.e. can support identical axis sizes with
    ! different domain decomposition)

    do i = 1, num_files_w
       filename = files_write(i)%name

       !--- check if any field in this file present domain.
       domain_present = .false.
       do j = 1, files_write(i)%nvar
          if (files_write(i)%var(j)%domain_present) then
              domain_present = .true.
              exit
          end if
       end do

       !--- get the unique axes for all the fields.
       num_x_axes = unique_axes(files_write(i), 1, id_x_axes, siz_x_axes, domain_x)
       num_y_axes = unique_axes(files_write(i), 2, id_y_axes, siz_y_axes, domain_y)
       num_z_axes = unique_axes(files_write(i), 3, id_z_axes, siz_z_axes          )

       if( domain_present ) then
          call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form, &
               is_root_pe=files_write(i)%is_root_pe, domain=array_domain(files_write(i)%var(j)%domain_idx))
       else  ! global data
          call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form,threading=MPP_SINGLE,&
               fileset=MPP_SINGLE, is_root_pe=files_write(i)%is_root_pe)
       end if

       write_on_this_pe = .false.
       if(domain_present) then
          io_domain => mpp_get_io_domain(array_domain(files_write(i)%var(j)%domain_idx))
          if(associated(io_domain)) then
             if(mpp_domain_is_tile_root_pe(io_domain)) write_on_this_pe = .true.
          endif
       endif
       !--- always write out from root pe
       if( files_write(i)%is_root_pe ) write_on_this_pe = .true.

       do j = 1, num_x_axes
         if (j < 10) then
             write(axisname,'(a,i1)') 'xaxis_',j
          else
             write(axisname,'(a,i2)') 'xaxis_',j
          endif
          if(id_x_axes(j) > 0) then
             call mpp_write_meta(unit,x_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),domain=domain_x(id_x_axes(j)),cartesian='X')
          else
             call mpp_write_meta(unit,x_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),cartesian='X')
          endif
       end do

       do j = 1, num_y_axes
         if (j < 10) then
             write(axisname,'(a,i1)') 'yaxis_',j
          else
             write(axisname,'(a,i2)') 'yaxis_',j
          endif
          if(id_y_axes(j) > 0) then
             call mpp_write_meta(unit,y_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),domain=domain_y(id_y_axes(j)),cartesian='Y')
          else
             call mpp_write_meta(unit,y_axes(j),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),cartesian='Y')
          endif
       end do

       do j = 1, num_z_axes
          if (j < 10) then
             write(axisname,'(a,i1)') 'zaxis_',j
          else
             write(axisname,'(a,i2)') 'zaxis_',j
          endif
          call mpp_write_meta(unit,z_axes(j),axisname,'none',axisname, &
               data=axisdata(1:siz_z_axes(j)),cartesian='Z')
       end do


       ! write time axis  (comment out if no time axis)
       call mpp_write_meta(unit,t_axes,&
            'Time','time level','Time',cartesian='T')

       ! write metadata for fields
       do j = 1, files_write(i)%nvar
          cur_var => files_write(i)%var(j)
          call mpp_write_meta(unit,cur_var%field, (/x_axes(cur_var%id_axes(1)), &
               y_axes(cur_var%id_axes(2)), z_axes(cur_var%id_axes(3)), t_axes/), cur_var%name, &
               'none',cur_var%name,pack=pack_size)
       enddo

       ! write values for ndim of spatial axes
       do j = 1, num_x_axes
          call mpp_write(unit,x_axes(j))
       enddo
       do j = 1, num_y_axes
          call mpp_write(unit,y_axes(j))
       enddo
       do j = 1, num_z_axes
          call mpp_write(unit,z_axes(j))
       enddo

       ! write data of each field
       do k = 1, files_write(i)%max_ntime
          do j = 1, files_write(i)%nvar
             cur_var => files_write(i)%var(j)
             tlev=k
             ! If some fields only have one time level, we do not need to write the second level, just keep
             ! the data missing.
             ! If some fields only have one time level, we just write out 0 to the other level
             if(k > cur_var%siz(4)) then
                cur_var%buffer(:,:,:,1) = 0.0
                kk = 1
             else
                kk = k
             end if
             if(cur_var%domain_present) then
                call mpp_write(unit, cur_var%field,array_domain(cur_var%domain_idx), cur_var%buffer(:,:,:,kk), tlev, &
                               default_data=cur_var%default_data)
             else if (write_on_this_pe) then
                call mpp_write(unit, cur_var%field, cur_var%buffer(:,:,:,kk), tlev)
             end if
          enddo ! end j loop
       enddo ! end k loop
       call mpp_close(unit)
    enddo ! end i loop

    !--- release the memory

    do i = 1,  num_files_w
       do j = 1, files_write(i)%nvar
          deallocate(files_write(i)%var(j)%buffer)
       end do
    end do

  cur_var=>NULL()
  module_is_initialized = .false.
  num_files_w = 0
  num_files_r = 0

end subroutine fms_io_exit
!.....................................................................
! </SUBROUTINE>

!---- return true if all fields in fileObj is read only
function all_field_read_only(fileObj)
  type(restart_file_type), intent(in) :: fileObj
  logical                             :: all_field_read_only
  integer :: j

  all_field_read_only = .TRUE.
  do j = 1, fileObj%nvar
     if( .not. fileObj%var(j)%read_only) then
        all_field_read_only = .FALSE.
        exit
     endif
  enddo

  return

end function all_field_read_only

!-------------------------------------------------------------------------------
!    This subroutine will calculate chksum and print out chksum information.
!
subroutine write_chksum(fileObj, action)
  type(restart_file_type), intent(inout) :: fileObj
  integer,                 intent(in)    :: action
  integer(LONG_KIND)                     :: data_chksum
  integer                                :: j, k, outunit
  integer                                :: isc, iec, jsc, jec
  integer                                :: isg, ieg, jsg, jeg
  integer                                :: ishift, jshift, iadd, jadd
  type(var_type), pointer, save          :: cur_var=>NULL()
  character(len=32)                      :: routine_name

  if(action == MPP_OVERWR) then
     routine_name = "save_restart"
  else if(action == MPP_RDONLY) then
     routine_name = "restore_state"
  else
     call mpp_error(FATAL, "fms_io_mod(write_chksum): action should be MPP_OVERWR or MPP_RDONLY")
  endif

  do j=1,fileObj%nvar
     cur_var => fileObj%var(j)

     if ( cur_var%domain_idx > 0) then
        call mpp_get_compute_domain(array_domain(cur_var%domain_idx), isc, iec, jsc, jec)
        call mpp_get_global_domain(array_domain(cur_var%domain_idx), isg, ieg, jsg, jeg)
        call mpp_get_domain_shift(array_domain(cur_var%domain_idx), ishift, jshift, cur_var%position)
     else if (ASSOCIATED(Current_domain)) then
        call mpp_get_compute_domain(Current_domain, isc, iec, jsc, jec)
        call mpp_get_global_domain(Current_domain, isg, ieg, jsg, jeg)
        call mpp_get_domain_shift(Current_domain, ishift, jshift, cur_var%position)
     else
        iec = cur_var%ie
        isc = cur_var%is
        ieg = cur_var%ie
        jec = cur_var%je
        jsc = cur_var%js
        jeg = cur_var%je
        ishift = 0
        jshift = 0
     endif
     iadd = iec-isc ! Size of the i-dimension on this processor (-1 as it is an increment)
     jadd = jec-jsc ! Size of the j-dimension on this processor
     if(iec == ieg) iadd = iadd + ishift
     if(jec == jeg) jadd = jadd + jshift

     if(action == MPP_OVERWR .OR. (action == MPP_RDONLY .AND. cur_var%initialized) ) then
        do k = 1, cur_var%siz(4)
           if ( Associated(fileObj%p0dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p0dr(k,j)%p, (/mpp_pe()/) )
           else if ( Associated(fileObj%p1dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p1dr(k,j)%p, (/mpp_pe()/) )
           else if ( Associated(fileObj%p2dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p2dr(k,j)%p(cur_var%is:cur_var%is+iadd,cur_var%js:cur_var%js+jadd) )
           else if ( Associated(fileObj%p3dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p3dr(k,j)%p(cur_var%is:cur_var%is+iadd,cur_var%js:cur_var%js+jadd, :) )
           else if ( Associated(fileObj%p4dr(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p4dr(k,j)%p(cur_var%is:cur_var%is+iadd, &
                                      & cur_var%js:cur_var%js+jadd, :, :) )
           else if ( Associated(fileObj%p0di(k,j)%p) ) then
              data_chksum = fileObj%p0di(k,j)%p
           else if ( Associated(fileObj%p1di(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p1di(k,j)%p, (/mpp_pe()/) )
           else if ( Associated(fileObj%p2di(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p2di(k,j)%p(cur_var%is:cur_var%is+iadd,cur_var%js:cur_var%js+jadd) )
           else if ( Associated(fileObj%p3di(k,j)%p) ) then
              data_chksum = mpp_chksum(fileObj%p3di(k,j)%p(cur_var%is:cur_var%is+iadd,cur_var%js:cur_var%js+jadd, :))
           else
              call mpp_error(FATAL, "fms_io(write_chksum): There is no pointer associated with the data of  field "// &
                   trim(cur_var%name)//" of file "//trim(fileObj%name) )
           end if
           outunit = stdout()
           write(outunit,'(a, I1, a, Z16)')'fms_io('//trim(routine_name)//'): At time level = ', k, ', chksum for "'//&
                trim(cur_var%name)// '" of "'// trim(fileObj%name)// '" = ', data_chksum

        enddo
     endif
  enddo
  cur_var =>NULL()

end subroutine write_chksum

function lookup_field_r(nfile,fieldname)
! Given fieldname, this function returns the field position in the model's fields list

  integer, intent(in)          :: nfile
  character(len=*), intent(in) :: fieldname
  integer                      :: lookup_field_r
  integer                      :: j

  lookup_field_r=-1
  do j = 1, files_read(nfile)%nvar
     if (trim(files_read(nfile)%var(j)%name) == trim(fieldname)) then
        lookup_field_r = j
        exit
     endif
  enddo
  return
end function lookup_field_r


!..........................................................

function lookup_domain(domain)
! given domain, this function returns the position of domain in array_domain or -1 if not found

  type(domain2d), intent(in) :: domain
  integer                    :: i, lookup_domain
  lookup_domain = -1
  do i =1, num_domains
     if(domain .EQ. array_domain(i)) then
        lookup_domain = i
        exit
     endif
  enddo
end function lookup_domain
!.........................................................
function lookup_axis(axis_sizes,siz,domains,dom)

! Given axis size (global), this function returns the axis id

  integer, intent(in)      :: axis_sizes(:), siz
  type(domain1d), optional :: domains(:)
  type(domain1d), optional :: dom
  integer :: lookup_axis
  integer :: j


  lookup_axis=-1
  do j=1,size(axis_sizes(:))
     if (siz == axis_sizes(j)) then
        if (PRESENT(domains)) then
           if (dom .EQ. domains(j)) then
              lookup_axis = j
              exit
           endif
        else
           lookup_axis = j
           exit
        endif
     endif
  enddo
  if (lookup_axis == -1) call mpp_error(FATAL,'fms_io(lookup_axis): could not find axis in set of axes')
end function lookup_axis

!=====================================================================================
!--- we assume any text data are at most 2-dimensional and level is for first dimension
subroutine read_data_text(filename,fieldname,data,level)
  character(len=*), intent(in)   :: filename, fieldname
  character(len=*), intent(out)  :: data
  integer, intent(in) , optional :: level
  logical                        :: file_opened, found_file, read_dist, io_domain_exist
  integer                        :: lev, unit, index_field
  integer                        :: file_index
  character(len=256)             :: fname

! Initialize files to default values
  if(.not.module_is_initialized) call mpp_error(FATAL,'fms_io(read_data_text):  module not initialized')

  file_opened=.false.
  if (PRESENT(level)) then
     lev = level
  else
     lev = 1
  endif

  found_file = get_file_name(filename, fname, read_dist, io_domain_exist, no_domain=.true. )
 if(.not.found_file) call mpp_error(FATAL, 'fms_io_mod(read_data_text): file ' //trim(filename)// &
          '(with the consideration of tile number) and corresponding distributed file are not found')
  call get_file_unit(fname, unit, file_index, read_dist, io_domain_exist )

! Get info of this file and field
  call get_field_id(unit, file_index, fieldname, index_field, .true., .true. )

  if ( lev < 1 .or. lev > files_read(file_index)%var(index_field)%siz(1) )  then
     write(error_msg,'(I5,"/",I5)') lev, files_read(file_index)%var(index_field)%siz(1)
     call mpp_error(FATAL,'fms_io(read_data_text): text level out of range, level/max_level=' &
          //trim(error_msg)//' in field/file: '//trim(fieldname)//'/'//trim(filename))
  endif

  call mpp_read(unit,files_read(file_index)%var(index_field)%field,data, level=level)
  return
end subroutine read_data_text
!..............................................................
! </SUBROUTINE>

function unique_axes(file, index, id_axes, siz_axes, dom)
  type(restart_file_type),   intent(inout)           :: file
  integer,                      intent(in)           :: index
  integer, dimension(:),       intent(out)           :: id_axes
  integer, dimension(:),       intent(out)           :: siz_axes
  type(domain1d), dimension(:), intent(in), optional :: dom
  integer                                            :: unique_axes
  type(var_type), pointer, save :: cur_var => NULL()
  integer :: i,j
  logical :: found

  unique_axes=0

  if(index <0 .OR. index > 4) call mpp_error(FATAL,"unique_axes(fms_io_mod): index should be 1, 2, 3 or 4")

  do i = 1, file%nvar
     cur_var => file%var(i)
     if(cur_var%read_only) cycle
     if(cur_var%ndim < index) cycle
     found = .false.
     do j = 1, unique_axes
        if(siz_axes(j) == cur_var%gsiz(index) ) then
           if(PRESENT(dom)) then
              if(cur_var%domain_idx == id_axes(j) ) then
                 found = .true.
                 exit
              else if(cur_var%domain_idx >0 .AND. id_axes(j) >0) then
                 if(dom(cur_var%domain_idx) .EQ. dom(id_axes(j)) ) then
                    found = .true.
                    exit
                 end if
              end if
           else
              found = .true.
              exit
           end if
        end if
     end do
     if(found) then
        cur_var%id_axes(index) = j
     else
        unique_axes = unique_axes+1
        if(unique_axes > max_axes) then
           write(error_msg,'(I3,"/",I3)') unique_axes, max_axes
           if(index == 1 ) then
              call mpp_error(FATAL,'# x axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           else if(index == 2 ) then
              call mpp_error(FATAL,'# y axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           else
              call mpp_error(FATAL,'# z axes exceeded max_axes in fms_io,num_axes/max_axes= '//trim(error_msg))
           end if
        endif
        id_axes(unique_axes)   = cur_var%domain_idx
        siz_axes(unique_axes) = cur_var%gsiz(index)
        if(siz_axes(unique_axes) > max_axis_size) then
           call mpp_error(FATAL, 'fms_io_mod(unique_axes): size_axes is greater than max_axis_size, '//&
              'increase fms_io_nml variable max_axis_size to at least ', siz_axes(unique_axes))
        endif
        cur_var%id_axes(index) = unique_axes
     end if
  end do

  cur_var => NULL()

  return

end function unique_axes

  !#######################################################################
  !#######################################################################
  !   --------- routines for reading distributed data ---------
  ! before calling these routines the domain decompostion must be set
  ! by calling "set_domain" with the appropriate domain2d data type
  !
  ! reading can be done either by all PEs (default) or by only the root PE
  ! this is controlled by namelist variable "read_all_pe".

  ! By default, array data is expected to be declared in data domain and no_halo
  !is NOT needed, however IF data is decalared in COMPUTE domain then optional NO_HALO should be .true.

  !#######################################################################

!#######################################################################
! private routines (read_eof,do_read)
! this routine is called when an EOF is found while
! reading a distributed data file using read_data

subroutine read_eof (end_found)
  logical, intent(out), optional :: end_found

  if (present(end_found))then
     end_found = .true.
  else
     call mpp_error(FATAL,'fms_io(read_eof): unexpected EOF')
  endif
end subroutine read_eof

!#######################################################################
! determines if current pe should read data
! checks namelist variable read_all_pe

function do_read ( )
  logical :: do_read
  do_read = mpp_pe() == mpp_root_pe() .or. read_all_pe
end function do_read

!#######################################################################
!   This function returns .true. if the field referred to by id has
! initialized from a restart file, and .false. otherwise.
!
! Arguments: id - A integer that is the index of the field in fileObj.
!  (in)  fileObj - The control structure returned by a previous call to
!                  register_restart_field
function query_initialized_id(fileObj, id)
  type(restart_file_type), intent(in) :: fileObj
  integer,                 intent(in) :: id

  logical :: query_initialized_id

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_id): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  if(id < 1 .OR. id > fileObj%nvar) call mpp_error(FATAL, "fms_io(query_initialized_id): " // &
      "argument id must be between 1 and nvar in the restart_file_type object")

  query_initialized_id = fileObj%var(id)%initialized

  return

end function query_initialized_id

!#########################################################################
!   This function returns .true. if the field referred to by name has
! initialized from a restart file, and .false. otherwise.
!
! Arguments: name - A pointer to the field that is being queried.
!  (in)  fileObj - The control structure returned by a previous call to
!                  register_restart_field
function query_initialized_name(fileObj, name)
  type(restart_file_type), intent(inout) :: fileObj
  character(len=*),           intent(in) :: name

  logical :: query_initialized_name

  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_name): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_name = .false.
  do m=1,fileObj%nvar
    if (trim(name) == fileObj%var(m)%name) then
      if (fileObj%var(m)%initialized) query_initialized_name = .true.
      exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if ((m>fileObj%nvar) .and. (mpp_pe() == mpp_root_pe())) then
    call mpp_error(NOTE,"fms_io(query_initialized_name): Unknown restart variable "//name// &
                        " queried for initialization.")
  end if

end function query_initialized_name

!#########################################################################
!   This function returns 1 if the field pointed to by f_ptr has
! initialized from a restart file, and 0 otherwise.  If f_ptr is
! NULL, it tests whether the entire restart file has been success-
! fully read.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
function query_initialized_r2d(fileObj, f_ptr, name)
  type(restart_file_type),   intent(inout) :: fileObj
  real, dimension(:,:), target, intent(in) :: f_ptr
  character(len=*),             intent(in) :: name

  logical :: query_initialized_r2d
  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_r2d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_r2d = .false.
  do m=1, fileObj%nvar
     if (ASSOCIATED(fileObj%p2dr(1,m)%p,f_ptr)) then
        if (fileObj%var(m)%initialized) query_initialized_r2d = .true.
        exit
     endif
  enddo
  ! Assume that you are going to initialize it now, so set flag to initialized if
  ! queried again.
  if (m>fileObj%nvar) then
     if (mpp_pe() == mpp_root_pe() ) call mpp_error(NOTE, "fms_io(query_initialized_r2d): Unable to find "// &
          trim(name)//" queried by pointer, "//"probably because of the suspect comparison of pointers by ASSOCIATED.")
     query_initialized_r2d = query_initialized_name(fileObj, name)
     if (mpp_pe() == mpp_root_pe() .AND. query_initialized_r2d) call mpp_error(NOTE, &
          "fms_io(query_initialized_r2d): "//trim(name)// " initialization confirmed by name.")
  endif

  return

end function query_initialized_r2d

!#########################################################################
!   This function returns 1 if the field pointed to by f_ptr has
! initialized from a restart file, and 0 otherwise.  If f_ptr is
! NULL, it tests whether the entire restart file has been success-
! fully read.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
function query_initialized_r3d(fileObj, f_ptr, name)
  type(restart_file_type),     intent(inout) :: fileObj
  real, dimension(:,:,:), target, intent(in) :: f_ptr
  character(len=*),               intent(in) :: name

  logical :: query_initialized_r3d
  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_r3d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_r3d = .false.
  do m=1, fileObj%nvar
     if (ASSOCIATED(fileObj%p3dr(1,m)%p,f_ptr)) then
        if (fileObj%var(m)%initialized) query_initialized_r3d = .true.
        exit
     endif
  enddo
  ! Assume that you are going to initialize it now, so set flag to initialized if
  ! queried again.
  if (m>fileObj%nvar) then
     if (mpp_pe() == mpp_root_pe() ) call mpp_error(NOTE, "fms_io(query_initialized_r3d): Unable to find "// &
          trim(name)//" queried by pointer, "//"probably because of the suspect comparison of pointers by ASSOCIATED.")
     query_initialized_r3d = query_initialized_name(fileObj, name)
     if (mpp_pe() == mpp_root_pe() .AND. query_initialized_r3d) call mpp_error(NOTE, &
          "fms_io(query_initialized_r3d): "//trim(name)// " initialization confirmed by name.")
  endif

  return

end function query_initialized_r3d


!#########################################################################
!   This function returns 1 if the field pointed to by f_ptr has
! initialized from a restart file, and 0 otherwise.  If f_ptr is
! NULL, it tests whether the entire restart file has been success-
! fully read.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
function query_initialized_r4d(fileObj, f_ptr, name)
  type(restart_file_type),       intent(inout) :: fileObj
  real, dimension(:,:,:,:), target, intent(in) :: f_ptr
  character(len=*),                 intent(in) :: name

  logical :: query_initialized_r4d
  integer :: m

  if (.not.associated(fileObj%var)) call mpp_error(FATAL, "fms_io(query_initialized_r4d): " // &
      "restart_file_type data must be initialized by calling register_restart_field before using it")

  query_initialized_r4d = .false.
  do m=1, fileObj%nvar
     if (ASSOCIATED(fileObj%p4dr(1,m)%p,f_ptr)) then
        if (fileObj%var(m)%initialized) query_initialized_r4d = .true.
        exit
     endif
  enddo
  ! Assume that you are going to initialize it now, so set flag to initialized if
  ! queried again.
  if (m>fileObj%nvar) then
     if (mpp_pe() == mpp_root_pe() ) call mpp_error(NOTE, "fms_io(query_initialized_r4d): Unable to find "// &
          trim(name)//" queried by pointer, "//"probably because of the suspect comparison of pointers by ASSOCIATED.")
     query_initialized_r4d = query_initialized_name(fileObj, name)
     if (mpp_pe() == mpp_root_pe() .AND. query_initialized_r4d) call mpp_error(NOTE, &
          "fms_io(query_initialized_r4d): "//trim(name)// " initialization confirmed by name.")
  endif

  return

end function query_initialized_r4d

!#######################################################################
!#######################################################################
!
! routines for opening specific types of files:
!
!                       form        action
! open_namelist_file  MPP_ASCII   MPP_RDONLY
! open restart_file   MPP_NATIVE
! open_ieee32_file    MPP_IEEE32
!
! all have: access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true.
! use the close_file interface to close these files
!
! if other types of files need to be opened the mpp_open and
! mpp_close routines in the mpp_io_mod should be used
!
!#######################################################################


! <FUNCTION NAME="open_namelist_file">
!   <DESCRIPTION>
! Opens single namelist file for reading only by all PEs
! the default file opened is called "input.nml".
!   </DESCRIPTION>
! <IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_namelist_file (file) result (unit)
  character(len=*), intent(in), optional :: file
  integer :: unit
! local variables necessary for nesting code and alternate input.nmls
  character(len=32) :: pelist_name
  character(len=128) :: filename

#ifdef INTERNAL_FILE_NML
  if(show_open_namelist_file_warning) call mpp_error(WARNING, &
     &  "fms_io_mod: open_namelist_file should not be called when INTERNAL_FILE_NML is defined")
#endif

  if (.not.module_is_initialized) call fms_io_init ( )
  if (present(file)) then
     call mpp_open ( unit, file, form=MPP_ASCII, action=MPP_RDONLY, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
  else
!  the following code is necessary for using alternate namelist files (nests, stretched grids, etc)
     pelist_name = mpp_get_current_pelist_name()
     if ( file_exist('input_'//trim(pelist_name)//'.nml', no_domain=.true.) ) then
        filename='input_'//trim(pelist_name)//'.nml'
     else
        filename='input.nml'
     endif
     call mpp_open ( unit, trim(filename), form=MPP_ASCII, action=MPP_RDONLY, &
          access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
  endif
end function open_namelist_file
! </FUNCTION>

! <FUNCTION NAME="open_restart_file">
!   <DESCRIPTION>
! Opens single restart file for reading by all PEs or
! writing by root PE only
! the file has native format and no mpp header records.
!   </DESCRIPTION>
!<IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
!<IN NAME="action" TYPE="character">
! action to be performed: can be 'read' or 'write'
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_restart_file (file, action) result (unit)
  character(len=*), intent(in) :: file, action
  integer :: unit
  integer :: mpp_action

  if (.not.module_is_initialized) call fms_io_init ( )

!   --- action (read,write) ---

  select case (lowercase(trim(action)))
  case ('read')
     mpp_action = MPP_RDONLY
  case ('write')
     mpp_action = MPP_OVERWR
  case default
     call mpp_error(FATAL,'fms_io(open_restart_file): action should be either read or write in file'//trim(file))
  end select

  call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
       access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

end function open_restart_file
! </FUNCTION>


! <FUNCTION NAME="open_direct_file">
!   <DESCRIPTION>
! Opens single direct access file for reading by all PEs or
! writing by root PE only
! the file has native format and no mpp header records.
!   </DESCRIPTION>

  function open_direct_file (file, action, recl) result (unit)
    character(len=*), intent(in) :: file, action
    integer,          intent(in) :: recl
    integer :: unit

    integer :: mpp_action

    if (.not.module_is_initialized) call fms_io_init ( )

    !   --- action (read,write) ---

    select case (lowercase(trim(action)))
    case ('read')
       mpp_action = MPP_RDONLY
    case ('write')
       mpp_action = MPP_OVERWR
    case default
       call mpp_error(FATAL,'invalid option for argument action')
    end select

    call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
         access=MPP_DIRECT, threading=MPP_SINGLE, nohdrs=.true., recl=recl )

  end function open_direct_file
! </FUNCTION>

! <FUNCTION NAME=" open_ieee32_file">
!   <DESCRIPTION>
! Opens single 32-bit ieee file for reading by all PEs or
! writing by root PE only (writing is not recommended)
! the file has no mpp header records.
!   </DESCRIPTION>
!<IN NAME="file" TYPE="character">
! name of the file to be opened
! </IN>
!<IN NAME="action" TYPE="character">
! action to be performed: can be 'read' or 'write'
! </IN>
! <OUT NAME="unit" TYPE="integer">
! unit number returned by this function
! </OUT>
function open_ieee32_file (file, action) result (unit)
  character(len=*), intent(in) :: file, action
  integer :: unit
  integer :: mpp_action

  if (.not.module_is_initialized) call fms_io_init ( )

!   --- action (read,write) ---
  select case (lowercase(trim(action)))
  case ('read')
     mpp_action = MPP_RDONLY
  case ('write')
     mpp_action = MPP_OVERWR
  case default
     call mpp_error (FATAL,'fms_io(open_ieee32_file): action should be either read or write in file'//trim(file))
  end select

  call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
       access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
       nohdrs=.true. )
end function open_ieee32_file
! </FUNCTION>

!#######################################################################
! <FUNCTION NAME=" close_file">
!   <DESCRIPTION>
!  Closes files that are opened by: open_namelist_file, open restart_file,
! and open_ieee32_file. Users should use mpp_close for other cases.
!   </DESCRIPTION>
!<IN NAME="unit" TYPE="integer">
! unit number of the file to be closed
! </IN>
!<IN NAME="status" TYPE="character, optional">
! action to be performed: can be 'delete'
! </IN>

subroutine close_file (unit, status, dist)
  integer,          intent(in)           :: unit
  character(len=*), intent(in), optional :: status
  logical,          intent(in), optional :: dist

  if (.not.module_is_initialized) call fms_io_init ( )
  if(PRESENT(dist))then
    ! If distributed, return if not I/O root
    if(dist)then
      if(.not. mpp_is_dist_ioroot(dr_set_size)) return
    endif
  endif

  if (unit == stdlog()) return
  if (present(status)) then
     if (lowercase(trim(status)) == 'delete') then
        call mpp_close (unit, action=MPP_DELETE)
     else
        call mpp_error(FATAL,'fms_io(close_file): status should be DELETE')
     endif
  else
     call mpp_close (unit)
  endif
end subroutine close_file
! </FUNCTION>

!#######################################################################


! The following function is here as a last resort.
! This is copied from what was utilities_mod in order that redundant code
! could be deleted.

function open_file(file, form, action, access, threading, recl, dist) result(unit)

 character(len=*), intent(in) :: file
 character(len=*), intent(in), optional :: form, action, access, threading
 integer         , intent(in), optional :: recl
 logical         , intent(in), optional :: dist  ! Distributed open?
 integer  :: unit

 character(len=32) :: form_local, action_local, access_local, thread_local
 character(len=32) :: action_ieee32
 logical :: open, no_headers, do_ieee32
 integer :: mpp_format, mpp_action, mpp_access, mpp_thread
!-----------------------------------------------------------------------

   if ( .not. module_is_initialized ) call fms_io_init ( )

   if (present(action)) then    ! must be present
      action_local = action
   else
      call mpp_error (FATAL, 'open_file in fms_mod : argument action not present')
   endif

   unit = 0  ! Initialize return value. Note that mpp_open will call mpi_abort on error
   if(PRESENT(dist))then
     if(lowercase(trim(action_local)) /= 'read') &
       call mpp_error(FATAL,'open_file in fms_mod: distributed'//lowercase(trim(action_local))// &
                              ' not currently supported')
     ! If distributed, return if not I/O root
     if(dist) then
       if(.not. mpp_is_dist_ioroot(dr_set_size)) return
     endif
   endif

!   ---- return stdlog if this is the logfile ----

    if (trim(file) == 'logfile.out') then
       unit = stdlog()
       return
    endif

!   ---- is this file open and connected to a unit ?? ----

   inquire (file=trim(file), opened=open, number=unit)

!  cannot open a file that is already open
!  except for the log file

   if ( open .and. unit >= 0 ) then
      call mpp_error (FATAL, 'open_file in fms_mod : '// &
                       'file '//trim(file)//' is already open')
   endif

!  --- defaults ---

   form_local   = 'formatted';  if (present(form))      form_local   = form
   access_local = 'sequential'; if (present(access))    access_local = access
   thread_local = 'single';     if (present(threading)) thread_local = threading
   no_headers   = .true.
   do_ieee32    = .false.

!   --- file format ---

    select case (lowercase(trim(form_local)))
       case ('formatted')
           mpp_format = MPP_ASCII
       case ('ascii')
           mpp_format = MPP_ASCII
       case ('unformatted')
           mpp_format = MPP_NATIVE
       case ('native')
           mpp_format = MPP_NATIVE
       case ('ieee32')
           do_ieee32 = .true.
       case ('netcdf')
           mpp_format = MPP_NETCDF
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument form')
    end select

!   --- action (read,write,append) ---

    select case (lowercase(trim(action_local)))
       case ('read')
           mpp_action = MPP_RDONLY
       case ('write')
           mpp_action = MPP_OVERWR
       case ('append')
           mpp_action = MPP_APPEND
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument action')
    end select

!   --- file access (sequential,direct) ---

    select case (lowercase(trim(access_local)))
       case ('sequential')
           mpp_access = MPP_SEQUENTIAL
       case ('direct')
           mpp_access = MPP_DIRECT
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument access')
    end select

!   --- threading (single,multi) ---

    select case (lowercase(trim(thread_local)))
       case ('single')
           mpp_thread = MPP_SINGLE
       case ('multi')
           mpp_thread = MPP_MULTI
       case default
           call mpp_error (FATAL, 'open_file in fms_mod : '// &
                            'invalid option for argument thread')
           if (trim(file) /= '_read_error.nml') no_headers = .false.
    end select

!   ---------------- open file -----------------------

    if ( .not.do_ieee32 ) then
       call mpp_open ( unit, file, form=mpp_format, action=mpp_action, &
                       access=mpp_access, threading=mpp_thread,        &
                       fileset=MPP_SINGLE,nohdrs=no_headers, recl=recl )
    else
     ! special open for ieee32 file
     ! pass local action flag to open changing append to write
       action_ieee32 = action_local
       if (lowercase(trim(action_ieee32)) == 'append') action_ieee32 = 'write'
       unit = open_ieee32_file ( file, action_ieee32 )
    endif

!-----------------------------------------------------------------------

 end function open_file

  !#######################################################################

  function string_from_integer(n)
    integer, intent(in) :: n
    character(len=16) :: string_from_integer

    if (mpp_pe() == mpp_root_pe() .and. warn_string_function ) &
            call mpp_error(WARNING, "The function named string has been moved "// &
            "from fms_io_mod to fms_mod.  Please update your call.")
    warn_string_function = .false.
    if(n<0) then
       call mpp_error(FATAL, 'fms_io_mod: n should be non-negative integer, contact developer')
    else if( n<10 ) then
       write(string_from_integer,'(i1)') n
    else if( n<100 ) then
       write(string_from_integer,'(i2)') n
    else if( n<1000 ) then
       write(string_from_integer,'(i3)') n
    else if( n<10000 ) then
       write(string_from_integer,'(i4)') n
    else if( n<100000 ) then
       write(string_from_integer,'(i5)') n
    else if( n<1000000 ) then
       write(string_from_integer,'(i6)') n
    else if( n<10000000 ) then
       write(string_from_integer,'(i7)') n
    else if( n<100000000 ) then
       write(string_from_integer,'(i8)') n
    else
       call mpp_error(FATAL, 'fms_io_mod: n is too big, contact developer')
    end if

    return

  end function string_from_integer

  !#######################################################################
  function string_from_real(a)
    real, intent(in) :: a
    character(len=32) :: string_from_real
    if (mpp_pe() == mpp_root_pe() .and. warn_string_function ) &
            call mpp_error(WARNING, "The function named string has been moved "// &
            "from fms_io_mod to fms_mod.  Please update your call.")
    warn_string_function = .false.

    write(string_from_real,*) a

    return

  end function string_from_real

  !#######################################################################
  subroutine get_mosaic_tile_file_sg(file_in, file_out, is_no_domain, domain, tile_count)
    character(len=*), intent(in)                   :: file_in
    character(len=*), intent(out)                  :: file_out
    logical,          intent(in)                   :: is_no_domain
    type(domain2D),   intent(in), optional, target :: domain
    integer,          intent(in), optional         :: tile_count
    character(len=256)                             :: basefile, tilename
    integer                                        :: lens, ntiles, ntileMe, tile, my_tile_id
    integer, dimension(:), allocatable             :: tile_id
    type(domain2d), pointer, save                  :: d_ptr =>NULL()
    logical                                        :: domain_exist

    if(index(file_in, '.nc', back=.true.)==0) then
       basefile = trim(file_in)
    else
       lens = len_trim(file_in)
       if(file_in(lens-2:lens) .NE. '.nc') call mpp_error(FATAL, &
            'fms_io_mod: .nc should be at the end of file '//trim(file_in))
       basefile = file_in(1:lens-3)
    end if

    !--- get the tile name
    ntiles = 1
    my_tile_id = 1
    domain_exist = .false.
    if(PRESENT(domain))then
       domain_exist = .true.
       ntiles = mpp_get_ntile_count(domain)
       d_ptr => domain
    elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
       domain_exist = .true.
       ntiles = mpp_get_ntile_count(Current_domain)
       d_ptr => Current_domain
    endif

    if(domain_exist) then
       ntileMe = mpp_get_current_ntile(d_ptr)
       allocate(tile_id(ntileMe))
       tile_id = mpp_get_tile_id(d_ptr)
       tile = 1
       if(present(tile_count)) tile = tile_count
       my_tile_id = tile_id(tile)
    endif

    if(ntiles > 1 .or. my_tile_id > 1 )then
       tilename = 'tile'//string(my_tile_id)
       if(index(basefile,'.'//trim(tilename),back=.true.) == 0)then
          basefile = trim(basefile)//'.'//trim(tilename);
       end if
    end if
    if(allocated(tile_id)) deallocate(tile_id)

    file_out = trim(basefile)//'.nc'

    d_ptr =>NULL()

  end subroutine get_mosaic_tile_file_sg

  subroutine get_var_att_value_text(file, varname, attname, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: attname
    character(len=*), intent(inout) :: attvalue
    integer                         :: unit

    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_att_value(unit, varname, attname, attvalue)
    call mpp_close(unit)

    return

  end subroutine get_var_att_value_text

  !#############################################################################
  ! return false if the attribute is not found in the file.
  function get_global_att_value_text(file, att, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: att
    character(len=*), intent(inout) :: attvalue
    logical                         :: get_global_att_value_text
    integer                         :: unit, ndim, nvar, natt, ntime, i
    type(atttype), allocatable      :: global_atts(:)

    get_global_att_value_text = .false.
    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)
    do i=1,natt
       if( trim(mpp_get_att_name(global_atts(i))) == trim(att) ) then
          attvalue = trim(mpp_get_att_char(global_atts(i)))
          get_global_att_value_text = .true.
          exit
       end if
    end do
    deallocate(global_atts)

    return

  end function get_global_att_value_text

  !#############################################################################
  ! return false if the attribute is not found in the file.
  function get_global_att_value_real(file, att, attvalue)
    character(len=*), intent(in)    :: file
    character(len=*), intent(in)    :: att
    real,             intent(inout) :: attvalue
    logical                         :: get_global_att_value_real
    integer                         :: unit, ndim, nvar, natt, ntime, i
    type(atttype), allocatable      :: global_atts(:)

    get_global_att_value_real = .false.
    call mpp_open(unit,trim(file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)
    do i=1,natt
       if( trim(mpp_get_att_name(global_atts(i))) == trim(att) ) then
          attvalue = mpp_get_att_real_scalar(global_atts(i))
          get_global_att_value_real = .true.
          exit
       end if
    end do
    deallocate(global_atts)

    return

  end function get_global_att_value_real

  !#############################################################################
  ! This routine will get the actual file name, as well as if read_dist is true or false.
  ! return true if such file exist and return false if not.
  function get_file_name(orig_file, actual_file, read_dist, io_domain_exist, no_domain, domain, &
                           tile_count)
    character(len=*),                 intent(in) :: orig_file
    character(len=*),                intent(out) :: actual_file
    logical,                         intent(out) :: read_dist
    logical,                         intent(out) :: io_domain_exist
    logical,                optional, intent(in) :: no_domain
    type(domain2D), target, optional, intent(in) :: domain
    integer,                optional, intent(in) :: tile_count
    logical                                      :: get_file_name

    type(domain2d), pointer, save :: d_ptr, io_domain
    logical                       :: fexist, is_no_domain
    integer                       :: tile_id(1)
    character(len=256)            :: fname
    character(len=512)            :: actual_file_tmp

    is_no_domain=.false.
    if(PRESENT(no_domain)) is_no_domain = no_domain


    fexist          = .false.
    read_dist       = .false.
    get_file_name   = .false.
    io_domain_exist = .false.

  !--- The file maybe not netcdf file, we just check the original file.
    if(index(orig_file, '.nc', back=.true.) == 0) then
       inquire (file=trim(orig_file), exist=fexist)
       if(fexist) then
          actual_file = orig_file
          get_file_name = .true.
          return
       endif
    endif

    if(present(domain)) then
       d_ptr => domain
    elseif (ASSOCIATED(Current_domain) .AND. .NOT. is_no_domain ) then
       d_ptr => Current_domain
    endif


    !JWD:  This is likely a temporary fix. Since fms_io needs to know tile_count,
    !JWD:  I just don't see how the physics can remain "tile neutral"
    call get_mosaic_tile_file(orig_file, actual_file, is_no_domain, domain, tile_count)

    !--- check if the file is group redistribution.
    if(ASSOCIATED(d_ptr)) then
       io_domain => mpp_get_io_domain(d_ptr)
       if(associated(io_domain)) then
          tile_id = mpp_get_tile_id(io_domain)
          write(fname, '(a,i4.4)' ) trim(actual_file)//'.', tile_id(1)
          inquire (file=trim(fname), exist=fexist)
          if(.not. fexist) then
             write(fname, '(a,i6.6)' ) trim(actual_file)//'.', tile_id(1)
             inquire (file=trim(fname), exist=fexist)
          endif
          if(fexist) io_domain_exist = .true.
       endif
       io_domain=>NULL()
    endif

    if(fexist) then
       read_dist = .true.
       d_ptr => NULL()
       get_file_name = .true.
       return
    endif

    inquire (file=trim(actual_file), exist=fexist)
    if(fexist) then
       d_ptr => NULL()
       get_file_name = .true.
       return
    endif

    !Perhaps the file has an ensemble instance appendix
    if(len_trim(filename_appendix) > 0) then
       call get_instance_filename(orig_file, actual_file)
       if(index(orig_file, '.nc', back=.true.) == 0) then
          inquire (file=trim(actual_file), exist=fexist)
          if(fexist) then
             d_ptr => NULL()
             get_file_name = .true.
             return
          endif
       endif

       ! Set actual_file to tmp for passing to get_mosaic_tile_file
       actual_file_tmp = actual_file
       call get_mosaic_tile_file(actual_file_tmp, actual_file, is_no_domain, domain, tile_count)

       !--- check if the file is group redistribution.
       if(ASSOCIATED(d_ptr)) then
          io_domain => mpp_get_io_domain(d_ptr)
          if(associated(io_domain)) then
             tile_id = mpp_get_tile_id(io_domain)
             if(mpp_npes()>10000) then
                write(fname, '(a,i6.6)' ) trim(actual_file)//'.', tile_id(1)
             else
                write(fname, '(a,i4.4)' ) trim(actual_file)//'.', tile_id(1)
             endif
             inquire (file=trim(fname), exist=fexist)
             if(fexist) io_domain_exist = .true.
          endif
          io_domain=>NULL()
       endif

       if(fexist) then
          read_dist = .true.
          d_ptr => NULL()
          get_file_name = .true.
          return
       endif

       inquire (file=trim(actual_file), exist=fexist)

       if(fexist) then
          d_ptr => NULL()
          get_file_name = .true.
          return
       endif
    endif

  end function get_file_name


  !#############################################################################
  subroutine get_file_unit(filename, unit, index_file, read_dist, io_domain_exist, domain )
    character(len=*),         intent(in) :: filename
    integer,                 intent(out) :: unit, index_file
    logical,                  intent(in) :: read_dist, io_domain_exist
    type(domain2d), optional, intent(in) :: domain

    logical  :: file_opened
    integer  :: i

    ! Need to check if filename has been opened or not
    file_opened=.false.
    do i=1,num_files_r
       if (files_read(i)%name == trim(filename))  then
          index_file = i
          unit = files_read(index_file)%unit
          return
       endif
    enddo

    ! need to open the file now
    ! Increase num_files_r and set file_type
    if(num_files_r == max_files_r) &  ! need to have bigger max_files_r
         call mpp_error(FATAL,'fms_io(get_file_unit): max_files_r exceeded, increase it via fms_io_nml')
    num_files_r=num_files_r + 1
    if(read_dist) then
       if(io_domain_exist) then
          if(present(domain)) then
             call mpp_open(unit,filename,form=form,action=MPP_RDONLY,threading=MPP_MULTI, &
                fileset=MPP_MULTI, domain=domain)
          else if(ASSOCIATED(current_domain) ) then
             call mpp_open(unit,filename,form=form,action=MPP_RDONLY,threading=MPP_MULTI, &
                fileset=MPP_MULTI, domain=current_domain)
          else
             call mpp_error(FATAL,'fms_io(get_file_unit): when io_domain_exsit = .true., '// &
                   'either domain is present or current_domain is associated')
          endif
       else
          call mpp_open(unit,trim(filename),form=form,action=MPP_RDONLY,threading=MPP_MULTI, &
            fileset=MPP_MULTI)
       endif
    else
       call mpp_open(unit,trim(filename),form=form,action=MPP_RDONLY,threading=MPP_MULTI, &
            fileset=MPP_SINGLE)
    end if
    files_read(num_files_r)%name = trim(filename)
    allocate(files_read(num_files_r)%var (max_fields) )
    files_read(num_files_r)%nvar = 0
    index_file = num_files_r
    files_read(index_file)%unit = unit

  end subroutine get_file_unit

  !#############################################################################
  subroutine get_field_id(unit, index_file, fieldname, index_field, is_no_domain, is_not_dim)
    integer,          intent(in) :: unit
    integer,          intent(in) :: index_file
    character(len=*), intent(in) :: fieldname
    integer,         intent(out) :: index_field
    logical,          intent(in) :: is_no_domain
    logical,          intent(in) :: is_not_dim

    character(len=128)                     :: name
    type(axistype),  dimension(max_axes)   :: axes
    type(fieldtype), dimension(max_fields) :: fields
    integer                                :: i, j, ndim, nvar, natt, var_dim
    integer                                :: siz_in(4)

    index_field = -1
    do j = 1, files_read(index_file)%nvar
       if (trim(files_read(index_file)%var(j)%name) == trim(fieldname)) then
          index_field = j
          return
       endif
    enddo

    !--- fieldname is not read, so need to get fieldname from file
    files_read(index_file)%nvar = files_read(index_file)%nvar + 1
    if(files_read(index_file)%nvar > max_fields) then
       write(error_msg,'(I3,"/",I3)') files_read(index_file)%nvar, max_fields
       call  mpp_error(FATAL,'fms_io(get_field_id): max_fields exceeded, needs increasing, nvar/max_fields=' &
            //trim(error_msg))
    endif
    call mpp_get_info(unit, ndim, nvar, natt, files_read(index_file)%max_ntime)
    if(files_read(index_file)%max_ntime < 1)  files_read(index_file)%max_ntime = 1
    if(nvar > max_fields) then
       write(error_msg,'(I3,"/",I3)') files_read(index_file)%nvar,max_fields
       call mpp_error(FATAL,'fms_io(get_field_id): max_fields too small needs increasing,nvar/max_fields=' &
            //trim(error_msg)//'in file'//trim(files_read(index_file)%name))
    endif
    call mpp_get_fields(unit, fields(1:nvar))
    siz_in = 1
    index_field = files_read(index_file)%nvar
    files_read(index_file)%var(index_field)%is_dimvar = .false.

    do i=1, nvar
       call mpp_get_atts(fields(i),name=name,ndim=var_dim,siz=siz_in)
       if(var_dim .GT. 4) call mpp_error(FATAL, 'fms_io(get_field_id): number of dimension of field '// &
                trim(name)//' in file '//trim(files_read(index_file)%name)//' should not be greater than 4')
       if (lowercase(trim(name)) == lowercase(trim(fieldname))) then ! found the variable
          if(var_dim .lt.3) then
             do j=var_dim+1,3
                siz_in(j)=1
             enddo
          endif
          files_read(index_file)%var(index_field)%name    = fieldname
          files_read(index_file)%var(index_field)%field   = fields(i)
          files_read(index_file)%var(index_field)%siz(1:4)  = siz_in(1:4)
          files_read(index_file)%var(index_field)%gsiz(1:3) = siz_in(1:3)
          return
       endif
    enddo

    !--- the fieldname may be a dimension variable.
    if( .not. is_not_dim) then
       if (ndim > max_axes) then
          write(error_msg,'(I3,"/",I3)') ndim, max_axes
          call  mpp_error(FATAL,'fms_io(get_field_id): max_axes exceeded, needs increasing, ndim/max_fields=' &
               //trim(error_msg)//' in file '//trim(files_read(index_file)%name))
       endif
       call mpp_get_axes(unit, axes(1:ndim))
       do i=1,ndim
          call mpp_get_atts(axes(i), name=name, len = siz_in(1))
          if (lowercase(trim(name)) == lowercase(trim(fieldname))) then
!             if(.not. is_no_domain) call mpp_error(FATAL, &
!                  'fms_io(get_field_id): the field is a dimension variable, no_domain should be true.')
             files_read(index_file)%var(index_field)%is_dimvar = .true.
             files_read(index_file)%var(index_field)%name      = fieldname
             files_read(index_file)%var(index_field)%axis      = axes(i)
             files_read(index_file)%var(index_field)%siz(1:4)    = siz_in(1:4)
             files_read(index_file)%var(index_field)%gsiz(1:3)   = siz_in(1:3)
             return
          endif
       enddo
    end if
    !--- the field is not in the file when reaching here.
    call mpp_error(FATAL, 'fms_io(get_field_id): field '//trim(fieldname)// &
                   ' NOT found in file '//trim(files_read(index_file)%name))

  end subroutine get_field_id

!#######################################################################
! check the existence of the given file name
! if the file_name string has zero length or the
! first character is blank return a false result
! <FUNCTION NAME="file_exist">

!   <OVERVIEW>
!     Checks the existence of a given file name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Checks the existence of the given file name.
!     If the file_name string has zero length or the
!     first character is blank return a false result.
!   </DESCRIPTION>
!   <TEMPLATE>
!     file_exist ( file_name )
!   </TEMPLATE>

!   <IN NAME="file_name"  TYPE="character" >
!     A file name (or path name) that is checked for existence.
!   </IN>
!   <OUT NAME=""  TYPE="logical" >
!     This function returns a logical result.  If file_name exists the result
!     is true, otherwise false is returned.
!     If the length of character string "file_name" is zero or the first
!     character is blank, then the returned value will be false.
!     When reading a file, this function is often used in conjunction with
!     routine open_file.
!   </OUT>
!   <ERROR MSG="set_domain not called" STATUS="FATAL">
!     Before calling write_data you must first call set_domain with domain2d data
!     type associated with the distributed data you are writing.
!   </ERROR>

 function file_exist (file_name, domain, no_domain)
  character(len=*), intent(in)         :: file_name
  type(domain2d), intent(in), optional :: domain
  logical,        intent(iN), optional :: no_domain

  logical                              :: file_exist, is_no_domain
  character(len=256)                   :: fname
  logical                              :: read_dist, io_domain_exist

  is_no_domain = .false.
  if(present(no_domain)) is_no_domain = no_domain
   !--- to deal with mosaic file, in this case, the file is assumed to be in netcdf format
   file_exist = get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=is_no_domain, domain=domain)
   if(is_no_domain) return
   if(.not.file_exist) file_exist=get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=.true.)

   return

 end function file_exist
! </FUNCTION>


!#######################################################################
! <FUNCTION NAME="field_exist">

!   <OVERVIEW>
!     check if a given field name exists in a given file name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     check if a given field name exists in a given file name.
!     If the field_name string has zero length or the
!     first character is blank return a false result.
!     if the file file_name don't exist, return a false result.
!   </DESCRIPTION>
!   <TEMPLATE>
!     field_exist ( file_name, field_name )
!   </TEMPLATE>

!   <IN NAME="file_name"  TYPE="character" >
!     A file name (or path name) that is checked for existence.
!   </IN>
!   <IN NAME="field_name"  TYPE="character" >
!     A field name that is checked for existence.
!   </IN>
!   <OUT NAME=""  TYPE="logical" >
!     This function returns a logical result.  If field exists in the
!     file file_name, the result is true, otherwise false is returned.
!     If the length of character string "field_name" is zero or the first
!     character is blank, then the returned value will be false.
!     if the file file_name don't exist, return a false result.
!   </OUT>

 function field_exist (file_name, field_name, domain, no_domain)
  character(len=*),                 intent(in) :: file_name
  character(len=*),                 intent(in) :: field_name
  type(domain2d), intent(in), optional, target :: domain
  logical,       intent(in),  optional         :: no_domain
  logical                      :: field_exist, is_no_domain
  integer                      :: unit, ndim, nvar, natt, ntime, i, nfile
  character(len=64)            :: name
  type(fieldtype), allocatable :: fields(:)
  logical                      :: file_exist, read_dist, io_domain_exist
  character(len=256)           :: fname

   field_exist = .false.
   if (len_trim(field_name) == 0) return
   if (field_name(1:1) == ' ')    return

   is_no_domain = .false.
   if(present(no_domain)) is_no_domain = no_domain

   file_exist=get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=is_no_domain, domain=domain)
   if(file_exist) then
      call get_file_unit(fname, unit, nfile, read_dist, io_domain_exist, domain=domain)
      call mpp_get_info(unit, ndim, nvar, natt, ntime)
      allocate(fields(nvar))
      call mpp_get_fields(unit,fields)

      do i=1, nvar
         call mpp_get_atts(fields(i),name=name)
         if(lowercase(trim(name)) == lowercase(trim(field_name))) field_exist = .true.
      enddo
      deallocate(fields)
    endif
    if(field_exist .or. is_no_domain) return
    file_exist =  get_file_name(file_name, fname, read_dist, io_domain_exist, no_domain=.true.)
    if(file_exist) then
       call get_file_unit(fname, unit, nfile, read_dist, io_domain_exist)
       call mpp_get_info(unit, ndim, nvar, natt, ntime)
       allocate(fields(nvar))
       call mpp_get_fields(unit,fields)
       do i=1, nvar
          call mpp_get_atts(fields(i),name=name)
          if(lowercase(trim(name)) == lowercase(trim(field_name))) field_exist = .true.
       enddo
       deallocate(fields)
    endif

    return

 end function field_exist
! </FUNCTION>


subroutine get_filename_appendix(string_out)
  character(len=*) , intent(out) :: string_out

  string_out = trim(filename_appendix)


end subroutine get_filename_appendix


subroutine nullify_filename_appendix()

  filename_appendix = ''

end subroutine nullify_filename_appendix


subroutine set_filename_appendix(string_in)
  character(len=*) , intent(in) :: string_in

  integer :: index_num

  ! Check if string has already been added
  index_num = index(filename_appendix, string_in)
  if ( index_num .le. 0 ) then
     filename_appendix = trim(filename_appendix)//trim(string_in)
  end if

end subroutine set_filename_appendix

subroutine get_instance_filename(name_in,name_out)
  character(len=*)  , intent(in)  :: name_in
  character(len=*), intent(inout) :: name_out
  integer :: length

  length = len_trim(name_in)
  name_out = name_in(1:length)

  if(len_trim(filename_appendix) > 0) then
     if(name_in(length-2:length) == '.nc') then
        name_out = name_in(1:length-3)//'.'//trim(filename_appendix)//'.nc'
     else
        name_out = name_in(1:length)  //'.'//trim(filename_appendix)
     end if
  end if

end subroutine get_instance_filename

!#######################################################################
subroutine parse_mask_table_2d(mask_table, maskmap, modelname)

  character(len=*), intent(in) :: mask_table
  logical,         intent(out) :: maskmap(:,:)
  character(len=*), intent(in) :: modelname
  integer                      :: nmask, layout(2)
  integer, allocatable         :: mask_list(:,:)
  integer                      :: unit, mystat, n, stdoutunit
  character(len=128)           :: record

  maskmap = .true.
  nmask = 0
  stdoutunit = stdout()
  if( mpp_pe() == mpp_root_pe() ) then
     call mpp_open(unit, mask_table, action=MPP_RDONLY)
     read(unit, FMT=*, IOSTAT=mystat) nmask
     if( mystat /= 0 ) call mpp_error(FATAL, &
          "fms_io(parse_mask_table_2d): Error reading nmask from file " //trim(mask_table))
     write(stdoutunit,*)"parse_mask_table: Number of domain regions masked in ", trim(modelname), " = ", nmask
     if( nmask > 0 ) then
        !--- read layout from mask_table and confirm it matches the shape of maskmap
        read(unit, FMT=*, IOSTAT=mystat) layout
        if( mystat /= 0 ) call mpp_error(FATAL, &
             "fms_io(parse_mask_talbe_2d): Error reading layout from file " //trim(mask_table))
        if( (layout(1) .NE. size(maskmap,1)) .OR. (layout(2) .NE. size(maskmap,2)) )then
           write(stdoutunit,*)"layout=", layout, ", size(maskmap) = ", size(maskmap,1), size(maskmap,2)
           call mpp_error(FATAL, "fms_io(parse_mask_table_2d): layout in file "//trim(mask_table)// &
                  "does not match size of maskmap for "//trim(modelname))
        endif
        !--- make sure mpp_npes() == layout(1)*layout(2) - nmask
        if( mpp_npes() .NE. layout(1)*layout(2) - nmask ) call mpp_error(FATAL, &
           "fms_io(parse_mask_table_2d): mpp_npes() .NE. layout(1)*layout(2) - nmask for "//trim(modelname))
      endif
   endif

   call mpp_broadcast(nmask, mpp_root_pe())

   if(nmask==0) then
      if( mpp_pe() == mpp_root_pe() ) call mpp_close(unit)
      return
   endif

   allocate(mask_list(nmask,2))

   if( mpp_pe() == mpp_root_pe() ) then
     n = 0
     do while( .true. )
        read(unit,'(a)',end=999) record
        if (record(1:1) == '#') cycle
        if (record(1:10) == '          ') cycle
        n = n + 1
        if( n > nmask ) then
           call mpp_error(FATAL, "fms_io(parse_mask_table_2d): number of mask_list entry "// &
                "is greater than nmask in file "//trim(mask_table) )
        endif
        read(record,*,err=888) mask_list(n,1), mask_list(n,2)
     enddo
888  call mpp_error(FATAL, "fms_io(parse_mask_table_2d):  Error in reading mask_list from file "//trim(mask_table))

999  continue
     !--- make sure the number of entry for mask_list is nmask
     if( n .NE. nmask) call mpp_error(FATAL, &
        "fms_io(parse_mask_table_2d): number of mask_list entry does not match nmask in file "//trim(mask_table))
     call mpp_close(unit)
  endif

  call mpp_broadcast(mask_list, 2*nmask, mpp_root_pe())
  do n = 1, nmask
     if(debug_mask_list) then
       write(stdoutunit,*) "==>NOTE from parse_mask_table_2d: ", trim(modelname), " mask_list = ", mask_list(n,1), &
                         & mask_list(n,2)
     endif
     maskmap(mask_list(n,1),mask_list(n,2)) = .false.
  enddo

  deallocate(mask_list)

end subroutine parse_mask_table_2d


!#######################################################################
subroutine parse_mask_table_3d(mask_table, maskmap, modelname)

  character(len=*), intent(in) :: mask_table
  logical,         intent(out) :: maskmap(:,:,:)
  character(len=*), intent(in) :: modelname
  integer                      :: nmask, layout(2)
  integer, allocatable         :: mask_list(:,:)
  integer                      :: unit, mystat, n, stdoutunit, ntiles
  character(len=128)           :: record

  maskmap = .true.
  nmask = 0
  stdoutunit = stdout()
  if( mpp_pe() == mpp_root_pe() ) then
     call mpp_open(unit, mask_table, action=MPP_RDONLY)
     read(unit, FMT=*, IOSTAT=mystat) nmask
     if( mystat /= 0 ) call mpp_error(FATAL, &
          "fms_io(parse_mask_table_3d): Error reading nmask from file " //trim(mask_table))
     write(stdoutunit,*)"parse_mask_table: Number of domain regions masked in ", trim(modelname), " = ", nmask
     if( nmask > 0 ) then
        !--- read layout from mask_table and confirm it matches the shape of maskmap
        read(unit, FMT=*, IOSTAT=mystat) layout(1), layout(2), ntiles
        if( mystat /= 0 ) call mpp_error(FATAL, &
             "fms_io(parse_mask_talbe_3d): Error reading layout from file " //trim(mask_table))
        if( (layout(1) .NE. size(maskmap,1)) .OR. (layout(2) .NE. size(maskmap,2)) )then
           write(stdoutunit,*)"layout=", layout, ", size(maskmap) = ", size(maskmap,1), size(maskmap,2)
           call mpp_error(FATAL, "fms_io(parse_mask_table_3d): layout in file "//trim(mask_table)// &
                  "does not match size of maskmap for "//trim(modelname))
        endif
        if( ntiles .NE. size(maskmap,3) ) then
           write(stdoutunit,*)"ntiles=", ntiles, ", size(maskmap,3) = ", size(maskmap,3)
           call mpp_error(FATAL, "fms_io(parse_mask_table_3d): ntiles in file "//trim(mask_table)// &
                  "does not match size of maskmap for "//trim(modelname))
        endif
        !--- make sure mpp_npes() == layout(1)*layout(2) - nmask
        if( mpp_npes() .NE. layout(1)*layout(2)*ntiles - nmask ) then
           print*, "layout=", layout, nmask, mpp_npes()
           call mpp_error(FATAL, &
              "fms_io(parse_mask_table_3d): mpp_npes() .NE. layout(1)*layout(2) - nmask for "//trim(modelname))
        endif
      endif
   endif

   call mpp_broadcast(nmask, mpp_root_pe())

   if(nmask==0) then
      if( mpp_pe() == mpp_root_pe() ) call mpp_close(unit)
      return
   endif

   allocate(mask_list(nmask,3))

   if( mpp_pe() == mpp_root_pe() ) then
     n = 0
     do while( .true. )
        read(unit,'(a)',end=999) record
        if (record(1:1) == '#') cycle
        if (record(1:10) == '          ') cycle
        n = n + 1
        if( n > nmask ) then
           call mpp_error(FATAL, "fms_io(parse_mask_table_3d): number of mask_list entry "// &
                "is greater than nmask in file "//trim(mask_table) )
        endif
        read(record,*,err=888) mask_list(n,1), mask_list(n,2), mask_list(n,3)
     enddo
888  call mpp_error(FATAL, "fms_io(parse_mask_table_3d):  Error in reading mask_list from file "//trim(mask_table))

999  continue
     !--- make sure the number of entry for mask_list is nmask
     if( n .NE. nmask) call mpp_error(FATAL, &
        "fms_io(parse_mask_table_3d): number of mask_list entry does not match nmask in file "//trim(mask_table))
     call mpp_close(unit)
  endif

  call mpp_broadcast(mask_list, 3*nmask, mpp_root_pe())
  do n = 1, nmask
     if(debug_mask_list) then
       write(stdoutunit,*) "==>NOTE from parse_mask_table_3d: ", trim(modelname), " mask_list = ", &
                           mask_list(n,1), mask_list(n,2), mask_list(n,3)
     endif
     maskmap(mask_list(n,1),mask_list(n,2),mask_list(n,3)) = .false.
  enddo

  deallocate(mask_list)

end subroutine parse_mask_table_3d


function get_great_circle_algorithm()
   logical :: get_great_circle_algorithm

   if(.NOT. module_is_initialized) call mpp_error(FATAL, &
        "fms_io(use_great_circle_algorithm): fms_io_init is not called yet")

   get_great_circle_algorithm = great_circle_algorithm

end function get_great_circle_algorithm

! </SUBROUTINE>

end module fms_io_mod
!> @}
! close documentation grouping
