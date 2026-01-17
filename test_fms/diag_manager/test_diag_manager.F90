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

! This program runs only one of many possible tests with each execution.
! Each test ends with an intentional fatal error.
! diag_manager_mod is not a stateless module, and there are situations
! where a fatal error leaves the module in a state that does not allow
! it to function properly if used again. Therefore, the program must
! be terminated after each intentional fatal error.

! Each test is dependent on the diag_table, and different diag_tables
! exist for each test. Depending on the test, an intentional fatal error
! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
! Because of this, the calls to all of those routines differ depending on the test.

! The diag_table for each test is included below.

!--------------------------------------------------------------------------------------------------
! diag_table for test 1

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 2

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 3

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 4

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 5

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 6

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 7

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 8

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 9

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 10

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 11

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 12

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
! # Test of the error check that duplicate field names do not appear in same file,
!  "test_mod",              "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 13

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
! # Test of WARNING message that no data is written when run length is less than output interval
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 14

! test_diag_manager
! 1990 1 29 0 0 0
! #output files
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
! # Test of check for invalid date. (Jan 29 1990 + one month = Feb 29 1990)
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 16

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
! # Test for output file name to be modified with appended string
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 17

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2_rms", "diag_test2", "all", "rms",  "none", 2,
!  "test_diag_manager_mod", "dat2", "dat2",     "diag_test2", "all", .true., "none", 2,
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
PROGRAM test
  ! This program runs only one of many possible tests with each execution.
  ! Each test ends with an intentional fatal error.
  ! diag_manager_mod is not a stateless module, and there are situations
  ! where a fatal error leaves the module in a state that does not allow
  ! it to function properly if used again. Therefore, the program must
  ! be terminated after each intentional fatal error.

  ! Each test is dependent on the diag_table, and different diag_tables
  ! exist for each test. Depending on the test, an intentional fatal error
  ! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
  ! Because of this, the calls to all of those routines differ depending on the test.

  USE mpp_mod, ONLY: mpp_pe, mpp_root_pe, mpp_debug, mpp_set_stack_size
  USE mpp_io_mod, ONLY: mpp_io_init
  USE mpp_domains_mod, ONLY: domain2d, mpp_define_domains, mpp_get_compute_domain
  USE mpp_domains_mod, ONLY: mpp_define_io_domain, mpp_define_layout
  USE mpp_domains_mod, ONLY: mpp_domains_init, mpp_domains_set_stack_size
  USE fms_mod, ONLY: fms_init, fms_end, mpp_npes, file_exist, check_nml_error
  USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdlog, stdout
  USE mpp_mod, ONLY: input_nml_file
  USE fms_io_mod, ONLY: fms_io_init
  USE fms_io_mod, ONLY: set_filename_appendix
  USE constants_mod, ONLY: constants_init, PI, RAD_TO_DEG

  USE time_manager_mod, ONLY: time_type, set_calendar_type, set_date, decrement_date, OPERATOR(+), set_time
  USE time_manager_mod, ONLY: NOLEAP, JULIAN, GREGORIAN, THIRTY_DAY_MONTHS, OPERATOR(*), assignment(=)
  use time_manager_mod, ONLY: OPERATOR(+), OPERATOR(-), OPERATOR(/), days_in_month

  USE diag_manager_mod, ONLY: diag_manager_init, send_data, diag_axis_init, diag_manager_end
  USE diag_manager_mod, ONLY: register_static_field, register_diag_field, diag_send_complete
  USE diag_manager_mod, ONLY: diag_manager_set_time_end, diag_field_add_attribute
  USE diag_manager_mod, ONLY: diag_field_add_cell_measures
  USE diag_manager_mod, ONLY: get_diag_field_id, DIAG_FIELD_NOT_FOUND
  USE diag_axis_mod, ONLY: get_axis_num
  USE platform_mod

  IMPLICIT NONE

  TYPE(domain2d) :: Domain1
  TYPE(domain2d) :: Domain2

  REAL, ALLOCATABLE, DIMENSION(:) :: lon_global1, lonb_global1
  REAL, ALLOCATABLE, DIMENSION(:) :: lat_global1, latb_global1
  REAL, ALLOCATABLE, DIMENSION(:) :: lon_global2, lonb_global2
  REAL, ALLOCATABLE, DIMENSION(:) :: lat_global2, latb_global2
  REAL, ALLOCATABLE, DIMENSION(:) :: pfull, bk, phalf
  REAL, ALLOCATABLE, DIMENSION(:) :: lon1, lat1, lonb1, latb1
  REAL, ALLOCATABLE, DIMENSION(:) :: lon2, lat2, lonb2, latb2
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dat1, dat1h
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dat2, dat2h
  REAL, ALLOCATABLE, DIMENSION(:,:) :: dat2_2d
  REAL :: solar_constant = 1600
  REAL :: surf_press = 1.e5
  REAL :: dp
  INTEGER :: id_phalf, id_pfull, id_bk
  INTEGER :: id_lon1, id_lonb1, id_latb1, id_lat1, id_dat1
  INTEGER :: id_lon2, id_lat2, id_dat2, id_dat2_2d, id_sol_con, id_dat2h, id_dat2h_2
  INTEGER :: id_dat2_got, id_none_got
  INTEGER :: i, j, k, is1, ie1, js1, je1, nml_unit, ierr, log_unit, out_unit, m
  INTEGER :: is_in, ie_in, js_in, je_in
  INTEGER :: is2, ie2, js2, je2, hi=1, hj=1
  INTEGER :: nlon1, nlat1, nlon2, nlat2
  INTEGER, DIMENSION(2) :: layout = (/1,1/)
  INTEGER :: test_number=1
  INTEGER :: nlon=10, nlat=10, nlev=10
  INTEGER :: io_layout(2) = (/1,1/)
  INTEGER :: nstep = 2
  TYPE(time_type) :: Time, Time_step, Time_end, Time_start, Run_length
  LOGICAL :: used, test_successful
  CHARACTER(len=256) :: err_msg
  INTEGER :: omp_get_num_threads

  INTEGER :: nyc1, n, jsw, jew, isw, iew
  INTEGER :: numthreads=1, ny_per_thread, idthread
  INTEGER :: months=0, days=1, dt_step=1

  ! Variables needed for test 22
  INTEGER :: id_nv, id_nv_init
  integer(kind=i4_kind)              :: npes             !<Total number of ranks in the current pelist.
!!!!!!



  NAMELIST /test_diag_manager_nml/ layout, test_number, nlon, nlat, nlev, io_layout, numthreads, &
                                   dt_step, months, days
  ! Initialize all id* vars to be -1
  id_nv = -1
  id_nv_init = -1
  id_phalf = -1
  id_pfull = -1
  id_bk = -1
  id_lon1 = -1
  id_lonb1 = -1
  id_latb1 = -1
  id_lat1 = -1
  id_dat1 = -1
  id_lon2 = -1
  id_lat2 = -1
  id_dat2 = -1
  id_dat2_2d = -1
  id_sol_con = -1
  id_dat2h = -1
  id_dat2h_2 = -1
  id_dat2_got = -1
  id_none_got = -1

  CALL fms_init
  log_unit = stdlog()
  out_unit = stdout()
  CALL constants_init
  CALL set_calendar_type(JULIAN)
  npes = mpp_npes()
  READ (input_nml_file, NML=test_diag_manager_nml, IOSTAT=ierr)
  ! Check the status of reading the diag_manager_nml
  IF ( check_nml_error(IOSTAT=ierr, NML_NAME='DIAG_MANAGER_NML') < 0 ) THEN
     IF ( mpp_pe() == mpp_root_pe() ) THEN
        CALL error_mesg('diag_manager_mod::diag_manager_init', &
               & 'TEST_DIAG_MANAGER_NML not found in input.nml.  Using defaults.', WARNING)
     END IF
  END IF
  WRITE (log_unit,test_diag_manager_nml)

SELECT CASE ( test_number ) ! Closes just before the CONTAINS block.
   ! If the test_number == 12, check for the correct error and skip everything else.
   CASE ( 12 )
     CALL diag_manager_init(err_msg=err_msg)
     IF ( err_msg /= '' ) THEN
        WRITE (out_unit,'(a)') 'test12 successful: err_msg='//TRIM(err_msg)
        CALL error_mesg('test_diag_manager','test12 successful.',NOTE)
     ELSE
        WRITE (out_unit,'(a)') 'test12 fails'
        CALL error_mesg('test_diag_manager','test12 fails',FATAL)
     END IF

  ! If the test number is not 12 or 23, run all other tests.
  CASE DEFAULT ! Contains all remaining code up to CONTAINS block.
     CALL diag_manager_init

  IF ( layout(1)*layout(2) .NE. mpp_npes() ) THEN
     CALL mpp_define_layout((/1,nlon,1,nlat/), mpp_npes(), layout )
  END IF

  nlon1 = nlon
  nlat1 = nlat
  nlon2 = nlon * 2
  nlat2 = nlat * 2

  CALL mpp_define_domains((/1,nlon1,1,nlat1/), layout, Domain1, name='test_diag_manager')
  CALL mpp_get_compute_domain(Domain1, is1, ie1, js1, je1)
  ALLOCATE(lon_global1(nlon1), lonb_global1(nlon1+1))
  ALLOCATE(lat_global1(nlat1), latb_global1(nlat1+1))
  ALLOCATE(lon_global2(nlon2), lonb_global2(nlon2+1))
  ALLOCATE(lat_global2(nlat2), latb_global2(nlat2+1))
  ALLOCATE(pfull(nlev), bk(nlev), phalf(nlev+1))

  ALLOCATE(lon1(is1:ie1), lat1(js1:je1), lonb1(is1:ie1+1), latb1(js1:je1+1))
  CALL compute_grid(nlon1, nlat1, is1, ie1, js1, je1, lon_global1, lat_global1, lonb_global1, latb_global1, lon1, &
                  & lat1, lonb1, latb1)
  CALL mpp_define_domains((/1,nlon2,1,nlat2/), layout, Domain2, name='test_diag_manager')
  CALL mpp_get_compute_domain(Domain2, is2, ie2, js2, je2)
  CALL mpp_define_io_domain(Domain1, io_layout)
  CALL mpp_define_io_domain(Domain2, io_layout)

  ALLOCATE(lon2(is2:ie2), lat2(js2:je2), lonb2(is2:ie2+1), latb2(js2:je2+1))
  CALL compute_grid(nlon2, nlat2, is2, ie2, js2, je2, lon_global2, lat_global2, lonb_global2, latb_global2, lon2, &
                  & lat2, lonb2, latb2)
  dp = surf_press/nlev
  DO k=1, nlev+1
     phalf(k) = dp*(k-1)
  END DO
  DO k=1, nlev
     pfull(k) = .5*(phalf(k) + phalf(k+1))
     bk(k) = pfull(k)/surf_press
  END DO

  ALLOCATE(dat1(is1:ie1,js1:je1,nlev))
  ALLOCATE(dat1h(is1-hi:ie1+hi,js1-hj:je1+hj,nlev))
  dat1h = 0.
  dat1 = 0.
  DO j=js1, je1
     DO i=is1, ie1
        dat1(i,j,1) = SIN(lon1(i))*COS(lat1(j))
     END DO
  END DO
  dat1h(is1:ie1,js1:je1,1) = dat1(:,:,1)
  dat1(:,:,2) = -dat1(:,:,1)
  dat1h(:,:,2) = -dat1h(:,:,1)

  ALLOCATE(dat2(is2:ie2,js2:je2,nlev))
  ALLOCATE(dat2_2d(is2:ie2,js2:je2))
  ALLOCATE(dat2h(is2-hi:ie2+hi,js2-hj:je2+hj,nlev))
  dat2h = 0.
  dat2 = 0.
  DO j=js2, je2
     DO i=is2, ie2
        dat2(i,j,1) = SIN(lon2(i))*COS(lat2(j))
     END DO
  END DO
  dat2h(is2:ie2,js2:je2,1) = dat2(:,:,1)
  dat2(:,:,2) = -dat2(:,:,1)
  dat2h(:,:,2) = -dat2h(:,:,1)
  dat2_2d = dat2(:,:,1)

  id_lonb1 = diag_axis_init('lonb1', RAD_TO_DEG*lonb_global1, 'degrees_E', 'x', long_name='longitude edges', &
                          & Domain2=Domain1)
  id_latb1 = diag_axis_init('latb1', RAD_TO_DEG*latb_global1, 'degrees_N', 'y', long_name='latitude edges',  &
                          & Domain2=Domain1)

  id_lon1  = diag_axis_init('lon1',  RAD_TO_DEG*lon_global1, 'degrees_E','x',long_name='longitude',Domain2=Domain1, &
                          & edges=id_lonb1)
  id_lat1  = diag_axis_init('lat1',  RAD_TO_DEG*lat_global1, 'degrees_N','y',long_name='latitude', Domain2=Domain1, &
                          & edges=id_latb1)

  id_phalf= diag_axis_init('phalf', phalf, 'Pa', 'z', long_name='half pressure level', direction=-1)
  id_pfull= diag_axis_init('pfull', pfull, 'Pa', 'z', long_name='full pressure level', direction=-1, edges=id_phalf)

  id_lon2 = diag_axis_init('lon2',  RAD_TO_DEG*lon_global2,  'degrees_E', 'x', long_name='longitude', Domain2=Domain2)
  id_lat2 = diag_axis_init('lat2',  RAD_TO_DEG*lat_global2,  'degrees_N', 'y', long_name='latitude',  Domain2=Domain2)

  IF ( test_number == 22 ) THEN
     ! Can we get the 'nv' axis ID?
     id_nv = get_axis_num('nv', 'nv')
     IF ( id_nv .GT. 0 ) THEN
        write (out_unit,'(a)') 'test22.1 Passes: id_nv has a positive value'
     ELSE
        write (out_unit,'(a)') 'test22.1 Failed: id_nv does not have a positive value'
     END IF

     ! Can I call diag_axis_init on 'nv' again, and get the same ID back?
     id_nv_init = diag_axis_init( 'nv',(/1.,2./),'none','N','vertex number', set_name='nv')
     IF ( id_nv_init .EQ. id_nv ) THEN
        write (out_unit,'(a)') 'test22.2 Passes: Can call diag_axis_init on "nv" and get same ID'
     ELSE
        write (out_unit,'(a)') 'test22.2 Failed: Cannot call diag_axis_init on "nv" and get same ID'
     END IF
  END IF

  IF ( test_number == 14 ) THEN
     Time = set_date(1990,1,29,0,0,0)
  ELSE
     Time = set_date(1990,1,1,0,0,0)
  END IF

  IF ( test_number == 16 ) THEN
     ! Test 16 tests the filename appendix
     CALL set_filename_appendix('g01')
  END IF
  id_dat1 = register_diag_field('test_diag_manager_mod', 'dat1', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data','K')
  IF ( test_number == 18 ) THEN
     CALL diag_field_add_attribute(id_dat1, 'real_att', 2.3)
     CALL diag_field_add_attribute(id_dat1, 'cell_methods', 'area: mean')
     CALL diag_field_add_attribute(id_dat1, 'cell_methods', 'lon: mean')
  END IF
  IF ( test_number == 18 .OR. test_number == 19 ) THEN
     id_dat2 = register_diag_field('test_diag_manager_mod', 'dat2', (/id_lon1,id_lat1,id_pfull/), Time, &
                                 & 'sample data', 'K')
     CALL diag_field_add_attribute(id_dat2, 'interp_method', 'none')
     CALL diag_field_add_attribute(id_dat2, 'int_att', (/ 1, 2 /) )
  ELSE
     id_dat2 = register_diag_field('test_diag_manager_mod', 'dat2', (/id_lon2,id_lat2,id_pfull/), Time, &
                                 & 'sample data', 'K')
  END IF
  id_sol_con = register_diag_field ('test_diag_manager_mod', 'solar_constant', Time, &
                  'solar constant', 'watts/m2')

  IF ( test_number == 20 ) THEN
     id_dat2_got = get_diag_field_id('test_diag_manager_mod', 'dat2')
     IF ( id_dat2_got == id_dat2 ) THEN
        WRITE (out_unit,'(a)') 'test20.1 Passes, id_dat2.EQ.id_dat2_got'
     ELSE
        WRITE (out_unit,'(a)') 'test20.1 Failed, id_dat2.NE.id_dat2_got'
     END IF

     id_none_got = get_diag_field_id('no_mod', 'no_var')
     IF ( id_none_got == DIAG_FIELD_NOT_FOUND ) THEN
        write (out_unit,'(a)') 'test20.2 Passes, id_none_got.EQ.DIAG_FIELD_NOT_FOUND'
     ELSE
        write (out_unit,'(a)') 'test20.2 Failed, id_none_got.NE.DIAG_FIELD_NOT_FOUND'
     END IF
  END IF

  IF ( dt_step == 0 ) CALL error_mesg ('test_diag_manager',&
       & 'dt_step is not set', FATAL)

  Time_step = set_time(dt_step,0)
  Time_start = Time
  Time_end = Time
  DO m = 1,months
     Time_end = Time_end + set_time(0,days_in_month(Time_end))
  END DO
  Time_end   = Time_end + set_time(0, days)
  Run_length = Time_end - Time_start
  nstep = Run_length / Time_step

  IF ( test_number == 18 ) THEN
     id_dat2h = register_diag_field('test_mod', 'dat2h', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K',&
          & volume=id_dat1, area=id_dat2, realm='myRealm', err_msg=err_msg)
     IF ( err_msg /= '' .OR. id_dat2h <= 0 ) THEN
        CALL error_mesg ('test_diag_manager',&
             & 'Unexpected error registering dat2h '//err_msg, FATAL)
     END IF
     id_dat2h_2 = register_diag_field('test_mod', 'dat2h_2', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K',&
          & err_msg=err_msg)
     CALL diag_field_add_cell_measures(id_dat2h_2, area=id_dat2, volume=id_dat1)
  ELSE IF ( test_number == 19 ) THEN
     id_dat2h = register_diag_field('test_mod', 'dat2h', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K',&
          & volume=id_dat1, area=id_dat1, err_msg=err_msg)
     IF ( err_msg /= '' .OR. id_dat2h <= 0 ) THEN
        CALL error_mesg ('test_diag_manager',&
             & 'Expected error registering dat2h '//err_msg, NOTE)
     END IF
  END IF

  IF(test_number == 16 .OR. test_number == 17 .OR. test_number == 18 .OR. test_number == 21 .OR. test_number == 22)THEN
     is_in = 1
     js_in = 1
     ie_in = nlon
     je_in = nlat

     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     IF ( id_dat2h > 0 ) used = send_data(id_dat2h, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, &
                                        & err_msg=err_msg)
     IF ( id_dat2h_2 > 0 ) used = send_data(id_dat2h_2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, &
                                          & je_in=je_in, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     IF ( id_dat2h > 0 ) used = send_data(id_dat2h, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, &
                                        & err_msg=err_msg)
     IF ( id_dat2h_2 > 0 ) used = send_data(id_dat2h_2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, &
                                          & je_in=je_in, err_msg=err_msg)
  END IF

  !-- The following is used to test openMP
  IF ( test_number == 15 ) THEN
!$      call omp_set_num_threads(numthreads)
      nyc1 = je1 - js1 + 1
      IF (MOD(nyc1, numthreads ) /= 0) THEN
         CALL error_mesg ('test_diag_manager',&
              & 'The number of OpenMP threads must be an integral multiple &
              &of the number of rows in the compute domain', FATAL)
     END IF
     ny_per_thread = nyc1/numthreads

     dat1 = 1
     CALL diag_manager_set_time_end(Time_end)
     DO n = 1, nstep

        Time = Time + Time_step
        !$OMP parallel do default(shared) private(isw, iew, jsw, jew )

        DO jsw = js1, je1, ny_per_thread
           jew = jsw + ny_per_thread -1
           isw = is1
           iew = ie1
           if(id_dat1>0) used = send_data(id_dat1, dat1(isw:iew, jsw:jew,:), Time, &
                                is_in=isw-is1+1, js_in=jsw-js1+1,err_msg=err_msg)
           if(id_sol_con>0) used = send_data(id_sol_con, solar_constant, Time, err_msg=err_msg )
        END DO
        !$OMP END parallel do
        !CALL diag_send_complete(Time_step)
     END DO
 END IF


  IF ( test_number == 14 ) THEN
     id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data','K',err_msg=err_msg)
     IF ( err_msg /= '' ) THEN
        WRITE (out_unit,'(a)') 'test14 successful. err_msg='//TRIM(err_msg)
     ELSE
        WRITE (out_unit,'(a)') 'test14 fails.'
     END IF
  ELSE
     id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data', 'K')
  END IF

  id_bk = register_static_field('test_diag_manager_mod', 'bk', (/id_pfull/), 'half level sigma', 'none')

  IF ( test_number == 13 ) THEN
     IF ( id_dat2_2d > 0 ) used=send_data(id_dat2_2d, dat2(:,:,1), Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') &
              & 'test13: successful if a WARNING message appears that refers to output interval greater than runlength'
     ELSE
        WRITE (out_unit,'(a)') 'test13 fails: err_msg='//TRIM(err_msg)
     END IF
  END IF

  ! Note: test12 involves diag_manager_init, it does not require a call to send_data.
  !       See call to diag_manager_init above.

  IF ( test_number == 11 ) THEN
     is_in = 1+hi
     js_in = 1+hj
     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj

     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, &
                                     & err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test11.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test11.1 fails. err_msg='//TRIM(err_msg)
     END IF

     ! intentional_error: je_in is missing
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test11.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test11.2 successful. err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 10 ) THEN
     !  1 window, no halos, static, 1 dimension, global data.

     IF ( id_bk > 0 ) used = send_data(id_bk, bk, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test10.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test10.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too large.
     IF ( id_bk > 0 ) used = send_data(id_bk, phalf, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE(out_unit,'(a)') 'test10.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test10.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 9 ) THEN
     !  1 window, no halos, static, 1 dimension, global data
     IF ( id_bk > 0 ) used = send_data(id_bk, bk, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test9.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test9.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small
     IF ( id_bk > 0 ) used = send_data(id_bk, bk(1:nlev-1), err_msg=err_msg) ! intentional_error
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test9.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test9.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 8 ) THEN
     !  1 window with halos
     is_in = 1+hi
     js_in = 1+hj

     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test8.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test8.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small in both x and y directions
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     ie_in = ie1-is1+1+hi
     je_in = je1-js1+1+hj
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in, &
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test8.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test8.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 7 ) THEN
     !  1 window with halos
     is_in = 1+hi
     js_in = 1+hj

     ie_in = ie1-is1+1+hi
     je_in = je1-js1+1+hj
     IF ( id_dat1 > 0 ) used=send_data(id_dat1, dat1h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test7.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test7.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too large in both x and y directions
     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj
     IF ( id_dat1 > 0 ) used=send_data(id_dat1, dat2h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test7.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test7.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 6 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within either do loop for test6.1
     test_successful = .TRUE.
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test6.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test6.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test6.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test6.1 fails.'
     END IF

     !  intentional_error: data array too small in y direction
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1)
     END DO
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test6.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test6.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 5 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within either do loop for test5.1
     test_successful = .TRUE.
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test5.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test5.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test5.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test5.1 fails.'
     END IF

     !  intentional_error: data array too small in x direction.
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1)
     END DO
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test5.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test5.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 4 ) THEN
     !  1 window, no halos
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test4.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test4.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small in both x and y directions
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test4.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test4.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 3 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within do loop for test3.1
     test_successful = .TRUE.
     DO i=is1, ie1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test3.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test3.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test3.1 fails.'
     END IF

     !  intentional_error: data array too large in y direction
     DO i=is1, ie1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test3.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test3.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 2 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within do loop for test2.1
     test_successful = .TRUE.
     DO j=js1, je1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test2.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test2.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test2.1 fails.'
     END IF

     !  intentional_error: data array too large in x direction
     DO j=js1, je1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test2.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test2.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 1 ) THEN
     !  1 window, no halos
     ! Here dat2 is too large in both x and y direction so you should get an error.
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test1.1 fails: Intentional error not detected'
     ELSE
        WRITE (out_unit,'(a)') 'test1.1 successful: '//TRIM(err_msg)
     END IF

     ! Here dat1 has the correct shape, so no error
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)

     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test1.2 successful'
     ELSE
        WRITE (out_unit,'(a)') 'test1.2 fails: '//TRIM(err_msg)
     END IF
  END IF
  CALL diag_manager_end(Time)
END SELECT ! End of case handling opened for test 12.

  CALL fms_end

CONTAINS

  SUBROUTINE compute_grid(nlon, nlat, is, ie, js, je, lon_global, lat_global, lonb_global, latb_global, lon, lat, &
                        & lonb, latb)
    INTEGER, INTENT(in) :: nlon, nlat, is, ie, js, je
    REAL, INTENT(out), DIMENSION(:) :: lon_global, lat_global, lonb_global, latb_global, lon, lat, lonb, latb

    REAL :: dlon, dlat
    INTEGER :: i, j

    dlon = 2*PI/nlon
    dlat = PI/nlat

    DO i=1, nlon+1
       lonb_global(i) = dlon*(i-1)
    END DO
    DO j=1,nlat+1
       latb_global(j) = dlat*(j-1) - .5*PI
    END DO
    DO i=1,nlon
       lon_global(i) = .5*(lonb_global(i) + lonb_global(i+1))
    END DO
    DO j=1,nlat
       lat_global(j) = .5*(latb_global(j) + latb_global(j+1))
    END DO
    lon  = lon_global(is:ie)
    lat  = lat_global(js:je)
    lonb = lonb_global(is:ie+1)
    latb = latb_global(js:je+1)
  END SUBROUTINE compute_grid

END PROGRAM test
