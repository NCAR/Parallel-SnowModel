!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nonetcdf_module.f90

module nonetcdf_module

  integer,parameter::max_dims=10,total_vars = 8,num_ctrl_vars = 3
  integer,parameter :: k15 = selected_int_kind(15)
  character(len=15), parameter :: precip_char = 'PREC_ACC_NC',tair_char = 'T2'
  character(len=15), parameter :: wdir_char = 'WDIR_EREL',wspd_char = 'WSPD'
  character(len=15), parameter :: time_char = 'Time',relh_char = 'RELH'
  character(len=15), parameter :: snow_char = 'SNOW_ACC',rain_char = 'RAIN_ACC'

  interface io_read
    module procedure read_check_double_1D,read_check_float_2D,read_check_float_3D,check
  end interface

  interface dimensions
    module procedure get_dimensions
  end interface

  interface calendar
    module procedure calndr, date_timestamp_arrays
  end interface

  interface stns
    module procedure get_nearest_stns_netcdf_1,model_netcdf_time,get_netcdf_data,netcdf_resolution
  end interface
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_check_float_2D(file,var_name,nx,ny,data1,dim_check)
! Reads 2D float variable
    implicit none
    character(len = *), intent(in) :: file
    character*(*),intent(in) :: var_name
    integer,intent(in) :: nx,ny,dim_check
    real, intent(out) ::  data1(:,:)
    integer :: ncid, varid, i, j
    integer :: xtype, ndims, dimids(max_dims), natts
    integer :: chunksizes(max_dims), deflate_level, endianness
    logical :: contiguous, shuffle, fletcher32

! Function called but nothing performed.

  end subroutine read_check_float_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_check_double_1D(file,var_name,nx,ny,data1,dim_check)

! Reads 2D float variable
    implicit none
    character(len = *), intent(in) :: file
    character*(*),intent(in) :: var_name
    integer,intent(in) :: nx,ny,dim_check
    double precision, intent(out) ::  data1(:)
    integer :: ncid, varid, i, j
    integer :: xtype, ndims, dimids(max_dims), natts
    integer :: chunksizes(max_dims), deflate_level, endianness
    logical :: contiguous, shuffle, fletcher32

! Function called but nothing performed.

  end subroutine read_check_double_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_check_float_3D(file,var_name,nx,ny,data1,dim_check,netcdf_index,&
                            &    met_time,var_flag,netcdf_nx,netcdf_ny,dt)
! Reads 2D float variable
    implicit none
    character(len = *), intent(in) :: file
    character(len = *),intent(in) :: var_name
    integer,intent(in) :: nx,ny,dim_check,netcdf_index,met_time,var_flag
    integer,intent(in) :: netcdf_nx,netcdf_ny
    real, intent(in) :: dt
    real, intent(out) ::  data1(:)
    integer :: ncid, varid, i, j, k, kk
    integer :: xtype, ndims, dimids(max_dims), natts
    integer :: chunksizes(max_dims), deflate_level, endianness
    logical :: contiguous, shuffle, fletcher32
    real :: val

! Function called but nothing performed.

  end subroutine read_check_float_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check(status)
   implicit none
   integer, intent ( in) :: status
   
! Function called but nothing performed.

  end subroutine check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_dimensions(file,namein,nameout,length)
    implicit none
    character*50,intent(in) :: namein
    character(len = *), intent(in) :: file
    character*50,intent(out) :: nameout
    integer, intent(out) :: length
    integer :: ncid, varid, dimid
    
! Function called but nothing performed.

  end subroutine get_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calndr(ioptn, iday, month, iyear, idayct)
! Subroutine calndr() performs calendar calculations using either
! the standard Gregorian calendar or the old Julian calendar.
! This subroutine extends the definitions of these calendar systems
! to any arbitrary year.  The algorithms in this subroutine
! will work with any date in the past or future,
! but overflows will occur if the numbers are sufficiently large.
! For a computer using a 32-bit integer, this routine can handle
! any date between roughly 5.8 million BC and 5.8 million AD
! without experiencing overflow during calculations.
   implicit none
   integer, intent (in) :: ioptn
   integer, intent(inout) :: iday,month,iyear,idayct
   integer :: jdref,jmonth,jyear,leap,n1yr,n4yr,n100yr,n400yr
   integer :: ndays,ndy400,ndy100,nyrs,yr400,yrref
   
! Function called but nothing performed.

  end subroutine calndr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine date_timestamp_arrays(run,julian_start,julian_end,max_maxiter,iyr,imo,&
                            &      idy,ihr)
    implicit none
    character*50,intent(in) :: run
    integer, intent(in) :: julian_start,julian_end,max_maxiter
    integer :: julian,ihrs_in_day,ioptn,k,kk,iday,imonth,iyear,last
    integer, intent(inout) :: iyr(max_maxiter),imo(max_maxiter),idy(max_maxiter)
    integer, intent(inout) :: ihr(max_maxiter)
    integer :: lastday(12)
    data lastday/31,28,31,30,31,30,31,31,30,31,30,31/
    
! Function called but nothing performed.

  end subroutine date_timestamp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_nearest_stns_netcdf_1(nx,ny,xmn,ymn,deltax,deltay,&
 &  n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line,netcdf_fpath)

    implicit none

    integer, parameter :: nvars_in_ctr = 3

    double precision xg_line(nx,ny),yg_line(nx,ny)
    real snowmodel_line_flag
    character*80 netcdf_fpath

    double precision xg,yg,xmn,ymn,dist_min,dist_min2
    real deltax,deltay

    integer i,j,k,kk,n_stns_used,nx,ny,i1,i2,i3,i4,new_j,test_length
    integer u,w,k_min,station_ny,station_nx,stn,station,v,num
    integer k_stn_lst_ct
    integer k_stn(nx,ny,9)
    integer cnt(9)
    character*150 fname_ctr
    character*50 var_name_ctr(nvars_in_ctr),nameout
    character*50 dimension_ctr(2)
    integer netcdf_nx,netcdf_ny,nstns

! Function called but nothing performed.
    print*,'YEAAAA'

  end subroutine get_nearest_stns_netcdf_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine model_netcdf_time(iyear,imonth,iday,xhour,julian_t)

    implicit none
    integer, intent(inout) :: iyear, imonth, iday
    real, intent(inout) :: xhour
    integer,intent(out) :: julian_t
    integer :: iyr_archive_start,imo_archive_start,idy_archive_start,ioffset_local
    integer :: julian_start,julian_end,ioffset_daily,ioffset_3hrly,ioffset_hrly,ioptn
    
! Function called but nothing performed.

  end subroutine model_netcdf_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_netcdf_data(nstns_orig,Tair_orig,rh_orig,xstn_orig,&
     &  ystn_orig,elev_orig,iyear,imonth,iday,xhour,undef,&
     &  windspd_orig,winddir_orig,prec_orig,isingle_stn_flag,&
     &  igrads_metfile,iter,netcdf_time_index,dt,netcdf_fpath,&
     &  rain_orig,snow_orig,nstns_max)

      implicit none

      integer iyr,imo,idy      ! year, month, and day of data
      real xhr,dt                 ! decimal hour
      integer idstn            ! station id number

      integer k,nstns_orig,isingle_stn_flag,igrads_metfile,iter
      integer iyear,imonth,iday,netcdf_nx,netcdf_ny,i,j,nstns

      real Tair_orig(nstns_max),rh_orig(nstns_max)
      real winddir_orig(nstns_max),windspd_orig(nstns_max)
      real rain_orig(nstns_max),snow_orig(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)
      real elev_orig(nstns_max),xhour,prec_orig(nstns_max)
      real undef               ! undefined value
      integer met_time,netcdf_time_index,netcdf_idx
      character(len=120) :: netcdf_dir
      character(len=150) :: file, fname_ctr
      character*50 met_vars(total_vars),nameout
      character*50 ctr_vars(num_ctrl_vars),dimension_ctr(2)
      character*80 netcdf_fpath
      integer nstns_max
      
! Function called but nothing performed.

    end subroutine get_netcdf_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine netcdf_resolution(netcdf_fpath,nstns_max)

      implicit none
      
      integer,intent(inout) :: nstns_max
      character*80,intent(in) :: netcdf_fpath
      
      integer netcdf_nx,netcdf_ny
      character(len=150) :: fname_ctr
      character*50 nameout
      character*50 dimension_ctr(2)
      
! Function called but nothing performed.

    end subroutine netcdf_resolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_netcdf_directory(imonth,iyear,netcdf_fpath,dt,netcdf_dir)

      implicit none
      integer, intent(in) :: iyear,imonth
      real, intent(in) :: dt
      character*80, intent(in) :: netcdf_fpath
      character(len=120), intent(out) :: netcdf_dir
              ! undefined value
      character*4 string,iyear_string,imonth_string,iyr_string
      character*2 imo_string
      character*9 fmt
      character*50 t_string
      
! Function called but nothing performed.

    end subroutine get_netcdf_directory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_1D_met_var(netcdf_dir,netcdf_var,time_var,nx,ny,netcdf_time_index,&
                            &   data_1D,var_flag,dt,imonth)

      implicit none
      integer, intent(in) :: nx,ny,netcdf_time_index,var_flag,imonth !imonth for debugging DELETE
      real, intent(in) :: dt
      character(len=120),intent(in) :: netcdf_dir
      character(len=50), intent(in) :: netcdf_var,time_var
      real, intent(out) :: data_1D(:)

      character(len=150) :: file
      character(len=50) :: nameout
      integer :: met_time,netcdf_idx,i
      
! Function called but nothing performed.

    end subroutine get_1D_met_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module nonetcdf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
