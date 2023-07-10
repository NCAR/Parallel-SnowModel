!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! netcdf_module.f90

module netcdf_module
  use netcdf

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
    use netcdf
    implicit none
    character(len = *), intent(in) :: file
    character*(*),intent(in) :: var_name
    integer,intent(in) :: nx,ny,dim_check
    real, intent(out) ::  data1(:,:)
    integer :: ncid, varid, i, j
    character (len = nf90_max_name) :: name
    integer :: xtype, ndims, dimids(max_dims), natts
    integer :: chunksizes(max_dims), deflate_level, endianness
    logical :: contiguous, shuffle, fletcher32


    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access
    ! to the file.
    call check( nf90_open(file, NF90_NOWRITE, ncid) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, var_name, varid) )

    ! Learn about the variable. This uses all optional parameters.
    call check( nf90_inquire_variable(ncid, varid, name, xtype, ndims, &
         dimids, natts, contiguous = contiguous, chunksizes = chunksizes, &
         deflate_level = deflate_level, shuffle = shuffle, &
         fletcher32 = fletcher32, endianness = endianness))
    if (name .ne. var_name .or. xtype .ne. NF90_FLOAT .or. ndims .ne. dim_check) stop 3
    call check( nf90_get_var(ncid, varid, data1) )

    call check( nf90_close(ncid)) ! close file

  end subroutine read_check_float_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_check_double_1D(file,var_name,nx,ny,data1,dim_check)

! Reads 2D float variable
    use netcdf
    implicit none
    character(len = *), intent(in) :: file
    character*(*),intent(in) :: var_name
    integer,intent(in) :: nx,ny,dim_check
    double precision, intent(out) ::  data1(:)
    integer :: ncid, varid, i, j
    character (len = nf90_max_name) :: name
    integer :: xtype, ndims, dimids(max_dims), natts
    integer :: chunksizes(max_dims), deflate_level, endianness
    logical :: contiguous, shuffle, fletcher32


    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access
    ! to the file.
    call check( nf90_open(file, NF90_NOWRITE, ncid) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, var_name, varid) )

    ! Learn about the variable. This uses all optional parameters.
    call check( nf90_inquire_variable(ncid, varid, name, xtype, ndims, &
         dimids, natts, contiguous = contiguous, chunksizes = chunksizes, &
         deflate_level = deflate_level, shuffle = shuffle, &
         fletcher32 = fletcher32, endianness = endianness))
    ! print*,'ndims',ndims
    ! print*,'name',name
    ! print*,'var_name',var_name
    ! print*,'xtype',xtype
    ! print*,'NF90_INT64',NF90_INT64
    if (name .ne. var_name .or. xtype .ne. NF90_DOUBLE .or. ndims .ne. dim_check) stop 3
    call check( nf90_get_var(ncid, varid, data1) )

    call check( nf90_close(ncid)) ! close file

  end subroutine read_check_double_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_check_float_3D(file,var_name,nx,ny,data1,dim_check,netcdf_index,&
                            &    met_time,var_flag,netcdf_nx,netcdf_ny,dt)
! Reads 2D float variable
    use netcdf
    implicit none
    character(len = *), intent(in) :: file
    character(len = *),intent(in) :: var_name
    integer,intent(in) :: nx,ny,dim_check,netcdf_index,met_time,var_flag
    integer,intent(in) :: netcdf_nx,netcdf_ny
    real, intent(in) :: dt
    real, intent(out) ::  data1(:)
    integer, allocatable :: start(:),cnt(:)
    real,allocatable :: data_full(:,:,:)
    integer :: ncid, varid, i, j, k, kk
    character (len = nf90_max_name) :: name
    integer :: xtype, ndims, dimids(max_dims), natts
    integer :: chunksizes(max_dims), deflate_level, endianness
    logical :: contiguous, shuffle, fletcher32
    real :: val


    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access
    ! to the file.

    ! allocate start and count arrays for each dimension
    allocate(start(3),cnt(3))

    ! populate start and count arrays
    start(1:2) = 1 ! start netcdf_nx and netcdf_ny at 1
    cnt(1) = netcdf_nx
    cnt(2) = netcdf_ny
    ! if (dt.lt.10800.0) then ! 1 hour

    start(3) = netcdf_index ! start time at netcdf_index
    cnt(3) = 1
    allocate(data_full(nx,ny,1))

    call check( nf90_open(file, NF90_NOWRITE, ncid) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, var_name, varid) )

    ! Learn about the variable. This uses all optional parameters.
    call check( nf90_inquire_variable(ncid, varid, name, xtype, ndims, &
         dimids, natts, contiguous = contiguous, chunksizes = chunksizes, &
         deflate_level = deflate_level, shuffle = shuffle, &
         fletcher32 = fletcher32, endianness = endianness))
    ! print*,'ndims',ndims
    ! print*,'name',name
    ! print*,'var_name',var_name
    ! print*,'xtype',xtype
    ! print*,'NF90_FLOAT',NF90_FLOAT
    if (name .ne. var_name .or. xtype .ne. NF90_FLOAT .or. ndims .ne. dim_check) stop 3
    call check( nf90_get_var(ncid, varid, data_full,start,cnt) )

! Flattening 3D into 1D array and perform perprocessing

    k = 1
    do j=1,ny
      do i=1,nx
        data1(k) = data_full(i,j,1)
        k = k + 1
      enddo
    enddo

    deallocate(data_full)
    call check( nf90_close(ncid)) ! close file

  end subroutine read_check_float_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check(status)
   integer, intent ( in) :: status

   if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop 2
   end if
  end subroutine check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_dimensions(file,namein,nameout,length)
    use netcdf
    character*50,intent(in) :: namein
    character(len = *), intent(in) :: file
    character*50,intent(out) :: nameout
    integer, intent(out) :: length
    integer :: ncid, varid, dimid

    call check( nf90_open(file, NF90_NOWRITE, ncid) ) ! gets ncid
    call check( nf90_inq_dimid(ncid,namein,dimid)) ! returns dimension id
    call check( nf90_inquire_dimension(ncid,dimid,nameout,length)) ! gets name and dimension
    call check( nf90_close(ncid)) ! close file

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

! Do preparation work.
!
! Look for out-of-range option values.
     if ((ioptn .eq. 0) .or. (abs(ioptn) .ge. 6)) then
        write(*,*)'For calndr(), you specified ioptn = ', ioptn
        write(*,*) &
    &   'Allowable values are 1 to 5 for the Gregorian calendar'
        write(*,*) &
    &   'and -1 to -5 for the Julian calendar.'
        stop
     endif

! Options 1-3 have "iyear" as an input value.
! Internally, we use variable "jyear" that does not have a jump
! from -1 (for 1 BC) to +1 (for 1 AD).
     if (abs(ioptn) .le. 3) then
        if (iyear .gt. 0) then
           jyear = iyear
        elseif (iyear .eq. 0) then
           write(*,*) &
    &      'For calndr(), you specified the nonexistent year 0'
           stop
        else
           jyear = iyear + 1
        endif
!
!        Set "leap" equal to 0 if "jyear" is not a leap year
!        and equal to 1 if it is a leap year.
        leap = 0
        if ((jyear/4)*4 .eq. jyear) then
           leap = 1
        endif
        if ((ioptn .gt. 0).and.((jyear/100)*100 .eq. jyear).and. &
    &       ((jyear/400)*400 .ne. jyear)) then
              leap = 0
        endif
     endif
!
! Options 3-5 involve Julian Day numbers, which need a reference year
! and the Julian Days that began at noon on 1 March of the reference
! year under the Gregorian and Julian calendars.  Any year for which
! "jyear" is divisible by 400 can be used as a reference year.
! We chose 1600 AD as the reference year because it is the closest
! multiple of 400 to the institution of the Gregorian calendar, making
! it relatively easy to compute the Julian Day for 1 March 1600
! given that, on 15 October 1582 under the Gregorian calendar,
! the Julian Day was 2299161.  Similarly, we need to do the same
! calculation for the Julian calendar.  We can compute this Julian
! Day knwoing that on 4 October 1582 under the Julian calendar,
! the Julian Day number was 2299160.  The details of these calculations
! is next.
!    From 15 October until 1 March, the number of days is the remainder
! of October plus the days in November, December, January, and February:
! 17+30+31+31+28 = 137, so 1 March 1583 under the Gregorian calendar
! was Julian Day 2,299,298.  Because of the 10 day jump ahead at the
! switch from the Julian calendar to the Gregorian calendar, 1 March
! 1583 under the Julian calendar was Julian Day 2,299,308.  Making use
! of the rules for the two calendar systems, 1 March 1600 was Julian
! Day 2,299,298 + (1600-1583)*365 + 5 (due to leap years) =
! 2,305,508 under the Gregorian calendar and day 2,305,518 under the
! Julian calendar.
!    We also set the number of days in 400 years and 100 years.
! For reference, 400 years is 146097 days under the Gregorian calendar
! and 146100 days under the Julian calendar.  100 years is 36524 days
! under the Gregorian calendar and 36525 days under the Julian calendar.
     if (abs(ioptn) .ge. 3) then
!
!        Julian calendar values.
        yrref  =    1600
        jdref  = 2305518
!               = Julian Day reference value for the day that begins
!                 at noon on 1 March of the reference year "yrref".
        ndy400 = 400*365 + 100
        ndy100 = 100*365 +  25
!
!        Adjust for Gregorian calendar values.
        if (ioptn .gt. 0) then
           jdref  = jdref  - 10
           ndy400 = ndy400 -  3
           ndy100 = ndy100 -  1
        endif
     endif
!
!----------------------------------------------------------------
! OPTIONS -1 and +1:
! Given a calendar date (iday,month,iyear), compute the day number
! of the year (idayct), where 1 January is day number 1 and 31 December
! is day number 365 or 366, depending on whether it is a leap year.
     if (abs(ioptn) .eq. 1) then
!
!     Compute the day number during the year.
     if (month .le. 2) then
        idayct = iday + (month-1)*31
     else
        idayct = iday + int(30.6001 * (month+1)) - 63 + leap
     endif
!
!----------------------------------------------------------------
! OPTIONS -2 and +2:
! Given the day number of the year (idayct) and the year (iyear),
! compute the day of the month (iday) and the month (month).
     elseif (abs(ioptn) .eq. 2) then
!
     if (idayct .lt. 60+leap) then
        month  = (idayct-1)/31
        iday   = idayct - month*31
        month  = month + 1
     else
        ndays  = idayct - (60+leap)
!               = number of days past 1 March of the current year.
        jmonth = (10*(ndays+31))/306 + 3
!               = month counter, =4 for March, =5 for April, etc.
        iday   = (ndays+123) - int(30.6001*jmonth)
        month  = jmonth - 1
     endif
!
!----------------------------------------------------------------
! OPTIONS -3 and +3:
! Given a calendar date (iday,month,iyear), compute the Julian Day
! number (idayct) that starts at noon.
     elseif (abs(ioptn) .eq. 3) then
!
!     Shift to a system where the year starts on 1 March, so January
!     and February belong to the preceding year.
!     Define jmonth=4 for March, =5 for April, ..., =15 for February.
     if (month .le. 2) then
       jyear  = jyear -  1
       jmonth = month + 13
     else
       jmonth = month +  1
     endif
!
!     Find the closest multiple of 400 years that is .le. jyear.
     yr400 = (jyear/400)*400
!           = multiple of 400 years at or less than jyear.
     if (jyear .lt. yr400) then
        yr400 = yr400 - 400
     endif
!
     n400yr = (yr400 - yrref)/400
!            = number of 400-year periods from yrref to yr400.
     nyrs   = jyear - yr400
!            = number of years from the beginning of yr400
!              to the beginning of jyear.
!
!     Compute the Julian Day number.
     idayct = iday + int(30.6001*jmonth) - 123 + 365*nyrs + nyrs/4 &
    &       + jdref + n400yr*ndy400
!
!     If we are using the Gregorian calendar, we must not count
!     every 100-th year as a leap year.  nyrs is less than 400 years,
!     so we do not need to consider the leap year that would occur if
!     nyrs were divisible by 400, i.e., we do not add nyrs/400.
     if (ioptn .gt. 0) then
        idayct = idayct - nyrs/100
     endif
!
!----------------------------------------------------------------
! OPTIONS -5, -4, +4, and +5:
! Given the Julian Day number (idayct) that starts at noon,
! compute the corresponding calendar date (iday,month,iyear)
! (abs(ioptn)=4) or day number during the year (abs(ioptn)=5).
       else
!
!     Create a new reference date which begins on the nearest
!     400-year cycle less than or equal to the Julian Day for 1 March
!     in the year in which the given Julian Day number (idayct) occurs.
       ndays  = idayct - jdref
       n400yr = ndays / ndy400
!            = integral number of 400-year periods separating
!              idayct and the reference date, jdref.
       jdref  = jdref + n400yr*ndy400
       if (jdref .gt. idayct) then
          n400yr = n400yr - 1
          jdref  = jdref  - ndy400
       endif
!
       ndays  = idayct - jdref
!            = number from the reference date to idayct.
!
       n100yr = min(ndays/ndy100, 3)
!            = number of complete 100-year periods
!              from the reference year to the current year.
!              The min() function is necessary to avoid n100yr=4
!              on 29 February of the last year in the 400-year cycle.
!
       ndays  = ndays - n100yr*ndy100
!            = remainder after removing an integral number of
!              100-year periods.
!
       n4yr   = ndays / 1461
!            = number of complete 4-year periods in the current century.
!              4 years consists of 4*365 + 1 = 1461 days.
!
       ndays  = ndays - n4yr*1461
!            = remainder after removing an integral number
!              of 4-year periods.
!
       n1yr   = min(ndays/365, 3)
!            = number of complete years since the last leap year.
!              The min() function is necessary to avoid n1yr=4
!              when the date is 29 February on a leap year,
!              in which case ndays=1460, and 1460/365 = 4.
!
       ndays  = ndays - 365*n1yr
!            = number of days so far in the current year,
!              where ndays=0 on 1 March.
!
       iyear  = n1yr + 4*n4yr + 100*n100yr + 400*n400yr + yrref
!            = year, as counted in the standard way,
!              but relative to 1 March.
!
! At this point, we need to separate ioptn=abs(4), which seeks a
! calendar date, and ioptn=abs(5), which seeks the day number during
! the year.  First compute the calendar date if desired (abs(ioptn)=4).
       if (abs(ioptn) .eq. 4) then
          jmonth = (10*(ndays+31))/306 + 3
!               = offset month counter.  jmonth=4 for March, =13 for
!                 December, =14 for January, =15 for February.
          iday   = (ndays+123) - int(30.6001*jmonth)
!               = day of the month, starting with 1 on the first day
!                 of the month.
!
!        Now adjust for the fact that the year actually begins
!        on 1 January.
          if (jmonth .le. 13) then
             month = jmonth - 1
          else
             month = jmonth - 13
             iyear = iyear + 1
          endif
!
! This code handles abs(ioptn)=5, finding the day number during the year.
       else
!        ioptn=5 always returns month=1, which we set now.
          month = 1
!
!        We need to determine whether this is a leap year.
          leap = 0
          if ((jyear/4)*4 .eq. jyear) then
             leap = 1
          endif
          if ((ioptn .gt. 0)               .and. &
      &       ((jyear/100)*100 .eq. jyear) .and. &
      &       ((jyear/400)*400 .ne. jyear)      ) then
                leap = 0
          endif
!
!        Now find the day number "iday".
!        ndays is the number of days since the most recent 1 March,
!        so ndays=0 on 1 March.
          if (ndays .le.305) then
             iday  = ndays + 60 + leap
          else
             iday  = ndays - 305
             iyear = iyear + 1
          endif
       endif
!
!     Adjust the year if it is .le. 0, and hence BC (Before Christ).
       if (iyear .le. 0) then
          iyear = iyear - 1
       endif
!
! End the code for the last option, ioptn.
     endif

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

    ihrs_in_day = 24

    if (run.eq.'3hrly') then

!This will create daily date-time stamps with hour = 0.0 to 23.0
!  for each day.
      ioptn = 4
      do julian=julian_start,julian_end
        call calndr (ioptn,iday,imonth,iyear,julian)
        do k=3,ihrs_in_day,3
          kk = (julian - julian_start) * ihrs_in_day/3 + k/3

!Write the date stamp out using the 0 hour at the end of the day
!  format, while taking account for leap years and the change-over
!  from one month and one year to the next.
          if (k.eq.24) then
            last = lastday(imonth)
            if (imonth.eq.2 .and. mod(iyear,4).eq.0 .and.&
   &          (mod(iyear,100).ne.0 .or. mod(iyear,1000).eq.0)) then
              last = last + 1
            endif
            if (iday.eq.last) then
              iday = 1
              imonth = imonth + 1
              if (imonth.gt.12) then
                imonth = 1
                iyear = iyear + 1
              endif
            else
              iday = iday + 1
            endif
            iyr(kk) = iyear
            imo(kk) = imonth
            idy(kk) = iday
            ihr(kk) = 0
          else
            iyr(kk) = iyear
            imo(kk) = imonth
            idy(kk) = iday
            ihr(kk) = k
          endif
!          print *,kk,iyr(kk),imo(kk),idy(kk),ihr(kk)
        enddo
      enddo

    elseif (run.eq.'daily') then

!This will create daily date stamps with hour = 12.0 (in the
!  middle of the day).
      ioptn = 4
      do julian=julian_start,julian_end
        call calndr (ioptn,iday,imonth,iyear,julian)
        kk = julian - julian_start + 1
        iyr(kk) = iyear
        imo(kk) = imonth
        idy(kk) = iday
        ihr(kk) = 12
!        print *,kk,iyr(kk),imo(kk),idy(kk),ihr(kk)
      enddo
    endif

  end subroutine date_timestamp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_nearest_stns_netcdf_1(nx,ny,xmn,ymn,deltax,deltay,&
 &  n_stns_used,k_stn,snowmodel_line_flag,xg_line,yg_line,netcdf_fpath)

#ifdef SERIAL
    use serial_module, only: global_map,global_ny,me
#elif CAF
    use caf_module, only: global_map,global_ny,me
#endif
    implicit none

    integer, parameter :: nvars_in_ctr = 3

    real, allocatable :: xstn(:,:), ystn(:,:)
    double precision, allocatable :: dsq(:)
    double precision xg_line(nx,ny),yg_line(nx,ny)
    real snowmodel_line_flag
    character*80 netcdf_fpath

    double precision xg,yg,xmn,ymn,dist_min,dist_min2
    real deltax,deltay

    integer i,j,k,kk,n_stns_used,nx,ny,i1,i2,i3,i4,new_j,test_length
    integer u,w,k_min,station_ny,station_nx,stn,station,v,num
    integer k_stn_lst_ct
    integer k_stn(nx,ny,9)
    integer, allocatable :: k_stn_lst(:)
    integer cnt(9)
    character*150 fname_ctr
    character*50 var_name_ctr(nvars_in_ctr),nameout
    character*50 dimension_ctr(2)
    integer netcdf_nx,netcdf_ny,nstns

! Variable assignment
      fname_ctr = trim(netcdf_fpath) // 'CTRL/wrf_ctrl.nc' ! file path for netcdf control

      var_name_ctr(1) = 'HGT' ! topo variable name
      var_name_ctr(2) = 'PROJ_XLONG' ! longitude variable name
      var_name_ctr(3) = 'PROJ_XLAT' ! latitude variable name

      dimension_ctr(1) = 'x' ! nx dimension of netcdf
      dimension_ctr(2) = 'y' ! ny dimension of netcdf


! Get wrf dimensions
      call get_dimensions(fname_ctr,dimension_ctr(1),nameout,netcdf_nx)
      call get_dimensions(fname_ctr,dimension_ctr(2),nameout,netcdf_ny)

! Calculate total number of netcdf grid points
      nstns = netcdf_nx*netcdf_ny

! Allocate arrays based on netcdf_nx and netcdf_ny shape
      allocate(dsq(nstns))
      allocate(xstn(netcdf_nx,netcdf_ny),ystn(netcdf_nx,netcdf_ny))
      allocate(k_stn_lst(nx*ny*n_stns_used))

! Pull in netcdf longitude data
      call read_check_float_2D(fname_ctr,var_name_ctr(2),netcdf_nx,netcdf_ny,xstn,2)

! Pull in netcdf latitude data
      call read_check_float_2D(fname_ctr,var_name_ctr(3),netcdf_nx,netcdf_ny,ystn,2)


      k_stn_lst_ct = 1
      do j=1,ny
        new_j = global_map(j)

        if (snowmodel_line_flag.eq.1.0) then
          xg = xg_line(i,j)
          yg = yg_line(i,j)
        else
          xg = xmn + deltax * (real(1) - 1.0)
          yg = ymn + deltay * (real(new_j) - 1.0)
        endif

! Loop through all of the stations, calculating the distance
!   between the current grid point and each of the stations.

        k = 1
        do u=1,netcdf_ny
          do w=1,netcdf_nx
            dsq(k) = (xg - xstn(w,u))**2 + (yg - ystn(w,u))**2
            k = k+1
          enddo
        enddo



! Loop through each of the station distances and find the
!   stations closest to the grid point in question.
        dist_min2 = 1.0e30
        do kk=1,n_stns_used
          dist_min = 1.0e30
          do k=1,nstns
            if (dsq(k).le.dist_min) then
              k_stn(1,j,kk) = k
              dist_min = dsq(k)
              if (dist_min.le.dist_min2) then
                k_min = k
                dist_min2 = dist_min
              endif
            endif
          enddo
          dsq(k_stn(1,j,kk)) = 1.0e30
        enddo

        do i=1,nx ! test to see if we can set this to 2 as opposed to re-doing the same grid cell

! calculate nearest stn from previous grid cell
          station_ny = k_min / netcdf_nx
          station_nx = k_min - (station_ny*netcdf_nx)
          ! added 1 for indexing to work
          station_ny = station_ny + 1


! xcoords of grid nodes at index i,j
! ycoords of grid nodes at index i,j
          if (snowmodel_line_flag.eq.1.0) then
            xg = xg_line(i,j)
            yg = yg_line(i,j)
          else
            xg = xmn + deltax * (real(i) - 1.0)
            yg = ymn + deltay * (real(new_j) - 1.0)
          endif

! Loop through all of the stations, calculating the distance
!   between the current grid point and each of the stations.

          k = 1
          do u=station_ny-1,station_ny+1
            do w=station_nx-1,station_nx+1
              dsq(k) = (xg - xstn(w,u))**2 + (yg - ystn(w,u))**2
              cnt(k) = ((u-1) * netcdf_nx) + w
              k = k+1
            enddo
          enddo



! Loop through each of the station distances and find the
!   stations closest to the grid point in question.
          dist_min2 = 1.0e30
          do kk=1,n_stns_used
            dist_min = 1.0e30
            do v=1,k-1
              if (dsq(v).le.dist_min) then
                k_stn(i,j,kk) = cnt(v)
                dist_min = dsq(v)
                num = v
                if (dist_min.le.dist_min2) then
                  k_min = cnt(v)
                  dist_min2 = dist_min
                endif
              endif
            enddo

! Eliminate the last found minimum from the next search by making
!   its distance a big number.
            dsq(num) = 1.0e30

          enddo
        enddo
      enddo

! Print k_stn for debugging

      ! open (unit=113,file='./outputs_wrf/wrf/kstn.gdat',form='unformatted',&
      ! &     access = 'direct',status = 'replace',recl = 4)
      ! k = 1
      ! do j=1,ny
      !   do i=1,nx
      !     do kk = 1,n_stns_used
      !       write(113,rec = k) k_stn(i,j,kk)
      !       k = k + 1
      !     enddo
      !   enddo
      ! enddo
      ! close(113)

      deallocate(xstn,ystn,dsq)



  end subroutine get_nearest_stns_netcdf_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine model_netcdf_time(iyear,imonth,iday,xhour,julian_t)

#ifdef SERIAL
    use serial_module
#elif CAF
    use caf_module
#endif
    implicit none
    integer, intent(inout) :: iyear, imonth, iday
    real, intent(inout) :: xhour
    integer,intent(out) :: julian_t
    integer :: iyr_archive_start,imo_archive_start,idy_archive_start,ioffset_local
    integer :: julian_start,julian_end,ioffset_daily,ioffset_3hrly,ioffset_hrly,ioptn


    iyr_archive_start = 1901 ! Corresponds variable of netcdf that defines hours since
    imo_archive_start = 1
    idy_archive_start = 1

    ! local time variable
  ! If you are running in local time (not UTC; i.e., UTC_flag = 0.0),
  !  account for that here.  For example, for Alaska time, shift
  !  this back by 9 hours (or 3, 3-hour records), for Colorado
  !  time shift this back by 6 hours (or 2, 3-hour records).  This
  !  is only used if you are running with 3-hour time steps, not
  !  daily time steps.  This offset must be set to zero if you are
  !  running with UTC_flag = 1.0.

  !Alaska:
  !    ioffset_local = -3

  !Colorado:
  !     ioffset_local = -2

  !If UTC_flag = 1.0:
     ioffset_local = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GET CALENDAR INFORMATION !!
     ioptn = 3
     ! Calculate offset from netcdf timing start
     call calndr (ioptn,idy_archive_start,imo_archive_start,&
    &  iyr_archive_start,julian_start)


     ! get julian end timing
     call calndr (ioptn,iday,imonth,iyear,julian_end)

     julian_end = julian_end - 1

     julian_t = ((julian_end - julian_start + 1) * 24) + INT(xhour)

     ! find the offset
     ! ioffset_daily = julian_end - julian_start + 1
     ! ioffset_3hrly = 8 * ioffset_daily + ioffset_local
     ! ioffset_hrly = ioffset_3hrly * 3

     ! Calculate the number of time steps in simulations
     ! ioptn = 3
     ! call calndr (ioptn,idy_start,imo_start,iyr_start,julian_start)
     !
     ! call calndr (ioptn,idy_end,imo_end,iyr_end,julian_end)
     !
     ! ! calculate max iterations
     ! maxiter_daily = julian_end - julian_start + 1
     ! maxiter_3hrly = 8 * maxiter_daily
     ! maxiter_hrly = maxiter_3hrly * 3
     !
     ! ! printing variables
     ! print*,'iyr_start', iyr_start
     ! print*,'imo_start', imo_start
     ! print*,'idy_start', idy_start
     ! print*,''
     !
     ! print*,'iyr_end', iyr_end
     ! print*,'imo_end', imo_end
     ! print*,'idy_end', idy_end
     ! print*,''
     !
     ! print*,'maxiter_3hrly', maxiter_3hrly
     ! print*,'ioffset_3hrly', ioffset_3hrly
     ! print*,''
     !
     ! print*,'ioffset_hrly', ioffset_hrly
     ! print*,'maxiter_hrly',maxiter_hrly
     ! print*,''
     !
     ! print*,'maxiter_daily', maxiter_daily
     ! print*,'ioffset_daily', ioffset_daily

  end subroutine model_netcdf_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_netcdf_data(nstns_orig,Tair_orig,rh_orig,xstn_orig,&
     &  ystn_orig,elev_orig,iyear,imonth,iday,xhour,undef,&
     &  windspd_orig,winddir_orig,prec_orig,isingle_stn_flag,&
     &  igrads_metfile,iter,netcdf_time_index,dt,netcdf_fpath,&
     &  rain_orig,snow_orig,nstns_max)

#ifdef SERIAL
      use serial_module, only: me
#elif CAF
      use caf_module, only: me
#endif
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
      double precision, allocatable :: time(:)
      real, allocatable :: t2(:,:,:)
      real, allocatable ::  xstn(:,:),ystn(:,:),elev(:,:)
      real elev_orig(nstns_max),xhour,prec_orig(nstns_max)
      real undef               ! undefined value
      integer met_time,netcdf_time_index,netcdf_idx
      character(len=120) :: netcdf_dir
      character(len=150) :: file, fname_ctr
      character*50 met_vars(total_vars),nameout
      character*50 ctr_vars(num_ctrl_vars),dimension_ctr(2)
      character*80 netcdf_fpath
      integer nstns_max

! User inputs
      fname_ctr = trim(netcdf_fpath) // 'CTRL/wrf_ctrl.nc' ! file path for netcdf control

      ctr_vars(1) = 'HGT' ! topo variable name
      ctr_vars(2) = 'PROJ_XLONG' ! longitude variable name
      ctr_vars(3) = 'PROJ_XLAT' ! latitude variable name

      dimension_ctr(1) = 'x' ! nx dimension of wrf
      dimension_ctr(2) = 'y' ! ny dimension of wrf

      met_vars(1) = time_char ! time variable
      met_vars(2) = precip_char ! prec variable name
      met_vars(3) = tair_char ! tair variable name
      met_vars(4) = wdir_char ! wdir variable name
      met_vars(5) = wspd_char ! wspd variable name
      met_vars(6) = relh_char ! relh variable name
      met_vars(7) = rain_char ! rain_acc variable
      met_vars(8) = snow_char ! snow_acc variable


! Identify directory where files are stored.
      call get_netcdf_directory(imonth,iyear,netcdf_fpath,dt,netcdf_dir)

! Get control variable dimensions
      call get_dimensions(fname_ctr,dimension_ctr(1),nameout,netcdf_nx)
      call get_dimensions(fname_ctr,dimension_ctr(2),nameout,netcdf_ny)
      ! print*,'2'

      nstns = netcdf_nx * netcdf_ny

! Assign nstns_orig variable to total stations
      nstns_orig = nstns

      allocate(xstn(netcdf_nx,netcdf_ny),ystn(netcdf_nx,netcdf_ny))
      allocate(elev(netcdf_nx,netcdf_ny))

! Pull station elevation height
      call read_check_float_2D(fname_ctr,ctr_vars(1),netcdf_nx,netcdf_ny,elev,3)

! Pull in wrf longitude data
      call read_check_float_2D(fname_ctr,ctr_vars(2),netcdf_nx,netcdf_ny,xstn,2)

! Pull in wrf latitude data
      call read_check_float_2D(fname_ctr,ctr_vars(3),netcdf_nx,netcdf_ny,ystn,2)
      ! print*,'3'

! Flatten xstn, ystn, and elev
      k = 1
      do j=1,netcdf_ny
        do i=1,netcdf_nx
          xstn_orig(k) = xstn(i,j)
          ystn_orig(k) = ystn(i,j)
          elev_orig(k) = elev(i,j)
          k = k + 1
        enddo
      enddo


! Get Tair 1D array
      call get_1D_met_var(netcdf_dir,met_vars(3),met_vars(1),netcdf_nx,netcdf_ny,netcdf_time_index,&
                            &   Tair_orig,0,dt,imonth)
! Get Precip 1D array
      call get_1D_met_var(netcdf_dir,met_vars(2),met_vars(1),netcdf_nx,netcdf_ny,netcdf_time_index,&
                            &   prec_orig,1,dt,imonth)
! Get wdir 1D array
      call get_1D_met_var(netcdf_dir,met_vars(4),met_vars(1),netcdf_nx,netcdf_ny,netcdf_time_index,&
                            &   winddir_orig,0,dt,imonth)
! Get wspd 1D array
      call get_1D_met_var(netcdf_dir,met_vars(5),met_vars(1),netcdf_nx,netcdf_ny,netcdf_time_index,&
                            &   windspd_orig,0,dt,imonth)
! Get relh 1D array
      call get_1D_met_var(netcdf_dir,met_vars(6),met_vars(1),netcdf_nx,netcdf_ny,netcdf_time_index,&
                            &   rh_orig,0,dt,imonth)
! Get rain 1D array
      call get_1D_met_var(netcdf_dir,met_vars(7),met_vars(1),netcdf_nx,netcdf_ny,netcdf_time_index,&
                            &   rain_orig,0,dt,imonth)
! Get snow 1D array
      call get_1D_met_var(netcdf_dir,met_vars(8),met_vars(1),netcdf_nx,netcdf_ny,netcdf_time_index,&
                            &   snow_orig,0,dt,imonth)



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

! User inputs
      fname_ctr = trim(netcdf_fpath) // 'CTRL/wrf_ctrl.nc' ! file path for netcdf control

      dimension_ctr(1) = 'x' ! nx dimension of wrf
      dimension_ctr(2) = 'y' ! ny dimension of wrf

! Get control variable dimensions
      call get_dimensions(fname_ctr,dimension_ctr(1),nameout,netcdf_nx)
      call get_dimensions(fname_ctr,dimension_ctr(2),nameout,netcdf_ny)
      ! print*,'2'

      nstns_max = (netcdf_nx * netcdf_ny)+1

    end subroutine netcdf_resolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_netcdf_directory(imonth,iyear,netcdf_fpath,dt,netcdf_dir)

#ifdef SERIAL
      use serial_module, only: me
#elif CAF
      use caf_module, only: me
#endif
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

! convert year and month to strings
      fmt = '(I4)'
      write(string,fmt) iyear
      iyear_string = string
      write(string,fmt) imonth
      imonth_string = string

! string adjust
      if (floor(real(imonth) / 10.0).eq.0) then
        imo_string = '0' // trim(adjustl(imonth_string))
        iyr_string = trim(adjustl(iyear_string))
      else
        imo_string = trim(adjustl(imonth_string))
        iyr_string = trim(adjustl(iyear_string))
      endif

! concatenate strings to get directory
      if (dt.lt.10700.0) then
        t_string = '1_hrly/year/'
      elseif (dt.gt.10801.0) then
        t_string = 'daily/year/'
      else
        t_string = '3_hrly/year/'
      endif
      netcdf_dir = trim(netcdf_fpath) // trim(t_string) // iyr_string // '/' // &
               &   'month' // '/' // imo_string // '/'

    end subroutine get_netcdf_directory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_1D_met_var(netcdf_dir,netcdf_var,time_var,nx,ny,netcdf_time_index,&
                            &   data_1D,var_flag,dt,imonth)

#ifdef SERIAL
      use serial_module, only: me
#elif CAF
      use caf_module, only: me
#endif
      implicit none
      integer, intent(in) :: nx,ny,netcdf_time_index,var_flag,imonth !imonth for debugging DELETE
      real, intent(in) :: dt
      character(len=120),intent(in) :: netcdf_dir
      character(len=50), intent(in) :: netcdf_var,time_var
      real, intent(out) :: data_1D(:)

      character(len=150) :: file
      character(len=50) :: nameout
      integer :: met_time,netcdf_idx,i
      double precision, allocatable :: time(:)

! file path string concatenation
      file = trim(netcdf_dir) //  trim(netcdf_var) // '.nc'
! Get time dimension
      call get_dimensions(file,time_var,nameout,met_time)

! allocate time array
      allocate(time(met_time))

! pull time data from netcdf
      call read_check_double_1D(file,time_var,nx,ny,time,1)

! identify index where snowmodel time equals netcdf time
      do i=1,met_time
        if (int(time(i)).eq. netcdf_time_index) then
          netcdf_idx = i
        endif
      enddo




      call read_check_float_3D(file,trim(netcdf_var),nx,ny,data_1D,3,&
                       &       netcdf_idx,met_time,var_flag,nx,ny,dt)




    end subroutine get_1D_met_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module netcdf_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
