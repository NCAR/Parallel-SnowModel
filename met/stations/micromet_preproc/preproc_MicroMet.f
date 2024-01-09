c preproc_MicroMet.f
c
c This is the MicroMet preprocessor.  It assumes a meteorological
c   station input file with the following 9 data columns:
c
c     year,month,day,hour,air_temp,rh,wind_speed,wind_dir,precip
c
c   Thus, it only processes the hourly variables,
c     tair, rh, spd, dir, prec,
c
c     with the units: Tair (deg C), rh (%), wind speed (m/s),
c     wind direction (0-360 True N), and precipitation (mm/hour).
c
c All inputs, outputs, and thresholds are defined through the
c   'preproc_MicroMet.par' file.
c
c The script 'compile.preproc_MicroMet' compiles the program,
c   linking preproc_MicroMet.f with arima.progs.f.  It creates an
c   executable called 'preproc_MicroMet' that is run by typing
c   preproc_MicroMet.
c
c The messages stating out-of-bounds limits that have been exceeded
c   and missing-data sections that were not filled in, are normally
c   printed to the screen.  If you want to save them to a file,
c   issue the command: preproc_MicroMet > preproc_summary.txt.
c
c The output format is such that the output file can be readily
c   used by the MicroMet program.
c
c The MicroMet preprocessor does the following:
c
c   1) Fills any missing dates with the missing-data value.
c   2) Performs a series of data quality assurance/quality control
c        QA/QC tests.
c   3) Fills any missing data with reasonable values.
c
c First, run the CALENDAR program, with an hourly addition, to fill
c   in any missing dates.  Any missing data values corresponding to
c   these missing dates are filled in with the undef value.  Note
c   that this must be done before the QA/QC tests because it would
c   make no sense to be comparing two adjacent values that differ
c   in time over some missing data span.
c
c Second, perform a series of QA/QC tests on the data based on:
c   Meek, D. W., and J. L. Hatfield, 1994: Data quality checking
c   for single station meteorological variables. Agric. For.
c   Meteorol., 69, 85-109.
c
c   For the QA/QC tests, consider three conditions/cases:
c
c   CASE 1.
c   Check for values outside acceptable ranges:
c     High/low ramge limits [LIM].  
c
c   CASE 2.
c   Check for consecutive values that exceed acceptable increments:
c     Rate-of-change limits [ROC].
c
c   CASE 3.
c   Check for constant consecutive values (for example, a constant
c     wind direction or consecutive zero wind speeds that indicate an
c     iced instrument:
c     Continuous no-observed-change with time limits (4 hours) [NOC].
c
c   The formulation assumes the inputs are hourly data.  The input
c     hours in a day can run from 1, 2, ... , 22, 23, 0; or they can
c     run from 1, 2, ... , 22, 23, 24.  The model uses the starting
c     and ending hour inputs to define iflag_1_24:
c       iflag_1_24 = 1 for 1-24
c       iflag_1_24 = 0 for 1-0
c     The reason I do it this way is that the met forcing is the
c     average of the previous hour, so hour 0 is really the last
c     hour of the day, not the first hour of the day.
c
c   Also note that partial days are not allowed with the input data.
c
c   The output values will always be in the 1-0 hour format.  Also
c     note that after the processing the final data arrays will be
c     filled with full days of data (no < 24-hour days are allowed at
c     the start and end of the data file), e.g., it will look like:
c       day 1  2002  9 28 22
c       day 1  2002  9 28 23
c       day 1  2002  9 29  0
c       day 2  2002  9 29  1
c       day 2  2002  9 29  2
c       ....................
c       day 2  2002  9 29 23
c       day 2  2002  9 30  0
c     Another way of saying all of this is that, on input
c       ihr_start == 1 and ihr_end == 0
c     or
c       ihr_start == 1 and ihr_end == 24 
c     And on output
c       ihr_start == 1 and ihr_end == 0
c
c Third, fill in any missing data with reasonable values.  This
c   model deals with tair, rh, spd, dir, prec by assuming there is a
c   diurnal cycle for these variables.  If it turns out that
c   there is no diurnal cycle, this is okay too.  The data is filled
c   differently for the following three cases:
c
c   1) for missing data = 1 hour, average the values from each side
c        of that hour.
c   2) for 2 hours <= missing data <= 24 hours, average the values
c        from 24 hours before and after each of the hours in that
c        period.
c   3) for missing data > 24 hours, use the arima model.  Here look
c        the missing data number of hours before and after the
c        missing period to fill the missing period.
c
c   For the case of missing data spanning more than 24 hours, this
c   formulation closely follows the ideas presented in Walton
c   (1996).  The time-series prediction is made using an
c   AutoRegressive Intergrated Moving Average (ARIMA) model (Box
c   and Jenkins 1976).
c
c   Walton, T. L., Jr., 1996: Fill-in of missing data in univariate
c   coastal data.  J. Applied Statistics, 23, 31-39.
c
c   Box, G. E. P., and G. M. Jenkins, 1976: Time Series Analysis,
c   Forecasting, and Control.  Holden-Day, San Francisco, California.
c
c   There are 4 cases where data will not be filled in by the ARIMA
c   programs, and will be left as missing:
c
c   Case A: The maximum allowed number of undefined days was
c     exceeded.
c   Case B: The ARIMA is required to look BEFORE available data; for
c     example, 5 days of missing data starting on day 3 requires the
c     model to use data before the start of day 1.
c   Case C: The ARIMA is required to look AFTER available data; for
c     example, 5 days of missing data starting two days before the
c     end of the data file equires the model to use data after the
c     end of the last day.
c   Case D: The ARIMA found missing values while searching for good
c     data to fill another missing segment.
c
c   Unless you override it with the "fill_full_flag", the program
c   will not fill missing regions over 24 hours that do not have
c   good data at least an equal size before and after the missing
c   region; it will send a message to the screen if these occur,
c   and leave the missing segments as missing.  A similar thing
c   occurs if there are two missing sections that are too close to
c   each other to get sufficient data to do the data filling.  If
c   fill_full_flag = 1.0 the program will look at that available
c   valid data before and after these missing sections, calculate
c   the average, and use that to fill in the missing hours.  Note
c   that for these special cases, this procedure will leave "flat
c   spots" in your filled data sections.
 
      implicit none

      integer maxiter2
      parameter (maxiter2=100000)

      integer ihrs_in_day,max_fill_days,max_fill_hrs,n_out_vars,
     &  iyr_start,imo_start,idy_start,iyr_end,imo_end,idy_end,ioptn,
     &  maxiter,julian_end,julian_start,julian,k,kk,iday,nrecs,
     &  imonth,iyear,iter,iflag_1_24,ihr_start,ihr_end,last,iheader,
     &  id_stn,tair_check_flag,rh_check_flag,wspd_check_flag,
     &  wdir_check_flag,prec_check_flag,tair_fill_flag,rh_fill_flag,
     &  wspd_fill_flag,wdir_fill_flag,prec_fill_flag,print_flag,
     &  max_arima_fill_days,fill_full_flag

      real undef,Tmin_hi,Tmin_lo,Tmax_hi,Tmax_lo,Tair_hrly_maxinc,
     &  RHmin,RHmax,RH_hrly_maxinc,WSPDmin,WSPDmax,WSPD_hrly_maxinc,
     &  WDIRmin,WDIRmax,WDIR_hrly_maxinc,PRECmin,PRECmax,
     &  PREC_hrly_maxinc,wdir_flag,x_loc,y_loc,elev,so_hem

      real var1(maxiter2)
      real var2(maxiter2)

      real predicted1(maxiter2)
      real predicted2(maxiter2)
      real y1(maxiter2) 
      real y2(maxiter2) 

      integer iyr(maxiter2)
      integer imo(maxiter2)
      integer idy(maxiter2)
      integer ihr(maxiter2)

      integer iiyr(maxiter2)
      integer iimo(maxiter2)
      integer iidy(maxiter2)
      integer iihr(maxiter2)

      real tair(maxiter2)
      real rh(maxiter2)
      real wspd(maxiter2)
      real wdir(maxiter2)
      real prec(maxiter2)

      real tair_orig(maxiter2)
      real rh_orig(maxiter2)
      real wspd_orig(maxiter2)
      real wdir_orig(maxiter2)
      real prec_orig(maxiter2)

      real xtair(maxiter2)
      real xrh(maxiter2)
      real xwspd(maxiter2)
      real xwdir(maxiter2)
      real xprec(maxiter2)

      real Tmin_daily(maxiter2)
      real Tmax_daily(maxiter2)

      character*80 infname
      character*80 outfname_txt
      character*80 outfname_grads

      integer lastday(12)
      data lastday/31,28,31,30,31,30,31,31,30,31,30,31/

c Read the input parameter information.
      call read_param(infname,outfname_txt,outfname_grads,
     &  iyr_start,imo_start,idy_start,ihr_start,iyr_end,imo_end,
     &  idy_end,ihr_end,undef,Tmin_hi,Tmin_lo,Tmax_hi,Tmax_lo,
     &  Tair_hrly_maxinc,RHmin,RHmax,RH_hrly_maxinc,WSPDmin,
     &  WSPDmax,WSPD_hrly_maxinc,WDIRmin,WDIRmax,WDIR_hrly_maxinc,
     &  PRECmin,PRECmax,PREC_hrly_maxinc,max_fill_days,iheader,
     &  id_stn,x_loc,y_loc,elev,tair_check_flag,rh_check_flag,
     &  wspd_check_flag,wdir_check_flag,prec_check_flag,
     &  tair_fill_flag,rh_fill_flag,wspd_fill_flag,wdir_fill_flag,
     &  prec_fill_flag,print_flag,max_arima_fill_days,so_hem,
     &  fill_full_flag)

c The number of hours in a day.
      ihrs_in_day = 24

c The maximum allowed length of continuous missing data that will
c   be filled.
      max_fill_hrs = max_fill_days * ihrs_in_day

c Open the met input file.
      open (55,file=infname,form='formatted')

c Open the output files.
      if (print_flag.eq.1 .or. print_flag.eq.3) then
        open (75,file=outfname_txt,form='formatted')
      endif

      if (print_flag.eq.2 .or. print_flag.eq.3) then
        n_out_vars = 16
        open (65,file=outfname_grads,
     &    form='unformatted',access='direct',recl=4*n_out_vars)
      endif

c Find the Julian day at the start of the data stream.
      ioptn = 3
      call calndr (ioptn,idy_start,imo_start,iyr_start,julian_start)

c Find the Julian day at the end of the data stream.
      call calndr (ioptn,idy_end,imo_end,iyr_end,julian_end)
      if (ihr_end.eq.0) then
        julian_end = julian_end - 1
      elseif (ihr_end.eq.24) then
c       julian_end = julian_end
      endif

c Calculate the total number of time slices in the data file.
      maxiter = (julian_end - julian_start + 1) * ihrs_in_day

      print *
      print *,'NUMBER OF TIME SLICES (HOURS, DAYS) = ',
     &  maxiter,real(maxiter)/real(ihrs_in_day)

c Fill the entire time domain with dates and stn info, and set all
c   variables to the undefined values.
      ioptn = 4
      do julian=julian_start,julian_end
        call calndr (ioptn,iday,imonth,iyear,julian)
        do k=1,ihrs_in_day
          kk = (julian - julian_start) * ihrs_in_day + k

c Write the date stamp out using the 0 hour at the end of the day
c   format, while taking account for leap years and the change-over
c   from one month and one year to the next.
          if (k.eq.24) then
            last = lastday(imonth)
            if (imonth.eq.2 .and. mod(iyear,4).eq.0 .and.
     &        (mod(iyear,100).ne.0 .or. mod(iyear,1000).eq.0)) then
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

            iiyr(kk) = iyear
            iimo(kk) = imonth
            iidy(kk) = iday
            iihr(kk) = 0
          else
            iiyr(kk) = iyear
            iimo(kk) = imonth
            iidy(kk) = iday
            iihr(kk) = k
          endif

          tair(kk) = undef
          rh(kk) = undef
          wspd(kk) = undef
          wdir(kk) = undef
          prec(kk) = undef
        enddo
      enddo

c Read through the met data and assign the obs to the correct record.
c   Read past the header.
      do iter=1,iheader
        read (55,*)
      end do

      ioptn = 3
      nrecs = 0
      do iter=1,maxiter2

        read (55,*,end=99) iyr(iter),imo(iter),idy(iter),ihr(iter),
     &    xtair(iter),xrh(iter),xwspd(iter),xwdir(iter),xprec(iter)

        call calndr (ioptn,idy(iter),imo(iter),iyr(iter),julian)
        kk = (julian - julian_start) * ihrs_in_day + ihr(iter)

        tair(kk) = xtair(iter)
        rh(kk) = xrh(iter)
        wspd(kk) = xwspd(iter)
        wdir(kk) = xwdir(iter)
        prec(kk) = xprec(iter)

        nrecs = nrecs + 1

      enddo

   99 continue

c Check to make sure the starting and ending hours are as assumed.
      print *
      print *,'CHECKING THE STARTING AND ENDING HOURS'

      if (ihr_start.ne.1) then
        print *,'you must define the starting hour to be 1'
        stop
      endif
c Define whether the daily input hours look like:
c   1, 2, ... , 22, 23, 0: iflag_1_24 = 0, or
c   1, 2, ... , 22, 23, 24: iflag_1_24 = 1
      if (ihr_end.eq.0) then
        iflag_1_24 = 0
      elseif (ihr_end.eq.24) then
        iflag_1_24 = 1
      else
        print *,'you must define the ending hour to be 0 or 24'
        stop
      endif

      if (ihr(1).ne.1) then
        print *,'the starting hour in the data must be 1'
        stop
      endif
      if (iflag_1_24.eq.1) then
        if (ihr(nrecs).ne.24) then
          print *,'the ending hour in the data must be 24'
          stop
        endif
      elseif (iflag_1_24.eq.0) then
        if (ihr(nrecs).ne.0) then
          print *,'the ending hour in the data must be 0'
          stop
        endif
      endif

c Save a copy of the unmodified data inputs.
      do k=1,maxiter
        tair_orig(k) = tair(k)
        rh_orig(k) = rh(k)
        wspd_orig(k) = wspd(k)
        wdir_orig(k) = wdir(k)
        prec_orig(k) = prec(k)
      enddo

c RUN THE AIR TEMPERATURE CHECK.
      if (tair_check_flag.eq.1) then
        print *
        print *,'PRINTING TAIR SUMMARY'
        call check_tair(Tmin_hi,Tmin_lo,Tmax_hi,Tmax_lo,
     &    Tair_hrly_maxinc,maxiter,iiyr,iimo,iidy,tair,Tmin_daily,
     &    Tmax_daily,undef,iihr,so_hem)
      endif

c RUN THE RELATIVE HUMIDITY CHECK.
      if (rh_check_flag.eq.1) then
        print *
        print *,'PRINTING RELH SUMMARY'
        call check_rh(RHmin,RHmax,RH_hrly_maxinc,maxiter,rh,undef,
     &    iiyr,iimo,iidy,iihr)
      endif

c RUN THE WIND SPEED CHECK.
      if (wspd_check_flag.eq.1) then
        print *
        print *,'PRINTING WSPD SUMMARY'
        call check_windspeed(WSPDmin,WSPDmax,
     &    WSPD_hrly_maxinc,maxiter,wspd,undef,iiyr,iimo,iidy,iihr)
      endif

c RUN THE WIND DIRECTION CHECK.
      if (wdir_check_flag.eq.1) then
        print *
        print *,'PRINTING WDIR SUMMARY'
        call check_winddir(WDIRmin,WDIRmax,
     &    WDIR_hrly_maxinc,maxiter,wdir,undef,iiyr,iimo,iidy,iihr)
      endif

c RUN THE PRECIPITATION CHECK.
      if (prec_check_flag.eq.1) then
        print *
        print *,'PRINTING PREC SUMMARY'
        call check_precip(PRECmin,PRECmax,
     &    PREC_hrly_maxinc,maxiter,prec,undef,iiyr,iimo,iidy,iihr)
      endif

c Check to make sure there is no missing data in the first and last
c   data records.  This is required for the missing-data filling
c   procedure.
      if (tair_fill_flag.eq.1 .or. rh_fill_flag.eq.1 .or.
     &  wspd_fill_flag.eq.1 .or. wdir_fill_flag.eq.1 .or.
     &  prec_fill_flag.eq.1) then

        print *
        print *,'CHECKING STARTING AND ENDING RECORDS'

        if (tair_fill_flag.eq.1) then
          if (tair(1).eq.undef) then
            print *,'missing tair at start time'
            print *,'  making tair start data'
            print *
            call make_start_data(maxiter,tair,undef)
          endif
          if (tair(maxiter).eq.undef) then
            print *,'missing tair at end time'
            print *,'  making tair end data'
            print *
            call make_end_data(maxiter,tair,undef)
          endif
        endif

        if (rh_fill_flag.eq.1) then
          if (rh(1).eq.undef) then
            print *,'missing rh at start time'
            print *,'  making rh start data'
            print *
            call make_start_data(maxiter,rh,undef)
          endif
          if (rh(maxiter).eq.undef) then
            print *,'missing rh at end time'
            print *,'  making rh end data'
            print *
            call make_end_data(maxiter,rh,undef)
          endif
        endif

        if (wspd_fill_flag.eq.1) then
          if (wspd(1).eq.undef) then
            print *,'missing wspd at start time'
            print *,'  making wspd start data'
            print *
            call make_start_data(maxiter,wspd,undef)
          endif
          if (wspd(maxiter).eq.undef) then
            print *,'missing wspd at end time'
            print *,'  making wspd end data'
            print *
            call make_end_data(maxiter,wspd,undef)
          endif
        endif

        if (wdir_fill_flag.eq.1) then
          if (wdir(1).eq.undef) then
            print *,'missing wdir at start time'
            print *,'  making wdir start data'
            print *
            call make_start_data(maxiter,wdir,undef)
          endif
          if (wdir(maxiter).eq.undef) then
            print *,'missing wdir at end time'
            print *,'  making wdir end data'
            print *
            call make_end_data(maxiter,wdir,undef)
          endif
        endif

        if (prec_fill_flag.eq.1) then
          if (prec(1).eq.undef) then
            print *,'missing prec at start time'
            print *,'  making prec start data'
            print *
            call make_start_data(maxiter,prec,undef)
          endif
          if (prec(maxiter).eq.undef) then
            print *,'missing prec at end time'
            print *,'  making prec end data'
            print *
            call make_end_data(maxiter,prec,undef)
          endif
        endif

c Run the data filling procedure to fill missing values with valid
c   numbers.  The wdir_flag=1.0 tells the code to do special wind
c   direction averaging to handle the 0-360 value range.
        print *
        print *,'RUNNING THE DATA FILL PROGRAMS'
        print *

      endif

c TAIR:
      if (tair_fill_flag.eq.1) then

c Use the ARIMA program to fill the missing values.
        print *,'RUNNING tair DATA FILL'
        call fill_missing_arima(maxiter,tair,undef,predicted1,
     &    predicted2,y1,y2,var1,var2,max_fill_hrs,ihrs_in_day,'TAIR',
     &    iiyr,iimo,iidy,max_fill_days,max_arima_fill_days)

        if (fill_full_flag.gt.0) then
          print *,'tair: FILLING ANY Case B, C, AND/OR D SEGMENTS'
          print *
          wdir_flag = 0.0

          if (fill_full_flag.eq.1) then
c Use the average of the values before and after any remaining
c   missing data sections (Case B, C, and D) to fill the missing
c   values.
            call fill_missing_diurnal(maxiter,tair,undef,max_fill_hrs,
     &        'TAIR',ihrs_in_day,max_fill_days)

          elseif (fill_full_flag.eq.2) then
c Note: If you don't like the diurnal cycle the above produces, the
c   subroutine below will produce flat lines across any remaining
c   missing data segments.
            call fill_missing_ave(maxiter,tair,undef,max_fill_hrs,
     &        wdir_flag,'TAIR',julian_start,ihrs_in_day)
          endif
        endif

      endif

c RH:
      if (rh_fill_flag.eq.1) then

c Use the ARIMA program to fill the missing values.
        print *,'RUNNING  rh  DATA FILL'
        call fill_missing_arima(maxiter,rh,undef,predicted1,
     &    predicted2,y1,y2,var1,var2,max_fill_hrs,ihrs_in_day,'RELH',
     &    iiyr,iimo,iidy,max_fill_days,max_arima_fill_days)

        if (fill_full_flag.gt.0) then
          print *,'relh: FILLING ANY Case B, C, AND/OR D SEGMENTS'
          print *
          wdir_flag = 0.0

          if (fill_full_flag.eq.1) then
c Use the average of the values before and after any remaining
c   missing data sections (Case B, C, and D) to fill the missing
c   values.
            call fill_missing_diurnal(maxiter,rh,undef,max_fill_hrs,
     &        'RELH',ihrs_in_day,max_fill_days)

          elseif (fill_full_flag.eq.2) then
c Note: If you don't like the diurnal cycle the above produces, the
c   subroutine below will produce flat lines across any remaining
c   missing data segments.
            call fill_missing_ave(maxiter,rh,undef,max_fill_hrs,
     &        wdir_flag,'RELH',julian_start,ihrs_in_day)
          endif
        endif

      endif

c WSPD:
      if (wspd_fill_flag.eq.1) then

c Use the ARIMA program to fill the missing values.
        print *,'RUNNING wspd DATA FILL'
        call fill_missing_arima(maxiter,wspd,undef,predicted1,
     &    predicted2,y1,y2,var1,var2,max_fill_hrs,ihrs_in_day,'WSPD',
     &    iiyr,iimo,iidy,max_fill_days,max_arima_fill_days)

        if (fill_full_flag.gt.0) then
          print *,'wspd: FILLING ANY Case B, C, AND/OR D SEGMENTS'
          print *
          wdir_flag = 0.0

          if (fill_full_flag.eq.1) then
c Use the average of the values before and after any remaining
c   missing data sections (Case B, C, and D) to fill the missing
c   values.
            call fill_missing_diurnal(maxiter,wspd,undef,max_fill_hrs,
     &        'WSPD',ihrs_in_day,max_fill_days)

          elseif (fill_full_flag.eq.2) then
c Note: If you don't like the diurnal cycle the above produces, the
c   subroutine below will produce flat lines across any remaining
c   missing data segments.
            call fill_missing_ave(maxiter,wspd,undef,max_fill_hrs,
     &        wdir_flag,'WSPD',julian_start,ihrs_in_day)
          endif
        endif

      endif

c WDIR:
      if (wdir_fill_flag.eq.1) then

c Use the ARIMA program to fill the missing values.
        print *,'RUNNING wdir DATA FILL'
        call fill_missing_arima(maxiter,wdir,undef,predicted1,
     &    predicted2,y1,y2,var1,var2,max_fill_hrs,ihrs_in_day,'WDIR',
     &    iiyr,iimo,iidy,max_fill_days,max_arima_fill_days)

        if (fill_full_flag.gt.0) then
          print *,'wdir: FILLING ANY Case B, C, AND/OR D SEGMENTS'
          print *
          wdir_flag = 1.0

          if (fill_full_flag.eq.1) then
c Use the average of the values before and after any remaining
c   missing data sections (Case B, C, and D) to fill the missing
c   values.
            call fill_missing_diurnal(maxiter,wdir,undef,max_fill_hrs,
     &        'WDIR',ihrs_in_day,max_fill_days)

          elseif (fill_full_flag.eq.2) then
c Note: If you don't like the diurnal cycle the above produces, the
c   subroutine below will produce flat lines across any remaining
c   missing data segments.
            call fill_missing_ave(maxiter,wdir,undef,max_fill_hrs,
     &        wdir_flag,'WDIR',julian_start,ihrs_in_day)
          endif
        endif

      endif

c PREC:
      if (prec_fill_flag.eq.1) then

c Use the ARIMA program to fill the missing values.
        print *,'RUNNING prec DATA FILL'
        call fill_missing_arima(maxiter,prec,undef,predicted1,
     &    predicted2,y1,y2,var1,var2,max_fill_hrs,ihrs_in_day,'PREC',
     &    iiyr,iimo,iidy,max_fill_days,max_arima_fill_days)

        if (fill_full_flag.gt.0) then
          print *,'prec: FILLING ANY Case B, C, AND/OR D SEGMENTS'
          print *
          wdir_flag = 0.0

          if (fill_full_flag.eq.1) then
c Use the average of the values before and after any remaining
c   missing data sections (Case B, C, and D) to fill the missing
c   values.
            call fill_missing_diurnal(maxiter,prec,undef,max_fill_hrs,
     &        'PREC',ihrs_in_day,max_fill_days)

          elseif (fill_full_flag.eq.2) then
c Note: If you don't like the diurnal cycle the above produces, the
c   subroutine below will produce flat lines across any remaining
c   missing data segments.
            call fill_missing_ave(maxiter,prec,undef,max_fill_hrs,
     &        wdir_flag,'PREC',julian_start,ihrs_in_day)
          endif
        endif

      endif

c Perform some data cleanup.
      do k=1,maxiter

c Constrain the relativel humidity.
        if (rh(k).ne.undef) then
          rh(k) = min(100.0,rh(k))
          rh(k) = max(0.0,rh(k))
        endif

c Clip negative wind speeds.
        if (wspd(k).ne.undef) wspd(k) = max(0.0,wspd(k))

c Make sure the wind directions range between 0 and 360 degrees.
        if (wdir(k).ge.360.0) wdir(k) = wdir(k) - 360.0
        if (wdir(k).lt.0.0 .and. wdir(k).ne.undef)
     &    wdir(k) = wdir(k) + 360.0

c Clip negative precipitation.
        if (prec(k).ne.undef) prec(k) = max(0.0,prec(k))

      enddo

      print *

c Save the data.
      do julian=julian_start,julian_end
        do k=1,ihrs_in_day
          kk = (julian - julian_start) * ihrs_in_day + k

c Write out a text file.
          if (print_flag.eq.1 .or. print_flag.eq.3) then
            write (75,80) iiyr(kk),iimo(kk),iidy(kk),real(iihr(kk)),
     &        id_stn,x_loc,y_loc,elev,tair(kk),rh(kk),wspd(kk),
     &        wdir(kk),prec(kk)
          endif

c Write out a grads file.
          if (print_flag.eq.2 .or. print_flag.eq.3) then
            write (65,rec=kk) iiyr(kk),iimo(kk),iidy(kk),iihr(kk),
     &        tair(kk),rh(kk),wspd(kk),wdir(kk),prec(kk),
     &        tair_orig(kk),rh_orig(kk),wspd_orig(kk),wdir_orig(kk),
     &        prec_orig(kk),Tmin_daily(kk),Tmax_daily(kk)
          endif

        enddo
      enddo

   80 format(i5,i3,i3,f6.2,i6,2f10.1,f8.1,5f9.2)

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine read_param(infname,outfname_txt,outfname_grads,
     &  iyr_start,imo_start,idy_start,ihr_start,iyr_end,imo_end,
     &  idy_end,ihr_end,undef,Tmin_hi,Tmin_lo,Tmax_hi,Tmax_lo,
     &  Tair_hrly_maxinc,RHmin,RHmax,RH_hrly_maxinc,WSPDmin,
     &  WSPDmax,WSPD_hrly_maxinc,WDIRmin,WDIRmax,WDIR_hrly_maxinc,
     &  PRECmin,PRECmax,PREC_hrly_maxinc,max_fill_days,iheader,
     &  id_stn,x_loc,y_loc,elev,tair_check_flag,rh_check_flag,
     &  wspd_check_flag,wdir_check_flag,prec_check_flag,
     &  tair_fill_flag,rh_fill_flag,wspd_fill_flag,wdir_fill_flag,
     &  prec_fill_flag,print_flag,max_arima_fill_days,so_hem,
     &  fill_full_flag)

c This program reads and processes the input parameter data.
c
c The following must be true:
c
c   All comment lines start with a ! in the first position.
c
c   Blank lines are permitted.
c
c   All parameter statements start with the parameter name, followed
c   by a space, followed by an =, followed by a space, followed by
c   the actual value, with nothing after that.  These statements can
c   have leading blanks, and must fit within 80 columns.
c
c   It is set up to deal with integer, real, and string input values.

      implicit none

c Input parameter names.
      integer max_fill_days,ihr_start,ihr_end,iyr_start,imo_start,
     &  idy_start,iyr_end,imo_end,idy_end,iheader,id_stn,
     &  tair_check_flag,rh_check_flag,wspd_check_flag,
     &  wdir_check_flag,prec_check_flag,tair_fill_flag,rh_fill_flag,
     &  wspd_fill_flag,wdir_fill_flag,prec_fill_flag,print_flag,
     &  max_arima_fill_days,fill_full_flag

      real undef,Tmin_hi,Tmin_lo,Tmax_hi,Tmax_lo,Tair_hrly_maxinc,
     &  RHmin,RHmax,RH_hrly_maxinc,WSPDmin,WSPDmax,WSPD_hrly_maxinc,
     &  WDIRmin,WDIRmax,WDIR_hrly_maxinc,PRECmin,PRECmax,
     &  PREC_hrly_maxinc,x_loc,y_loc,elev,so_hem

      character*80 infname
      character*80 outfname_txt
      character*80 outfname_grads

c Working parameter names.
      character*80 input_string,c_param,c_value
      integer k,max_par_lines,i_param_chars,i_value_chars,
     &  icomment_flag

      open (40,file='preproc_MicroMet.par',status='old')

      max_par_lines = 1000

      do k=1,max_par_lines
        read (40,'(a80)',end=99) input_string

        call get_param_data(input_string,c_param,c_value,
     &    i_param_chars,i_value_chars,icomment_flag)

c Process the string if it is not a comment.
        if (icomment_flag.eq.0) then

c Integer example.
c         if (c_param(1:i_param_chars).eq.'nx')
c    &      call char2int(nx,i_value_chars,c_value,
c    &        c_param(1:i_param_chars))

c Real example.
c         if (c_param(1:i_param_chars).eq.'deltax')
c    &      call char2real(deltax,i_value_chars,c_value,
c    &        c_param(1:i_param_chars))

c Character example.
c         if (c_param(1:i_param_chars).eq.'fname_out')
c    &      call char2char(fname_out,c_value,i_value_chars,
c    &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'infname')
     &      call char2char(infname,c_value,i_value_chars,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'outfname_txt')
     &      call char2char(outfname_txt,c_value,i_value_chars,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'outfname_grads')
     &      call char2char(outfname_grads,c_value,i_value_chars,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'print_flag')
     &      call char2int(print_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'iyr_start')
     &      call char2int(iyr_start,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'imo_start')
     &      call char2int(imo_start,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'idy_start')
     &      call char2int(idy_start,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'ihr_start')
     &      call char2int(ihr_start,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'iyr_end')
     &      call char2int(iyr_end,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'imo_end')
     &      call char2int(imo_end,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'idy_end')
     &      call char2int(idy_end,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'ihr_end')
     &      call char2int(ihr_end,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'iheader')
     &      call char2int(iheader,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'undef')
     &      call char2real(undef,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'so_hem')
     &      call char2real(so_hem,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'Tmin_hi')
     &      call char2real(Tmin_hi,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'Tmin_lo')
     &      call char2real(Tmin_lo,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'Tmax_hi')
     &      call char2real(Tmax_hi,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'Tmax_lo')
     &      call char2real(Tmax_lo,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'Tair_hrly_maxinc')
     &      call char2real(Tair_hrly_maxinc,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'tair_check_flag')
     &      call char2int(tair_check_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'RHmin')
     &      call char2real(RHmin,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'RHmax')
     &      call char2real(RHmax,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'RH_hrly_maxinc')
     &      call char2real(RH_hrly_maxinc,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'rh_check_flag')
     &      call char2int(rh_check_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'WSPDmin')
     &      call char2real(WSPDmin,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'WSPDmax')
     &      call char2real(WSPDmax,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'WSPD_hrly_maxinc')
     &      call char2real(WSPD_hrly_maxinc,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'wspd_check_flag')
     &      call char2int(wspd_check_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'WDIRmin')
     &      call char2real(WDIRmin,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'WDIRmax')
     &      call char2real(WDIRmax,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'WDIR_hrly_maxinc')
     &      call char2real(WDIR_hrly_maxinc,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'wdir_check_flag')
     &      call char2int(wdir_check_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'PRECmin')
     &      call char2real(PRECmin,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'PRECmax')
     &      call char2real(PRECmax,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'PREC_hrly_maxinc')
     &      call char2real(PREC_hrly_maxinc,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'prec_check_flag')
     &      call char2int(prec_check_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'max_fill_days')
     &      call char2int(max_fill_days,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'max_arima_fill_days')
     &      call char2int(max_arima_fill_days,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'fill_full_flag')
     &      call char2int(fill_full_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'tair_fill_flag')
     &      call char2int(tair_fill_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'rh_fill_flag')
     &      call char2int(rh_fill_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'wspd_fill_flag')
     &      call char2int(wspd_fill_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'wdir_fill_flag')
     &      call char2int(wdir_fill_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'prec_fill_flag')
     &      call char2int(prec_fill_flag,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'id_stn')
     &      call char2int(id_stn,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'x_loc')
     &      call char2real(x_loc,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'y_loc')
     &      call char2real(y_loc,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

          if (c_param(1:i_param_chars).eq.'elev')
     &      call char2real(elev,i_value_chars,c_value,
     &        c_param(1:i_param_chars))

        endif
      enddo

  99  continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine make_start_data(maxiter,var,undef)

      implicit none

      integer maxiter,icount,k
      real undef,sum
      real var(maxiter)

c If the starting variable is undefined, calculate the average of
c   the closest (in time) 24 valid values and use that to fill in the
c   missing value.
      icount = 0
      sum = 0.0
      do k=2,maxiter
        if (var(k).ne.undef) then
          icount = icount + 1
          sum = sum + var(k)
          if (icount.eq.24) go to 77
        endif
      enddo
      print *,'DIDN`T FIND 24 VALID VALUES IN make_start_data'
      stop

  77  var(1) = sum / real(icount)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine make_end_data(maxiter,var,undef)

      implicit none

      integer maxiter,icount,k
      real undef,sum
      real var(maxiter)

c If the ending variable is undefined, calculate the average of
c   the closest (in time) 24 valid values and use that to fill in the
c   missing value.
      icount = 0
      sum = 0.0
      do k=maxiter-1,1,-1
        if (var(k).ne.undef) then
          icount = icount + 1
          sum = sum + var(k)
          if (icount.eq.24) go to 77
        endif
      enddo
      print *,'DIDN`T FIND 24 VALID VALUES IN make_end_data'
      stop

  77  var(maxiter) = sum / real(icount)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fill_missing_diurnal(maxiter,var,undef,max_fill_hrs,
     &  vname,ihrs_in_day,max_fill_days)

      implicit none

      character*4 vname

      integer maxiter,k,kk,location,n_undefs,kkk,max_fill_days,
     &  location1,location2,max_fill_hrs,ihrs_in_day

      real undef,var1_tmp,var2_tmp,var_u_ave,var_v_ave,pi,
     &  deg2rad,rad2deg
      real var(maxiter)

      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0
      rad2deg = 180.0 / pi

c Read through the data record and fill any missing data.  Do this
c   by averaging the first valid data found at 24 hour increments
c   before and after the missing data point.

c Loop though the data, searching for groupings of consecutive
c   undefined values.
      do k=2,maxiter-1
        if (var(k).eq.undef .and. var(k-1).ne.undef) then
          n_undefs = 1
          do kk=1,maxiter-k+1
            location = k + kk
            if (var(location).eq.undef) then
              n_undefs = n_undefs + 1
            else
c As soon as you find a valid value, exit this counting.
              goto 77
            endif
          enddo
  77      continue

c Do the data filling.  Increment forward and backwards in 24-hour
c   steps, up to max_fill_days, until you find valid values to be
c   used in the averages.
          do kk=1,n_undefs
            do kkk=1,max_fill_days
              location1 = k - kkk * ihrs_in_day + kk - 1
              if (location1.lt.1) then
                var1_tmp = undef
                goto 78
              endif
              if (var(location1).ne.undef) then
                var1_tmp = var(location1)
                goto 78
              endif
            enddo
  78        continue
            do kkk=1,max_fill_days
              location2 = k + kkk * ihrs_in_day + kk - 1
              if (location2.gt.maxiter) then
                var2_tmp = undef
                goto 79
              endif
              if (var(location2).ne.undef) then
                var2_tmp = var(location2)
                goto 79
              endif
            enddo
  79        continue
            if (var1_tmp.eq.undef .and. var2_tmp.ne.undef) then
              var(k+kk-1) = var2_tmp
            elseif (var1_tmp.ne.undef .and. var2_tmp.eq.undef) then
              var(k+kk-1) = var1_tmp
            else
              if (vname.eq.'WDIR') then
                var_u_ave = 0.5 * (- sin(deg2rad*var1_tmp) -
     &            sin(deg2rad*var2_tmp))
                var_v_ave = 0.5 * (- cos(deg2rad*var1_tmp) -
     &            cos(deg2rad*var2_tmp))
                var(k+kk-1) = 270.0 -
     &            rad2deg*atan2(var_v_ave,var_u_ave)
                if (var(k+kk-1).ge.360.0) var(k+kk-1) =
     &            var(k+kk-1) - 360.0
              else
                var(k+kk-1) = 0.5 * (var1_tmp + var2_tmp)
              endif
            endif
          enddo
        endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fill_missing_ave(maxiter,var,undef,max_fill_hrs,
     &  wdir_flag,vname,julian_start,ihrs_in_day)

      implicit none

      character*4 vname

      integer maxiter,k,kk,location,n_undefs,kkk,julian_start,
     &  location1,location2,max_fill_hrs,n_counts,julian,ioptn,
     &  iday,imonth,iyear,ihrs_in_day

      real undef,ave,sum,wdir_flag,pi,deg2rad,rad2deg,sum_u,sum_v,
     &  ave_u,ave_v,ave_dir
      real var(maxiter)

      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0
      rad2deg = 180.0 / pi

c Read through the data record and fill any missing data.  Count the
c   number of consecutive missing values, average that many values
c   before and after the missing segment (if a missing value is
c   found during this sweep, only use the available values up to
c   that point), and replace all missing values in that segment with
c   the average of the available before and after values.
c
c When searching for missing segments in the data stream, first
c   look for and fix single missing values, then two missing values,
c   etc., up to some maximum allowed missing segment length.  This
c   way each time through it makes the data cleaner for the next
c   cycle.

c Loop though the data, searching for groupings of consecutive
c   undefined values.
      do kkk=1,max_fill_hrs
        do k=2,maxiter-1
          if (var(k).eq.undef .and. var(k-1).ne.undef) then
            n_undefs = 1
            do kk=1,maxiter-k+1
              location = k + kk
              if (var(location).eq.undef) then
                n_undefs = n_undefs + 1
              else
c As soon as you find a valid value, exit this counting.
                goto 77
              endif
            enddo
  77        continue

c Do the data filling.
            if (n_undefs.eq.kkk) then
              sum = 0.0
              sum_u = 0.0
              sum_v = 0.0
              n_counts = 0

c Search BEFORE the missing section.  Escape if you find a missing
c   value before coming to the end of the missing length.
              do kk=1,n_undefs
                location1 = k - kk
                if (var(location1).eq.undef .or. location1.lt.1)
     &            goto 91
                n_counts = n_counts + 1
                if (wdir_flag.eq.1.0) then
                  sum_u = sum_u - sin(deg2rad*var(location1))
                  sum_v = sum_v - cos(deg2rad*var(location1))
                else
                  sum = sum + var(location1)
                endif
              enddo

  91  continue

c Search AFTER the missing section.  Escape if you find a missing
c   value before coming to the end of the missing length.
              do kk=1,n_undefs
                location2 = k + n_undefs - 1 + kk
                if (var(location2).eq.undef .or. location2.gt.maxiter)
     &            goto 92
                n_counts = n_counts + 1
                if (wdir_flag.eq.1.0) then
                  sum_u = sum_u - sin(deg2rad*var(location2))
                  sum_v = sum_v - cos(deg2rad*var(location2))
                else
                  sum = sum + var(location2)
                endif
              enddo

  92  continue

              if (wdir_flag.eq.1.0) then
                ave_u = sum_u / (2.0 * real(n_undefs))
                ave_v = sum_v / (2.0 * real(n_undefs))
                ave_dir = 270.0 - rad2deg*atan2(ave_v,ave_u)
                if (ave_dir.ge.360.0) ave_dir = ave_dir - 360.0
                do kk=1,n_undefs
                  var(k+kk-1) = ave_dir
                enddo
              else
                ave = sum / real(n_counts)
                do kk=1,n_undefs
                  var(k+kk-1) = ave
                enddo
              endif

            elseif (kkk.eq.max_fill_hrs .and.
     &        n_undefs.gt.max_fill_hrs) then

c Use k to find the date that this missing segment starts.
              julian =  julian_start + k /ihrs_in_day
              ioptn = 4
              call calndr (ioptn,iday,imonth,iyear,julian)
              write (*,101) vname,iyear,imonth,iday,
     &          real(n_undefs)/24.0

            endif
          endif
        enddo
      enddo

c Exceeded max undef values.
 101  format (a5,' cannot define missing data [case A]: starting at',
     &  i5,i3,i3,', for',f5.1,' days')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fill_missing_arima(maxiter,var,undef,predicted1,
     &  predicted2,y1,y2,var1,var2,max_fill_hrs,ihrs_in_day,vname,
     &  iiyr,iimo,iidy,max_fill_days,max_arima_fill_days)

      implicit none

      character*4 vname

      integer maxiter,k,kk,ihrs_in_day,location,n_undefs,n,i1,i2,
     &  n_undef_start,location1,location2,nfcst,max_fill_hrs,i,ii,
     &  i1_min,i2_max,max_fill_days,kkk,ii1,ii2,max_arima_fill_days

      integer iiyr(maxiter)
      integer iimo(maxiter)
      integer iidy(maxiter)

      real undef,w1,w2,var1_tmp,var2_tmp
      real var(maxiter)
      real var1(maxiter)
      real var2(maxiter)
      real predicted1(maxiter)
      real predicted2(maxiter)
      real y1(maxiter) 
      real y2(maxiter) 
      real y_sum_u,y_sum_v,y_ave_u,y_ave_v,y_ave,var_u_ave,var_v_ave
      real pi,deg2rad,rad2deg

c Build some temporary variable arrays.
      do k=1,maxiter
        var1(k) = var(k)
        var2(k) = var(k)
      enddo

c Read through the data record and fill any missing data.  Consider
c   three cases:
c   1) for missing data = 1 hour, average the values from each side
c        of that hour.
c   2) for 2 hours <= missing data <= 24 hours, average the values
c        from 24 hours before and after each of the hours in that
c        period.
c   3) for missing data > 24 hours, use the arima model.  Here look
c        the missing data number of hours before and after the
c        missing period to fill the missing period.

c Search the entire file three times in a row.  The first time take
c   care of Case 1, the second time Case 2, and the third time
c   Case 3.  This way each time through it makes things cleaner for
c   the next cycle.

c CASE 1.
c Read through the entire data string, identifying any missing data.
c   When a missing point is found, look ahead and behind to see
c   whether the adjacent points are missing.  If not, fill it in.
      if (var(1).eq.undef .or. var(maxiter).eq.undef)
     &  print *, 'THE FIRST AND LAST DATA POINTS CANNOT BE UNDEFINED'
      do k=2,maxiter-1
        if (var(k-1).ne.undef .and. var(k).eq.undef .and.
     &    var(k+1).ne.undef) then
          if (vname.eq.'WDIR') then
            var_u_ave = 0.5 * (- sin(deg2rad*var(k-1)) -
     &        sin(deg2rad*var(k+1)))
            var_v_ave = 0.5 * (- cos(deg2rad*var(k-1)) -
     &        cos(deg2rad*var(k+1)))
            var(k) = 270.0 - rad2deg*atan2(var_v_ave,var_u_ave)
            if (var(k).ge.360.0) var(k) = var(k) - 360.0
          else
            var(k) = 0.5 * (var(k-1) + var(k+1))
          endif
        endif
      enddo

c CASE 2.
c Loop though the data, searching for groupings of less than 25
c   consecutive undefined values.
      do k=2,maxiter-1
        if (var(k).eq.undef .and. var(k-1).ne.undef) then
          n_undefs = 1
          do kk=1,ihrs_in_day
            location = k + kk
            if (var(location).eq.undef) then
              n_undefs = n_undefs + 1
            else
c As soon as you find a valid value, exit this counting.
              goto 77
            endif
          enddo
  77      continue
c Do the data filling.  Increment forward and backwards in 24-hour
c   steps, up to max_fill_days, until you find valid values to be
c   used in the averages.
          if (n_undefs.le.ihrs_in_day) then
            do kk=1,n_undefs
              do kkk=1,max_fill_days
                location1 = k - kkk * ihrs_in_day + kk - 1
                if (location1.lt.1) then
                  var1_tmp = undef
                  goto 78
                endif
                if (var(location1).ne.undef) then
                  var1_tmp = var(location1)
                  goto 78
                endif
              enddo
  78          continue
              do kkk=1,max_fill_days
                location2 = k + kkk * ihrs_in_day + kk - 1
                if (location2.gt.maxiter) then
                  var2_tmp = undef
                  goto 79
                endif
                if (var(location2).ne.undef) then
                  var2_tmp = var(location2)
                  goto 79
                endif
              enddo
  79          continue
              if (var1_tmp.eq.undef .and. var2_tmp.ne.undef) then
                var(k+kk-1) = var2_tmp
              elseif (var1_tmp.ne.undef .and. var2_tmp.eq.undef) then
                var(k+kk-1) = var1_tmp
              else
                if (vname.eq.'WDIR') then
                  var_u_ave = 0.5 * (- sin(deg2rad*var1_tmp) -
     &              sin(deg2rad*var2_tmp))
                  var_v_ave = 0.5 * (- cos(deg2rad*var1_tmp) -
     &              cos(deg2rad*var2_tmp))
                  var(k+kk-1) = 270.0 -
     &              rad2deg*atan2(var_v_ave,var_u_ave)
                  if (var(k+kk-1).ge.360.0) var(k+kk-1) =
     &              var(k+kk-1) - 360.0
                else
                  var(k+kk-1) = 0.5 * (var1_tmp + var2_tmp)
                endif
              endif
            enddo
          endif
        endif
      enddo

c CASE 3.
c Search for groupings of greater than 24 missing values.  Any
c   groupings of less than 25 have already been taken care of, so
c   we don't have to worry about the small stuff anymore.
      do k=2,maxiter-1
        if (var(k).eq.undef .and. var(k-1).ne.undef) then
          n_undef_start = k
          n_undefs = 1
          do kk=1,maxiter-k+1
            location = k + kk
            if (var(location).eq.undef) then
              n_undefs = n_undefs + 1
            else
c As soon as you find a valid value, exit this counting.
              goto 88
            endif
          enddo
  88      continue

c Do the data filling using the arima model.

c First figure out whether the model will be looking outside the
c   range of the available data set.
          i1_min = n_undef_start - 1 - n_undefs
          i2_max = n_undefs - 1 + n_undef_start + n_undefs

c Print a message stating that the data cannot be filled for the
c   case where the forecast was looking outside the range of
c   available data at the begining and/or end of the data record.
          if (i1_min.lt.1) then
            write (*,102) vname,iiyr(k),iimo(k),iidy(k),
     &        real(n_undefs)/24.0

          elseif (i2_max.gt.maxiter) then
            write (*,103) vname,iiyr(k),iimo(k),iidy(k),
     &        real(n_undefs)/24.0

          elseif (n_undefs.le.max_fill_hrs .and. i1_min.ge.1 .and.
     &      i2_max.le.maxiter) then

c Extract the obs to be used to generate the prediction.
c   1 indicates pre-missing, 2 indicates post-missing.  Make
c   sure the 2 data are running backwards.  Note that the arima
c   program will produce garbage if nfcst is greater then n.
c Also perform a check to make sure that the obs used as part of
c   the prediction calculation contain no undefined values.  If
c   this happens do not perform the data fill for the missing
c    section.
            nfcst = n_undefs
            n = n_undefs
            do i=1,n
              i1 = i - 1 + n_undef_start - n_undefs
              y1(i) = var(i1)
              if (y1(i).eq.undef) then
                write (*,104) vname,iiyr(k),iimo(k),iidy(k),
     &            real(n_undefs)/24.0
                goto 98
              endif
              i2 = n - i + n_undef_start + n_undefs
              y2(i) = var(i2)
              if (y2(i).eq.undef) then
                write (*,103) vname,iiyr(k),iimo(k),iidy(k),
     &            real(n_undefs)/24.0
                goto 98
              endif
            enddo

c To handle the problems with the 360-0 line, average the wind
c   direction data and, if needed, add a constant value to get
c   most of the data away from the 360-0 line.
            if (vname.eq.'WDIR') then
              pi = 2.0 * acos(0.0)
              deg2rad = pi / 180.0
              rad2deg = 180.0 / pi

              y_sum_u = 0.0
              y_sum_v = 0.0
              do i=1,n
                y_sum_u = y_sum_u - sin(deg2rad*y1(i)) -
     &            sin(deg2rad*y2(i))
                y_sum_v = y_sum_v - cos(deg2rad*y1(i)) -
     &            cos(deg2rad*y2(i))
              enddo
              y_ave_u = y_sum_u / (2.0 * real(n))
              y_ave_v = y_sum_v / (2.0 * real(n))

              y_ave = 270.0 - rad2deg*atan2(y_ave_v,y_ave_u)
              if (y_ave.ge.360.0) y_ave = y_ave - 360.0

c Here I am making the discontinuity occur 180 degree off of the
c   data average.
              if (y_ave.lt.90.0) then
                do i=1,n
                  if (y1(i).lt.y_ave+180.0) then
                    y1(i) = y1(i) + 360.0
                  endif
                  if (y2(i).lt.y_ave+180.0) then
                    y2(i) = y2(i) + 360.0
                  endif
                enddo
              elseif (y_ave.lt.180.0) then
                do i=1,n
                  if (y1(i).lt.y_ave+180.0) then
                    y1(i) = y1(i) + 360.0
                  endif
                  if (y2(i).lt.y_ave+180.0) then
                    y2(i) = y2(i) + 360.0
                  endif
                enddo
              elseif (y_ave.lt.270.0) then
                do i=1,n
                  if (y1(i).lt.y_ave-180.0) then
                    y1(i) = y1(i) + 360.0
                  endif
                  if (y2(i).lt.y_ave-180.0) then
                    y2(i) = y2(i) + 360.0
                  endif
                enddo
              elseif (y_ave.le.360.0) then
                do i=1,n
                  if (y1(i).lt.y_ave-180.0) then
                    y1(i) = y1(i) + 360.0
                  endif
                  if (y2(i).lt.y_ave-180.0) then
                    y2(i) = y2(i) + 360.0
                  endif
                enddo
              endif
            endif

c Perform the arima forecasting analysis.
c   Forward.
            call forecast (y1,n,nfcst,predicted1)

c   Backward.
            call forecast (y2,n,nfcst,predicted2)

c Fill the missing section with the forecast data.
            do i=n_undef_start,n_undef_start+n_undefs-1
              if (n_undefs.le.max_arima_fill_days*ihrs_in_day) then
                if (vname.eq.'WDIR') then
c Do the filling by using only the forecast diurnal cycle, not the
c   trend over the missing period.
                  ii1 = i + 1 - n_undef_start
                  i1 = mod(ii1,24)
                  if (i1.eq.0) i1 = 24
                  var1(i) = predicted1(i1)

                  ii2 = n - (i - n_undef_start)
                  i2 = mod(ii2,24)
                  if (i2.eq.0) i2 = 24
                  var2(i) = predicted2(i2)
                else
                  i1 = i + 1 - n_undef_start
                  var1(i) = predicted1(i1)
                  i2 = n - (i - n_undef_start)
                  var2(i) = predicted2(i2)
                endif
              else
c Do the filling by using only the forecast diurnal cycle, not the
c   trend over the missing period.
                ii1 = i + 1 - n_undef_start
                i1 = mod(ii1,24)
                if (i1.eq.0) i1 = 24
                var1(i) = predicted1(i1)

                ii2 = n - (i - n_undef_start)
                i2 = mod(ii2,24)
                if (i2.eq.0) i2 = 24
                var2(i) = predicted2(i2)
              endif
            enddo

c Merge the two predictions together using the weighting functions.
c   Provide a linear weighting to the predicted values.
            do i=1,n
              w1 = real(n-i+1) / real(n+1)
              w2 = real(i) / real(n+1)
              ii = i - 1 + n_undef_start
              var(ii) = w1 * var1(ii) + w2 * var2(ii)
            enddo

   98       continue

          else

c Print a message stating that the data cannot be filled for the
c   case where the maximum permitted number of undefined values has
c   been exceeded.
            write (*,101) vname,iiyr(k),iimo(k),iidy(k),
     &        real(n_undefs)/24.0

          endif
        endif
      enddo

c Exceeded max undef values.
 101  format (a5,' cannot define missing data [case A]: starting at',
     &  i5,i3,i3,', for',f5.1,' days')

c Arima looking BEFORE available data.
 102  format (a5,' cannot define missing data [case B]: starting at',
     &  i5,i3,i3,', for',f5.1,' days')

c Arima looking AFTER available data.
 103  format (a5,' cannot define missing data [case C]: starting at',
     &  i5,i3,i3,', for',f5.1,' days')

c Arima found missing values.
 104  format (a5,' cannot define missing data [case D]: starting at',
     &  i5,i3,i3,', for',f5.1,' days')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine forecast (y,n,nfcst,predicted)

c This subroutine sets up the call the ARIMA programs (about 85
c   pages of code provided by the Standards Time Series and
c   Regression Package (STARPAC).  ftp://ftp.ucar/edu/starpac

      parameter (nfcsto=1)

      real predicted(n)
      real y(n) 
      real fcst(nfcst,nfcsto)
      integer mspec(4,2)
      real par(3)

      integer ifcsto(nfcsto)

      double precision dstak(5000)
      common /cstak/ dstak
      ldstak = 5000

c Define the input parameters.
c   nfac = number of factors in the model.
      data nfac/2/
c   mspec = orders of p, d, q, and s for each factor in the model.
c   Needs to look like ((mspec(i,j), i=1,4), j=1,nfac).  Note that
c   to handle periodicity, the 24 is for hourly data, 12 would be
c   for monthly data, etc.
c     data mspec/0,1,1,1,0,1,1,12/
      data mspec/0,1,1,1,0,1,1,24/
c   npar = number of parameters in the model.
      data npar/3/
c   par = parameter values used in the model.
      data par/0.0,0.4,0.6/

c Printing control to arima.output.notes.
c     nprt = 22222
c     nprt = 11112
      nprt = 00000

c     CALL IPRINT(IPRT)
c     OPEN (UNIT=IPRT, FILE='arima.output.notes')

      ifcsto(nfcsto) = n
      ifcst = nfcst

c Perform the arima forecasting analysis.
      call AIMFS(Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK,
     &   NFCST, NFCSTO, IFCSTO, NPRT, FCST, IFCST, FCSTSD)

c Extract the predicted data values into a 1-d array.
      do i=1,nfcst
        predicted(i) = fcst(i,nfcsto)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_precip(PRECmin,PRECmax,
     &  PREC_hrly_maxinc,maxiter,prec,undef,iiyr,iimo,iidy,iihr)

      implicit none

      integer iter,maxiter

      integer iiyr(maxiter)
      integer iimo(maxiter)
      integer iidy(maxiter)
      integer iihr(maxiter)

      real PRECmin,PRECmax,PREC_hrly_maxinc,PREC_diff,undef
      real prec(maxiter)

      do iter=1,maxiter

        if (prec(iter).ne.undef) then

c CASE 1: LIM
          if (prec(iter).gt.PRECmax) write (*,101)
     &      prec(iter),PRECmax,iiyr(iter),iimo(iter),iidy(iter),
     &      iihr(iter)
          if (prec(iter).lt.PRECmin) write (*,102)
     &      prec(iter),PRECmin,iiyr(iter),iimo(iter),iidy(iter),
     &      iihr(iter)

c CASE 2: ROC
          if (iter.gt.1) then
            if (prec(iter-1).ne.undef) then
              if (prec(iter).gt.0.0) then
                PREC_diff = abs(prec(iter) - prec(iter-1))
                if (PREC_diff.gt.PREC_hrly_maxinc) write (*,103)
     &            PREC_diff,PREC_hrly_maxinc,iiyr(iter),iimo(iter),
     &            iidy(iter),iihr(iter)
              endif
            endif
          endif

c CASE 3: NOC
          if (iter.gt.3) then
            if (prec(iter-1).ne.undef .and. prec(iter-2).ne. undef
     &        .and. prec(iter-3).ne.undef) then
              if (prec(iter).gt.0.0) then
                if (prec(iter).eq.prec(iter-1) .and.
     &            prec(iter).eq.prec(iter-2) .and. 
     &            prec(iter).eq.prec(iter-3)) write (*,104)
     &            prec(iter),iiyr(iter),iimo(iter),iidy(iter),
     &            iihr(iter)
              endif
            endif
          endif
        endif

      enddo

 101  format ('PREC [',f5.2,'] above bounds [',f5.2,'] at',i5,i3,i3,i3)
 102  format ('PREC [',f5.2,'] below bounds [',f5.2,'] at',i5,i3,i3,i3)
 103  format ('PREC hourly inc [',f5.2,'] above bounds [',f5.2,
     &  '] at',i5,i3,i3,i3)
 104  format ('Constant PREC [',f5.2,'] found at',i5,i3,i3,i3)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_winddir(WDIRmin,WDIRmax,
     &  WDIR_hrly_maxinc,maxiter,wdir,undef,iiyr,iimo,iidy,iihr)

      implicit none

      integer iter,maxiter

      integer iiyr(maxiter)
      integer iimo(maxiter)
      integer iidy(maxiter)
      integer iihr(maxiter)

      real WDIRmin,WDIRmax,WDIR_hrly_maxinc,WDIR_diff,undef
      real wdir(maxiter)

      do iter=1,maxiter

        if (wdir(iter).ne.undef) then

c CASE 1: LIM
          if (wdir(iter).gt.WDIRmax) write (*,101)
     &      wdir(iter),WDIRmax,iiyr(iter),iimo(iter),iidy(iter),
     &      iihr(iter)
          if (wdir(iter).lt.WDIRmin) write (*,102)
     &      wdir(iter),WDIRmin,iiyr(iter),iimo(iter),iidy(iter),
     &      iihr(iter)
  
c CASE 2: ROC
          if (iter.gt.1) then
            if (wdir(iter-1).ne.undef) then
              WDIR_diff = abs(wdir(iter) - wdir(iter-1))
              WDIR_diff = min(WDIR_diff,360.0-WDIR_diff)
              if (WDIR_diff.gt.WDIR_hrly_maxinc) write (*,103)
     &          WDIR_diff,WDIR_hrly_maxinc,iiyr(iter),iimo(iter),
     &          iidy(iter),iihr(iter)
            endif
          endif

c CASE 3: NOC
          if (iter.gt.3) then
            if (wdir(iter-1).ne.undef .and. wdir(iter-2).ne. undef
     &        .and. wdir(iter-3).ne.undef) then
              if (wdir(iter).eq.wdir(iter-1) .and.
     &          wdir(iter).eq.wdir(iter-2) .and. 
     &          wdir(iter).eq.wdir(iter-3)) write (*,104)
     &          wdir(iter),iiyr(iter),iimo(iter),iidy(iter),
     &          iihr(iter)
            endif
          endif
        endif

      enddo

 101  format ('WDIR [',f5.1,'] above bounds [',f5.1,'] at',i5,i3,i3,i3)
 102  format ('WDIR [',f5.1,'] below bounds [',f5.1,'] at',i5,i3,i3,i3)
 103  format ('WDIR hourly inc [',f5.1,'] above bounds [',f5.1,
     &  '] at',i5,i3,i3,i3)
 104  format ('Constant WDIR [',f5.1,'] found at',i5,i3,i3,i3)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_windspeed(WSPDmin,WSPDmax,
     &  WSPD_hrly_maxinc,maxiter,wspd,undef,iiyr,iimo,iidy,iihr)

      implicit none

      integer iter,maxiter

      integer iiyr(maxiter)
      integer iimo(maxiter)
      integer iidy(maxiter)
      integer iihr(maxiter)

      real WSPDmin,WSPDmax,WSPD_hrly_maxinc,WSPD_diff,undef
      real wspd(maxiter)

      do iter=1,maxiter

        if (wspd(iter).ne.undef) then

c CASE 1: LIM
          if (wspd(iter).gt.WSPDmax) write (*,101)
     &      wspd(iter),WSPDmax,iiyr(iter),iimo(iter),iidy(iter),
     &      iihr(iter)
          if (wspd(iter).lt.WSPDmin) write (*,102)
     &      wspd(iter),WSPDmin,iiyr(iter),iimo(iter),iidy(iter),
     &      iihr(iter)

c CASE 2: ROC
          if (iter.gt.1) then
            if (wspd(iter-1).ne.undef) then
              WSPD_diff = abs(wspd(iter) - wspd(iter-1))
              if (WSPD_diff.gt.WSPD_hrly_maxinc) write (*,103)
     &          WSPD_diff,WSPD_hrly_maxinc,iiyr(iter),iimo(iter),
     &          iidy(iter),iihr(iter)
            endif
          endif

c CASE 3: NOC
          if (iter.gt.3) then
            if (wspd(iter-1).ne.undef .and. wspd(iter-2).ne. undef
     &        .and. wspd(iter-3).ne.undef) then
              if (wspd(iter).eq.wspd(iter-1) .and.
     &          wspd(iter).eq.wspd(iter-2) .and. 
     &          wspd(iter).eq.wspd(iter-3)) write (*,104)
     &          wspd(iter),iiyr(iter),iimo(iter),iidy(iter),
     &          iihr(iter)
            endif
          endif
        endif

      enddo

 101  format ('WSPD [',f5.1,'] above bounds [',f5.1,'] at',i5,i3,i3,i3)
 102  format ('WSPD [',f5.1,'] below bounds [',f5.1,'] at',i5,i3,i3,i3)
 103  format ('WSPD hourly inc [',f5.1,'] above bounds [',f5.1,
     &  '] at',i5,i3,i3,i3)
 104  format ('Constant WSPD [',f5.1,'] found at',i5,i3,i3,i3)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_tair(Tmin_hi,Tmin_lo,Tmax_hi,Tmax_lo,
     &  Tair_hrly_maxinc,maxiter,iiyr,iimo,iidy,tair,Tmin_daily,
     &  Tmax_daily,undef,iihr,so_hem)

      implicit none

      integer ioptn,iter,maxiter,idayct
      integer iiyr(maxiter)
      integer iimo(maxiter)
      integer iidy(maxiter)
      integer iihr(maxiter)

      real d,Tmin_hi,Tmin_lo,Tmax_hi,Tmax_lo,Tair_hrly_maxinc,pi,
     &  temp_diff,undef,so_hem,xoffset
      real tair(maxiter)
      real Tmin_daily(maxiter)
      real Tmax_daily(maxiter)

      pi = 2.0 * acos(0.0)

c Have the calendar program calculate the day-of-year.
      ioptn = 1

      do iter=1,maxiter

c Calculate the day_of_year.
        call calndr (ioptn,iidy(iter),iimo(iter),iiyr(iter),idayct)
        d = real(idayct)

c Calculate the air temperature limits.
        if (so_hem.eq.1.0) then
          xoffset = 182.0
        else
          xoffset = 0.0
        endif

        Tmin_daily(iter) = (Tmin_hi+Tmin_lo)/2.0 +
     &    (Tmin_hi-Tmin_lo)/2.0 * cos(2.0*pi*(d-(200.0-xoffset))/365.0)

        Tmax_daily(iter) = (Tmax_hi+Tmax_lo)/2.0 +
     &    (Tmax_hi-Tmax_lo)/2.0 * cos(2.0*pi*(d-(190.0-xoffset))/365.0)

c Run the checks.
        if (tair(iter).ne.undef) then

c CASE 1: LIM
          if (tair(iter).gt.Tmax_daily(iter)) write (*,101)
     &      tair(iter),Tmax_daily(iter),iiyr(iter),iimo(iter),
     &      iidy(iter),iihr(iter)

          if (tair(iter).lt.Tmin_daily(iter)) write (*,102)
     &      tair(iter),Tmin_daily(iter),iiyr(iter),iimo(iter),
     &      iidy(iter),iihr(iter)

c CASE 2: ROC
          if (iter.gt.1) then
            if (tair(iter-1).ne.undef) then
              TEMP_diff = abs(tair(iter) - tair(iter-1))
              if (TEMP_diff.gt.Tair_hrly_maxinc) write (*,103)
     &          TEMP_diff,Tair_hrly_maxinc,iiyr(iter),iimo(iter),
     &          iidy(iter),iihr(iter)
            endif
          endif

c CASE 3: NOC
          if (iter.gt.3) then
            if (tair(iter-1).ne.undef .and. tair(iter-2).ne. undef
     &        .and. tair(iter-3).ne.undef) then
              if (tair(iter).eq.tair(iter-1) .and.
     &          tair(iter).eq.tair(iter-2) .and. 
     &          tair(iter).eq.tair(iter-3)) write (*,104)
     &          tair(iter),iiyr(iter),iimo(iter),iidy(iter),
     &          iihr(iter)
            endif
          endif

        endif
      enddo

 101  format ('TAIR [',f6.1,'] above bounds [',f6.1,'] at',i5,i3,i3,i3)
 102  format ('TAIR [',f6.1,'] below bounds [',f6.1,'] at',i5,i3,i3,i3)
 103  format ('TAIR hourly inc [',f5.1,'] above bounds [',f5.1,
     &  '] at',i5,i3,i3,i3)
 104  format ('Constant TAIR [',f6.1,'] found at',i5,i3,i3,i3)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_rh(RHmin,RHmax,RH_hrly_maxinc,maxiter,rh,undef,
     &  iiyr,iimo,iidy,iihr)

      implicit none

      integer iter,maxiter

      integer iiyr(maxiter)
      integer iimo(maxiter)
      integer iidy(maxiter)
      integer iihr(maxiter)

      real RHmin,RHmax,RH_hrly_maxinc,RH_diff,undef
      real rh(maxiter)

      do iter=1,maxiter

        if (rh(iter).ne.undef) then

c CASE 1: LIM
          if (rh(iter).gt.RHmax) write (*,101)
     &      rh(iter),RHmax,iiyr(iter),iimo(iter),iidy(iter),iihr(iter)
          if (rh(iter).lt.RHmin) write (*,102)
     &      rh(iter),RHmin,iiyr(iter),iimo(iter),iidy(iter),iihr(iter)

c CASE 2: ROC
          if (iter.gt.1) then
            if (rh(iter-1).ne.undef) then
              RH_diff = abs(rh(iter) - rh(iter-1))
              if (RH_diff.gt.RH_hrly_maxinc) write (*,103)
     &          RH_diff,RH_hrly_maxinc,iiyr(iter),iimo(iter),
     &          iidy(iter),iihr(iter)
            endif
          endif

c CASE 3: NOC
          if (iter.gt.3) then
            if (rh(iter-1).ne.undef .and. rh(iter-2).ne. undef
     &        .and. rh(iter-3).ne.undef) then
              if (rh(iter).eq.rh(iter-1) .and.
     &          rh(iter).eq.rh(iter-2) .and. 
     &          rh(iter).eq.rh(iter-3)) write (*,104)
     &          rh(iter),iiyr(iter),iimo(iter),iidy(iter),iihr(iter)
            endif
          endif
        endif

      enddo

 101  format ('RELH [',f6.1,'] above bounds [',f6.1,'] at',i5,i3,i3,i3)
 102  format ('RELH [',f6.1,'] below bounds [',f6.1,'] at',i5,i3,i3,i3)
 103  format ('RELH hourly inc [',f6.1,'] above bounds [',f6.1,
     &  '] at',i5,i3,i3,i3)
 104  format ('Constant RELH [',f6.1,'] found at',i5,i3,i3,i3)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine char2char(cvar,c_value,i_value_chars,
     &  c_param)

      implicit none

      character*(*) c_value,c_param,cvar
      integer i_value_chars

      cvar = c_value(1:i_value_chars)

      print *,c_param,' = ',cvar(1:i_value_chars)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine char2real(rvar,i_value_chars,c_value,
     &  c_param)

      implicit none

      character*(*) c_value,c_param
      integer i_value_chars
      real rvar
      character*8 form

c Read an real value (rvar) from the character string (c_value).
      write (form,90) i_value_chars
   90 format ('(f',i2,'.0)')
      read (c_value,form) rvar

      print *,c_param,' =',rvar

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine char2int(ivar,i_value_chars,c_value,
     &  c_param)

      implicit none

      character*(*) c_value,c_param
      integer i_value_chars,ivar
      character*8 form

c Read an integer value (ivar) from the character string (c_value).
      write (form,90) i_value_chars
   90 format ('(i',i2,')')
      read (c_value,form) ivar

      print *,c_param,' =',ivar

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_param_data(input_string,c_param,c_value,
     &  i_param_chars,i_value_chars,icomment_flag)

      implicit none

      character*(*) input_string,c_param,c_value

      integer leading_blanks,trailing_blanks
      integer i,icomment_flag
      integer i_param_start,i_equals_position,i_param_end,
     &  i_value_start,i_value_end,i_param_chars,i_value_chars,
     &  i_loc,i_leading_blanks,i_trailing_blanks

c If the input string is not a comment line, process the data.
      if (input_string(1:1).ne.'!') then

c First count the number of leading and trailing blanks.
        i_leading_blanks = leading_blanks(input_string)
        i_trailing_blanks = trailing_blanks(input_string)

c If the input string is not completely blank, process the data.
        if (i_leading_blanks.ne.len(input_string)) then
          icomment_flag = 0

c Define the starting and ending points of the parameter name and
c   parameter value.
          i_param_start = i_leading_blanks + 1
          i_equals_position = index(input_string,'=')
          i_param_end = i_equals_position - 2
          i_value_start = i_equals_position + 2
          i_value_end =  len(input_string) - i_trailing_blanks
          i_param_chars = i_param_end - i_param_start + 1
          i_value_chars = i_value_end - i_value_start + 1

c Pull out the parameter name and value.
          do i=1,i_param_chars
            i_loc = i + i_param_start - 1
            c_param(i:i) = input_string(i_loc:i_loc)
          enddo

          do i=1,i_value_chars
            i_loc = i + i_value_start - 1
            c_value(i:i) = input_string(i_loc:i_loc)
          enddo

        else
          icomment_flag = 1
        endif

      else
        icomment_flag = 1
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function leading_blanks(input_string)

      implicit none

      integer k
      character*(*) input_string

c Count the number of blanks preceeding the first non-blank
c   character.
      leading_blanks = 0
      do k=1,len(input_string)
        if (input_string(k:k).eq.' ') then
          leading_blanks = leading_blanks + 1
        else
          return
        endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function trailing_blanks(input_string)

      implicit none

      integer k
      character*(*) input_string

c Count the number of blanks following the last non-blank
c   character.
      trailing_blanks = 0
      do k=len(input_string),1,-1
        if (input_string(k:k).eq.' ') then
          trailing_blanks = trailing_blanks + 1
        else
          return
        endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c file = calendar.f  version 1.0
c
c Note from Glen: This is the version that should be used as
c   part of my computer programs!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c This program performs various date conversion calculations
c using subroutine calndr().
c On a 32-bit computer, it can handle any date between roughly
c 5 million BC and 5 million AD.  This limitation is due to
c the range of integers that can be expressed with 32 bits.
c The algorithm has no limitation.
c
c Using function idaywk(), the day of the week is computed
c along with the answer to the user's calendar calculation.
c
c External routines called:
c calndr  calendar conversions
c idaywk  day of the week determination
c
c Portability
c This routine is coded to Fortran 77 standards except that
c lower case is used.
c
c Copyright (C) 1999 Jon Ahlquist.
c Issued under the second GNU General Public License.
c See www.gnu.org for details.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c If you find any errors, please notify:
c Jon Ahlquist <ahlquist@met.fsu.edu>   
c Dept of Meteorology
c Florida State University
c Tallahassee, FL 32306-4520
c 15 March 1999.
c
c----------
c
c Declare variables.
c     implicit     none
c     integer      day, daynum, ioptn, julian, month, year
c     character*9  daynam(0:6)
c
c Declare the integer function used to compute the day of the week.
c     integer  idaywk
c
c Define the day names.
c     data  daynam /'Sunday',   'Monday', 'Tuesday', 'Wednesday',
c    &              'Thursday', 'Friday', 'Saturday'/
c
c Variables and their meanings
c day     day of the month.
c daynam  array of day names.  (daynam(0)='Sunday', daynam(1)='Monday',
c            ..., daynam(6)='Saturday')
c daynum  day number during the year.  (1 for 1 January, 2 for
c            2 January, 32 for 1 February, etc.)
c idaywk  integer function that returns an integer counter indicating
c            the day of the week, where 0 refers to Sunday, 1 to Monday,
c            up to 6 for Saturday.
c ioptn   option indicator where 0 < abs(ioptn) < 6.
c            See below and especially subroutine calndr for details.
c julian  Julian Day number.
c month   month counter (1=January, 2=February, ..., 12=December)
c year    year expressed with ALL digits.  DO NOT abbreviate years
c            by using only the last two digits.
c
c----------
c
      subroutine calndr (ioptn, iday, month, iyear, idayct)
c
c----------
c
c CALNDR = CALeNDaR conversions, version 1.0
c
c Input variable specifying the desired calendar conversion option.
      integer ioptn
c
c Input/Output variables (sometimes input, sometimes output,
c depending on the value of the desired option, ioptn.)
      integer  iday, month, iyear, idayct
c
c----------
c
c Subroutine calndr() performs calendar calculations using either
c the standard Gregorian calendar or the old Julian calendar.
c This subroutine extends the definitions of these calendar systems
c to any arbitrary year.  The algorithms in this subroutine
c will work with any date in the past or future,
c but overflows will occur if the numbers are sufficiently large.
c For a computer using a 32-bit integer, this routine can handle
c any date between roughly 5.8 million BC and 5.8 million AD
c without experiencing overflow during calculations.
c
c No external functions or subroutines are called.
c
c----------
c
c INPUT/OUTPUT ARGUMENTS FOR SUBROUTINE CALNDR()
c
c "ioptn" is the desired calendar conversion option explained below.
c Positive option values use the standard modern Gregorian calendar.
c Negative option values use the old Julian calendar which was the
c standard in Europe from its institution by Julius Caesar in 45 BC
c until at least 4 October 1582.  The Gregorian and Julian calendars
c are explained further below.
c
c (iday,month,iyear) is a calendar date where "iday" is the day of
c the month, "month" is 1 for January, 2 for February, etc.,
c and "iyear" is the year.  If the year is 1968 AD, enter iyear=1968,
c since iyear=68 would refer to 68 AD.
c For BC years, iyear should be negative, so 45 BC would be iyear=-45.
c By convention, there is no year 0 under the BC/AD year numbering
c scheme.  That is, years proceed as 2 BC, 1 BC, 1 AD, 2 AD, etc.,
c without including 0.  Subroutine calndr() will print an error message
c and stop if you specify iyear=0.
c
c "idayct" is a day count.  It is either the day number during the
c specified year or the Julian Day number, depending on the value
c of ioptn.  By day number during the specified year, we mean
c idayct=1 on 1 January, idayct=32 on 1 February, etc., to idayct=365
c or 366 on 31 December, depending on whether the specified year
c is a leap year.
c
c The values of input variables are not changed by this subroutine.
c
c
c ALLOWABLE VALUES FOR "IOPTN" and the conversions they invoke.
c Positive option values ( 1 to  5) use the standard Gregorian calendar.
c Negative option values (-1 to -5) use the old      Julian    calendar.
c
c Absolute
c  value
c of ioptn   Input variable(s)     Output variable(s)
c
c    1       iday,month,iyear      idayct
c Given a calendar date (iday,month,iyear), compute the day number
c (idayct) during the year, where 1 January is day number 1 and
c 31 December is day number 365 or 366, depending on whether it is
c a leap year.
c
c    2       idayct,iyear          iday,month
c Given the day number of the year (idayct) and the year (iyear),
c compute the day of the month (iday) and the month (month).
c
c    3       iday,month,iyear      idayct
c Given a calendar date (iday,month,iyear), compute the Julian Day
c number (idayct) that starts at noon of the calendar date specified.
c
c    4       idayct                iday,month,iyear
c Given the Julian Day number (idayct) that starts at noon,
c compute the corresponding calendar date (iday,month,iyear).
c
c    5       idayct                iday,month,iyear
c Given the Julian Day number (idayct) that starts at noon,
c compute the corresponding day number for the year (iday)
c and year (iyear).  On return from calndr(), "month" will always
c be set equal to 1 when ioptn=5.
c
c No inverse function is needed for ioptn=5 because it is
c available through option 3.  One simply calls calndr() with:
c ioptn = 3,
c iday  = day number of the year instead of day of the month,
c month = 1, and
c iyear = whatever the desired year is.
c
c----------
c
c EXAMPLES
c The first 6 examples are for the standard Gregorian calendar.
c All the examples deal with 15 October 1582, which was the first day
c of the Gregorian calendar.  15 October is the 288-th day of the year.
c Julian Day number 2299161 began at noon on 15 October 1582.
c
c Find the day number during the year on 15 October 1582
c     ioptn = 1
c     call calndr (ioptn, 15, 10, 1582,  idayct)
c calndr() should return idayct=288
c
c Find the day of the month and month for day 288 in year 1582.
c     ioptn = 2
c     call calndr (ioptn, iday, month, 1582, 288)
c calndr() should return iday=15 and month=10.
c
c Find the Julian Day number for 15 October 1582.
c     ioptn = 3
c     call calndr (ioptn, 15, 10, 1582, julian)
c calndr() should return julian=2299161
c
c Find the Julian Day number for day 288 during 1582 AD.
c When the input is day number of the year, one should specify month=1
c     ioptn = 3
c     call calndr (ioptn, 288, 1, 1582, julian)
c calndr() should return dayct=2299161
c
c Find the date for Julian Day number 2299161.
c     ioptn = 4
c     call calndr (ioptn, iday, month, iyear, 2299161)
c calndr() should return iday=15, month=10, and iyear=1582
c 
c Find the day number during the year (iday) and year
c for Julian Day number 2299161.
c     ioptn = 5
c     call calndr (ioptn, iday, month, iyear, 2299161)
c calndr() should return iday=288, month=1, iyear=1582
c
c Given 15 October 1582 under the Gregorian calendar,
c find the date (idayJ,imonthJ,iyearJ) under the Julian calendar.
c To do this, we call calndr() twice, using the Julian Day number
c as the intermediate value.
c     call calndr ( 3, 15,        10, 1582,    julian)
c     call calndr (-4, idayJ, monthJ, iyearJ,  julian)
c The first call to calndr() should return julian=2299161, and
c the second should return idayJ=5, monthJ=10, iyearJ=1582
c
c----------
c
c BASIC CALENDAR INFORMATION
c
c The Julian calendar was instituted by Julius Caesar in 45 BC.
c Every fourth year is a leap year in which February has 29 days.
c That is, the Julian calendar assumes that the year is exactly
c 365.25 days long.  Actually, the year is not quite this long.
c The modern Gregorian calendar remedies this by omitting leap years
c in years divisible by 100 except when the year is divisible by 400.
c Thus, 1700, 1800, and 1900 are leap years under the Julian calendar
c but not under the Gregorian calendar.  The years 1600 and 2000 are
c leap years under both the Julian and the Gregorian calendars.
c Other years divisible by 4 are leap years under both calendars,
c such as 1992, 1996, 2004, 2008, 2012, etc.  For BC years, we recall
c that year 0 was omitted, so 1 BC, 5 BC, 9 BC, 13 BC, etc., and 401 BC,
c 801 BC, 1201 BC, etc., are leap years under both calendars, while
c 101 BC, 201 BC, 301 BC, 501 BC, 601 BC, 701 BC, 901 BC, 1001 BC,
c 1101 BC, etc., are leap years under the Julian calendar but not
c the Gregorian calendar.
c
c The Gregorian calendar is named after Pope Gregory XIII.  He declared
c that the last day of the old Julian calendar would be Thursday,
c 4 October 1582 and that the following day, Friday, would be reckoned
c under the new calendar as 15 October 1582.  The jump of 10 days was
c included to make 21 March closer to the spring equinox.
c
c Only a few Catholic countries (Italy, Poland, Portugal, and Spain)
c switched to the Gregorian calendar on the day after 4 October 1582.
c It took other countries months to centuries to change to the
c Gregorian calendar.  For example, England's first day under the
c Gregorian calendar was 14 September 1752.  The same date applied to
c the entire British empire, including America.  Japan, Russia, and many
c eastern European countries did not change to the Gregorian calendar
c until the 20th century.  The last country to change was Turkey,
c which began using the Gregorian calendar on 1 January 1927.
c
c Therefore, between the years 1582 and 1926 AD, you must know
c the country in which an event was dated to interpret the date
c correctly.  In Sweden, there was even a year (1712) when February
c had 30 days.  Consult a book on calendars for more details
c about when various countries changed their calendars.
c
c DAY NUMBER DURING THE YEAR
c The day number during the year is simply a counter equal to 1 on
c 1 January, 32 on 1 February, etc., thorugh 365 or 366 on 31 December,
c depending on whether the year is a leap year.  Sometimes this is
c called the Julian Day, but that term is better reserved for the
c day counter explained below.
c
c JULIAN DAY NUMBER
c The Julian Day numbering system was designed by Joseph Scaliger
c in 1582 to remove ambiguity caused by varying calendar systems.
c The name "Julian Day" was chosen to honor Scaliger's father,
c Julius Caesar Scaliger (1484-1558), an Italian scholar and physician
c who lived in France.  Because Julian Day numbering was especially
c designed for astronomers, Julian Days begin at noon so that the day
c counter does not change in the middle of an astronmer's observing
c period.  Julian Day 0 began at noon on 1 January 4713 BC under the
c Julian calendar.  A modern reference point is that 23 May 1968
c (Gregorian calendar) was Julian Day 2,440,000.
c
c JULIAN DAY NUMBER EXAMPLES
c
c The table below shows a few Julian Day numbers and their corresponding
c dates, depending on which calendar is used.  A negative 'iyear' refers
c to BC (Before Christ).
c
c                     Julian Day under calendar:
c iday  month   iyear     Gregorian   Julian
c  24     11   -4714            0        -38
c   1      1   -4713           38          0
c   1      1       1      1721426    1721424
c   4     10    1582      2299150    2299160
c  15     10    1582      2299161    2299171
c   1      3    1600      2305508    2305518
c  23      5    1968      2440000    2440013
c   5      7    1998      2451000    2451013
c   1      3    2000      2451605    2451618
c   1      1    2001      2451911    2451924
c
c From this table, we can see that the 10 day difference between the
c two calendars in 1582 grew to 13 days by 1 March 1900, since 1900 was
c a leap year under the Julian calendar but not under the Gregorian
c calendar.  The gap will widen to 14 days after 1 March 2100 for the
c same reason.
c 
c----------
c
c PORTABILITY
c
c This subroutine is written in standard FORTRAN 77.
c It calls no external functions or subroutines and should run
c without problem on any computer having a 32-bit word or longer.
c 
c----------
c
c ALGORITHM
c
c The goal in coding calndr() was clear, clean code, not efficiency.
c Calendar calculations usually take a trivial fraction of the time
c in any program in which dates conversions are involved.
c Data analysis usually takes the most time.
c
c Standard algorithms are followed in this subroutine.  Internal to
c this subroutine, we use a year counter "jyear" such that
c  jyear=iyear   when iyear is positive
c       =iyear+1 when iyear is negative.
c Thus, jyear does not experience a 1 year jump like iyear does
c when going from BC to AD.  Specifically, jyear=0 when iyear=-1,
c i.e., when the year is 1 BC.
c
c For simplicity in dealing with February, inside this subroutine,
c we let the year begin on 1 March so that the adjustable month,
c February is the last month of the year.
c It is clear that the calendar used to work this way because the
c months September, October, November, and December refer to
c 7, 8, 9, and 10.  For consistency, jyear is incremented on 1 March
c rather than on 1 January.  Of course, everything is adjusted back to
c standard practice of years beginning on 1 January before answers
c are returned to the routine that calls calndr().
c
c Lastly, we use a trick to calculate the number of days from 1 March
c until the end of the month that precedes the specified month.
c That number of days is int(30.6001*(month+1))-122,
c where 30.6001 is used to avoid the possibility of round-off and
c truncation error.  For example, if 30.6 were used instead,
c 30.6*5 should be 153, but round-off error could make it 152.99999,
c which would then truncated to 152, causing an error of 1 day.
c
c Algorithm reference:
c Dershowitz, Nachum and Edward M. Reingold, 1990: Calendrical
c Calculations.  Software-Practice and Experience, vol. 20, number 9
c (September 1990), pp. 899-928.
c
c Copyright (C) 1999 Jon Ahlquist.
c Issued under the second GNU General Public License.
c See www.gnu.org for details.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c If you find any errors, please notify:
c Jon Ahlquist <ahlquist@met.fsu.edu>   
c Dept of Meteorology
c Florida State University
c Tallahassee, FL 32306-4520
c 15 March 1999.
c
c-----
c Declare internal variables.
      integer  jdref,  jmonth, jyear, leap,
     &         n1yr, n4yr, n100yr, n400yr,
     &         ndays, ndy400, ndy100, nyrs,
     &         yr400, yrref
c
c Explanation of all internal variables.
c jdref   Julian Day on which 1 March begins in the reference year.
c jmonth  Month counter which equals month+1 if month .gt. 2
c          or month+13 if month .le. 2.
c jyear   Year index,  jyear=iyear if iyear .gt. 0, jyear=iyear+1
c            if iyear .lt. 0.  Thus, jyear does not skip year 0
c            like iyear does between BC and AD years.
c leap    =1 if the year is a leap year, =0 if not.
c n1yr    Number of complete individual years between iyear and
c            the reference year after all 4, 100,
c            and 400 year periods have been removed.
c n4yr    Number of complete 4 year cycles between iyear and
c            the reference year after all 100 and 400 year periods
c            have been removed.
c n100yr  Number of complete 100 year periods between iyear and
c            the reference year after all 400 year periods
c            have been removed.
c n400yr  Number of complete 400 year periods between iyear and
c            the reference year.
c ndays   Number of days since 1 March during iyear.  (In intermediate
c            steps, it holds other day counts as well.)
c ndy400  Number of days in 400 years.  Under the Gregorian calendar,
c            this is 400*365 + 100 - 3 = 146097.  Under the Julian
c            calendar, this is 400*365 + 100 = 146100.
c ndy100  Number of days in 100 years,  Under the Gregorian calendar,
c            this is 100*365 + 24 = 36524.   Under the Julian calendar,
c            this is 100*365 + 25 = 36525.
c nyrs    Number of years from the beginning of yr400
c              to the beginning of jyear.  (Used for option +/-3).
c yr400   The largest multiple of 400 years that is .le. jyear.
c
c
c----------------------------------------------------------------
c Do preparation work.
c
c Look for out-of-range option values.
      if ((ioptn .eq. 0) .or. (abs(ioptn) .ge. 6)) then
         write(*,*)'For calndr(), you specified ioptn = ', ioptn
         write(*,*)
     &   'Allowable values are 1 to 5 for the Gregorian calendar'
         write(*,*)
     &   'and -1 to -5 for the Julian calendar.'
         stop
      endif
c
c Options 1-3 have "iyear" as an input value.
c Internally, we use variable "jyear" that does not have a jump
c from -1 (for 1 BC) to +1 (for 1 AD).
      if (abs(ioptn) .le. 3) then
         if (iyear .gt. 0) then
            jyear = iyear
         elseif (iyear .eq. 0) then
            write(*,*)
     &      'For calndr(), you specified the nonexistent year 0'
            stop
         else
            jyear = iyear + 1
         endif
c
c        Set "leap" equal to 0 if "jyear" is not a leap year
c        and equal to 1 if it is a leap year.
         leap = 0
         if ((jyear/4)*4 .eq. jyear) then
            leap = 1
         endif
         if ((ioptn .gt. 0)               .and.
     &       ((jyear/100)*100 .eq. jyear) .and.
     &       ((jyear/400)*400 .ne. jyear)      ) then
               leap = 0
         endif
      endif
c
c Options 3-5 involve Julian Day numbers, which need a reference year
c and the Julian Days that began at noon on 1 March of the reference
c year under the Gregorian and Julian calendars.  Any year for which
c "jyear" is divisible by 400 can be used as a reference year.
c We chose 1600 AD as the reference year because it is the closest
c multiple of 400 to the institution of the Gregorian calendar, making
c it relatively easy to compute the Julian Day for 1 March 1600
c given that, on 15 October 1582 under the Gregorian calendar,
c the Julian Day was 2299161.  Similarly, we need to do the same
c calculation for the Julian calendar.  We can compute this Julian
c Day knwoing that on 4 October 1582 under the Julian calendar,
c the Julian Day number was 2299160.  The details of these calculations
c is next. 
c    From 15 October until 1 March, the number of days is the remainder
c of October plus the days in November, December, January, and February:
c 17+30+31+31+28 = 137, so 1 March 1583 under the Gregorian calendar
c was Julian Day 2,299,298.  Because of the 10 day jump ahead at the
c switch from the Julian calendar to the Gregorian calendar, 1 March
c 1583 under the Julian calendar was Julian Day 2,299,308.  Making use
c of the rules for the two calendar systems, 1 March 1600 was Julian
c Day 2,299,298 + (1600-1583)*365 + 5 (due to leap years) =
c 2,305,508 under the Gregorian calendar and day 2,305,518 under the
c Julian calendar.
c    We also set the number of days in 400 years and 100 years.
c For reference, 400 years is 146097 days under the Gregorian calendar
c and 146100 days under the Julian calendar.  100 years is 36524 days
c under the Gregorian calendar and 36525 days under the Julian calendar.
      if (abs(ioptn) .ge. 3) then
c
c        Julian calendar values.
         yrref  =    1600
         jdref  = 2305518
c               = Julian Day reference value for the day that begins
c                 at noon on 1 March of the reference year "yrref".
         ndy400 = 400*365 + 100
         ndy100 = 100*365 +  25
c
c        Adjust for Gregorian calendar values.
         if (ioptn .gt. 0) then
            jdref  = jdref  - 10
            ndy400 = ndy400 -  3
            ndy100 = ndy100 -  1
         endif
      endif
c
c----------------------------------------------------------------
c OPTIONS -1 and +1:
c Given a calendar date (iday,month,iyear), compute the day number
c of the year (idayct), where 1 January is day number 1 and 31 December
c is day number 365 or 366, depending on whether it is a leap year.
      if (abs(ioptn) .eq. 1) then
c
c     Compute the day number during the year.
      if (month .le. 2) then
         idayct = iday + (month-1)*31
      else
         idayct = iday + int(30.6001 * (month+1)) - 63 + leap
      endif
c
c----------------------------------------------------------------
c OPTIONS -2 and +2:
c Given the day number of the year (idayct) and the year (iyear),
c compute the day of the month (iday) and the month (month).
      elseif (abs(ioptn) .eq. 2) then
c
      if (idayct .lt. 60+leap) then
         month  = (idayct-1)/31
         iday   = idayct - month*31
         month  = month + 1
      else
         ndays  = idayct - (60+leap)
c               = number of days past 1 March of the current year.
         jmonth = (10*(ndays+31))/306 + 3
c               = month counter, =4 for March, =5 for April, etc.
         iday   = (ndays+123) - int(30.6001*jmonth) 
         month  = jmonth - 1
      endif
c
c----------------------------------------------------------------
c OPTIONS -3 and +3:
c Given a calendar date (iday,month,iyear), compute the Julian Day
c number (idayct) that starts at noon.
      elseif (abs(ioptn) .eq. 3) then
c
c     Shift to a system where the year starts on 1 March, so January
c     and February belong to the preceding year.
c     Define jmonth=4 for March, =5 for April, ..., =15 for February.
      if (month .le. 2) then
        jyear  = jyear -  1
        jmonth = month + 13
      else
        jmonth = month +  1
      endif
c
c     Find the closest multiple of 400 years that is .le. jyear.
      yr400 = (jyear/400)*400
c           = multiple of 400 years at or less than jyear.
      if (jyear .lt. yr400) then
         yr400 = yr400 - 400
      endif
c
      n400yr = (yr400 - yrref)/400
c            = number of 400-year periods from yrref to yr400.
      nyrs   = jyear - yr400
c            = number of years from the beginning of yr400
c              to the beginning of jyear.
c
c     Compute the Julian Day number.
      idayct = iday + int(30.6001*jmonth) - 123 + 365*nyrs + nyrs/4
     &       + jdref + n400yr*ndy400
c
c     If we are using the Gregorian calendar, we must not count
c     every 100-th year as a leap year.  nyrs is less than 400 years,
c     so we do not need to consider the leap year that would occur if
c     nyrs were divisible by 400, i.e., we do not add nyrs/400.
      if (ioptn .gt. 0) then
         idayct = idayct - nyrs/100
      endif
c
c----------------------------------------------------------------
c OPTIONS -5, -4, +4, and +5:
c Given the Julian Day number (idayct) that starts at noon,
c compute the corresponding calendar date (iday,month,iyear)
c (abs(ioptn)=4) or day number during the year (abs(ioptn)=5).
      else
c
c     Create a new reference date which begins on the nearest
c     400-year cycle less than or equal to the Julian Day for 1 March
c     in the year in which the given Julian Day number (idayct) occurs.
      ndays  = idayct - jdref
      n400yr = ndays / ndy400
c            = integral number of 400-year periods separating
c              idayct and the reference date, jdref.
      jdref  = jdref + n400yr*ndy400
      if (jdref .gt. idayct) then
         n400yr = n400yr - 1
         jdref  = jdref  - ndy400
      endif
c
      ndays  = idayct - jdref
c            = number from the reference date to idayct.
c
      n100yr = min(ndays/ndy100, 3)
c            = number of complete 100-year periods
c              from the reference year to the current year.
c              The min() function is necessary to avoid n100yr=4
c              on 29 February of the last year in the 400-year cycle.
c
      ndays  = ndays - n100yr*ndy100
c            = remainder after removing an integral number of
c              100-year periods.
c
      n4yr   = ndays / 1461
c            = number of complete 4-year periods in the current century.
c              4 years consists of 4*365 + 1 = 1461 days.
c
      ndays  = ndays - n4yr*1461
c            = remainder after removing an integral number
c              of 4-year periods.
c
      n1yr   = min(ndays/365, 3)
c            = number of complete years since the last leap year.
c              The min() function is necessary to avoid n1yr=4
c              when the date is 29 February on a leap year,
c              in which case ndays=1460, and 1460/365 = 4.
c
      ndays  = ndays - 365*n1yr
c            = number of days so far in the current year,
c              where ndays=0 on 1 March.
c
      iyear  = n1yr + 4*n4yr + 100*n100yr + 400*n400yr + yrref 
c            = year, as counted in the standard way,
c              but relative to 1 March.
c
c At this point, we need to separate ioptn=abs(4), which seeks a
c calendar date, and ioptn=abs(5), which seeks the day number during
c the year.  First compute the calendar date if desired (abs(ioptn)=4).
      if (abs(ioptn) .eq. 4) then
         jmonth = (10*(ndays+31))/306 + 3
c               = offset month counter.  jmonth=4 for March, =13 for
c                 December, =14 for January, =15 for February.
         iday   = (ndays+123) - int(30.6001*jmonth)
c               = day of the month, starting with 1 on the first day
c                 of the month.
c
c        Now adjust for the fact that the year actually begins
c        on 1 January.
         if (jmonth .le. 13) then
            month = jmonth - 1
         else
            month = jmonth - 13
            iyear = iyear + 1
         endif
c
c This code handles abs(ioptn)=5, finding the day number during the year.
      else
c        ioptn=5 always returns month=1, which we set now.
         month = 1
c
c        We need to determine whether this is a leap year.
         leap = 0
         if ((jyear/4)*4 .eq. jyear) then
            leap = 1
         endif
         if ((ioptn .gt. 0)               .and.
     &       ((jyear/100)*100 .eq. jyear) .and.
     &       ((jyear/400)*400 .ne. jyear)      ) then
               leap = 0
         endif
c
c        Now find the day number "iday".
c        ndays is the number of days since the most recent 1 March,
c        so ndays=0 on 1 March.
         if (ndays .le.305) then
            iday  = ndays + 60 + leap
         else
            iday  = ndays - 305
            iyear = iyear + 1
         endif
      endif
c
c     Adjust the year if it is .le. 0, and hence BC (Before Christ).
      if (iyear .le. 0) then
         iyear = iyear - 1
      endif
c
c End the code for the last option, ioptn.
      endif
c
      return
      end


      integer function idaywk(jdayno)
c
c IDAYWK = compute the DAY of the WeeK given the Julian Day number,
c          version 1.0.
c
c Input variable
      integer  jdayno
c jdayno = Julian Day number starting at noon of the day in question.
c
c Output variable:
c idaywk = day of the week, where 0=Sunday, 1=Monday, ..., 6=Saturday.
c
c----------
c Compute the day of the week given the Julian Day number.
c You can find the Julian Day number given (day,month,year)
c using subroutine calndr.f.
c Example: For the first day of the Gregorian calendar,
c 15 October 1582, compute the Julian day number (option 3 of
c subroutine calndr) and compute the day of the week.
c     call calndr (3, 15, 10, 1582, jdayno) 
c     write(*,*) jdayno, idaywk(jdayno)
c The numbers printed should be 2299161 and 5,
c where 6 refers to Friday.
c
c Copyright (C) 1999 Jon Ahlquist.
c Issued under the second GNU General Public License.
c See www.gnu.org for details.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c If you find any errors, please notify:
c Jon Ahlquist <ahlquist@met.fsu.edu>   
c Dept of Meteorology
c Florida State University
c Tallahassee, FL 32306-4520
c 15 March 1999.
c
c-----
c Declare internal variable.
c jdSun is the Julian Day number starting at noon on any Sunday.
c I arbitrarily chose the first Sunday after Julian Day 1,
c which is Julian Day 6.
      integer  jdSun
	data     jdSun /6/
      idaywk = mod(jdayno-jdSun,7)
c If jdayno-jdSun < 0, then we are taking the modulus of a negative
c number. Fortran's built-in mod function returns a negative value
c when the argument is negative.  In that case, we adjust the result
c to a positive value.
      if (idaywk .lt. 0) idaywk = idaywk + 7
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

