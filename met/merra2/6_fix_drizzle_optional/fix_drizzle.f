c fix_drizzle.f

c This program reads the .gdat preciptation archive for this
c   SnowModel run and determines what the "clipping" level
c   needs to be to eliminate (re)analysis "drizzle".  This
c   code also puts the drizzle volume that was removed back in
c   the remaining precipitation inputs so the (re)analysis
c   precipitation volume is conserved.  It also does the
c   analysis year-by-year (each year is different).

c I want to just keep the main snowfall 'storms'.  So, I don't
c   want to just clip the low precipitation values, because
c   there may be a 'storm' that has low snowfall in one corner
c   of the simulation domain; I want to keep that information.
c   So, here I am creating a 'storm' time series mask that can
c   be used to clip and scale entire 'storm' events, and toss
c   everything else.

c Also note that this likely won't work for large domains
c   where one part of the domain has a different, or no storm
c   than the storm in another part of the domain.  It is not
c   clear to me how to deal with that situation.  For now I
c   will assume that the multi-layer model will not be used
c   for continental-scale domains.

      implicit none

      integer max_nx,max_ny,max_maxiter,max_clips,max_maxiter_yr

c This provides enough array space to read and process the
c   reanalysis spatial arrays.  These array dimensions must be
c   equal or greater than the reanalysis dimensions.
      parameter (max_nx=1000,max_ny=1000)

c This provides enough time array space.  300,000 is enough space
c   for a 100-year simulation at a 3hrly time step.  The 366*8 is
c   provides space for 1 year at the 3hrly time step.
      parameter (max_maxiter=300000)
      parameter (max_maxiter_yr=366*8)

c This is the number of clipping tests that are going to be
c   performed.
      parameter (max_clips=100)

      integer nx,ny,i,maxiter,iter,k_clip_index,irecs_in_day,
     &  nlayers_target,nyears

      real spre_clip(max_maxiter_yr,max_clips)
      real clip_mask(max_maxiter_yr)

      real clip_increment,clip_amount(max_clips),ratio(max_clips)

      character*3 time_increment

      integer  iyr_start,imo_start,idy_start,nyear, irec_start,
     &  ndays,irec,iimo_start,iidy_start,iimo_end,iidy_end,idays,
     &  iday_start

      real clip_amount_all(max_clips,max_maxiter)
      real ratio_all(max_clips,max_maxiter)
      real clip_mask_all(max_maxiter)
      real spre_clip_all(max_maxiter,max_clips)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Find the clipping amount that produces approximately 12 snow
c   accumulation layers (or the closest value to 12 layers), where
c   this '12' layers was defined by nlayers_target.
      nlayers_target = 12

c Define the 'winter' period when you want nlayers_target to be
c   calculated for.  Here you provide the months and days that
c   define the start and end of your 'winter' period.  The reason
c   this is an option is so that summer snowfalls are not included
c   in your winter snowpack-layer/snowstorm count.  If you want
c   this period to span the entire simulation year, define these
c   dates to be the start and end of the simulation year.  Note
c   that in the code I have assumed that the 'end year' is the
c   'start year + 1'; if the start and end are the same year, you
c   will have to modify the code to deal with this case.  Also
c   note that how the code is now written, if you get snow storms
c   and snow layers outside of this time window you have defined,
c   they will produce layers in addition to the nlayers_target
c   you defined above.
      iimo_start = 10
      iidy_start = 1
      iimo_end = 5
      iidy_end = 31

c The clipping amount is the clip_increment (in mm/day) times
c   k_clipped (see the subroutine for additional details).
      clip_increment = 0.5

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Read in the domain space and time dimensions for this run.
      open (21,file='../5_mk_gdat_optional/sm_met_inputs.ctl')

      read (21,*)
      read (21,*)
      read (21,*)
      read (21,91) nx
      read (21,91) ny
      read (21,*)
      read (21,91) maxiter,time_increment

      print *
      print *,nx,ny,maxiter
      print *

   91 format (5x,i8,22x,a3)

c Read in the start date for this run.
      open (22,file=
     &  '../4_maxiter_offset/start_end_dates_maxiter_ioffset.dat')
      read (22,92) iyr_start
      read (22,92) imo_start
      read (22,92) idy_start

   92 format (14x,i10)

c Convert the time_increment to an 'irecs_in_day' value, and
c   calculate the number of years you want to process (the number
c   of years in the simulation).  Here I have assumed that this
c   simulation runs for full year(s), not some fraction of a
c   year.
      if (time_increment.eq.'1dy') then
        irecs_in_day = 1
        nyears = int(real(maxiter) / (real(irecs_in_day) * 364.0))
      elseif (time_increment.eq.'3hr') then
        irecs_in_day = 8
        nyears = int(real(maxiter) / (real(irecs_in_day) * 364.0))
      endif

      print *
      print *,'nyears =',nyears
      print *

c Adjust the clipping increment for 3-hour time steps.
      clip_increment = clip_increment / real(irecs_in_day)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Open a general information text file.
      open (51,file='drizzle_info.dat')

c Open the (re)analysis input file.
      open (31,file='../5_mk_gdat_optional/sm_met_inputs.gdat',
     &  form='unformatted',access='direct',recl=4*5*nx*ny)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Loop through the years in the simulation.
      do nyear=1,nyears

        print *,'simulation year =',nyear
        print *,'simulation year =',nyear
        print *,'simulation year =',nyear
        print *

c Calculate the starting record and number of records in this year.
        call get_irecs (iyr_start,imo_start,idy_start,nyear,
     &    irec_start,ndays,irecs_in_day,iimo_start,iidy_start,
     &    iimo_end,iidy_end,idays,iday_start)

        call clip_drizzle (max_nx,max_ny,max_maxiter_yr,max_clips,
     &    nx,ny,ndays,irec_start,spre_clip,clip_amount,clip_mask,
     &    nlayers_target,clip_increment,k_clip_index,ratio,
     &    irecs_in_day,idays,iday_start)

c Build the year-specific arrays.
        clip_amount_all(k_clip_index,nyear) = clip_amount(k_clip_index)
        ratio_all(k_clip_index,nyear) = ratio(k_clip_index)

c Save the general information to a simple text file.
        write (51,93) nyear,k_clip_index,
     &    clip_amount_all(k_clip_index,nyear),
     &    ratio_all(k_clip_index,nyear)

c Build the entire time series (stitch the years together).
        do iter=1,ndays*irecs_in_day
          irec = iter + irec_start - 1
          clip_mask_all(irec) = clip_mask(iter)
          do i=1,max_clips
            spre_clip_all(irec,i) = spre_clip(iter,i)
          enddo
        enddo

      enddo

   93 format (2i8,2f12.6)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Save the clipping results.  Here k_clipped is going into the "x"
c   position in the .gdat file.
      open (41,file='sprec_clipped.gdat',
     &  form='unformatted',access='direct',recl=4*max_clips,
     &  status='replace')

      do iter=1,maxiter
        write (41,rec=iter) (spre_clip_all(iter,i),i=1,max_clips)
      enddo

c Create the GrADS .ctl (control) file to go with the GrADS
c   .gdat output file that was generated.
      call mk_spre_clipped_ctl(max_clips,iyr_start,imo_start,
     &  idy_start,time_increment,maxiter)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Save the information required to build the precipitation-clipped
c   MicroMet file.
      open (61,file='drizzle_storm_clip.gdat',
     &  form='unformatted',access='direct',recl=4*1,
     &  status='replace')

      do iter=1,maxiter
        write (61,rec=iter) clip_mask_all(iter)
      enddo

c Create the GrADS .ctl (control) file to go with the GrADS
c   .gdat output file that was generated.
      call mk_clip_mask_ctl(iyr_start,imo_start,idy_start,
     &  time_increment,maxiter)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine clip_drizzle (max_nx,max_ny,max_maxiter_yr,max_clips,
     &  nx,ny,ndays,irec_start,spre_clip,clip_amount,clip_mask,
     &  nlayers_target,clip_increment,k_clip_index,ratio,
     &  irecs_in_day,idays,iday_start)

      implicit none

      integer max_nx,max_ny,max_maxiter_yr,max_clips,ndays,irec_start

      real tair(max_nx,max_ny),relh(max_nx,max_ny),
     &  wspd(max_nx,max_ny),wdir(max_nx,max_ny),prec(max_nx,max_ny)

      real spre(max_nx,max_ny)

      integer nx,ny,i,j,iter,k_clipped,k_clip,k_clip_index,
     &  nlayers_target,irec,irecs_in_day,idays,iday_start

      real spre_clip(max_maxiter_yr,max_clips)
      real spre_ave(max_maxiter_yr),clip_mask(max_maxiter_yr)
      real spre_sum(max_clips)

      real clip_increment,clip_amount(max_clips),ratio(max_clips)
      integer icount(max_clips),idiff(max_clips)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iter=1,ndays*irecs_in_day

        irec = iter + irec_start - 1

        read (31,rec=irec)
     &    ((tair(i,j),i=1,nx),j=1,ny),
     &    ((relh(i,j),i=1,nx),j=1,ny),
     &    ((wspd(i,j),i=1,nx),j=1,ny),
     &    ((wdir(i,j),i=1,nx),j=1,ny),
     &    ((prec(i,j),i=1,nx),j=1,ny)

c Create a simple snowfall array.
        do j=1,ny
          do i=1,nx
            if (tair(i,j).le.2.0) then
              spre(i,j) = prec(i,j)
            else
              spre(i,j) = 0.0
            endif
          enddo
        enddo

c Create a domain-averaged snowfall time series.
        spre_ave(iter) = 0.0
        do j=1,ny
          do i=1,nx
            spre_ave(iter) = spre_ave(iter) + spre(i,j)
          enddo
        enddo
        spre_ave(iter) = spre_ave(iter) / real(nx * ny)

      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Gradually increase the drizzle clipping level and count the
c   resulting number of continuous (in time) snowfall precipitation 
c   events in a year.

c The clipping amount is the clip_increment (in mm) times k_clipped.
c   Note that the first "clip" (k_clipped=1) does no clipping.
      do k_clipped=1,max_clips
        clip_amount(k_clipped) = clip_increment * real(k_clipped - 1)
        do iter=1,ndays*irecs_in_day
          if (spre_ave(iter).ge.clip_amount(k_clipped)) then
            spre_clip(iter,k_clipped) = spre_ave(iter)
          else
            spre_clip(iter,k_clipped) = 0.0
          endif
        enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Count the number of time-continuous precipitation events there
c   are per year in this time series.
      do k_clipped=1,max_clips
        icount(k_clipped) = 0
c This processes all days in the year.
c       do iter=1,ndays*irecs_in_day-1
c This processes the 'winter' days in the year.
        do iter=iday_start,idays*irecs_in_day-1
          if (spre_clip(iter+1,k_clipped).eq.0.0 .and.
     &      spre_clip(iter,k_clipped).ne.0.0) then
            icount(k_clipped) = icount(k_clipped) + 1
          endif
        enddo
      enddo

c Print the number of storms found for each clipping level.
      do k_clipped=1,max_clips
        print *,k_clipped,icount(k_clipped)
      enddo
      print *

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Sum all of the time series.  Note that k_clipped=1 is the
c   unclipped time series
      do k_clipped=1,max_clips
        spre_sum(k_clipped) = 0.0
        do iter=1,ndays*irecs_in_day
          spre_sum(k_clipped) = spre_sum(k_clipped) +
     &      spre_clip(iter,k_clipped)
        enddo
      enddo

c Add the resulting clipped volume back in the remaining non-zero
c   precipitation events.
      do k_clipped=1,max_clips
        ratio(k_clipped) = spre_sum(1) / spre_sum(k_clipped)
        do iter=1,ndays*irecs_in_day
          spre_clip(iter,k_clipped) = ratio(k_clipped) *
     &      spre_clip(iter,k_clipped)
        enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Find the clipping amount that produces the closest value to
c   nlayers_target.

c Find the clipping value that came closest to that target.
      do k_clipped=1,max_clips
        idiff(k_clipped) = abs(nlayers_target - icount(k_clipped))
c       print *,k_clipped,idiff(k_clipped),icount(k_clipped)
      enddo

c Here I have eliminated the first (unclipped) case from
c   consideration.  This is because often the unclipped case has
c   about 12 'storms' in it, some of which last about a month.
      k_clip = 1000
      k_clip_index = -9999
c     do k_clipped=1,max_clips
      do k_clipped=2,max_clips
        if (idiff(k_clipped).lt.k_clip) then
          k_clip = idiff(k_clipped)
          k_clip_index = k_clipped
        endif
      enddo

      print *,'the clip with closest to 12 storms is',k_clip_index
      print *

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c This is important: I want to just keep the main snowfall
c   'storms'.  So, I don't want to just clip the low precipitation
c   values, because there may be a 'storm' that has low snowfall
c   in one corner of the simulation domain; I want to keep that
c   information.  So, here I am creating a 'storm' time series
c   mask that can be used to clip and scale entire 'storm' events,
c   and toss everything else.

c Note that this is set up so you can just multiply the prec on
c   a given date (iter) by clip_mask(iter), in the mk_mm file.
      do iter=1,ndays*irecs_in_day
        if (spre_clip(iter,k_clip_index).eq.0.0) then
          clip_mask(iter) = 0.0
        else
          clip_mask(iter) = ratio(k_clip_index)
        endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mk_clip_mask_ctl (iyear_init,imonth_init,iday_init,
     &  time_increment,maxiter)

      implicit none

      integer iyear_init,imonth_init,iday_init,maxiter,igrads_dt

      real xhour_init,undef

      character*23 output_fname
      character*40 filename
      character*3 cmo(12)
      character*2 cdt
      character*3 time_increment

      data cmo /'jan','feb','mar','apr','may','jun',
     &          'jul','aug','sep','oct','nov','dec'/

      undef = -9999.0

      filename = 'drizzle_storm_clip.ctl'
      output_fname = 'drizzle_storm_clip.gdat'

      if (time_increment.eq.'3hr') then
        xhour_init = 3.0
        igrads_dt = 3
        cdt = 'hr'
      elseif (time_increment.eq.'1dy') then
        xhour_init = 12.0
        igrads_dt = 1
        cdt = 'dy'
      else
        print *, 'the mk_ctl program cannot deal with this dt/run value'
        stop
      endif

      open (71,file=filename)

      write (71,51) output_fname

      write (71,52)
      write (71,53) undef
c (i,j) indexing.
      write (71,54) 1,1.0,1.0
      write (71,55) 1,1.0,1.0

      write (71,56)
      write (71,57) maxiter,nint(xhour_init),iday_init,
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,68)

      close (71)

c This "a" by itself clips the trailing blanks in the a80 string.
   51 format ('DSET ^',a)
   52 format ('TITLE Precipitation clipping data for SnowModel domain')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF        1 LINEAR 1 1')
c This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     1')
   59 format ('pclip  0  0 precip clipping mask (0 and scaling factor)')
   68 format ('ENDVARS')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mk_spre_clipped_ctl (max_clips,iyear_init,imonth_init,
     &  iday_init,time_increment,maxiter)

      implicit none

      integer max_clips,iyear_init,imonth_init,iday_init,maxiter,
     &  igrads_dt

      real xhour_init,undef

      character*18 output_fname
      character*40 filename
      character*3 cmo(12)
      character*2 cdt
      character*3 time_increment

      data cmo /'jan','feb','mar','apr','may','jun',
     &          'jul','aug','sep','oct','nov','dec'/

      undef = -9999.0

      filename = 'sprec_clipped.ctl'
      output_fname = 'sprec_clipped.gdat'

      if (time_increment.eq.'3hr') then
        xhour_init = 3.0
        igrads_dt = 3
        cdt = 'hr'
      elseif (time_increment.eq.'1dy') then
        xhour_init = 12.0
        igrads_dt = 1
        cdt = 'dy'
      else
        print *, 'the mk_ctl program cannot deal with this dt/run value'
        stop
      endif

      open (71,file=filename)

      write (71,51) output_fname

      write (71,52)
      write (71,53) undef
c (i,j) indexing.
      write (71,54) max_clips,1.0,1.0
      write (71,55) 1,1.0,1.0

      write (71,56)
      write (71,57) maxiter,nint(xhour_init),iday_init,
     &  cmo(imonth_init),iyear_init,igrads_dt,cdt
      write (71,58)
      write (71,59)
      write (71,68)

      close (71)

c This "a" by itself clips the trailing blanks in the a80 string.
   51 format ('DSET ^',a)
   52 format ('TITLE Precipitation clipping data for SnowModel domain')
   53 format ('UNDEF ',f10.1)

   54 format ('XDEF ',i8,' LINEAR ',2f20.8)
   55 format ('YDEF ',i8,' LINEAR ',2f20.8)

   56 format ('ZDEF        1 LINEAR 1 1')
c This i2.2 puts a zero in front of single digit numbers like 1.
   57 format ('TDEF ',i8,' LINEAR ',i2.2,'Z',i2.2,a3,i4,' ',i2,a2)
   58 format ('VARS     1')
   59 format ('spre  0  0 precipitation (m/time_step)')
   68 format ('ENDVARS')

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_irecs (iyr_start,imo_start,idy_start,nyear,
     &  irec_start,ndays,irecs_in_day,iimo_start,iidy_start,
     &  iimo_end,iidy_end,idays,iday_start)

c Calculate the starting and ending record for this year.

      implicit none

      integer iyr_start,imo_start,idy_start,nyear,iyr_beg,iyr_end,
     &  irec_start,ndays,ioptn,julian_start,julian_beg,julian_end,
     &  irec_end,irecs_in_day,iimo_start,iidy_start,iimo_end,
     &  iidy_end,idays,iday_start

c 'start' is the begining of the simulation.
c 'beg' is the begining of the simulation year.
c 'end' is the end of the simulation year.
      iyr_beg = iyr_start + nyear - 1
      iyr_end = iyr_start + nyear

      ioptn = 3
      call calndr (ioptn,idy_start,imo_start,iyr_start,julian_start)
      call calndr (ioptn,idy_start,imo_start,iyr_beg,julian_beg)
      call calndr (ioptn,idy_start,imo_start,iyr_end,julian_end)

c The end day is one less than this start of the next year.
      julian_end = julian_end - 1

      irec_start = (julian_beg - julian_start) * irecs_in_day + 1
      irec_end = (julian_end - julian_start) * irecs_in_day + 1
      ndays = julian_end - julian_beg + 1

c     print *,nyear,iyr_beg,iyr_end,irec_start,irec_end,ndays

c Find the start and end and number of days in the 'winter' period.
      call calndr (ioptn,iidy_start,iimo_start,iyr_beg,julian_start)
      call calndr (ioptn,iidy_end,iimo_end,iyr_end,julian_end)

      iday_start = julian_start - julian_beg + 1
      idays = julian_end - julian_start + 1

      print *,nyear,iyr_beg,iyr_end,iday_start,idays

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

