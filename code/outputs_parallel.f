c outputs_parallel.f

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine OUTPUTS_PARALLEL(nx,ny,iter,Tair_grid,rh_grid,
     &  uwind_grid,vwind_grid,windspd_grid,winddir_grid,
     &  Qsi_grid,Qli_grid,prec_grid,Tsfc,Qle,Qh,Qe,Qc,Qm,Qf,
     &  e_balance,snow_depth,xro_snow,swe_depth,ro_nsnow,
     &  runoff,rain,sprec,sum_prec,sum_runoff,w_balance,
     &  snow_d,topo_land,wbal_qsubl,sum_sprec,wbal_salt,
     &  wbal_susp,ro_snow_grid,sum_Qcs,canopy_int,Qcs,
     &  iyear,imonth,iday,xhour,undef,deltax,xmn,ymn,
     &  wbal_subgrid,canopy_unload,sum_qsubl,sum_trans,
     &  sum_unload,sum_glacmelt,glacier_melt,swemelt,
     &  sfc_pressure,sum_swemelt,albedo,nrecs_max,
     &  icorr_factor_loop,swesublim,vegtype,iter_start,
     &  seaice_run,print_inc,cloud_frac_grid,
     &  output_path_wo_assim,output_path_wi_assim,print_var,
     &  print_outvars,Qsubl_depth,Qsalt,Qsusp)

c This subroutine is available to provide user-defined outputs.
c   These might be special-case situations, like just writing out
c   data at the end of every day, writing out a few grid cells,
c   saving each data arrays to individual files, etc.

      use caf_module,only: me,l_ny,max_l_ny
      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,iter,max_poly,num_poly,iyear,imonth,iday,
     &  icorr_factor_loop,iter_start,icorr_loop_new,
     &  individual_files,k,irec

      real Tair_grid(nx,max_l_ny),rh_grid(nx,max_l_ny),
     &  uwind_grid(nx,max_l_ny),vwind_grid(nx,max_l_ny),
     &  windspd_grid(nx,max_l_ny),winddir_grid(nx,max_l_ny),
     &  Qsi_grid(nx,max_l_ny),Qli_grid(nx,max_l_ny),
     &  prec_grid(nx,max_l_ny),Tsfc(nx,max_l_ny),
     &  Qle(nx,max_l_ny),Qh(nx,max_l_ny),Qe(nx,max_l_ny),
     &  Qc(nx,max_l_ny),Qm(nx,max_l_ny),Qf(nx,max_l_ny),
     &  e_balance(nx,max_l_ny),snow_depth(nx,max_l_ny),
     &  xro_snow(nx,max_l_ny),swe_depth(nx,max_l_ny),
     &  ro_nsnow(nx,max_l_ny),runoff(nx,max_l_ny),
     &  rain(nx,max_l_ny),sprec(nx,max_l_ny),
     &  sum_prec(nx,max_l_ny),sum_runoff(nx,max_l_ny),
     &  w_balance(nx,max_l_ny),snow_d(nx,max_l_ny),
     &  topo_land(nx,ny),wbal_qsubl(nx,max_l_ny),
     &  sum_sprec(nx,max_l_ny),wbal_salt(nx,max_l_ny),
     &  wbal_susp(nx,max_l_ny),ro_snow_grid(nx,max_l_ny),
     &  sum_Qcs(nx,max_l_ny),canopy_int(nx,max_l_ny),
     &  Qcs(nx,max_l_ny),wbal_subgrid(nx,max_l_ny),
     &  canopy_unload(nx,max_l_ny),sum_qsubl(nx,max_l_ny),
     &  sum_trans(nx,max_l_ny),glacier_melt(nx,max_l_ny),
     &  sum_unload(nx,max_l_ny),sum_glacmelt(nx,max_l_ny),
     &  swemelt(nx,max_l_ny),sfc_pressure(nx,max_l_ny),
     &  sum_swemelt(nx,max_l_ny),swesublim(nx,max_l_ny),
     &  vegtype(nx,max_l_ny),albedo(nx,max_l_ny),
     &  cloud_frac_grid(nx,max_l_ny),Qsubl_depth(nx,max_l_ny),
     &  Qsalt(nx,max_l_ny),Qsusp(nx,max_l_ny)

      real undef,xhour,deltax,pi,rad2deg,seaice_run,print_inc
      double precision xmn,ymn
      double precision nrecs_max,nrecs

      real uwnd(nx,max_l_ny)
      real vwnd(nx,max_l_ny)

c Define the output variable data block.  Note that the print_var
c   "yes/no" array was generated in readparam_code.f.
      real vars(nx_max,ny_max,n_print_vars)
      character*80 output_path_wo_assim,output_path_wi_assim
      character*1 print_var(n_print_vars)
      character*4 print_outvars(n_print_vars)
      character*4 string,me_string
      character*9 fmt

c This was defined in preprocess_code.f, in subroutine mk_ctl_files.
c     data print_outvars /'tair','relh','wspd','qsin','qlin',
c    &                    'qlem','albd','wdir','prec','rpre',
c    &                    'spre','smlt','ssub','roff','glmt',
c    &                    'snod','sden','swed','sspr','ssmt',
c    &                    'cldf','var1','var2','var3','var4',
c    &                    'var5','var6','var7','var8','var9'/

c These now come in from the .par file.
c     character path1*(*)
c     character path2*(*)

      integer i_trailing_blanks,trailing_blanks,i_len_wo,i_len_wi

c Calculate how long the paths are.
      i_len_wo = 80 - trailing_blanks(output_path_wo_assim)
      i_len_wi = 80 - trailing_blanks(output_path_wi_assim)
      fmt = '(I4)'
      write(string,fmt) me
      me_string = string

c     print *, i_len_wo,i_len_wi
c     print*, output_path_wo_assim(1:i_len_wo)
c     print *, output_path_wi_assim(1:i_len_wi)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c BEGIN USER EDIT SECTION.

c Define which variables you want to save, whether they will be
c   output at every time step or some time-step increment, which
c   directories the data will be put in, etc.

c Define the output file locations (paths).
c These now come in from the .par file.
c     parameter (path1 =
c    &  'outputs/wo_assim/')
c     parameter (path2 =
c    &  'outputs/wi_assim/')

c Write a seperate file for each variable (individual_files = 1).
c   No other option has been implemented here.
      individual_files = 1

c Define the number of time steps you are going to sum or average
c   over.  If you want to output data at every model time step, set
c   print_inc = 1.0.  For run with an hourly time step and data
c   writes once a day, print_inc = 24.0.  For a run with 3-hourly
c   time steps and data writes once a day, print_inc = 8.0.
c This now comes in from the .par file.
c     print_inc = 24.0
c     print_inc = 1.0
c     print_inc = 8.0

c Define the variables you want to save.  The following are the
c   variables this subroutine is currently set up to output.
c   Listed are the output variable name and the corresponding model
c   variable name.
c Note that these "yes/no" flags are now defined in snowmodel.par.

c VALUES AVERAGED OVER THE PERIOD.
c    1   tair(i,j) = Tair_grid(i,j) - 273.15
c    2   relh(i,j) = rh_grid(i,j)
c    3   wspd(i,j) = windspd_grid(i,j)
c    4   qsin(i,j) = Qsi_grid(i,j)
c    5   qlin(i,j) = Qli_grid(i,j)
c    6   qlem(i,j) = Qle(i,j)
c    7   albd(i,j) = albedo(i,j)
c    8   wdir(i,j) = from uwind_grid(i,j) and vwind_grid(i,j)

c VALUES SUMMED OVER THE PERIOD.
c    9   prec(i,j) = prec_grid(i,j)
c   10   rpre(i,j) = rain(i,j)
c   11   spre(i,j) = sprec(i,j)
c   12   smlt(i,j) = swemelt(i,j)
c   13   ssub(i,j) = swesublim(i,j)
c   14   roff(i,j) = runoff(i,j)
c   15   glmt(i,j) = glacier_melt(i,j)

c VALUES SAVED AT THE END OF THE PERIOD.
c   16   snod(i,j) = snow_depth(i,j)
c   17   sden(i,j) = xro_snow(i,j)
c   18   swed(i,j) = swe_depth(i,j)
c   19   sspr(i,j) = sum_sprec(i,j)
c   20   ssmt(i,j) = sum_swemelt(i,j)

c NEW VARIABLES ADDED TO THE OUTPUT DATA BLOCK.
c   (you have to modify the code below to
c    do the processing you want for these)
c   21   cldf(i,j) = cloud_fraction_grid(i,j)

c Define which variables you want to save by placing a yes = 'y' or
c   no = 'n' in front of the variable number.

c VALUES AVERAGED OVER THE PERIOD.
c VALUES AVERAGED OVER THE PERIOD.
c 1 = tair(i,j) = Tair_grid(i,j) - 273.15
c     print_var(1)  = 'y'

c 2 = relh(i,j) = rh_grid(i,j)
c     print_var(2)  = 'n'

c 3 = wspd(i,j) = windspd_grid(i,j)
c     print_var(3)  = 'n'

c 4 = qsin(i,j) = Qsi_grid(i,j)
c     print_var(4)  = 'n'

c 5 = qlin(i,j) = Qli_grid(i,j)
c     print_var(5)  = 'n'

c 6 = qlem(i,j) = Qle(i,j)
c     print_var(6)  = 'n'

c 7 = albd(i,j) = albedo(i,j)
c     print_var(7)  = 'n'

c 8 = wdir(i,j) = from uwind_grid(i,j) and vwind_grid(i,j)
c     print_var(8)  = 'n'

c VALUES SUMMED OVER THE PERIOD.
c VALUES SUMMED OVER THE PERIOD.
c  9 = prec(i,j) = prec_grid(i,j)
c     print_var(9)  = 'y'

c 10 = rpre(i,j) = rain(i,j)
c     print_var(10) = 'y'

c 11 = spre(i,j) = sprec(i,j)
c     print_var(11) = 'y'

c 12 = smlt(i,j) = swemelt(i,j)
c     print_var(12) = 'n'

c 13 = ssub(i,j) = swesublim(i,j)
c     print_var(13) = 'n'

c 14 = roff(i,j) = runoff(i,j)
c     print_var(14) = 'n'

c 15 = glmt(i,j) = glacier_melt(i,j)
c     print_var(15) = 'n'

c VALUES SAVED AT THE END OF THE PERIOD.
c VALUES SAVED AT THE END OF THE PERIOD.
c 16 = snod(i,j) = snow_depth(i,j)
c     print_var(16) = 'y'

c 17 = sden(i,j) = xro_snow(i,j)
c     print_var(17) = 'y'

c 18 = swed(i,j) = swe_depth(i,j)
c     print_var(18) = 'y'

c 19 = sspr(i,j) = sum_sprec(i,j)
c     print_var(19) = 'y'

c 20 = ssmt(i,j) = sum_swemelt(i,j)
c     print_var(20) = 'y'

c NEW VARIABLES ADDED TO THE OUTPUT DATA BLOCK.
c NEW VARIABLES ADDED TO THE OUTPUT DATA BLOCK.
c   (you have to modify the code below to
c    do the processing you want for these)
c 21 = cldf(i,j) = cloud_fraction_grid(i,j)
c     print_var(21) = 'n'

c Extra variables.
c     print_var(22) = 'n'
c     print_var(23) = 'n'
c     print_var(24) = 'n'
c     print_var(25) = 'n'
c     print_var(26) = 'n'
c     print_var(27) = 'n'
c     print_var(28) = 'n'
c     print_var(29) = 'n'
c     print_var(30) = 'n'

c Note that this data output implementation is currently configured
c   to mask out the ocean points (vegtype.eq.24.0) if this is a
c   land run (seaice_run = 0.0); mask out all land points (vegtype.
c   ne.24.0) if this is an ocean/sea ice run (seaice_run = 1.0, 3.0,
c   or 4.0); and to not mask out anything if this is a combined land
c   and sea ice run (seaice_run = 2.0).

c END USER EDIT SECTION.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Define the constants used in the wind-direction averaging.
      pi = 2.0 * acos(0.0)
      rad2deg = 180.0 / pi

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      uwnd = 0.0
      vwnd = 0.0

c Use individual output files for each variable.
      if (individual_files.eq.1) then

        if (iter.eq.iter_start) then
          nrecs = nx * ny
          if (nrecs.gt.nrecs_max) then
            print *,'Your simulation domain has too many grid cells'
            print *,'to print the .gdat files like the write statements'
            print *,'are currently configured.  You will have to change'
            print *,'them to look like:'
            print *,'    do j=1,ny'
            print *,'      write (51,rec=j) (var(i,j),i=1,nx)'
            print *,'    enddo'
            stop
          endif
        endif

c Open individual output files for each variable.
        if (iter.eq.iter_start) then
          if (icorr_factor_loop.eq.1) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                open (220+k+me,
     &file=output_path_wo_assim(1:i_len_wo)//'/'//print_outvars(k)//
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',
     &            form='unformatted',access='direct',recl=4*nx*ny,
     &            status='replace')
              endif
            enddo
          endif

          if (icorr_factor_loop.eq.2) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                open (320+k+me,
     &file=output_path_wi_assim(1:i_len_wi)//'/'//print_outvars(k)//
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',
     &            form='unformatted',access='direct',recl=4*nx*ny,
     &            status='replace')
              endif
            enddo
          endif
        endif

        if (iter.eq.iter_start) then
c Initialize the averaging and summing arrays.
          do j=1,ny
            do i=1,nx
              do k=1,n_print_vars
                vars(i,j,k) = 0.0
              enddo
              uwnd(i,j) = 0.0
              vwnd(i,j) = 0.0
            enddo
          enddo
        endif

c Perform the avaraging, summing, etc.
        do j=1,ny
          do i=1,nx
c Values averaged over the period.
            vars(i,j,1) = vars(i,j,1) + (Tair_grid(i,j) - 273.15) /
     &        print_inc
            vars(i,j,2) = vars(i,j,2) + rh_grid(i,j) / print_inc
            vars(i,j,3) = vars(i,j,3) + windspd_grid(i,j) / print_inc
            vars(i,j,4) = vars(i,j,4) + Qsi_grid(i,j) / print_inc
            vars(i,j,5) = vars(i,j,5) + Qli_grid(i,j) / print_inc
            vars(i,j,6) = vars(i,j,6) + Qle(i,j) / print_inc
            vars(i,j,7) = vars(i,j,7) + albedo(i,j) / print_inc

            uwnd(i,j) = uwnd(i,j) + uwind_grid(i,j) / print_inc
            vwnd(i,j) = vwnd(i,j) + vwind_grid(i,j) / print_inc

c Some compilers do not allow both u and v to be 0.0 in
c   the atan2 computation.
            if (abs(uwnd(i,j)).lt.1e-10) uwnd(i,j) = 1e-10

            vars(i,j,8) = rad2deg * atan2(uwnd(i,j),vwnd(i,j))
            if (vars(i,j,8).ge.180.0) then
              vars(i,j,8) = vars(i,j,8) - 180.0
            else
              vars(i,j,8) = vars(i,j,8) + 180.0
            endif

c Values summed over the period.
            vars(i,j,9) = vars(i,j,9) + prec_grid(i,j)
            vars(i,j,10) = vars(i,j,10) + rain(i,j)
            vars(i,j,11) = vars(i,j,11) + sprec(i,j)
            vars(i,j,12) = vars(i,j,12) + swemelt(i,j)
            vars(i,j,13) = vars(i,j,13) + swesublim(i,j)
            vars(i,j,14) = vars(i,j,14) + runoff(i,j)
            vars(i,j,15) = vars(i,j,15) + glacier_melt(i,j)

c Values saved at the end of the day.
            vars(i,j,16) = snow_depth(i,j)
            vars(i,j,17) = xro_snow(i,j)
            vars(i,j,18) = swe_depth(i,j)
            vars(i,j,19) = sum_sprec(i,j)
            vars(i,j,20) = sum_swemelt(i,j)

c Values averaged over the period.
          vars(i,j,21) = vars(i,j,21) + cloud_frac_grid(i,j) / print_inc

c New variables.
            vars(i,j,22) = vars(i,j,22) + Qsubl_depth(i,j)
            vars(i,j,23) = sum_trans(i,j)
            vars(i,j,24) = vars(i,j,24) + Qsalt(i,j)
            vars(i,j,25) = vars(i,j,25) + Qsusp(i,j)

          enddo
        enddo

c Check to see whether this is the data-write time step.
        if (mod(iter,nint(print_inc)).eq.0) then

c Mask out the ocean points (vegtype.eq.24.0) if this is
c   a land run (seaice_run = 0.0).  Mask out all land points
c   (vegtype.ne.24.0) if this is an ocean/sea ice run
c   (seaice_run = 1.0, 3.0, or 4.0).  Do not mask out anything
c   if this this is a combined land and sea ice run (seaice_run
c   = 2.0).
          if (seaice_run.eq.0.0) then
            do j=1,ny
              do i=1,nx
                if (vegtype(i,j).eq.24.0) then
                  do k=1,n_print_vars
                    vars(i,j,k) = undef
                  enddo
                endif
              enddo
            enddo
          elseif (seaice_run.eq.1.0 .or. seaice_run.eq.3.0 .or.
     &      seaice_run.eq.4.0) then
            do j=1,ny
              do i=1,nx
                if (vegtype(i,j).ne.24.0) then
                  do k=1,n_print_vars
                    vars(i,j,k) = undef
                  enddo
                endif
              enddo
            enddo
          endif

c Write out the data.
          irec = iter / nint(print_inc)
          if (icorr_factor_loop.eq.1) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                write (220+k+me,rec=irec) ((vars(i,j,k),i=1,nx),j=1,ny)
              endif
            enddo
          elseif (icorr_factor_loop.eq.2) then
            do k=1,n_print_vars
              if (print_var(k).eq.'y') then
                write (320+k+me,rec=irec) ((vars(i,j,k),i=1,nx),j=1,ny)
              endif
            enddo
          endif

c Reinitialize the averaging and summing arrays.
          do j=1,ny
            do i=1,nx
              do k=1,n_print_vars
                vars(i,j,k) = 0.0
              enddo
              uwnd(i,j) = 0.0
              vwnd(i,j) = 0.0
            enddo
          enddo
        endif

c Use more than one variable in an output file.
      else

        print *,'Use more than one variable in an output file:'
        print *,'  THIS HAS NOT BEEN IMPLEMENTED YET'

      endif

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c THIS IS AN EXAMPLE OF SAVING DATA IN ASCII/TEXT FORMAT.

c I have completely removed this example; I now do this with
c   improved codes as part of post-processing steps.  It is
c   just too slow to do it as part of the model simulation.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c THE CODE BELOW WAS USED TO SAVE AVERAGES OVER POLYGONS.

c I have completely removed this example; if I were to do this
c   again I would do it as a post-processing step.  If you really
c   want to see what I did here, you can look in one of the pre-
c   2018 code distributions.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
