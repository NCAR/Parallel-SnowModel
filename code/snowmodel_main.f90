program snowmodel_main
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c snowmodel_main.f
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !c All units are in m, kg, s, K.
  !c
  !c The author of this code is:
  !c
  ! c   Dr. Glen E. Liston
  ! c   Cooperative Institute for Research
  ! c     in the Atmosphere (CIRA)         |
  ! c   Colorado State University
  ! c   Fort Collins, Colorado 80523-1375
  ! c
  ! c   Voice: (970) 491-8220
  ! c   FAX: (970) 491-8241
  ! c
  ! c   glen.liston@colostate.edu

  use snowmodel_inc
  use snowmodel_vars
  use caf_module

  implicit none
  real :: start,finish,xlat_main
  integer*8 :: readpar_start,readpar_end, COUNT_RATE
  integer*8 :: preproc_start,preproc_end
  integer*8 :: master_start,master_end
  integer*8 :: micromet_start,micromet_end,total_micromet_time
  integer*8 :: enbal_start,enbal_end,total_enbal_time
  integer*8 :: snowpack_start,snowpack_end,total_snowpack_time
  integer*8 :: snowtran_start,snowtran_end,total_snowtran_time
  integer*8 :: output_start,output_end,total_output_time
  integer*8 :: distribute_start,distribute_end,total_distribute_time
  integer*8 :: gather_start,gather_end,total_gather_time
  integer*8 :: noninit_start, noninit_end
  total_micromet_time = 0
  total_enbal_time = 0
  total_snowpack_time = 0
  total_snowtran_time = 0
  total_output_time = 0
  total_distribute_time = 0
  total_gather_time = 0

  !include 'snowmodel.inc'
  ! include 'snowmodel_vars.inc'

  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! cccccccccccccccccccccc  INITIALIZE THE MODEL  cccccccccccccccccccccccc
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! c The following provides the ability to read in a 'snowmodel.par'
  ! c   file that is called something other than 'snowmodel.par' and/or
  ! c   located somewhere different than the directory where the
  ! c   'snowmodel' executable file is located.  If you run:
  ! c     snowmodel
  ! c   it assumes the .par file is called 'snowmodel.par' and is
  ! c   located in the default location.  You can also run the code as
  ! c   follows:
  ! c     snowmodel parpath/parname.par
  ! c   and it will look for the .par file in the different location
  ! c   and/or the different name.

  if (me == 0) then
    call System_Clock(master_start, COUNT_RATE)
  endif

  CALL GET_DOT_PAR_PATH_LOCATION(snowmodel_dot_par_fname)

  !c Read the input parameters.
  if (me == 0) then
    call System_Clock(readpar_start)
  endif
  CALL READPARAM_CODE(dt,deltax,deltay,Utau_t_flag,&
       &  subgrid_flag,twolayer_flag,snowmodel_dot_par_fname,&
       &  bc_flag,curve_len_scale,slopewt,curvewt,ht_windobs,&
       &  ht_rhobs,ro_snow,snow_d_init_const,const_veg_flag,&
       &  vegsnowdepth,nx,ny,max_iter,met_input_fname,xmn,ymn,&
       &  iyear_init,imonth_init,iday_init,xhour_init,undef,ifill,&
       &  iobsint,dn,xlat,i_tair_flag,i_rh_flag,i_wind_flag,&
       &  i_solar_flag,i_prec_flag,isingle_stn_flag,igrads_metfile,&
       &  windspd_min,icond_flag,run_micromet,run_enbal,run_snowpack,&
       &  run_snowtran,topoflag,topoveg_fname,snowtran_output_fname,&
       &  micromet_output_fname,enbal_output_fname,Utau_t_const,&
       &  snowpack_output_fname,print_micromet,print_enbal,&
       &  print_snowpack,print_snowtran,i_longwave_flag,print_user,&
       &  ascii_topoveg,topo_ascii_fname,veg_ascii_fname,&
       &  irun_data_assim,lapse_rate_user_flag,&
       &  iprecip_lapse_rate_user_flag,use_shortwave_obs,&
       &  use_longwave_obs,use_sfc_pressure_obs,calc_subcanopy_met,&
       &  sfc_sublim_flag,gap_frac,cloud_frac_factor,&
       &  albedo_snow_forest,albedo_snow_clearing,albedo_glacier,&
       &  barnes_lg_domain,n_stns_used,tabler_dir,slope_adjust,&
       &  lat_solar_flag,UTC_flag,iveg_ht_flag,ihrestart_flag,&
       &  ihrestart_inc,i_dataassim_loop,tsls_threshold,dz_snow_min,&
       &  print_multilayer,multilayer_snowpack,max_layers,&
       &  multilayer_output_fname,izero_snow_date,curve_lg_scale_flag,&
       &  check_met_data,seaice_run,snowmodel_line_flag,wind_lapse_rate,&
       &  iprecip_scheme,cf_precip_flag,snowfall_frac,print_inc,&
       &  output_path_wo_assim,output_path_wi_assim,Tabler_1_flag,&
       &  Tabler_2_flag,tabler_sfc_path_name,print_var,print_parallel)

  if(nx /= 0 .and. ny /= 0) then
    call allocate_arrays(nx,ny,me,lat_solar_flag,seaice_run,&
                        & UTC_flag,snowmodel_line_flag,&
                        & print_parallel,print_multilayer)
  endif
  if (me == 0) then
    call System_Clock(readpar_end)
  endif

  ! c This loop runs the correction/data assimilation adjustment
  ! c   iterations.
  if (ihrestart_flag.ge.0) then
     if (i_dataassim_loop.lt.0.0) then
        i_corr_start = 2
     else
        i_corr_start = 1
     endif
  else
     i_corr_start = 1
  endif

  do icorr_factor_loop=i_corr_start,irun_data_assim+1

     ! c Perform the correction (precipitation and melt) factor
     ! c   calculations.
     if (irun_data_assim.eq.1 .and. icorr_factor_loop.eq.2) then
        CALL DATAASSIM_USER(nx,ny,icorr_factor_index,&
             &      corr_factor,max_iter,deltax,deltay,xmn,ymn,nobs_dates,&
             &      print_inc,iday_init,imonth_init,iyear_init,dt,&
             &      output_path_wo_assim,xhour_init)
        if (ihrestart_flag.ge.-1) then
           CALL HRESTART_SAVE_DA(nx,ny,max_iter,corr_factor,&
                &        icorr_factor_index,nobs_dates)
        endif
     endif

     ! c Perform a variety of preprocessing and model setup steps, like
     ! c   read in topography and vegetation arrays, open input and output
     ! c   files, etc.
     if (me == 0) then
       call System_Clock(preproc_start)
     endif

     call caf_partitioning()
     sync all

     CALL PREPROCESS_CODE(topoveg_fname,const_veg_flag,&
          &    l_vegtype,l_veg_z0,vegsnowdepth,fetch,xmu,C_z,h_const,&
          &    wind_min,Up_const,dz_susp,ztop_susp,fall_vel,Ur_const,&
          &    ro_water,ro_air,gravity,vonKarman,pi,twopio360,snow_z0,&
          &    nx,ny,l_sum_sprec,l_sum_qsubl,l_sum_trans,l_sum_unload,l_topo_tmp,&
          &    l_topo_land,l_snow_d,topoflag,l_snow_d_init,snow_d_init_const,&
          &    l_soft_snow_d,met_input_fname,igrads_metfile,deltax,deltay,&
          &    snowtran_output_fname,micromet_output_fname,&
          &    enbal_output_fname,snowpack_output_fname,print_micromet,&
          &    print_enbal,print_snowpack,print_snowtran,run_micromet,&
          &    run_enbal,run_snowpack,run_snowtran,l_ro_snow_grid,l_swe_depth,&
          &    l_sum_runoff,l_sum_prec,ro_snow,twolayer_flag,l_sum_Qcs,&
          &    l_canopy_int,ascii_topoveg,topo_ascii_fname,icorr_factor_loop,&
          &    veg_ascii_fname,undef,isingle_stn_flag,max_iter,&
          &    i_tair_flag,i_rh_flag,i_wind_flag,i_prec_flag,l_sum_glacmelt,&
          &    l_snow_depth,l_sum_d_canopy_int,corr_factor,icorr_factor_index,&
          &    l_sum_sfcsublim,barnes_lg_domain,n_stns_used,l_k_stn,xmn,ymn,&
          &    l_ro_soft_snow_old,l_sum_swemelt,xlat,lat_solar_flag,l_xlat_grid,&
          &    l_xlon_grid,UTC_flag,dt,l_swe_depth_old,l_canopy_int_old,&
          &    l_vegsnowd_xy,iveg_ht_flag,ihrestart_flag,i_dataassim_loop,&
          &    multilayer_snowpack,max_layers,multilayer_output_fname,&
          &    print_multilayer,l_KK,l_tslsnowfall,tsls_threshold,&
          &    irun_data_assim,izero_snow_date,iclear_mn,iclear_dy,&
          &    xclear_hr,l_snod_layer,l_swed_layer,l_ro_layer,l_T_old,l_gamma,&
          &    icond_flag,curve_lg_scale_flag,l_curve_wt_lg,check_met_data,&
          &    seaice_run,snowmodel_line_flag,xg_line,yg_line,print_user,&
          &    cf_precip_flag,l_cf_precip,print_inc,xhour_init,Tabler_1_flag,&
          &    Tabler_2_flag,iyear_init,imonth_init,iday_init,print_var,&
          &    output_path_wo_assim,output_path_wi_assim,nrecs_max,&
          &    tabler_sfc_path_name,print_outvars,l_diam_layer)

      call distribute2(snowmodel_line_flag)
      call single_gather(topo,l_topo_tmp)
      call single_gather(topo_land,l_topo_land)
      xlat_main = l_xlat_grid(1,1)[1] !Distribute for micromet get_lapse_rates()
      sync all
      if (me == 0) then
        call System_Clock(preproc_end)
      endif

     ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! ccccccccccccccccccccccccc  RUN THE MODEL  cccccccccccccccccccccccccccc
     ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     if (me.eq.0) call System_Clock(noninit_start)
     ! c Loop through the number of time steps (iterations) in the
     ! c   simulation.
     if (ihrestart_flag.ge.0) then
        iter_start = ihrestart_flag + 1
     else
        iter_start = 1
     endif

     !c Distribute data to images

     do iter=iter_start,max_iter

        !c Distribute the meteorological station data.
        if (run_micromet.eq.1.0) then
           if (me == 0) then
             call System_Clock(micromet_start)
           endif
           CALL MICROMET_CODE(nx,ny,xmn,ymn,deltax,deltay,&
                &        iyear_init,imonth_init,iday_init,xhour_init,dt,undef,&
                &        ifill,iobsint,dn,iter,curve_len_scale,slopewt,curvewt,&
                &        topo,l_curvature,l_terrain_slope,l_slope_az,l_Tair_grid,&
                &        l_rh_grid,l_uwind_grid,l_vwind_grid,l_Qsi_grid,l_prec_grid,&
                &        i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,&
                &        i_prec_flag,isingle_stn_flag,igrads_metfile,&
                &        l_windspd_grid,l_winddir_grid,windspd_flag,winddir_flag,&
                &        l_sprec,windspd_min,l_Qli_grid,i_longwave_flag,l_vegtype,&
                &        forest_LAI,iyear,imonth,iday,xhour,corr_factor,&
                &        icorr_factor_index,lapse_rate_user_flag,&
                &        iprecip_lapse_rate_user_flag,use_shortwave_obs,&
                &        use_longwave_obs,use_sfc_pressure_obs,l_sfc_pressure,&
                &        run_enbal,run_snowpack,calc_subcanopy_met,l_vegsnowd_xy,&
                &        gap_frac,cloud_frac_factor,barnes_lg_domain,n_stns_used,&
                &        l_k_stn,l_xlat_grid,l_xlon_grid,UTC_flag,icorr_factor_loop,&
                &        snowmodel_line_flag,l_xg_line,l_yg_line,irun_data_assim,&
                &        wind_lapse_rate,iprecip_scheme,cf_precip_flag,l_cf_precip,&
                &        l_cloud_frac_grid,snowfall_frac,seaice_run,xlat_main)

          if (me == 0) then
            call System_Clock(micromet_end)
            total_micromet_time = total_micromet_time + (micromet_end - micromet_start)
          endif

!           call manual_bcast()
        endif

        !c Perform a surface energy balance over the domain.
        if (run_enbal.eq.1.0) then
           if (me == 0) then
             call System_Clock(enbal_start)
           endif
           CALL ENBAL_CODE(nx,ny,l_Tair_grid,l_uwind_grid,l_sfc_pressure,&
                &        l_vwind_grid,l_rh_grid,l_Tsfc,l_Qsi_grid,l_Qli_grid,l_Qle,l_Qh,l_Qe,&
                &        l_Qc,l_Qm,l_e_balance,l_Qf,l_snow_d,ht_windobs,icond_flag,&
                &        l_albedo,snow_z0,l_veg_z0,l_vegtype,undef,albedo_snow_forest,&
                &        albedo_snow_clearing,albedo_glacier,l_snod_layer,l_T_old,&
                &        l_gamma,l_KK)
           if (me == 0) then
              call System_Clock(enbal_end)
              total_enbal_time = total_enbal_time + (enbal_end - enbal_start)
           endif

 !          call manual_bcast()
        endif

        ! c Evolve the snowpack according to the defined melt and
        ! c   precipitation inputs.
        if (run_snowpack .eq. 1.0) then
           if (me == 0) then
             call System_Clock(snowpack_start)
           endif
           CALL SNOWPACK_CODE(nx,ny,l_Tair_grid,l_rh_grid,l_ro_nsnow, &
                &        dt,l_swe_depth,l_Tsfc,l_snow_d,l_prec_grid,l_runoff,l_Qm,l_rain,&
                &        l_sprec,iter,l_w_balance,l_sum_prec,l_sum_runoff,l_xro_snow,&
                &        undef,ro_snow,l_ro_snow_grid,l_soft_snow_d,l_sum_sprec,&
                &        l_snow_depth,l_windspd_grid,l_Qsi_grid,l_sum_Qcs,l_canopy_int,&
                &        l_Qcs,l_vegtype,forest_LAI,l_albedo,l_glacier_melt,&
                &        l_canopy_unload,l_sum_unload,l_sum_glacmelt,run_snowtran,&
                &        l_swemelt,l_d_canopy_int,l_sum_d_canopy_int,l_snow_d_init,&
                &        l_sfc_pressure,l_Qe,sfc_sublim_flag,l_sum_sfcsublim,&
                &        l_sum_swemelt,corr_factor,icorr_factor_index,l_swesublim,&
                &        l_swe_depth_old,l_canopy_int_old,l_KK,max_layers,l_melt_flag,&
                &        ro_snowmax,tsls_threshold,dz_snow_min,l_tslsnowfall,&
                &        l_change_layer,l_snod_layer,l_swed_layer,l_ro_layer,l_T_old,l_gamma,&
                &        multilayer_snowpack,seaice_run,seaice_conc,ht_windobs,&
                &        l_windspd_2m_grid,l_diam_layer,l_flux_layer, l_sum_trans)
           if (me == 0) then
              call System_Clock(snowpack_end)
              total_snowpack_time = total_snowpack_time + (snowpack_end - snowpack_start)
           endif
        endif

        !call manual_bcast()!Only for testing purposes!
        call single_gather(vwind_grid,l_vwind_grid)

        !c Run the blowing-snow model.
        if (run_snowtran.eq.1.0) then
           if (me == 0) then
             call System_Clock(snowtran_start)
           endif
           call single_gather(vwind_grid,l_vwind_grid)
           CALL SNOWTRAN_CODE(bc_flag,bs_flag,C_z,&
                &        l_conc_salt,deltax,deltay,l_dh_salt,l_dh_salt_u,l_dh_salt_v,&
                &        l_dh_susp,l_dh_susp_u,l_dh_susp_v,dt,dz_susp,fall_vel,fetch,&
                &        gravity,h_const,h_star,ht_rhobs,ht_windobs,index_ue,&
                &        index_uw,index_vn,index_vs,iter,nx,ny,pi,l_Qsalt,l_Qsalt_max,&
                &        l_Qsalt_maxu,l_Qsalt_maxv,l_Qsalt_u,l_Qsalt_v,l_Qsubl,l_Qsusp,&
                &        l_Qsusp_u,l_Qsusp_v,l_rh_grid,ro_air,ro_snow,ro_water,l_snow_d,&
                &        snow_d_init,snow_z0,l_soft_snow_d,l_sprec,l_sum_glacmelt,&
                &        subgrid_flag,l_wbal_salt,l_wbal_susp,l_wbal_qsubl,l_sum_sprec,&
                &        l_tabler_ee,l_tabler_ne,l_tabler_nn,l_tabler_nw,l_tabler_se,&
                &        l_tabler_ss,l_tabler_sw,l_tabler_ww,l_tair_grid,topo,topo_land,&
                &        topoflag,twolayer_flag,Up_const,Ur_const,Utau,&
                &        l_Utau_t,l_uwind_grid,l_veg_z0,l_vegsnowd_xy,vonKarman,&
                &        vwind_grid,wind_min,winddir_flag,l_winddir_grid,&
                &        windspd_flag,l_windspd_grid,xmu,z_0,ztop_susp,max_iter,&
                &        run_enbal,run_snowpack,l_wbal_subgrid,l_sum_qsubl,l_sum_trans,&
                &        l_swe_depth,l_snow_depth,l_ro_snow_grid,l_sum_prec,l_sum_runoff,&
                &        l_sum_Qcs,l_canopy_int,l_w_balance,l_sum_sfcsublim,tabler_dir,&
                &        slope_adjust,Utau_t_const,Utau_t_flag,l_ro_soft_snow_old,&
                &        l_ro_soft_snow,l_ro_nsnow,l_prec_grid,l_Qcs,l_runoff,l_d_canopy_int,&
                &        l_glacier_melt,l_swe_depth_old,l_swesublim,l_canopy_unload,&
                &        l_canopy_int_old,iter_start,multilayer_snowpack,l_swed_layer,&
                &        l_KK,l_snod_layer,l_ro_layer,curve_lg_scale_flag,l_curve_wt_lg,&
                &        seaice_run,seaice_conc,l_tslsnowfall,l_T_old,tsls_threshold,&
                &        curve_len_scale,Tabler_1_flag,Tabler_2_flag,undef,&
                &        tabler_sfc_path_name,output_path_wo_assim,&
                &        output_path_wi_assim,icorr_factor_loop,l_windspd_2m_grid,&
                &        l_Qsubl_depth)
            if (me == 0) then
                call System_Clock(snowtran_end)
                total_snowtran_time = total_snowtran_time + (snowtran_end - snowtran_start)
            endif
        endif
        if (print_parallel.eq.0.0) then
          if (me == 0) then
            call System_Clock(gather_start)
          endif
          call gather(print_multilayer)
          sync all
          if (me == 0) then
            call System_Clock(gather_end)
            total_gather_time = total_gather_time + (gather_end - gather_start)
            call System_Clock(output_start)
          endif
        else
          if (me==0) then
            call System_Clock(gather_start)
            call System_Clock(gather_end)
            total_gather_time = total_gather_time + (gather_end - gather_start)
            call System_Clock(output_start)
          endif
        endif

        !Required by the presence of "distribute" at the beginning of the main loop
!        call single_gather(curve_wt_lg,l_curve_wt_lg)

        ! c If this is a sea ice run with incremental remapping, perform the
        ! c   remapping here, before any data is written out for this time step.
        ! c   Note that the incremental remapping programs must be compiled with
        ! c   pgf90 or gfortran.
        ! c         if (seaice_run.eq.3.0) then
        ! c           call remapper_main(iter,nx,ny,seaice_conc,swe_depth,dt,
        ! c    &        deltax,deltay)
        ! c         endif

        ! c Note that here we write out the entire potential vertical domain
        ! c   (nz_max), because it may have been filled with snow at some point
        ! c   during the simulation.
        if (run_snowpack .eq. 1.0 .and. multilayer_snowpack .eq. 1 .and. &
             &      print_multilayer.eq.1.0 .and. me == 0 .and. &
             &      print_parallel.eq.0.0) then
           write(401,rec=iter)&
                &        ((real(KK(i,j)),i=1,nx),j=1,ny),&
                &        ((snow_depth(i,j),i=1,nx),j=1,ny),&
                &        ((xro_snow(i,j),i=1,nx),j=1,ny),&
                &        ((swe_depth(i,j),i=1,nx),j=1,ny),&
                &        (((snod_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),&
                &        (((ro_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),&
                &        (((swed_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max),&
                &        (((diam_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
        elseif (run_snowpack .eq. 1.0 .and. multilayer_snowpack .eq. 1 &
             &      .and. print_multilayer .eq. 2.0 .and. me == 0 .and. &
             &      print_parallel.eq.0.0) then
           write(401,rec=iter) &
                &        ((real(KK(i,j)),i=1,nx),j=1,ny),&
                &        ((snow_depth(i,j),i=1,nx),j=1,ny),&
                &        ((xro_snow(i,j),i=1,nx),j=1,ny),&
                &        ((swe_depth(i,j),i=1,nx),j=1,ny)
           write(402,rec=iter) &
                &        (((snod_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
           write(403,rec=iter)&
                &        (((ro_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
           write(404,rec=iter)&
                &        (((swed_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
           write(405,rec=iter)&
                &        (((diam_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
           write(406,rec=iter)&
                &        (((flux_layer(i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
           write(407,rec=iter)&
                &        (((T_old(i,j,k)-273.15,i=1,nx),j=1,ny),k=1,nz_max)
           write(408,rec=iter)&
                &        (((gamma(i,j,k),i=1,nx),j=1,ny),k=1,nz_max)
        endif

        if (print_micromet.eq.1.0 .and. me == 0 .and. print_parallel.eq.0.0) then
           write(81,rec=iter)&
                &          ((Tair_grid(i,j)-273.15,i=1,nx),j=1,ny),&
                &          ((rh_grid(i,j),i=1,nx),j=1,ny),&
                &          ((uwind_grid(i,j),i=1,nx),j=1,ny),&
                &          ((vwind_grid(i,j),i=1,nx),j=1,ny),&
                &          ((windspd_grid(i,j),i=1,nx),j=1,ny),&
                &          ((winddir_grid(i,j),i=1,nx),j=1,ny),&
                &          ((Qsi_grid(i,j),i=1,nx),j=1,ny),&
                &          ((Qli_grid(i,j),i=1,nx),j=1,ny),&
                &          ((prec_grid(i,j),i=1,nx),j=1,ny)
        endif

        if (print_enbal.eq.1.0 .and. me == 0 .and. print_parallel.eq.0.0) then
           write(82,rec=iter) &
                &          ((Tair_grid(i,j)-273.15,i=1,nx),j=1,ny),&
                &          ((Tsfc(i,j)-273.15,i=1,nx),j=1,ny),&
                &          ((Qsi_grid(i,j),i=1,nx),j=1,ny),&
                &          ((Qli_grid(i,j),i=1,nx),j=1,ny),&
                &          ((Qle(i,j),i=1,nx),j=1,ny),&
                &          ((Qh(i,j),i=1,nx),j=1,ny),&
                &          ((Qe(i,j),i=1,nx),j=1,ny),&
                &          ((Qc(i,j),i=1,nx),j=1,ny),&
                &          ((Qm(i,j),i=1,nx),j=1,ny),&
                &          ((albedo(i,j),i=1,nx),j=1,ny),&
                &          ((e_balance(i,j),i=1,nx),j=1,ny)
        endif

        !c Save the outputs from the SNOWPACK and SNOWTRAN routines.
        if (run_snowpack.eq.1.0 .and. print_snowpack.eq.1.0 .and. me == 0 .and. &
             &    print_parallel.eq.0.0) then
           write(83,rec=iter)&
                &        ((snow_depth(i,j),i=1,nx),j=1,ny),&
                &        ((xro_snow(i,j),i=1,nx),j=1,ny),&
                &        ((swe_depth(i,j),i=1,nx),j=1,ny),&
                &        ((runoff(i,j),i=1,nx),j=1,ny),&
                &        ((rain(i,j),i=1,nx),j=1,ny),&
                &        ((sprec(i,j),i=1,nx),j=1,ny),&
                &        ((Qcs(i,j),i=1,nx),j=1,ny),&
                &        ((canopy_int(i,j),i=1,nx),j=1,ny),&
                &        ((sum_Qcs(i,j),i=1,nx),j=1,ny),&
                &        ((sum_prec(i,j),i=1,nx),j=1,ny),&
                &        ((sum_sprec(i,j),i=1,nx),j=1,ny),&
                &        ((sum_unload(i,j),i=1,nx),j=1,ny),&
                &        ((sum_runoff(i,j),i=1,nx),j=1,ny),&
                &        ((sum_swemelt(i,j),i=1,nx),j=1,ny),&
                &        ((sum_sfcsublim(i,j),i=1,nx),j=1,ny),&
                &        ((w_balance(i,j),i=1,nx),j=1,ny)
        endif

        if (run_snowtran.eq.1.0 .and. print_snowtran.eq.1.0 .and. me == 0 .and. &
            &     print_parallel.eq.0.0) then
           write(84,rec=iter)&
                &        ((snow_d(i,j),i=1,nx),j=1,ny),&
                &        ((wbal_qsubl(i,j),i=1,nx),j=1,ny),&
                &        ((wbal_salt(i,j),i=1,nx),j=1,ny),&
                &        ((wbal_susp(i,j),i=1,nx),j=1,ny),&
                &        ((wbal_subgrid(i,j),i=1,nx),j=1,ny),&
                &        ((sum_qsubl(i,j),i=1,nx),j=1,ny),&
                &        ((sum_trans(i,j),i=1,nx),j=1,ny)
        endif


        ! c The call to outputs_user is available to provide user-defined
        ! c   outputs.  These might be special-case situations, like just
        ! c   writing out data at the end of every day, writing out a few
        ! c   grid cells, saving each data arrays to individual files, etc.
        if (print_user .eq. 1.0 .and. me == 0 .and. print_parallel .eq. 0.0) then
           CALL OUTPUTS_USER(nx,ny,iter,Tair_grid,rh_grid,&
                &        uwind_grid,vwind_grid,windspd_grid,winddir_grid,&
                &        Qsi_grid,Qli_grid,prec_grid,Tsfc,Qle,Qh,Qe,Qc,Qm,Qf,&
                &        e_balance,snow_depth,xro_snow,swe_depth,ro_nsnow,&
                &        runoff,rain,sprec,sum_prec,sum_runoff,w_balance,&
                &        snow_d,topo_land,wbal_qsubl,sum_sprec,wbal_salt,&
                &        wbal_susp,ro_snow_grid,sum_Qcs,canopy_int,Qcs,&
                &        iyear,imonth,iday,xhour,undef,deltax,xmn,ymn,&
                &        wbal_subgrid,canopy_unload,sum_qsubl,sum_trans,&
                &        sum_unload,sum_glacmelt,glacier_melt,swemelt,&
                &        sfc_pressure,sum_swemelt,albedo,nrecs_max,&
                &        icorr_factor_loop,swesublim,vegtype,iter_start,&
                &        seaice_run,print_inc,cloud_frac_grid,&
                &        output_path_wo_assim,output_path_wi_assim,print_var,&
                &        print_outvars,Qsubl_depth,Qsalt,Qsusp)
        endif

        if (print_user .eq. 1.0 .and. print_parallel .eq. 1.0) then
           CALL OUTPUTS_PARALLEL(nx,l_ny,iter,l_Tair_grid,l_rh_grid,&
                &        l_uwind_grid,l_vwind_grid,l_windspd_grid,l_winddir_grid,&
                &        l_Qsi_grid,l_Qli_grid,l_prec_grid,l_Tsfc,l_Qle,l_Qh,l_Qe,l_Qc,l_Qm,l_Qf,&
                &        l_e_balance,l_snow_depth,l_xro_snow,l_swe_depth,l_ro_nsnow,&
                &        l_runoff,l_rain,l_sprec,l_sum_prec,l_sum_runoff,l_w_balance,&
                &        l_snow_d,topo_land,l_wbal_qsubl,l_sum_sprec,l_wbal_salt,&
                &        l_wbal_susp,l_ro_snow_grid,l_sum_Qcs,l_canopy_int,l_Qcs,&
                &        iyear,imonth,iday,xhour,undef,deltax,xmn,ymn,&
                &        l_wbal_subgrid,l_canopy_unload,l_sum_qsubl,l_sum_trans,&
                &        l_sum_unload,l_sum_glacmelt,l_glacier_melt,l_swemelt,&
                &        l_sfc_pressure,l_sum_swemelt,l_albedo,nrecs_max,&
                &        icorr_factor_loop,l_swesublim,l_vegtype,iter_start,&
                &        seaice_run,print_inc,l_cloud_frac_grid,&
                &        output_path_wo_assim,output_path_wi_assim,print_var,&
                &        print_outvars,l_Qsubl_depth,l_Qsalt,l_Qsusp)
        endif

        ! c For multi-year simulations, sometimes it is desirable to zero
        ! c   out the snow cover arrays on a certain summer date, to prevent
        ! c   glaciers from forming.
        if (imonth .eq. iclear_mn .and. iday .eq. iclear_dy .and. &
             &      xhour.eq.xclear_hr) then

           if (seaice_run .eq. 4.0) then
              CALL ZERO_SNOW_SEAICE_4(nx,ny,l_snow_depth,l_ro_snow_grid,&
                   &          ro_snow,l_swe_depth,l_swe_depth_old,l_canopy_int_old,l_KK,&
                   &          l_sum_swemelt,l_tslsnowfall,l_snod_layer,l_swed_layer,&
                   &          l_ro_layer,l_T_old,l_sum_sprec,multilayer_snowpack,&
                   &          tsls_threshold,iyear,l_diam_layer,output_path_wo_assim)
           else
              CALL ZERO_SNOW(nx,ny,l_snow_depth,l_ro_snow_grid,ro_snow,&
                   &          l_swe_depth,l_swe_depth_old,l_canopy_int_old,l_KK,l_sum_swemelt,&
                   &          l_tslsnowfall,l_snod_layer,l_swed_layer,l_ro_layer,l_T_old,&
                   &          l_sum_sprec,multilayer_snowpack,tsls_threshold,&
                   &          l_sum_trans)
           endif

        endif

        !c Save the history restart information.
        if (ihrestart_flag .ge. -1 .and. me == 0) then
           if (mod(iter,ihrestart_inc) .eq. 0 &
                &        .or. iter .eq. max_iter) then

              !Gather data on image 0

              CALL HRESTART_SAVE(nx,ny,iter,snow_d,snow_depth,&
                   &          canopy_int,soft_snow_d,ro_snow_grid,swe_depth,&
                   &          ro_soft_snow_old,snow_d_init,swe_depth_old,&
                   &          canopy_int_old,topo,sum_sprec,icorr_factor_loop,&
                   &          max_iter)
           endif
        endif
        if (me == 0) then
            call System_Clock(output_end)
            total_output_time = total_output_time + (output_end - output_start)
        endif

        sync all
     enddo

  enddo
  if (me == 0) then
    call System_Clock(master_end)
    call System_Clock(noninit_end)
  endif

  !c Print a banner when the model run is finished.
  if (me == 0) then
    print*, "ReadParam Time : ",  ((readpar_end-readpar_start) / real(COUNT_RATE))
    print*, "Preproc Time : ",  ((preproc_end-preproc_start) / real(COUNT_RATE))
    print*, "Distribute Time : ",  ((total_distribute_time) / real(COUNT_RATE))
    print*, "Micromet Time : ",  ((total_micromet_time) / real(COUNT_RATE))
    print*, "Enbal Time : ",  ((total_enbal_time) / real(COUNT_RATE))
    print*, "Snowpack Time : ",  ((total_snowpack_time) / real(COUNT_RATE))
    print*, "SnowTran Time : ",  ((total_snowtran_time) / real(COUNT_RATE))
    print*, "Gather Time : ",  ((total_gather_time) / real(COUNT_RATE))
    print*, "Output Time : ",  ((total_output_time) / real(COUNT_RATE))
    print*, "Master Noninit Time : ", ((noninit_end - noninit_start) / real(COUNT_RATE))
    print*, "Master Total Time : ", ((master_end - master_start) / real(COUNT_RATE))
  endif
  print *
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *,&
       & '                 The SnowModel Run Has Finished                '
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *

  !stop
end program snowmodel_main

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
