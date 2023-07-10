module snowmodel_vars
  use snowmodel_inc
  implicit none

  !include 'snowmodel.inc'
!c Misc variables.
  character*100 snowmodel_dot_par_fname


!c SnowTran-3D variables.
  integer max_iter

  integer i,j,iter,nx,ny

  real ro_snow,ro_water,ro_air,gravity,vonKarman,snow_z0
  real deltax,deltay,dt,curvewt,utau_t_flag
  real fetch,xmu,C_z,h_const,wind_min,windspd_flag
  real Up_const,dz_susp,ztop_susp,fall_vel,Ur_const
  real pi,twopio360,bc_flag,topoflag,Utau_t_const
  real ht_windobs,ht_rhobs,slopewt,bs_flag,twolayer_flag
  real subgrid_flag,curve_len_scale
  real snow_d_init_const,const_veg_flag,winddir_flag
  real windspd_min,tabler_dir,slope_adjust

  real, allocatable :: topo_land(:,:)
!  real, allocatable :: tabler_nn(:,:)
  real, allocatable :: l_tabler_nn(:,:)[:]
  real, allocatable :: l_tabler_nn_orig(:,:)[:]
!  real, allocatable ::  tabler_ss(:,:)
  real, allocatable ::  l_tabler_ss(:,:)[:]
  real, allocatable ::  l_tabler_ss_orig(:,:)[:]
!  real, allocatable ::  tabler_ee(:,:)
  real, allocatable ::  l_tabler_ee(:,:)[:]
  real, allocatable ::  l_tabler_ee_orig(:,:)[:]
!  real, allocatable ::  tabler_ww(:,:)
  real, allocatable ::  l_tabler_ww(:,:)[:]
  real, allocatable ::  l_tabler_ww_orig(:,:)[:]
!  real, allocatable ::  tabler_ne(:,:)
  real, allocatable ::  l_tabler_ne(:,:)[:]
  real, allocatable ::  l_tabler_ne_orig(:,:)[:]
!  real, allocatable ::  tabler_se(:,:)
  real, allocatable ::  l_tabler_se(:,:)[:]
  real, allocatable ::  l_tabler_se_orig(:,:)[:]
!  real, allocatable ::  tabler_sw(:,:)
  real, allocatable ::  l_tabler_sw(:,:)[:]
  real, allocatable ::  l_tabler_sw_orig(:,:)[:]
!  real, allocatable ::  tabler_nw(:,:)
  real, allocatable ::  l_tabler_nw(:,:)[:]
  real, allocatable ::  l_tabler_nw_orig(:,:)[:]
  real, allocatable ::  topo(:,:)
  real, allocatable :: l_topo_tmp(:,:)[:]
  real, allocatable :: l_vegtype(:,:)[:]
  real, allocatable ::  vegtype(:,:)
  real, allocatable ::  l_snow_d_tabler(:,:)[:]
  real, allocatable ::  l_dh_dep(:,:)[:]
  real, allocatable ::  l_snow_d_tmp(:,:)[:]

  real, allocatable :: l_uwind_grid(:,:)[:],l_vwind_grid(:,:)[:]
  real, allocatable ::  uwind_grid(:,:),vwind_grid(:,:)
  real, allocatable :: l_windspd_grid(:,:)[:], l_winddir_grid(:,:)[:]
  real, allocatable ::  windspd_grid(:,:),winddir_grid(:,:)
  real, allocatable ::  tair_grid(:,:),sprec(:,:)
  real, allocatable :: l_tair_grid(:,:)[:]
  real, allocatable :: l_sprec(:,:)[:]
  real, allocatable ::  rh_grid(:,:), curvature(:,:)
  real, allocatable :: l_rh_grid(:,:)[:],l_curvature(:,:)[:]

  real, allocatable :: l_slope_az(:,:)[:],l_terrain_slope(:,:)[:]

  real, allocatable ::  slope_az(:,:)
  real, allocatable ::  terrain_slope(:,:)

  ! integer,allocatable :: index_ue(:,:),index_uw(:,:)
  ! integer,allocatable :: index_vn(:,:),index_vs(:,:)
  integer index_ue(ny_max,2*nx_max+1),index_uw(ny_max,2*nx_max+1)
  integer index_vn(nx_max,2*ny_max+1),index_vs(nx_max,2*ny_max+1)

  real, allocatable :: l_snow_d(:,:)[:]
  real, allocatable ::  snow_d(:,:)
  real, allocatable :: l_snow_d_init(:,:)[:]
  real, allocatable ::  snow_d_init(:,:)
  real, allocatable ::  Utau(:,:)
  real, allocatable ::  Utau_t(:,:)
  real, allocatable ::  l_Utau(:,:)[:]
  real, allocatable ::  l_Utau_t(:,:)[:]
  real, allocatable ::  z_0(:,:)
  real, allocatable ::  h_star(:,:)
  real, allocatable ::  conc_salt(:,:)
  real, allocatable ::  l_conc_salt(:,:)

!  real, allocatable ::  Qsalt_max(:,:)
  real, allocatable ::  l_Qsalt_max(:,:)[:]
!  real, allocatable ::  Qsalt_maxu(:,:),Qsalt_maxv(:,:)
  real, allocatable ::  l_Qsalt_maxu(:,:)[:],l_Qsalt_maxv(:,:)[:]
  real, allocatable ::  Qsalt(:,:)
  real, allocatable ::  l_Qsalt(:,:)[:]
!  real, allocatable ::  Qsalt_u(:,:)
  real, allocatable ::  Qsalt_v(:,:)
  real, allocatable ::  l_Qsalt_v(:,:)[:]
  real, allocatable ::  l_Qsalt_u(:,:)[:]
!  real, allocatable ::  dh_salt(:,:)
  real, allocatable ::  l_dh_salt(:,:)[:]
!  real, allocatable ::  dh_salt_u(:,:),dh_salt_v(:,:)
  real, allocatable ::  l_dh_salt_u(:,:)[:],l_dh_salt_v(:,:)[:]

  real, allocatable ::  Qsusp(:,:)
  real, allocatable ::  l_Qsusp(:,:)[:]
!  real, allocatable ::  Qsusp_u(:,:)
  real, allocatable ::  Qsusp_v(:,:)
  real, allocatable ::  l_Qsusp_u(:,:)[:],l_Qsusp_v(:,:)[:]
!  real, allocatable ::  dh_susp(:,:)
  real, allocatable ::  l_dh_susp(:,:)[:]
!  real, allocatable ::  dh_susp_u(:,:),dh_susp_v(:,:)
  real, allocatable ::  l_dh_susp_u(:,:)[:],l_dh_susp_v(:,:)[:]

!  real, allocatable ::  Qsubl(:,:)
  real, allocatable ::  l_Qsubl(:,:)[:]
  real, allocatable ::  Qsubl_depth(:,:)
  real, allocatable ::  l_Qsubl_depth(:,:)[:]

  real, allocatable :: l_sum_sprec(:,:)[:]
  real, allocatable ::  sum_sprec(:,:)

  real, allocatable ::  wbal_qsubl(:,:)
  real, allocatable ::  l_wbal_qsubl(:,:)[:]
  real, allocatable ::  wbal_salt(:,:)
  real, allocatable ::  l_wbal_salt(:,:)[:]
  real, allocatable ::  wbal_susp(:,:)
  real, allocatable ::  l_wbal_susp(:,:)[:]
  real, allocatable ::  wbal_subgrid(:,:)
  real, allocatable ::  l_wbal_subgrid(:,:)[:]
  real, allocatable ::  sum_qsubl(:,:)
  real, allocatable ::  l_sum_qsubl(:,:)[:]
  real, allocatable ::  sum_trans(:,:)
  real, allocatable ::  l_sum_trans(:,:)[:]
  real, allocatable ::  l_soft_snow_d(:,:)[:]
  real, allocatable ::  soft_snow_d(:,:)
  real, allocatable ::  ro_soft_snow(:,:)
  real, allocatable ::  ro_soft_snow_old(:,:)
  real, allocatable ::  l_ro_soft_snow(:,:)[:]
  real, allocatable ::  l_ro_soft_snow_old(:,:)[:]
  real, allocatable ::  l_dh_subgrid(:,:)[:]
  real, allocatable ::  l_snow_d_dep(:,:)[:]
  real, allocatable ::  l_soft_snow_d_dep(:,:)[:]

  real vegsnowdepth(nvegtypes)
  real, allocatable ::  l_veg_z0(:,:)[:]
  real, allocatable ::  veg_z0(:,:)
  real, allocatable ::  l_vegsnowd_xy(:,:)[:]
  real, allocatable ::  vegsnowd_xy(:,:)
  integer iveg_ht_flag

  real run_micromet,run_enbal,run_snowpack,run_snowtran

  character*80 topoveg_fname,met_input_fname,topo_ascii_fname,&
  &  veg_ascii_fname

  real curve_lg_scale_flag
  real, allocatable ::  curve_wt_lg(:,:)
  real, allocatable ::  l_curve_wt_lg(:,:)[:]

  character*80 tabler_sfc_path_name
  real Tabler_1_flag,Tabler_2_flag
!c End SnowTran-3D variables.

!c MicroMet variables.
  double precision xmn  ! center x coords of lower left grid cell
  double precision ymn  ! center y coords of lower left grid cell
  !double precision xg_line(nx_max,ny_max),yg_line(nx_max,ny_max)
  double precision, allocatable :: xg_line(:,:),yg_line(:,:)

  real dn                  ! average observation spacing

  real, allocatable :: l_Qsi_grid(:,:)[:]
  real, allocatable ::  Qsi_grid(:,:)    ! output
  real, allocatable :: l_Qli_grid(:,:)[:]
  real, allocatable ::  Qli_grid(:,:)    ! output
  real, allocatable :: l_prec_grid(:,:)[:]
  real, allocatable ::  prec_grid(:,:)   ! output
  real, allocatable ::  xlat_grid(:,:)   ! lat (dec deg) of cell centers
  real, allocatable ::  xlon_grid(:,:)   ! lon (dec deg) of cell centers

  integer iyear_init     ! model start year
  integer imonth_init    ! model start month
  integer iday_init      ! model start day
  real xhour_init        ! model start hour
  real xlat      ! approx. latitude of domain center, decimal deg

  real undef       ! undefined value
  integer ifill    ! flag (=1) forces a value in every cell
  integer iobsint  ! flag (=1) use dn value from .par file

  integer i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,&
       &  i_prec_flag,i_longwave_flag,isingle_stn_flag,igrads_metfile,&
       &  lapse_rate_user_flag,iprecip_lapse_rate_user_flag,n_stns_used,&
       &  lat_solar_flag,ihrestart_inc,iter_start,ihrestart_flag,&
       &  iprecip_scheme

  real xhour,ascii_topoveg,use_shortwave_obs,gap_frac,&
       &  use_longwave_obs,use_sfc_pressure_obs,calc_subcanopy_met,&
       &  cloud_frac_factor,barnes_lg_domain,UTC_flag,check_met_data,&
       &  snowmodel_line_flag,wind_lapse_rate

  integer, parameter :: nftypes = 5
  real forest_LAI(nftypes)

  integer iyear,imonth,iday
  !integer k_stn(nx_max,ny_max,9)
  integer, allocatable :: k_stn(:,:,:)!Last dim is 9

  real, allocatable ::  l_cf_precip(:,:)[:]
  real, allocatable ::  cf_precip(:,:)
  real cf_precip_flag

  real, allocatable ::  l_cloud_frac_grid(:,:)[:]
  real, allocatable ::  cloud_frac_grid(:,:)
  real snowfall_frac
!c End MicroMet variables.

!c EnBal variables.
  integer icond_flag

  real, allocatable :: l_albedo(:,:)[:]
  real, allocatable :: l_swe_depth(:,:)[:]
  real, allocatable :: l_sfc_pressure(:,:)[:]
  real, allocatable :: l_Qm(:,:)[:]
  real, allocatable :: l_Qc(:,:)[:]
  real, allocatable :: l_Qe(:,:)[:]
  real, allocatable :: l_Qf(:,:)[:]
  real, allocatable :: l_Tsfc(:,:)[:]
  real, allocatable :: l_Qle(:,:)[:]
  real, allocatable :: l_Qh(:,:)[:]
  real, allocatable :: l_e_balance(:,:)[:]

  real, allocatable :: Tsfc(:,:),Qle(:,:),&
       &  Qh(:,:),Qe(:,:),Qc(:,:),&
       &  Qm(:,:),e_balance(:,:),Qf(:,:),&
       &  swe_depth(:,:),sfc_pressure(:,:),&
       &  albedo(:,:)

  real albedo_snow_forest,albedo_snow_clearing,albedo_glacier
!c End EnBal variables.

  !c SnowPack variables.
  real, allocatable :: l_ro_nsnow(:,:)[:]
  real, allocatable :: l_snow_depth(:,:)[:]
  real, allocatable :: l_rain(:,:)[:]
  real, allocatable :: l_sum_runoff(:,:)[:]
  real, allocatable :: l_sum_prec(:,:)[:]
  real, allocatable :: l_runoff(:,:)[:]
  real, allocatable :: l_xro_snow(:,:)[:]
  real, allocatable :: l_canopy_int(:,:)[:]
  real, allocatable :: l_canopy_int_old(:,:)[:]
  real, allocatable :: l_canopy_unload(:,:)[:]
  real, allocatable :: l_d_canopy_int(:,:)[:]
  real, allocatable :: l_sum_d_canopy_int(:,:)[:]
  real, allocatable :: l_glacier_melt(:,:)[:]
  real, allocatable :: l_sum_glacmelt(:,:)[:]
  real, allocatable :: l_swemelt(:,:)[:]
  real, allocatable :: l_sum_swemelt(:,:)[:]
  real, allocatable :: l_swesublim(:,:)[:]
  real, allocatable :: l_swe_depth_old(:,:)[:]
  real, allocatable :: l_w_balance(:,:)[:]
  real, allocatable :: l_ro_snow_grid(:,:)[:]
  real, allocatable :: l_sum_unload(:,:)[:]
  real, allocatable :: l_Qcs(:,:)[:]
  real, allocatable :: l_sum_Qcs(:,:)[:]
  real, allocatable :: l_sum_sfcsublim(:,:)[:]
  real, allocatable :: l_windspd_2m_grid(:,:)[:]

  real, allocatable :: ro_nsnow(:,:),w_balance(:,:),&
       &  runoff(:,:),rain(:,:),&
       &  sum_prec(:,:),sum_runoff(:,:),&
       &  xro_snow(:,:),ro_snow_grid(:,:),&
       &  sum_Qcs(:,:),canopy_int(:,:),&
       &  Qcs(:,:),canopy_unload(:,:),&
       &  snow_depth(:,:),glacier_melt(:,:),&
       &  sum_unload(:,:),sum_glacmelt(:,:),&
       &  swemelt(:,:),&
       &  sum_d_canopy_int(:,:),sum_sfcsublim(:,:),&
       &  sum_swemelt(:,:),swesublim(:,:),&
       &  swe_depth_old(:,:),canopy_int_old(:,:),&
       &  windspd_2m_grid(:,:)

  real sfc_sublim_flag

  integer max_layers,multilayer_snowpack,k
  integer, allocatable :: l_KK(:,:)[:]
  integer, allocatable :: KK(:,:)
  integer, allocatable :: melt_flag(:,:,:)!nz_max


  real ro_snowmax,tsls_threshold,dz_snow_min
  real, allocatable :: l_tslsnowfall(:,:)[:]
  real, allocatable ::  tslsnowfall(:,:)
  real, allocatable :: l_change_layer(:,:)[:]
  real, allocatable ::  change_layer(:,:)
  real, allocatable ::  snod_layer(:,:,:)
  real, allocatable ::  l_snod_layer(:,:,:)[:]
  real, allocatable ::  swed_layer(:,:,:)
  real, allocatable ::  l_swed_layer(:,:,:)[:]
  real, allocatable ::  ro_layer(:,:,:)
  real, allocatable ::  l_ro_layer(:,:,:)[:]
  real, allocatable ::  T_old(:,:,:)
  real, allocatable ::  l_T_old(:,:,:)[:]
  real, allocatable ::  gamma(:,:,:)
  real, allocatable ::  l_gamma(:,:,:)[:]
  real, allocatable ::  diam_layer(:,:,:)
  real, allocatable ::  l_diam_layer(:,:,:)[:]
  real, allocatable ::  flux_layer(:,:,:)
  real, allocatable ::  l_flux_layer(:,:,:)[:]

  integer iclear_mn,iclear_dy,izero_snow_date
  real xclear_hr
!c End SnowPack variables.

!c SeaIce variables.
  real seaice_run
  real, allocatable :: seaice_conc(:,:,:)
!c End SeaIce variables.

!c Data assimilaion (precipitation and melt) factor variables.
  real corr_factor(nx_max,ny_max,max_obs_dates+1)
  integer icorr_factor_index(max_time_steps)

  integer icorr_factor_loop,irun_data_assim,nobs_dates,&
  &  i_dataassim_loop,i_corr_start

!c End data assimilaion (precipitation and melt) factor variables.

!c Print output variables.
  real print_micromet,print_enbal,print_snowpack,print_snowtran,&
       &  print_multilayer
  real print_user,print_inc
  double precision nrecs_max

  character*80 micromet_output_fname
  character*80 enbal_output_fname
  character*80 snowtran_output_fname
  character*80 snowpack_output_fname
  character*80 multilayer_output_fname
  character*80 output_path_wo_assim,output_path_wi_assim

  character*1 print_var(n_print_vars)
  character*4 print_outvars(n_print_vars)
!c End print output variables.
contains

  subroutine allocate_arrays(nx,ny,me)
    implicit none

    integer, intent(in) :: nx,ny,me
    integer, save :: first = 1

    real,dimension(nx,ny) :: zero_real

    zero_real = 0.0

    if(first == 0) return
    first = 0

    allocate(topo_land(nx,ny),source=zero_real)
!    allocate(tabler_nn(nx,ny),source=zero_real)
!    allocate( tabler_ss(nx,ny),source=zero_real)
!    allocate( tabler_ee(nx,ny),source=zero_real)
!    allocate( tabler_ww(nx,ny),source=zero_real)
!    allocate( tabler_ne(nx,ny),source=zero_real)
!    allocate( tabler_se(nx,ny),source=zero_real)
!    allocate( tabler_sw(nx,ny),source=zero_real)
!    allocate( tabler_nw(nx,ny),source=zero_real)
    allocate( topo(nx,ny),source=zero_real)
    allocate( vegtype(nx,ny),source=zero_real)
!    allocate( uwind_grid(nx,ny),source=zero_real)
    allocate( vwind_grid(nx,ny),source=zero_real)
!    allocate( windspd_grid(nx,ny),winddir_grid(nx,ny),source=zero_real)
    allocate( tair_grid(nx,ny),source=zero_real)
!    allocate( sprec(nx,ny),source=zero_real)
!    allocate( rh_grid(nx,ny),source=zero_real)
    allocate(curvature(nx,ny),source=zero_real)

    allocate( slope_az(nx,ny),source=zero_real)
    allocate( terrain_slope(nx,ny),source=zero_real)

    ! The following 2 variables have nx and ny swapped.
    ! Nx (cols) is also doubled
    !integer,allocatable :: index_ue(:,:),index_uw(:,:)
    ! The following 2 variables have ny (cols) doubled
    !integer,allocatable :: index_vn(:,:),index_vs(:,:)

    allocate( snow_d(nx,ny),source=zero_real)

    allocate( snow_d_init(nx,ny),source=zero_real)
    allocate( Utau(nx,ny),source=zero_real)
    allocate( Utau_t(nx,ny),source=zero_real)
    !allocate( z_0(nx,ny),source=zero_real)
    !allocate( h_star(nx,ny),source=zero_real)
    allocate( conc_salt(nx,ny),source=zero_real)

!    allocate( Qsalt_max(nx,ny),source=zero_real)
!    allocate( Qsalt_maxu(nx,ny),source=zero_real)
!    allocate( Qsalt_maxv(nx,ny),source=zero_real)
!    allocate( Qsalt(nx,ny),source=zero_real)
!    allocate( Qsalt_u(nx,ny),source=zero_real)
    allocate( Qsalt_v(nx,ny),source=zero_real)
!    allocate( dh_salt(nx,ny),source=zero_real)
!    allocate( dh_salt_u(nx,ny),source=zero_real)
!    allocate( dh_salt_v(nx,ny),source=zero_real)

!    allocate( Qsusp(nx,ny),source=zero_real)
!    allocate( Qsusp_u(nx,ny),source=zero_real)
    allocate( Qsusp_v(nx,ny),source=zero_real)
!    allocate( dh_susp(nx,ny),source=zero_real)
!    allocate( dh_susp_u(nx,ny),dh_susp_v(nx,ny),source=zero_real)

!    allocate( Qsubl(nx,ny),source=zero_real)
!    allocate( Qsubl_depth(nx,ny),source=zero_real)

    allocate( sum_sprec(nx,ny),source=zero_real)

!    allocate( wbal_qsubl(nx,ny),source=zero_real)
!    allocate( wbal_salt(nx,ny),source=zero_real)
!    allocate( wbal_susp(nx,ny),source=zero_real)
!    allocate( wbal_subgrid(nx,ny),source=zero_real)
    allocate( sum_qsubl(nx,ny),source=zero_real)
    allocate( sum_trans(nx,ny),source=zero_real)

    allocate( soft_snow_d(nx,ny),source=zero_real)
    allocate( ro_soft_snow(nx,ny),source=zero_real)
    allocate( ro_soft_snow_old(nx,ny),source=zero_real)
    allocate( curve_wt_lg(nx,ny),source=zero_real)

    !allocate( vegsnowdepth(nvegtypes))!vegsnowdepth(nvegtypes)
    allocate( veg_z0(nx,ny),source=zero_real)
    allocate( vegsnowd_xy(nx,ny),source=zero_real)

!   ADDED LINE TO ALLOCATE GLOBAL VARIABLE ONLY TO MAIN PROCESSOR
    if (me == 0) then
      allocate( Qsalt(nx,ny),source=zero_real)
      allocate( Qsusp(nx,ny),source=zero_real)
      allocate( rh_grid(nx,ny),source=zero_real)
      allocate( sprec(nx,ny),source=zero_real)
      allocate( wbal_qsubl(nx,ny),source=zero_real)
      allocate( wbal_salt(nx,ny),source=zero_real)
      allocate( wbal_susp(nx,ny),source=zero_real)
      allocate( wbal_subgrid(nx,ny),source=zero_real)
      allocate( uwind_grid(nx,ny),source=zero_real)
      allocate( windspd_grid(nx,ny),winddir_grid(nx,ny),source=zero_real)
      allocate( w_balance(nx,ny),source=zero_real)
      allocate( ro_nsnow(nx,ny),source=zero_real)
      allocate( prec_grid(nx,ny),source=zero_real)
      allocate( Qcs(nx,ny),source=zero_real)
      allocate( runoff(nx,ny),source=zero_real)
      allocate( glacier_melt(nx,ny),source=zero_real)
      allocate( swesublim(nx,ny),source=zero_real)
      allocate( canopy_unload(nx,ny),source=zero_real)
      allocate( Qsubl_depth(nx,ny),source=zero_real)
    endif

    allocate (xg_line(nx,ny),yg_line(nx,ny))
    xg_line = 0.0
    yg_line = 0.0
    allocate( Qsi_grid(nx,ny),source=zero_real)    ! output
    allocate( Qli_grid(nx,ny),source=zero_real)    ! output

!    allocate( prec_grid(nx,ny),source=zero_real)   ! output
    allocate( xlat_grid(nx,ny),source=zero_real)   ! lat (dec deg) of cell centers
    allocate( xlon_grid(nx,ny),source=zero_real)   ! lon (dec deg) of cell centers

    allocate( cf_precip(nx,ny),source=zero_real)

    allocate(k_stn(nx,ny,9))!Last dim is 9
    k_stn = 0

    allocate( cloud_frac_grid(nx,ny),source=zero_real)

    allocate(Tsfc(nx,ny),Qle(nx,ny),&
         &  Qh(nx,ny),Qe(nx,ny),Qc(nx,ny),&
         &  Qm(nx,ny),e_balance(nx,ny),Qf(nx,ny),&
         &  swe_depth(nx,ny),sfc_pressure(nx,ny),&
         &  albedo(nx,ny),source=zero_real)

    allocate(rain(nx,ny),&
         &  sum_prec(nx,ny),sum_runoff(nx,ny),&
         &  xro_snow(nx,ny),ro_snow_grid(nx,ny),&
         &  sum_Qcs(nx,ny),canopy_int(nx,ny),&
         &  snow_depth(nx,ny),&
         &  sum_unload(nx,ny),sum_glacmelt(nx,ny),&
         &  swemelt(nx,ny),&
         &  sum_d_canopy_int(nx,ny),sum_sfcsublim(nx,ny),&
         &  sum_swemelt(nx,ny),&
         &  swe_depth_old(nx,ny),canopy_int_old(nx,ny),&
         &  windspd_2m_grid(nx,ny),source=zero_real)

    allocate(KK(nx,ny))
    KK = 0
    allocate(melt_flag(nx,ny,nz_max))!nz_max
    melt_flag = 0

    allocate( tslsnowfall(nx,ny),source=zero_real)

    allocate( change_layer(nx,ny),source=zero_real)
    allocate( snod_layer(nx,ny,nz_max))
    snod_layer = 0.0
    allocate( swed_layer(nx,ny,nz_max))
    swed_layer = 0.0
    allocate( ro_layer(nx,ny,nz_max))
    ro_layer = 0.0
    allocate( T_old(nx,ny,nz_max))
    T_old = 0.0
    allocate( gamma(nx,ny,nz_max))
    gamma = 0.0
    allocate( diam_layer(nx,ny,nz_max))
    diam_layer = 0.0
    allocate( flux_layer(nx,ny,nz_max))
    flux_layer = 0.0
    allocate(seaice_conc(nx,ny,nz_max))
    seaice_conc = 0.0

  end subroutine allocate_arrays

end module snowmodel_vars
