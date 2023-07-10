module caf_module
  implicit none

  public

  integer :: me, np
  integer :: rem, l_ny, max_l_ny, max_l_nx
  integer, allocatable :: partition_ny(:)
  integer, allocatable :: prefix_sum(:)

  !Used for halo exchange
  real, allocatable :: left_rcv(:)[:], right_rcv(:)[:]

  interface single_gather

     module procedure single_gather_integer,single_gather_real

  end interface single_gather

  interface single_scatter

     module procedure single_scatter_integer,single_scatter_real

  end interface single_scatter

contains

  subroutine caf_partitioning()
    use snowmodel_inc
    use snowmodel_vars
    implicit none

    me = this_image() - 1
    np = num_images()

    allocate(partition_ny(0:np-1))
    allocate(prefix_sum(0:np-1))
    partition_ny = 0
    prefix_sum = 0

    rem = mod(ny,np)
    l_ny = (ny-rem)/np
    max_l_ny = l_ny

    if(rem > 0) max_l_ny = max_l_ny + 1
    if(me < rem) l_ny = l_ny + 1

    call allocation()
    call partition()
!    call distribute()

    sync all

  end subroutine caf_partitioning

  subroutine partition()
    use snowmodel_vars
    implicit none
    integer :: ii,idx

    do idx = 0, np-1
       partition_ny(idx) = (ny-rem)/np
       if(idx < rem) then
          partition_ny(idx) = partition_ny(idx) + 1
       end if
    end do

    do ii = 1, np-1
       prefix_sum(ii) = partition_ny(ii-1) + prefix_sum(ii-1)
    end do

  end subroutine partition

  subroutine allocation()
    use snowmodel_inc
    use snowmodel_vars
    implicit none

    real,dimension(nx,max_l_ny) :: zero_real
    real,dimension(nx,max_l_ny,nz_max) :: zero_real_three_d

    max_l_nx = nx

    zero_real = 0.0
    zero_real_three_d = 0.0

    !Not coarrays but used in the parallel version
    allocate(l_conc_salt(nx,max_l_ny),source=zero_real)
    allocate( z_0(nx,max_l_ny),source=zero_real)
    allocate( h_star(nx,max_l_ny),source=zero_real)

    allocate( l_wbal_subgrid(nx,max_l_ny)[*],source=zero_real)
    allocate( l_wbal_susp(nx,max_l_ny)[*],source=zero_real)
    allocate( l_wbal_salt(nx,max_l_ny)[*],source=zero_real)
    allocate( l_wbal_qsubl(nx,max_l_ny)[*],source=zero_real)
    allocate( l_ro_soft_snow(nx,max_l_ny)[*],source=zero_real)
    allocate( l_ro_soft_snow_old(nx,max_l_ny)[*],source=zero_real)
    allocate( l_curve_wt_lg(nx,max_l_ny)[*],source=zero_real)
    allocate( l_dh_salt_u(nx,max_l_ny)[*],source=zero_real)
    allocate( l_dh_salt_v(nx,max_l_ny)[*],source=zero_real)
    allocate( l_dh_salt(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Utau(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Utau_t(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsusp(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsusp_u(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsusp_v(nx,max_l_ny)[*],source=zero_real)
    allocate(l_dh_susp(nx,max_l_ny)[*],source=zero_real)
    allocate(l_dh_susp_u(nx,max_l_ny)[*],source=zero_real)
    allocate(l_dh_susp_v(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsubl(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsubl_depth(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsalt(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsalt_u(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsalt_v(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsalt_max(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsalt_maxu(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsalt_maxv(nx,max_l_ny)[*],source=zero_real)
    allocate(l_vegtype(nx,max_l_ny)[*],source=zero_real)
    allocate(l_vegsnowd_xy(nx,max_l_ny)[*],source=zero_real)
    allocate(l_veg_z0(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tair_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_cloud_frac_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_cf_precip(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sprec(nx,max_l_ny)[*],source=zero_real)
    allocate(l_rh_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_windspd_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_winddir_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_soft_snow_d(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qsi_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qli_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_prec_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_albedo(nx,max_l_ny)[*],source=zero_real)
    allocate(l_curvature(nx,max_l_ny)[*],source=zero_real)
    allocate(l_slope_az(nx,max_l_ny)[*],source=zero_real)
    allocate(l_terrain_slope(nx,max_l_ny)[*],source=zero_real)
    allocate(l_swe_depth(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sfc_pressure(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qm(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qh(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qc(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qe(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qf(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qle(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Tsfc(nx,max_l_ny)[*],source=zero_real)
    allocate(l_snow_d(nx,max_l_ny)[*],source=zero_real)
    allocate(l_snow_d_init(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tslsnowfall(nx,max_l_ny)[*],source=zero_real)
    allocate(l_change_layer(nx,max_l_ny)[*],source=zero_real)
    allocate(l_KK(nx,max_l_ny)[*])
    l_kk = 0
    allocate(l_ro_nsnow(nx,max_l_ny)[*],source=zero_real)
    allocate(l_snow_depth(nx,max_l_ny)[*],source=zero_real)
    allocate(l_rain(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_runoff(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_prec(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_qsubl(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_trans(nx,max_l_ny)[*],source=zero_real)
    allocate(l_runoff(nx,max_l_ny)[*],source=zero_real)
    allocate(l_xro_snow(nx,max_l_ny)[*],source=zero_real)
    allocate(l_canopy_int(nx,max_l_ny)[*],source=zero_real)
    allocate(l_canopy_int_old(nx,max_l_ny)[*],source=zero_real)
    allocate(l_canopy_unload(nx,max_l_ny)[*],source=zero_real)
    allocate(l_d_canopy_int(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_d_canopy_int(nx,max_l_ny)[*],source=zero_real)
    allocate(l_glacier_melt(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_glacmelt(nx,max_l_ny)[*],source=zero_real)
    allocate(l_swemelt(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_swemelt(nx,max_l_ny)[*],source=zero_real)
    allocate(l_swesublim(nx,max_l_ny)[*],source=zero_real)
    allocate(l_swe_depth_old(nx,max_l_ny)[*],source=zero_real)
    allocate(l_e_balance(nx,max_l_ny)[*],source=zero_real)
    allocate(l_w_balance(nx,max_l_ny)[*],source=zero_real)
    allocate(l_ro_snow_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_unload(nx,max_l_ny)[*],source=zero_real)
    allocate(l_Qcs(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_Qcs(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_sfcsublim(nx,max_l_ny)[*],source=zero_real)
    allocate(l_windspd_2m_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_sum_sprec(nx,max_l_ny)[*],source=zero_real)
    allocate(l_uwind_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_vwind_grid(nx,max_l_ny)[*],source=zero_real)

    allocate(l_dh_subgrid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_snow_d_dep(nx,max_l_ny)[*],source=zero_real)
    allocate(l_soft_snow_d_dep(nx,max_l_ny)[*],source=zero_real)

    allocate(l_tabler_nn(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_nn_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ee(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ee_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ss(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ss_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ww(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ww_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ne(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_ne_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_se(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_se_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_nw(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_nw_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_sw(nx,max_l_ny)[*],source=zero_real)
    allocate(l_tabler_sw_orig(nx,max_l_ny)[*],source=zero_real)
    allocate(l_topo_tmp(nx,max_l_ny)[*],source=zero_real)
    allocate(l_topo_land(nx,max_l_ny)[*],source=zero_real)
    allocate(l_snow_d_tabler(nx,max_l_ny)[*],source=zero_real)
    allocate(l_dh_dep(nx,max_l_ny)[*],source=zero_real)
    allocate(l_snow_d_tmp(nx,max_l_ny)[*],source=zero_real)
    allocate(l_xlat_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_xlon_grid(nx,max_l_ny)[*],source=zero_real)
    allocate(l_xg_line(nx,max_l_ny)[*])
    allocate(l_yg_line(nx,max_l_ny)[*])
    l_xg_line = 0.0
    l_yg_line = 0.0


    allocate(l_swed_layer(nx,max_l_ny,nz_max)[*],source=zero_real_three_d)
    allocate(l_snod_layer(nx,max_l_ny,nz_max)[*],source=zero_real_three_d)
    allocate(l_ro_layer(nx,max_l_ny,nz_max)[*],source=zero_real_three_d)
    allocate(l_T_old(nx,max_l_ny,nz_max)[*],source=zero_real_three_d)
    allocate(l_gamma(nx,max_l_ny,nz_max)[*],source=zero_real_three_d)
    allocate(l_diam_layer(nx,max_l_ny,nz_max)[*],source=zero_real_three_d)
    allocate(l_flux_layer(nx,max_l_ny,nz_max)[*],source=zero_real_three_d)
    allocate(l_k_stn(nx,max_l_ny,9)[*])
    l_k_stn = 0
    allocate(l_melt_flag(nx,max_l_ny,nz_max)[*])
    l_melt_flag = 0


    !Halo exchange buffer allocation
    allocate(left_rcv(nx)[*],right_rcv(nx)[*])
    left_rcv = 0.0
    right_rcv = 0.0

  end subroutine allocation

  subroutine distribute()
    use snowmodel_inc
    use snowmodel_vars
    implicit none

    integer :: ny_start, ny_end, idx

    if(me == 0) then
       do idx = 0, np-1

          ny_start = prefix_sum(idx) + 1
          ny_end = ny_start + partition_ny(idx) - 1

          l_wbal_subgrid(1:nx, 1:partition_ny(idx))[idx+1] = wbal_subgrid(1:nx, ny_start:ny_end)
          l_wbal_susp(1:nx, 1:partition_ny(idx))[idx+1] = wbal_susp(1:nx, ny_start:ny_end)
          l_wbal_salt(1:nx, 1:partition_ny(idx))[idx+1] = wbal_salt(1:nx, ny_start:ny_end)
          l_wbal_qsubl(1:nx, 1:partition_ny(idx))[idx+1] = wbal_qsubl(1:nx, ny_start:ny_end)
!          l_ro_soft_snow(1:nx, 1:partition_ny(idx))[idx+1] = ro_soft_snow(1:nx, ny_start:ny_end)
          l_ro_soft_snow_old(1:nx, 1:partition_ny(idx))[idx+1] = ro_soft_snow_old(1:nx, ny_start:ny_end)
!          l_curve_wt_lg(1:nx, 1:partition_ny(idx))[idx+1] = curve_wt_lg(1:nx, ny_start:ny_end)
!          l_vegtype(1:nx, 1:partition_ny(idx))[idx+1] = vegtype(1:nx, ny_start:ny_end)
!          l_vegsnowd_xy(1:nx, 1:partition_ny(idx))[idx+1] = vegsnowd_xy(1:nx, ny_start:ny_end)
!          l_veg_z0(1:nx, 1:partition_ny(idx))[idx+1] = veg_z0(1:nx, ny_start:ny_end)
!          l_tair_grid(1:nx, 1:partition_ny(idx))[idx+1] = tair_grid(1:nx, ny_start:ny_end)
          l_cloud_frac_grid(1:nx, 1:partition_ny(idx))[idx+1] = cloud_frac_grid(1:nx, ny_start:ny_end)
          l_sprec(1:nx, 1:partition_ny(idx))[idx+1] = sprec(1:nx, ny_start:ny_end)
!          l_cf_precip(1:nx, 1:partition_ny(idx))[idx+1] = cf_precip(1:nx, ny_start:ny_end)
          l_rh_grid(1:nx, 1:partition_ny(idx))[idx+1] = rh_grid(1:nx, ny_start:ny_end)
          l_windspd_grid(1:nx, 1:partition_ny(idx))[idx+1] = windspd_grid(1:nx, ny_start:ny_end)
          l_winddir_grid(1:nx, 1:partition_ny(idx))[idx+1] = winddir_grid(1:nx, ny_start:ny_end)
          l_soft_snow_d(1:nx, 1:partition_ny(idx))[idx+1] = soft_snow_d(1:nx, ny_start:ny_end)
          l_Qsi_grid(1:nx, 1:partition_ny(idx))[idx+1] = Qsi_grid(1:nx, ny_start:ny_end)
          l_Qli_grid(1:nx, 1:partition_ny(idx))[idx+1] = Qli_grid(1:nx, ny_start:ny_end)
          l_prec_grid(1:nx, 1:partition_ny(idx))[idx+1] = prec_grid(1:nx, ny_start:ny_end)
          l_albedo(1:nx, 1:partition_ny(idx))[idx+1] = albedo(1:nx, ny_start:ny_end)
!          l_curvature(1:nx, 1:partition_ny(idx))[idx+1] = curvature(1:nx, ny_start:ny_end)
!          l_slope_az(1:nx, 1:partition_ny(idx))[idx+1] = slope_az(1:nx, ny_start:ny_end)
!          l_terrain_slope(1:nx, 1:partition_ny(idx))[idx+1] = terrain_slope(1:nx, ny_start:ny_end)
          l_swe_depth(1:nx, 1:partition_ny(idx))[idx+1] = swe_depth(1:nx, ny_start:ny_end)
          l_sfc_pressure(1:nx, 1:partition_ny(idx))[idx+1] = sfc_pressure(1:nx, ny_start:ny_end)
          l_Qm(1:nx, 1:partition_ny(idx))[idx+1] = Qm(1:nx, ny_start:ny_end)
          l_Qh(1:nx, 1:partition_ny(idx))[idx+1] = Qh(1:nx, ny_start:ny_end)
          l_Qc(1:nx, 1:partition_ny(idx))[idx+1] = Qc(1:nx, ny_start:ny_end)
          l_Qe(1:nx, 1:partition_ny(idx))[idx+1] = Qe(1:nx, ny_start:ny_end)
          l_Qf(1:nx, 1:partition_ny(idx))[idx+1] = Qf(1:nx, ny_start:ny_end)
          l_Qle(1:nx, 1:partition_ny(idx))[idx+1] = Qle(1:nx, ny_start:ny_end)
          l_Tsfc(1:nx, 1:partition_ny(idx))[idx+1] = Tsfc(1:nx, ny_start:ny_end)
!          l_snow_d(1:nx, 1:partition_ny(idx))[idx+1] = snow_d(1:nx, ny_start:ny_end)
          l_snow_d_init(1:nx, 1:partition_ny(idx))[idx+1] = snow_d_init(1:nx, ny_start:ny_end)
!          l_tslsnowfall(1:nx, 1:partition_ny(idx))[idx+1] = tslsnowfall(1:nx, ny_start:ny_end)
!          l_change_layer(1:nx, 1:partition_ny(idx))[idx+1] = change_layer(1:nx, ny_start:ny_end)
          l_KK(1:nx, 1:partition_ny(idx))[idx+1] = KK(1:nx, ny_start:ny_end)
          l_ro_nsnow(1:nx, 1:partition_ny(idx))[idx+1] = ro_nsnow(1:nx, ny_start:ny_end)
          l_snow_depth(1:nx, 1:partition_ny(idx))[idx+1] = snow_depth(1:nx, ny_start:ny_end)
          l_rain(1:nx, 1:partition_ny(idx))[idx+1] = rain(1:nx, ny_start:ny_end)
          l_sum_runoff(1:nx, 1:partition_ny(idx))[idx+1] = sum_runoff(1:nx, ny_start:ny_end)
          l_sum_sprec(1:nx, 1:partition_ny(idx))[idx+1] = sum_sprec(1:nx, ny_start:ny_end)
          l_sum_prec(1:nx, 1:partition_ny(idx))[idx+1] = sum_prec(1:nx, ny_start:ny_end)
          l_runoff(1:nx, 1:partition_ny(idx))[idx+1] = runoff(1:nx, ny_start:ny_end)
          l_xro_snow(1:nx, 1:partition_ny(idx))[idx+1] = xro_snow(1:nx, ny_start:ny_end)
          l_canopy_int(1:nx, 1:partition_ny(idx))[idx+1] = canopy_int(1:nx, ny_start:ny_end)
          l_canopy_int_old(1:nx, 1:partition_ny(idx))[idx+1] = canopy_int_old(1:nx, ny_start:ny_end)
          l_canopy_unload(1:nx, 1:partition_ny(idx))[idx+1] = canopy_unload(1:nx, ny_start:ny_end)
!          l_d_canopy_int(1:nx, 1:partition_ny(idx))[idx+1] = d_canopy_int(1:nx, ny_start:ny_end)
          l_sum_d_canopy_int(1:nx, 1:partition_ny(idx))[idx+1] = sum_d_canopy_int(1:nx, ny_start:ny_end)
          l_glacier_melt(1:nx, 1:partition_ny(idx))[idx+1] = glacier_melt(1:nx, ny_start:ny_end)
          l_sum_glacmelt(1:nx, 1:partition_ny(idx))[idx+1] = sum_glacmelt(1:nx, ny_start:ny_end)
          l_swemelt(1:nx, 1:partition_ny(idx))[idx+1] = swemelt(1:nx, ny_start:ny_end)
          l_sum_swemelt(1:nx, 1:partition_ny(idx))[idx+1] = sum_swemelt(1:nx, ny_start:ny_end)
          l_swesublim(1:nx, 1:partition_ny(idx))[idx+1] = swesublim(1:nx, ny_start:ny_end)
          l_swe_depth_old(1:nx, 1:partition_ny(idx))[idx+1] = swe_depth_old(1:nx, ny_start:ny_end)
          l_e_balance(1:nx, 1:partition_ny(idx))[idx+1] = e_balance(1:nx, ny_start:ny_end)
          l_w_balance(1:nx, 1:partition_ny(idx))[idx+1] = w_balance(1:nx, ny_start:ny_end)
          l_ro_snow_grid(1:nx, 1:partition_ny(idx))[idx+1] = ro_snow_grid(1:nx, ny_start:ny_end)
          l_sum_unload(1:nx, 1:partition_ny(idx))[idx+1] = sum_unload(1:nx, ny_start:ny_end)
          l_Qcs(1:nx, 1:partition_ny(idx))[idx+1] = Qcs(1:nx, ny_start:ny_end)
          l_sum_Qcs(1:nx, 1:partition_ny(idx))[idx+1] = sum_Qcs(1:nx, ny_start:ny_end)
          l_sum_sfcsublim(1:nx, 1:partition_ny(idx))[idx+1] = sum_sfcsublim(1:nx, ny_start:ny_end)
!          l_windspd_2m_grid(1:nx, 1:partition_ny(idx))[idx+1] = windspd_2m_grid(1:nx, ny_start:ny_end)
          l_uwind_grid(1:nx, 1:partition_ny(idx))[idx+1] = uwind_grid(1:nx, ny_start:ny_end)
          l_vwind_grid(1:nx, 1:partition_ny(idx))[idx+1] = vwind_grid(1:nx, ny_start:ny_end)
       end do
    end if

  end subroutine distribute

  subroutine distribute2(snowmodel_line_fg)
    use snowmodel_inc
    use snowmodel_vars
    implicit none
    real, intent(in) :: snowmodel_line_fg

    integer :: ny_start, ny_end, idx

    if(me == 0) then
       do idx = 0, np-1

          ny_start = prefix_sum(idx) + 1
          ny_end = ny_start + partition_ny(idx) - 1

!          l_curve_wt_lg(1:nx, 1:partition_ny(idx))[idx+1] = curve_wt_lg(1:nx, ny_start:ny_end)
!          l_vegtype(1:nx, 1:partition_ny(idx))[idx+1] = vegtype(1:nx, ny_start:ny_end)
!          l_vegsnowd_xy(1:nx, 1:partition_ny(idx))[idx+1] = vegsnowd_xy(1:nx, ny_start:ny_end)
!          l_veg_z0(1:nx, 1:partition_ny(idx))[idx+1] = veg_z0(1:nx, ny_start:ny_end)
!          l_cf_precip(1:nx, 1:partition_ny(idx))[idx+1] = cf_precip(1:nx, ny_start:ny_end)
          l_snow_d(1:nx, 1:partition_ny(idx))[idx+1] = snow_d(1:nx, ny_start:ny_end)
          l_vwind_grid(1:nx, 1:partition_ny(idx))[idx+1] = vwind_grid(1:nx, ny_start:ny_end)

          if (snowmodel_line_fg.eq.1.0) then
            l_xg_line(1:nx, 1:partition_ny(idx))[idx+1] = xg_line(1:nx, ny_start:ny_end)
            l_yg_line(1:nx, 1:partition_ny(idx))[idx+1] = yg_line(1:nx, ny_start:ny_end)
         endif
       end do
    end if

  end subroutine distribute2


  subroutine gather(print_multilayer_param)
    use snowmodel_inc
    use snowmodel_vars
    implicit none

    real, intent(in) :: print_multilayer_param

    integer :: ny_start, ny_end, idx

    sync all
    if(me == 0) then
       do idx = 0, np-1
          if(idx == 0) then
             ny_start = 1
             ny_end = partition_ny(idx)
          else
             ny_start = prefix_sum(idx) + 1
             ny_end = ny_start + partition_ny(idx) - 1
          end if

          vegtype(1:nx, ny_start:ny_end) = l_vegtype(1:nx, 1:partition_ny(idx))[idx+1]
!          vegsnowd_xy(1:nx, ny_start:ny_end) = l_vegsnowd_xy(1:nx, 1:partition_ny(idx))[idx+1]
!          veg_z0(1:nx, ny_start:ny_end) = l_veg_z0(1:nx, 1:partition_ny(idx))[idx+1]
          tair_grid(1:nx, ny_start:ny_end) = l_tair_grid(1:nx, 1:partition_ny(idx))[idx+1]
          cloud_frac_grid(1:nx, ny_start:ny_end) = l_cloud_frac_grid(1:nx, 1:partition_ny(idx))[idx+1]
          sprec(1:nx, ny_start:ny_end) = l_sprec(1:nx, 1:partition_ny(idx))[idx+1]
!          cf_precip(1:nx, ny_start:ny_end) = l_cf_precip(1:nx, 1:partition_ny(idx))[idx+1]
          rh_grid(1:nx, ny_start:ny_end) = l_rh_grid(1:nx, 1:partition_ny(idx))[idx+1]
          windspd_grid(1:nx, ny_start:ny_end) = l_windspd_grid(1:nx, 1:partition_ny(idx))[idx+1]
          winddir_grid(1:nx, ny_start:ny_end) = l_winddir_grid(1:nx, 1:partition_ny(idx))[idx+1]
          soft_snow_d(1:nx, ny_start:ny_end) = l_soft_snow_d(1:nx, 1:partition_ny(idx))[idx+1]
          Qsi_grid(1:nx, ny_start:ny_end) = l_Qsi_grid(1:nx, 1:partition_ny(idx))[idx+1]
          Qli_grid(1:nx, ny_start:ny_end) = l_Qli_grid(1:nx, 1:partition_ny(idx))[idx+1]
          prec_grid(1:nx, ny_start:ny_end) = l_prec_grid(1:nx, 1:partition_ny(idx))[idx+1]
          albedo(1:nx, ny_start:ny_end) = l_albedo(1:nx, 1:partition_ny(idx))[idx+1]
!          curvature(1:nx, ny_start:ny_end) = l_curvature(1:nx, 1:partition_ny(idx))[idx+1]
!          slope_az(1:nx, ny_start:ny_end) = l_slope_az(1:nx, 1:partition_ny(idx))[idx+1]
!          terrain_slope(1:nx, ny_start:ny_end) = l_terrain_slope(1:nx, 1:partition_ny(idx))[idx+1]
          swe_depth(1:nx, ny_start:ny_end) = l_swe_depth(1:nx, 1:partition_ny(idx))[idx+1]
          sfc_pressure(1:nx, ny_start:ny_end) = l_sfc_pressure(1:nx, 1:partition_ny(idx))[idx+1]
          Qm(1:nx, ny_start:ny_end) = l_Qm(1:nx, 1:partition_ny(idx))[idx+1]
          Qh(1:nx, ny_start:ny_end) = l_Qh(1:nx, 1:partition_ny(idx))[idx+1]
          Qc(1:nx, ny_start:ny_end) = l_Qc(1:nx, 1:partition_ny(idx))[idx+1]
          Qe(1:nx, ny_start:ny_end) = l_Qe(1:nx, 1:partition_ny(idx))[idx+1]
          Qf(1:nx, ny_start:ny_end) = l_Qf(1:nx, 1:partition_ny(idx))[idx+1]
          Qle(1:nx, ny_start:ny_end) = l_Qle(1:nx, 1:partition_ny(idx))[idx+1]
          Tsfc(1:nx, ny_start:ny_end) = l_Tsfc(1:nx, 1:partition_ny(idx))[idx+1]
          snow_d(1:nx, ny_start:ny_end) = l_snow_d(1:nx, 1:partition_ny(idx))[idx+1]
          snow_d_init(1:nx, ny_start:ny_end) = l_snow_d_init(1:nx, 1:partition_ny(idx))[idx+1]
!          tslsnowfall(1:nx, ny_start:ny_end) = l_tslsnowfall(1:nx, 1:partition_ny(idx))[idx+1]
!          change_layer(1:nx, ny_start:ny_end) = l_change_layer(1:nx, 1:partition_ny(idx))[idx+1]
          KK(1:nx, ny_start:ny_end) = l_KK(1:nx, 1:partition_ny(idx))[idx+1]
          ro_nsnow(1:nx, ny_start:ny_end) = l_ro_nsnow(1:nx, 1:partition_ny(idx))[idx+1]
          snow_depth(1:nx, ny_start:ny_end) = l_snow_depth(1:nx, 1:partition_ny(idx))[idx+1]
          rain(1:nx, ny_start:ny_end) = l_rain(1:nx, 1:partition_ny(idx))[idx+1]
          sum_runoff(1:nx, ny_start:ny_end) = l_sum_runoff(1:nx, 1:partition_ny(idx))[idx+1]
          sum_sprec(1:nx, ny_start:ny_end) = l_sum_sprec(1:nx, 1:partition_ny(idx))[idx+1]
          sum_prec(1:nx, ny_start:ny_end) = l_sum_prec(1:nx, 1:partition_ny(idx))[idx+1]
          runoff(1:nx, ny_start:ny_end) = l_runoff(1:nx, 1:partition_ny(idx))[idx+1]
          xro_snow(1:nx, ny_start:ny_end) = l_xro_snow(1:nx, 1:partition_ny(idx))[idx+1]
          canopy_int(1:nx, ny_start:ny_end) = l_canopy_int(1:nx, 1:partition_ny(idx))[idx+1]
          canopy_int_old(1:nx, ny_start:ny_end) = l_canopy_int_old(1:nx, 1:partition_ny(idx))[idx+1]
          canopy_unload(1:nx, ny_start:ny_end) = l_canopy_unload(1:nx, 1:partition_ny(idx))[idx+1]
!          d_canopy_int(1:nx, ny_start:ny_end) = l_d_canopy_int(1:nx, 1:partition_ny(idx))[idx+1]
          sum_d_canopy_int(1:nx, ny_start:ny_end) = l_sum_d_canopy_int(1:nx, 1:partition_ny(idx))[idx+1]
          glacier_melt(1:nx, ny_start:ny_end) = l_glacier_melt(1:nx, 1:partition_ny(idx))[idx+1]
          sum_glacmelt(1:nx, ny_start:ny_end) = l_sum_glacmelt(1:nx, 1:partition_ny(idx))[idx+1]
          swemelt(1:nx, ny_start:ny_end) = l_swemelt(1:nx, 1:partition_ny(idx))[idx+1]
          sum_swemelt(1:nx, ny_start:ny_end) = l_sum_swemelt(1:nx, 1:partition_ny(idx))[idx+1]
          swesublim(1:nx, ny_start:ny_end) = l_swesublim(1:nx, 1:partition_ny(idx))[idx+1]
          swe_depth_old(1:nx, ny_start:ny_end) = l_swe_depth_old(1:nx, 1:partition_ny(idx))[idx+1]
          e_balance(1:nx, ny_start:ny_end) = l_e_balance(1:nx, 1:partition_ny(idx))[idx+1]
          w_balance(1:nx, ny_start:ny_end) = l_w_balance(1:nx, 1:partition_ny(idx))[idx+1]
          ro_snow_grid(1:nx, ny_start:ny_end) = l_ro_snow_grid(1:nx, 1:partition_ny(idx))[idx+1]
          ro_soft_snow_old(1:nx, ny_start:ny_end) = l_ro_soft_snow_old(1:nx, 1:partition_ny(idx))[idx+1]
          sum_unload(1:nx, ny_start:ny_end) = l_sum_unload(1:nx, 1:partition_ny(idx))[idx+1]
          Qcs(1:nx, ny_start:ny_end) = l_Qcs(1:nx, 1:partition_ny(idx))[idx+1]
          sum_Qcs(1:nx, ny_start:ny_end) = l_sum_Qcs(1:nx, 1:partition_ny(idx))[idx+1]
          sum_sfcsublim(1:nx, ny_start:ny_end) = l_sum_sfcsublim(1:nx, 1:partition_ny(idx))[idx+1]
!          windspd_2m_grid(1:nx, ny_start:ny_end) = l_windspd_2m_grid(1:nx, 1:partition_ny(idx))[idx+1]
          uwind_grid(1:nx, ny_start:ny_end) = l_uwind_grid(1:nx, 1:partition_ny(idx))[idx+1]
          vwind_grid(1:nx, ny_start:ny_end) = l_vwind_grid(1:nx, 1:partition_ny(idx))[idx+1]
          if (print_multilayer_param.eq.1.0) then
            snod_layer(1:nx, ny_start:ny_end,:) = l_snod_layer(1:nx, 1:partition_ny(idx),:)[idx+1]
            swed_layer(1:nx, ny_start:ny_end,:) = l_swed_layer(1:nx, 1:partition_ny(idx),:)[idx+1]
            ro_layer(1:nx, ny_start:ny_end,:) = l_ro_layer(1:nx, 1:partition_ny(idx),:)[idx+1]
            T_old(1:nx, ny_start:ny_end,:) = l_T_old(1:nx, 1:partition_ny(idx),:)[idx+1]
            gamma(1:nx, ny_start:ny_end,:) = l_gamma(1:nx, 1:partition_ny(idx),:)[idx+1]
            diam_layer(1:nx, ny_start:ny_end,:) = l_diam_layer(1:nx, 1:partition_ny(idx),:)[idx+1]
            flux_layer(1:nx, ny_start:ny_end,:) = l_flux_layer(1:nx, 1:partition_ny(idx),:)[idx+1]
          endif

       end do
    end if

    sync all

  end subroutine gather

  subroutine single_gather_integer(global_array, local_array)
    implicit none

    integer, intent(out) :: global_array(:,:)
    integer :: local_array(:,:)[*]
    integer :: idx,ny_start,ny_end

    sync all

    do idx = 0,np-1
       ny_start = prefix_sum(idx) + 1
       ny_end = ny_start + partition_ny(idx) - 1

       global_array(1:max_l_nx,ny_start:ny_end) = local_array(1:max_l_nx,1:partition_ny(idx))[idx+1]

    end do

    sync all

  end subroutine single_gather_integer

  subroutine single_gather_3D_integer(global_array, local_array)
    implicit none

    integer, intent(out) :: global_array(:,:,:)
    integer :: local_array(:,:,:)[*]
    integer :: idx,ny_start,ny_end

    sync all

    do idx = 0,np-1
       ny_start = prefix_sum(idx) + 1
       ny_end = ny_start + partition_ny(idx) - 1

       global_array(1:max_l_nx,ny_start:ny_end,:) = local_array(1:max_l_nx,1:partition_ny(idx),:)[idx+1]

    end do

    sync all

  end subroutine single_gather_3D_integer

  subroutine single_gather_real(global_array, local_array)
    implicit none

    real, intent(out) :: global_array(:,:)
    real :: local_array(:,:)[*]
    integer :: idx,ny_start,ny_end

    sync all

    do idx = 0,np-1
       ny_start = prefix_sum(idx) + 1
       ny_end = ny_start + partition_ny(idx) - 1

       global_array(1:max_l_nx,ny_start:ny_end) = local_array(1:max_l_nx,1:partition_ny(idx))[idx+1]

    end do

    sync all

  end subroutine single_gather_real


    subroutine single_scatter_integer(global_array, local_array)
    implicit none

    integer :: global_array(:,:)
    integer :: local_array(:,:)[*]
    integer :: idx,ny_start,ny_end,nx

    nx = max_l_nx

    sync all

    if(me == 0) then
       do idx = 0, np-1

          ny_start = prefix_sum(idx) + 1
          ny_end = ny_start + partition_ny(idx) - 1

          local_array(1:nx, 1:partition_ny(idx))[idx+1] = global_array(1:nx, ny_start:ny_end)
       end do
    end if

    sync all

  end subroutine single_scatter_integer

  subroutine single_scatter_real(global_array, local_array)
    implicit none

    real :: global_array(:,:)
    real :: local_array(:,:)[*]
    integer :: idx,ny_start,ny_end,nx

    nx = max_l_nx

    sync all

    if(me == 0) then
       do idx = 0, np-1

          ny_start = prefix_sum(idx) + 1
          ny_end = ny_start + partition_ny(idx) - 1

          local_array(1:nx, 1:partition_ny(idx))[idx+1] = global_array(1:nx, ny_start:ny_end)
       end do
    end if

    sync all

  end subroutine single_scatter_real

  subroutine single_scatter_double(global_array, local_array)
    implicit none

    double precision :: global_array(:,:)
    double precision :: local_array(:,:)[*]
    integer :: idx,ny_start,ny_end,nx

    nx = max_l_nx

    sync all

    if(me == 0) then
       do idx = 0, np-1

          ny_start = prefix_sum(idx) + 1
          ny_end = ny_start + partition_ny(idx) - 1

          local_array(1:nx, 1:partition_ny(idx))[idx+1] = global_array(1:nx, ny_start:ny_end)
       end do
    end if

    sync all

  end subroutine single_scatter_double

  subroutine manual_bcast()
    use snowmodel_inc
    use snowmodel_vars
    implicit none

    integer :: ny_start, ny_end, idx

    sync all

    do idx = 0,np-1

       ny_start = prefix_sum(idx) + 1
       ny_end = ny_start + partition_ny(idx) - 1

       vegtype(1:nx, ny_start:ny_end) = l_vegtype(1:nx, 1:partition_ny(idx))[idx+1]
       vegsnowd_xy(1:nx, ny_start:ny_end) = l_vegsnowd_xy(1:nx, 1:partition_ny(idx))[idx+1]
!       veg_z0(1:nx, ny_start:ny_end) = l_veg_z0(1:nx, 1:partition_ny(idx))[idx+1]
       tair_grid(1:nx, ny_start:ny_end) = l_tair_grid(1:nx, 1:partition_ny(idx))[idx+1]
       cloud_frac_grid(1:nx, ny_start:ny_end) = l_cloud_frac_grid(1:nx, 1:partition_ny(idx))[idx+1]
       sprec(1:nx, ny_start:ny_end) = l_sprec(1:nx, 1:partition_ny(idx))[idx+1]
       cf_precip(1:nx, ny_start:ny_end) = l_cf_precip(1:nx, 1:partition_ny(idx))[idx+1]
       rh_grid(1:nx, ny_start:ny_end) = l_rh_grid(1:nx, 1:partition_ny(idx))[idx+1]
       windspd_grid(1:nx, ny_start:ny_end) = l_windspd_grid(1:nx, 1:partition_ny(idx))[idx+1]
       winddir_grid(1:nx, ny_start:ny_end) = l_winddir_grid(1:nx, 1:partition_ny(idx))[idx+1]
       soft_snow_d(1:nx, ny_start:ny_end) = l_soft_snow_d(1:nx, 1:partition_ny(idx))[idx+1]
       Qsi_grid(1:nx, ny_start:ny_end) = l_Qsi_grid(1:nx, 1:partition_ny(idx))[idx+1]
       Qli_grid(1:nx, ny_start:ny_end) = l_Qli_grid(1:nx, 1:partition_ny(idx))[idx+1]
       prec_grid(1:nx, ny_start:ny_end) = l_prec_grid(1:nx, 1:partition_ny(idx))[idx+1]
       albedo(1:nx, ny_start:ny_end) = l_albedo(1:nx, 1:partition_ny(idx))[idx+1]
!       curvature(1:nx, ny_start:ny_end) = l_curvature(1:nx, 1:partition_ny(idx))[idx+1]
!       slope_az(1:nx, ny_start:ny_end) = l_slope_az(1:nx, 1:partition_ny(idx))[idx+1]
!       terrain_slope(1:nx, ny_start:ny_end) = l_terrain_slope(1:nx, 1:partition_ny(idx))[idx+1]
       swe_depth(1:nx, ny_start:ny_end) = l_swe_depth(1:nx, 1:partition_ny(idx))[idx+1]
       sfc_pressure(1:nx, ny_start:ny_end) = l_sfc_pressure(1:nx, 1:partition_ny(idx))[idx+1]
       Qm(1:nx, ny_start:ny_end) = l_Qm(1:nx, 1:partition_ny(idx))[idx+1]
       Qh(1:nx, ny_start:ny_end) = l_Qh(1:nx, 1:partition_ny(idx))[idx+1]
       Qc(1:nx, ny_start:ny_end) = l_Qc(1:nx, 1:partition_ny(idx))[idx+1]
       Qe(1:nx, ny_start:ny_end) = l_Qe(1:nx, 1:partition_ny(idx))[idx+1]
       Qf(1:nx, ny_start:ny_end) = l_Qf(1:nx, 1:partition_ny(idx))[idx+1]
       Qle(1:nx, ny_start:ny_end) = l_Qle(1:nx, 1:partition_ny(idx))[idx+1]
       Tsfc(1:nx, ny_start:ny_end) = l_Tsfc(1:nx, 1:partition_ny(idx))[idx+1]
       snow_d(1:nx, ny_start:ny_end) = l_snow_d(1:nx, 1:partition_ny(idx))[idx+1]
       snow_d_init(1:nx, ny_start:ny_end) = l_snow_d_init(1:nx, 1:partition_ny(idx))[idx+1]
!       tslsnowfall(1:nx, ny_start:ny_end) = l_tslsnowfall(1:nx, 1:partition_ny(idx))[idx+1]
!       change_layer(1:nx, ny_start:ny_end) = l_change_layer(1:nx, 1:partition_ny(idx))[idx+1]
       KK(1:nx, ny_start:ny_end) = l_KK(1:nx, 1:partition_ny(idx))[idx+1]
       ro_nsnow(1:nx, ny_start:ny_end) = l_ro_nsnow(1:nx, 1:partition_ny(idx))[idx+1]
       snow_depth(1:nx, ny_start:ny_end) = l_snow_depth(1:nx, 1:partition_ny(idx))[idx+1]
       rain(1:nx, ny_start:ny_end) = l_rain(1:nx, 1:partition_ny(idx))[idx+1]
       sum_runoff(1:nx, ny_start:ny_end) = l_sum_runoff(1:nx, 1:partition_ny(idx))[idx+1]
       sum_sprec(1:nx, ny_start:ny_end) = l_sum_sprec(1:nx, 1:partition_ny(idx))[idx+1]
       sum_prec(1:nx, ny_start:ny_end) = l_sum_prec(1:nx, 1:partition_ny(idx))[idx+1]
       runoff(1:nx, ny_start:ny_end) = l_runoff(1:nx, 1:partition_ny(idx))[idx+1]
       xro_snow(1:nx, ny_start:ny_end) = l_xro_snow(1:nx, 1:partition_ny(idx))[idx+1]
       canopy_int(1:nx, ny_start:ny_end) = l_canopy_int(1:nx, 1:partition_ny(idx))[idx+1]
       canopy_int_old(1:nx, ny_start:ny_end) = l_canopy_int_old(1:nx, 1:partition_ny(idx))[idx+1]
       canopy_unload(1:nx, ny_start:ny_end) = l_canopy_unload(1:nx, 1:partition_ny(idx))[idx+1]
!       d_canopy_int(1:nx, ny_start:ny_end) = l_d_canopy_int(1:nx, 1:partition_ny(idx))[idx+1]
       sum_d_canopy_int(1:nx, ny_start:ny_end) = l_sum_d_canopy_int(1:nx, 1:partition_ny(idx))[idx+1]
       glacier_melt(1:nx, ny_start:ny_end) = l_glacier_melt(1:nx, 1:partition_ny(idx))[idx+1]
       sum_glacmelt(1:nx, ny_start:ny_end) = l_sum_glacmelt(1:nx, 1:partition_ny(idx))[idx+1]
       swemelt(1:nx, ny_start:ny_end) = l_swemelt(1:nx, 1:partition_ny(idx))[idx+1]
       sum_swemelt(1:nx, ny_start:ny_end) = l_sum_swemelt(1:nx, 1:partition_ny(idx))[idx+1]
       swesublim(1:nx, ny_start:ny_end) = l_swesublim(1:nx, 1:partition_ny(idx))[idx+1]
       swe_depth_old(1:nx, ny_start:ny_end) = l_swe_depth_old(1:nx, 1:partition_ny(idx))[idx+1]
       e_balance(1:nx, ny_start:ny_end) = l_e_balance(1:nx, 1:partition_ny(idx))[idx+1]
       w_balance(1:nx, ny_start:ny_end) = l_w_balance(1:nx, 1:partition_ny(idx))[idx+1]
       ro_snow_grid(1:nx, ny_start:ny_end) = l_ro_snow_grid(1:nx, 1:partition_ny(idx))[idx+1]
       sum_unload(1:nx, ny_start:ny_end) = l_sum_unload(1:nx, 1:partition_ny(idx))[idx+1]
       Qcs(1:nx, ny_start:ny_end) = l_Qcs(1:nx, 1:partition_ny(idx))[idx+1]
       sum_Qcs(1:nx, ny_start:ny_end) = l_sum_Qcs(1:nx, 1:partition_ny(idx))[idx+1]
       sum_sfcsublim(1:nx, ny_start:ny_end) = l_sum_sfcsublim(1:nx, 1:partition_ny(idx))[idx+1]
!       windspd_2m_grid(1:nx, ny_start:ny_end) = l_windspd_2m_grid(1:nx, 1:partition_ny(idx))[idx+1]
       uwind_grid(1:nx, ny_start:ny_end) = l_uwind_grid(1:nx, 1:partition_ny(idx))[idx+1]
       vwind_grid(1:nx, ny_start:ny_end) = l_vwind_grid(1:nx, 1:partition_ny(idx))[idx+1]
    end do
    sync all
  end subroutine manual_bcast

  function global_map(j) result(new_j)
    implicit none
    integer, intent(in) :: j
    integer :: new_j

    new_j = j + prefix_sum(me)

  end function global_map

  subroutine new_j_start_end(new_j_start, new_j_end)
    implicit none

    integer, intent(out) :: new_j_start, new_j_end

    new_j_start = 1
    new_j_end = l_ny

    if(me == 0) new_j_start = 2
    if(me == np-1) new_j_end = l_ny - 1

  end subroutine new_j_start_end

  ! The halo exchange is done using "puts" instead of "gets"
  subroutine halo_exchange(array_2d)
    implicit none

    real, intent(in) :: array_2d(max_l_nx,max_l_ny)
    integer :: nx

    nx = max_l_nx

    sync all

    !only right exchange
    if(me < np-1) then
       ! +2 because CAF counts from 1 and we transfer on our right
       left_rcv(1:nx)[me+2] = array_2d(1:nx,l_ny)
    end if
    !only left exchange
    if(me > 0) then
       right_rcv(1:nx)[me] = array_2d(1:nx,1)
    end if

    sync all

  end subroutine halo_exchange

end module caf_module
