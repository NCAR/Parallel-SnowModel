! c snowtran_code.f

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccc     Snow-Transport Modeling System - 3D (SnowTran-3D)    cccccc
! ccccc                    Copyright (C) 1998                    cccccc
! ccccc          by Glen E. Liston, InterWorks Consulting        cccccc
! ccccc                    All Rights Reserved                   cccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! c This FORTRAN code receives inputs of wind speed, wind direction,
! c   air temperature, relative humidity, vegetation type, topography,
! c   and precipitation, and it outputs snow depth, saltation flux,
! c   suspended flux, sublimation of blowing snow, and the snow depth
! c   changes resulting from these processes.
! c
! c All units are in m, kg, s, K.
! c
! c This model is described in the paper:
! c   A Snow-Transport Model for Complex Terrain, by Glen E. Liston
! c   and Matthew Sturm, Journal of Glaciology, 1998, Vol. 44,
! c   No. 148, pages 498-516.
! c
! c The author of this code is:
! c   Dr. Glen E. Liston
! c   InterWorks Consulting
! c   15621 SnowMan Road
! c   Loveland, Colorado 80538
! c
! c To run in 2-D mode, set nx = 3 and look at the data at i = 2.
! c   This is required because of the boundary conditions imposed
! c   along i = 1 and 3.

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine SNOWTRAN_CODE(bc_flag,bs_flag,C_z,&
     &  conc_salt,deltax,deltay,dh_salt,dh_salt_u,dh_salt_v,&
     &  dh_susp,dh_susp_u,dh_susp_v,dt,dz_susp,fall_vel,fetch,&
     &  gravity,h_const,h_star,ht_rhobs,ht_windobs,index_ue,&
     &  index_uw,index_vn,index_vs,iter,nx,ny,pi,Qsalt,Qsalt_max,&
     &  Qsalt_maxu,Qsalt_maxv,Qsalt_u,Qsalt_v,Qsubl,Qsusp,&
     &  Qsusp_u,Qsusp_v,rh_grid,ro_air,ro_snow,ro_water,snow_d,&
     &  snow_d_init,snow_z0,soft_snow_d,sprec,sum_glacmelt,&
     &  subgrid_flag,wbal_salt,wbal_susp,wbal_qsubl,sum_sprec,&
     &  tabler_ee,tabler_ne,tabler_nn,tabler_nw,tabler_se,&
     &  tabler_ss,tabler_sw,tabler_ww,tair_grid,topo,topo_land,&
     &  topoflag,twolayer_flag,Up_const,Ur_const,Utau,&
     &  Utau_t,uwind_grid,veg_z0,vegsnowd_xy,vonKarman,&
     &  vwind_grid,wind_min,winddir_flag,winddir_grid,&
     &  windspd_flag,windspd_grid,xmu,z_0,ztop_susp,max_iter,&
     &  run_enbal,run_snowpack,wbal_subgrid,sum_qsubl,sum_trans,&
     &  swe_depth,snow_depth,ro_snow_grid,sum_prec,sum_runoff,&
     &  sum_Qcs,canopy_int,w_balance,sum_sfcsublim,tabler_dir,&
     &  slope_adjust,Utau_t_const,Utau_t_flag,ro_soft_snow_old,&
     &  ro_soft_snow,ro_nsnow,prec,Qcs,runoff,d_canopy_int,&
     &  glacier_melt,swe_depth_old,swesublim,canopy_unload,&
     &  canopy_int_old,iter_start,multilayer_snowpack,swed_layer,&
     &  KK,snod_layer,ro_layer,curve_lg_scale_flag,curve_wt_lg,&
     &  seaice_run,seaice_conc,tslsnowfall,T_old,tsls_threshold,&
     &  curve_len_scale,Tabler_1_flag,Tabler_2_flag,undef,&
     &  tabler_sfc_path_name,output_path_wo_assim,&
     &  output_path_wi_assim,icorr_factor_loop,windspd_2m_grid,&
     &  Qsubl_depth)

  use caf_module, only: l_ny, max_l_ny, single_gather, single_scatter,&
       global_map, me
  use snowmodel_vars, only: l_Utau,l_snow_d,l_topo_tmp, &
                          & l_dh_subgrid,l_Qsalt_v,l_vwind_grid,&
                          & l_Qsusp_v,l_ro_snow_grid,l_tabler_nn_orig,&
                          & l_tabler_ne_orig,l_tabler_nw_orig,l_tabler_ss_orig,&
                          & l_tabler_se_orig,l_tabler_sw_orig,l_tabler_ee_orig,&
                          & l_tabler_ww_orig
  use snowmodel_inc
  implicit none

  integer iter,nx,ny,i,j,iter_start,max_iter,new_j

  real ro_snow,ro_water,ro_air,gravity,vonKarman,snow_z0
  real deltax,deltay,dt,undef
  real fetch,xmu,C_z,h_const,wind_min,windspd_flag
  real Up_const,dz_susp,ztop_susp,fall_vel,Ur_const
  real Utau_t_const,pi,bc_flag,topoflag,Utau_t_flag
  real ht_windobs,ht_rhobs,bs_flag,twolayer_flag
  real subgrid_flag,winddir_flag,curve_len_scale
  real run_enbal,run_snowpack,tabler_dir,slope_adjust

  real topo_land(nx,ny)
  real tabler_nn(nx,max_l_ny)
  real tabler_ss(nx,max_l_ny)
  real tabler_ee(nx,max_l_ny)
  real tabler_ww(nx,max_l_ny)
  real tabler_ne(nx,max_l_ny)
  real tabler_se(nx,max_l_ny)
  real tabler_sw(nx,max_l_ny)
  real tabler_nw(nx,max_l_ny)
  real topo(nx,ny)

  real topo_tmp(nx,ny)

  real uwind_grid(nx,max_l_ny),vwind_grid(nx,ny)
  real windspd_grid(nx,max_l_ny),winddir_grid(nx,max_l_ny)
  real tair_grid(nx,max_l_ny),sprec(nx,max_l_ny)
  real rh_grid(nx,max_l_ny),windspd_2m_grid(nx,max_l_ny)

  integer index_ue(ny_max,2*nx_max+1),index_uw(ny_max,2*nx_max+1)
  integer index_vn(nx_max,2*ny_max+1),index_vs(nx_max,2*ny_max+1)

  real snow_d(nx,max_l_ny)
  real snow_depth(nx,max_l_ny)
  real swe_depth(nx,max_l_ny)
  real ro_snow_grid(nx,max_l_ny)
  real ro_soft_snow(nx,max_l_ny)
  real ro_soft_snow_old(nx,max_l_ny)
  real ro_nsnow(nx,max_l_ny)
  real snow_d_init(nx,ny)
  real Utau(nx,ny)
  real Utau_t(nx,max_l_ny)
  real z_0(nx,max_l_ny)
  real h_star(nx,max_l_ny)
  real conc_salt(nx,max_l_ny)

  real Qsalt_max(nx,max_l_ny)
  real Qsalt_maxu(nx,max_l_ny),Qsalt_maxv(nx,max_l_ny)
  real Qsalt(nx,max_l_ny)
  real Qsalt_u(nx,max_l_ny),Qsalt_v(nx,max_l_ny)
  real dh_salt(nx,max_l_ny)
  real dh_salt_u(nx,max_l_ny),dh_salt_v(nx,max_l_ny)


  real Qsusp(nx,max_l_ny)
  real Qsusp_u(nx,max_l_ny),Qsusp_v(nx,max_l_ny)
  real dh_susp(nx,max_l_ny)
  real dh_susp_u(nx,max_l_ny),dh_susp_v(nx,max_l_ny)

  real Qsubl(nx,max_l_ny)
  real Qsubl_depth(nx,max_l_ny)


  real sum_sprec(nx,max_l_ny)
  real wbal_qsubl(nx,max_l_ny)
  real wbal_salt(nx,max_l_ny)
  real wbal_susp(nx,max_l_ny)
  real wbal_subgrid(nx,max_l_ny)
  real sum_qsubl(nx,max_l_ny)
  real sum_trans(nx,max_l_ny)
  real soft_snow_d(nx,max_l_ny)

  real prec(nx,max_l_ny)
  real Qcs(nx,max_l_ny)
  real runoff(nx,max_l_ny)
  real d_canopy_int(nx,max_l_ny)
  real glacier_melt(nx,max_l_ny)
  real swe_depth_old(nx,max_l_ny)
  real swesublim(nx,max_l_ny)
  real canopy_unload(nx,max_l_ny)
  real canopy_int_old(nx,max_l_ny)

  real vegsnowd_xy(nx,max_l_ny)
  real veg_z0(nx,max_l_ny)

  real sum_glacmelt(nx,max_l_ny),w_balance(nx,max_l_ny),&
       &  sum_prec(nx,max_l_ny),sum_runoff(nx,max_l_ny),&
       &  sum_Qcs(nx,max_l_ny),canopy_int(nx,max_l_ny),&
       &  sum_sfcsublim(nx,max_l_ny)


  integer multilayer_snowpack,k
  integer KK(nx,max_l_ny)
  real swe_change_tmp,swe_change,tsls_threshold
  real swed_layer_z(nz_max)
  real swed_layer(nx,max_l_ny,nz_max)
  real snod_layer(nx,max_l_ny,nz_max)
  real ro_layer(nx,max_l_ny,nz_max)
  real T_old(nx,max_l_ny,nz_max)
  real tslsnowfall(nx,max_l_ny)

  real curve_lg_scale_flag
  real curve_wt_lg(nx,max_l_ny)

  real seaice_run
  real seaice_conc(nx,ny)

  real Tabler_1_flag,Tabler_2_flag
  character*80 output_path_wo_assim,output_path_wi_assim,&
       &  tabler_sfc_path_name
  integer i_len_tabler,trailing_blanks,icorr_factor_loop,&
       &  i_len_wo,i_len_wi
  character*4 string,me_string
  character*9 fmt

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! c Perform some intialization steps that are unique to SnowTran-3D.
  if (iter.eq.iter_start) then

! c This is now done in preprocess_code.f
! c       print *,
! c    & 'cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
! c       print *,
! c    & 'c     Snow-Transport Modeling System - 3D (SnowTran-3D)    c'
! c       print *,
! c    & 'c                    Copyright (C) 1998                    c'
! c       print *,
! c    & 'c          by Glen E. Liston, InterWorks Consulting        c'
! c       print *,
! c    & 'c                    All Rights Reserved                   c'
! c       print *,
! c    & 'cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'

! c SnowTran-3D can only be run with SnowPack and EnBal turned on.
     if (run_enbal.ne.1.0 .and. run_snowpack.ne.1.0) then
        print *,'You cannot run SnowTran-3D without'
        print *,'SnowPack and EnBal anymore.'
        stop
     endif

     if (subgrid_flag.eq.1.0) then

!c Check to make sure topoflag = 1.0.
        if (subgrid_flag.eq.1.0 .and. topoflag.eq.0.0) then
           print *,'If subgrid_flag=1.0, then topoflag must = 1.0'
           print *,'  Correct this in snowmodel.par to continue.'
           stop
        endif

! c The Tabler surfaces were originally developed assuming the grid
! c   increment would never be less than 1.0 m.  You probably should
! c   not run it with deltax and deltay less than 1.0 without some
! c   further testing.
        if (deltax.lt.1.0 .or. deltay.lt.1.0) then
           print *,'The Tabler subgrid algorithm has not been'
           print *,'tested for deltax and/or deltay less than'
           print *,'1.0 m.  You should probably do some testing'
           print *,'before running the model at less than 1.0-m'
           print *,'resolution.  Acually I am pretty sure it will'
           print *,'run, but it will not generate the correct'
           print *,'snow-depth profiles.'
           stop
        endif

! c If this is the first time through, generate the Tabler snow
! c   accumulation surfaces for the land topography.

        call single_scatter(topo_land, l_topo_tmp)
        call tabler_3d(nx,l_ny,l_topo_tmp,deltax,deltay,&
             &      l_tabler_ww_orig,l_tabler_ee_orig,l_tabler_ss_orig,&
             &      l_tabler_nn_orig,l_tabler_ne_orig,l_tabler_se_orig,&
             &      l_tabler_sw_orig,l_tabler_nw_orig,slope_adjust)

! c Include the vegetation snow-holding depths in the Tabler
! c   surfaces.
        do i=1,nx
           do j=1,l_ny
              l_tabler_nn_orig(i,j) =&
                   &          max(l_tabler_nn_orig(i,j),vegsnowd_xy(i,j))
              l_tabler_ne_orig(i,j) =&
                   &          max(l_tabler_ne_orig(i,j),vegsnowd_xy(i,j))
              l_tabler_ee_orig(i,j) =&
                   &          max(l_tabler_ee_orig(i,j),vegsnowd_xy(i,j))
              l_tabler_se_orig(i,j) =&
                   &          max(l_tabler_se_orig(i,j),vegsnowd_xy(i,j))
              l_tabler_ss_orig(i,j) =&
                   &          max(l_tabler_ss_orig(i,j),vegsnowd_xy(i,j))
              l_tabler_sw_orig(i,j) =&
                   &          max(l_tabler_sw_orig(i,j),vegsnowd_xy(i,j))
              l_tabler_ww_orig(i,j) =&
                   &          max(l_tabler_ww_orig(i,j),vegsnowd_xy(i,j))
              l_tabler_nw_orig(i,j) =&
                   &          max(l_tabler_nw_orig(i,j),vegsnowd_xy(i,j))
           enddo
        enddo

! c Save the Tabler equilibrium surfaces for each of the 8 main wind
! c   directions.  Also save the land DEM distribution.
        if (Tabler_1_flag.eq.1.0) then
           fmt = '(I4)'
           write(string,fmt) me
           me_string = string

           i_len_tabler = 80 - trailing_blanks(tabler_sfc_path_name)
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_nn'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_nn_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !tabler_ne
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_ne'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_ne_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !tabler_nw
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_nw'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_nw_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !tabler_ss
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_ss'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_ss_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !tabler_sw
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_sw'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_sw_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !tabler_se
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_se'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_se_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !tabler_ee
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_ee'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_ee_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !tabler_ww
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'tabler_ww'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_tabler_ww_orig(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           !topo
           open (51+me,&
     &file=tabler_sfc_path_name(1:i_len_tabler)//'topo'//'/me'//&
     &            '_'//TRIM(ADJUSTL(me_string))//'.gdat',&
     &            form='unformatted',access='direct',recl=4*nx*l_ny,&
     &            status='replace')
           write(51+me,rec=1) ((l_topo_tmp(i,j),i=1,nx),j=1,l_ny)
           close(51)
           sync all
           close (51)
        endif

! c This available if you want to save the Tabler surface that was
! c   used at each time step.
        if (Tabler_2_flag.eq.1.0) then
           if (icorr_factor_loop.eq.1) then
              i_len_wo = 80 - trailing_blanks(output_path_wo_assim)
              open(52,&
                   &  file=output_path_wo_assim(1:i_len_wo)//'tabler_sfcs_iter.gdat',&
                   &        form='unformatted',access='direct',recl=4*nx*ny)
           endif
           if (icorr_factor_loop.eq.2) then
              i_len_wi = 80 - trailing_blanks(output_path_wi_assim)
              open(52,&
                   &  file=output_path_wi_assim(1:i_len_wi)//'tabler_sfcs_iter.gdat',&
                   &        form='unformatted',access='direct',recl=4*nx*ny)
           endif
        endif
     endif

  endif

! c Print out some basic run information to the screen.
! c     print 102, windspd_flag,winddir_flag
! c 102 format(25x,'    wind spd = ',f5.2,'   wind dir = ',f4.0)

! c In SnowTran-3D, the summed snow precipitation must be in units
! c   of snow-depth.  The rest of the routines assume that it is in
! c   swe units.
  do j=1,l_ny
     do i=1,nx
        sum_sprec(i,j) = sum_sprec(i,j) * ro_water / ro_snow
     enddo
  enddo

  call single_scatter(Utau,l_Utau)

!c Update the threshold friction velocity.
  if (Utau_t_flag.eq.0.0) then
     if (curve_lg_scale_flag.eq.1.0) then
        do j=1,l_ny
           do i=1,nx
              Utau_t(i,j) = curve_wt_lg(i,j) * Utau_t_const
           enddo
        enddo
     else
        do j=1,l_ny
           do i=1,nx
              Utau_t(i,j) = Utau_t_const
           enddo
        enddo
     endif
  elseif (Utau_t_flag.eq.1.0) then
     do j=1,l_ny
        do i=1,nx
           call surface_snow_1(tair_grid(i,j),windspd_2m_grid(i,j),&
                &        sprec(i,j),ro_soft_snow(i,j),Utau_t(i,j),&
                &        ro_soft_snow_old(i,j),dt,ro_nsnow(i,j))
        enddo
     enddo
     if (curve_lg_scale_flag.eq.1.0) then
        do j=1,l_ny
           do i=1,nx
              Utau_t(i,j) = curve_wt_lg(i,j) * Utau_t(i,j)
           enddo
        enddo
     endif
  endif


! c Set the blowing snow flag to zero until it is clear that we will
! c   have blowing snow.
  bs_flag = 0.0

! c If the wind speed is lower that some threshold, then don't
! c   need to to any of the snow transport computations.
  if (windspd_flag.ge.4.0) then

! c Get the wind direction indexing arrays for this particular
! c   wind event (time step).
!     call getdirection(nx,ny,uwind_grid,vwind_grid,index_ue,&
!          &    index_uw,index_vn,index_vs)

! c Solve for Utau and z_0 if snow is saltating, else solve assuming
! c   z_0 is known from snow depth and/or veg type, and solve for
! c   Utau.
     call solveUtau(l_Utau,ht_windobs,windspd_grid,C_z,vonKarman,&
          &    gravity,z_0,h_star,h_const,vegsnowd_xy,snow_d,&
          &    snow_z0,veg_z0,bs_flag,nx,ny,Utau_t,soft_snow_d)

! Gather local variables to global after solveUtau
     call single_gather(Utau,l_Utau)
! Scatter global variables to local before saltation
!     call single_scatter(Qsalt_v,l_Qsalt_v)
     call single_scatter(vwind_grid,l_vwind_grid)

! c If the blowing snow flag indicates wind transported snow
! c   somewhere within the domain (bs_flag = 1.0), run the saltation
     ! c   and suspension models.
     if (bs_flag.eq.1.0) then

! c Solve for the saltation flux.
! c         print *,'         Saltation'
        call saltation(Qsalt,deltax,fetch,Utau,Utau_t,nx,ny,&
             &      ro_air,gravity,vegsnowd_xy,snow_d,&
             &      Qsalt_max,Qsalt_maxu,Qsalt_maxv,deltay,Qsalt_u,Qsalt_v,&
             &      uwind_grid,vwind_grid,xmu,soft_snow_d,bc_flag)


! c Solve for the suspension flux.
! c         print *,'         Suspension'
        call suspension(l_Utau,vonKarman,nx,ny,conc_salt,&
             &      Qsalt,Qsusp,z_0,h_star,dz_susp,ztop_susp,pi,&
             &      fall_vel,Ur_const,Up_const,Utau_t,Qsubl,ht_rhobs,&
             &      tair_grid,rh_grid,Qsusp_u,Qsusp_v,uwind_grid,&
             &      l_vwind_grid)

     elseif (bs_flag.eq.0.0) then


        call noblowsnow(nx,ny,Qsalt_max,Qsalt_maxu,&
             &      Qsalt_maxv,Qsalt,Qsalt_u,Qsalt_v,dh_salt,dh_salt_u,&
             &      dh_salt_v,conc_salt,Qsusp,Qsusp_u,Qsusp_v,dh_susp,&
             &      dh_susp_u,dh_susp_v,Qsubl,l_dh_subgrid)

     endif

  else

! c This 'noblowsnow' call zeros out data from a previous time step
! c   that had blowing snow.
     call noblowsnow(nx,ny,Qsalt_max,Qsalt_maxu,&
          &    Qsalt_maxv,Qsalt,Qsalt_u,Qsalt_v,dh_salt,dh_salt_u,&
          &    dh_salt_v,conc_salt,Qsusp,Qsusp_u,Qsusp_v,dh_susp,&
          &    dh_susp_u,dh_susp_v,Qsubl,l_dh_subgrid)

  endif

!  call single_gather(Qsusp_v, l_Qsusp_v)
!  call single_gather(Qsalt_v,l_Qsalt_v)


! c If this is a normal SnowTran-3D run, adjust the accumulations
! c   and erosions in reponse to varying transport fluxes across
! c   the simulation domain.  If it is a sea ice run with 25-km
! c   grid cells, just adjust the snow depth in response to the
! c   blowing snow sublimation fluxes only.

  if (seaice_run.eq.1.0 .and. deltax.eq.25000.0) then

! c Here Qsubl goes in as a flux, and comes out in snow depth units.
! c   And Qsubl_depth comes out in swe depth units.
     call bs_sublimation_only(nx,ny,dt,bs_flag,ro_water,&
     &    snow_depth,swe_depth,ro_snow_grid,soft_snow_d,Qsubl,&
     &    vegsnowd_xy,Qsubl_depth)


  elseif (seaice_run.eq.4.0 .and. deltax.eq.25000.0) then
! c Here Qsubl goes in as a flux, and comes out in snow depth units.
! c   And Qsubl_depth comes out in swe depth units.
     call bs_sublimation_only(nx,ny,dt,bs_flag,ro_water,&
          &    snow_depth,swe_depth,ro_snow_grid,soft_snow_d,Qsubl,&
          &    vegsnowd_xy,Qsubl_depth)

  else

! c Compute the new snow depth due to accumulation from precipitation,
! c   saltation, and suspension, and the mass loss due to
! c   sublimation.


     call accum(snow_d,nx,ny,ro_snow,dt,ro_water,&
          &    deltax,deltay,vegsnowd_xy,Tabler_2_flag,&
          &    index_ue,index_uw,index_vn,index_vs,undef,&
          &    Qsalt_u,Qsalt_v,Qsusp_u,Qsusp_v,Qsubl,dh_salt,&
          &    dh_salt_u,dh_salt_v,dh_susp,dh_susp_u,dh_susp_v,&
          &    wbal_qsubl,wbal_salt,wbal_susp,bs_flag,&
          &    soft_snow_d,topo,topo_land,topoflag,subgrid_flag,&
          &    tabler_nn,tabler_ss,tabler_ee,tabler_ww,&
          &    tabler_ne,tabler_se,tabler_sw,tabler_nw,&
          &    uwind_grid,vwind_grid,wbal_subgrid,sum_qsubl,&
          &    sum_trans,swe_depth,snow_depth,ro_snow_grid,&
          &    l_dh_subgrid,tabler_dir,iter,slope_adjust,&
          &    curve_len_scale)


  endif

! c Use the changes in swe due to saltation, suspension, and
! c   blowing snow sublimation to adjust the multilayer snowpack
! c   layers.
  if (multilayer_snowpack.eq.1) then

     if (seaice_run.ne.0.0) then
        do j=1,l_ny
           do i=1,nx
! c This is the sublimation in terms of swe depth.
              wbal_qsubl(i,j) = Qsubl(i,j) * ro_snow_grid(i,j) / &
                   &          ro_water
! c Because these simulations are done over grid cells that are
! c   around 25-km by 25-km, subgrid blowing snow processes are not
! c   simulated.
              wbal_salt(i,j) = 0.0
              wbal_susp(i,j) = 0.0
              wbal_subgrid(i,j) = 0.0
           enddo
        enddo
     endif

     do j=1,l_ny
!        new_j = global_map(j)
        do i=1,nx
           swe_change = wbal_qsubl(i,j) + wbal_salt(i,j) + &
                &        wbal_susp(i,j) + wbal_subgrid(i,j)

! c Net mass loss for this grid cell at this time step.
           if (swe_change.lt.0.0) then
              swe_change_tmp = -swe_change

! c Extract the vertical column for this i,j point, and send it
! c   to the subroutine. *** Note that I should use f95, then I would
! c   not have to do this (I could pass in subsections of the arrays).
              do k=1,nz_max
                 swed_layer_z(k) = swed_layer(i,j,k)
              enddo

! c Check to see whether a layer reduction is required.
              CALL REDUCE_LAYERS(swe_change_tmp,swed_layer_z,KK(i,j))

! c Re-build the 3-D array.  See note above about using f95 to avoid this.
              do k=1,nz_max
                 swed_layer(i,j,k) = swed_layer_z(k)
              enddo

! c Update the snow layer thicknesses, and recalculate the total
! c   snow and swe depths.  Assume this swe change does not change
! c   the snow density and does not change the soft snow depth.  It
! c   only reduces the snow depth and the associated swe depth.
              snow_depth(i,j) = 0.0
              swe_depth(i,j) = 0.0
              do k=1,KK(i,j)
                 snod_layer(i,j,k) = swed_layer(i,j,k) * ro_water / &
                      &            ro_layer(i,j,k)
! c               ro_layer(i,j,k) = ro_layer(i,j,k)
                 snow_depth(i,j) = snow_depth(i,j) + snod_layer(i,j,k)
                 swe_depth(i,j) = swe_depth(i,j) + swed_layer(i,j,k)
              enddo

!c Net mass gain for this grid cell at this time step.
           elseif (swe_change.gt.0.0) then

!c Add to the existing top layer.
              swed_layer(i,j,KK(i,j)) = swed_layer(i,j,KK(i,j)) + &
                   &          swe_change
!c             ro_layer(i,j,k) = ro_layer(i,j,k)
              snod_layer(i,j,KK(i,j)) = swed_layer(i,j,KK(i,j)) * &
                   &          ro_water / ro_layer(i,j,KK(i,j))

! c Update the snow layer thicknesses, and recalculate the total
! c   snow and swe depths.  Assume this swe change does not change
! c   the snow density and does not change the soft snow depth.  It
! c   only reduces the snow depth and the associated swe depth.
              snow_depth(i,j) = 0.0
              swe_depth(i,j) = 0.0
              do k=1,KK(i,j)
                 snow_depth(i,j) = snow_depth(i,j) + snod_layer(i,j,k)
                 swe_depth(i,j) = swe_depth(i,j) + swed_layer(i,j,k)
              enddo

           else

              snow_depth(i,j) = 0.0
              swe_depth(i,j) = 0.0
              do k=1,KK(i,j)
                 snow_depth(i,j) = snow_depth(i,j) + snod_layer(i,j,k)
                 swe_depth(i,j) = swe_depth(i,j) + swed_layer(i,j,k)
              enddo

           endif

        enddo
     enddo
  endif

!c Perform a water balance check (see notes in this subroutine).
  if (seaice_run.eq.0.0) then

! c Don't do this calculation for subgrid_flag = 1.0.  This can be
! c   turned back on and the required adjustments made to the code,
! c   if needed.
       if (subgrid_flag.eq.0.0) then
         call waterbal_snowtran(w_balance,prec,Qcs,&
    &      runoff,d_canopy_int,swe_depth,glacier_melt,iter,&
    &      wbal_qsubl,wbal_salt,wbal_susp,wbal_subgrid,nx,l_ny,&
    &      swe_depth_old,swesublim,canopy_unload,canopy_int,&
    &      canopy_int_old)
       endif
  endif


! c If this is a sea ice run, zero out the ocean grid cells that
! c   have no sea ice in them.
  if (seaice_run.ne.0.0) then
     CALL ZERO_SEAICE_SNOW(nx,ny,snow_depth,ro_snow_grid,&
          &    ro_snow,swe_depth,swe_depth_old,canopy_int_old,KK,&
          &    tslsnowfall,snod_layer,swed_layer,ro_layer,T_old,&
          &    multilayer_snowpack,tsls_threshold,seaice_conc,&
          &    sum_sprec,sum_trans)
  endif

!c Save the mass balance variables from this time step.
  do j=1,l_ny
     do i=1,nx
        swe_depth_old(i,j) = swe_depth(i,j)
        canopy_int_old(i,j) = canopy_int(i,j)
     enddo
  enddo

! c In SnowTran-3D, the summed snow precipitation were in units
! c   of snow-depth.  The rest of the routines assume that it is in
! c   swe units.
  do j=1,l_ny
     do i=1,nx
        sum_sprec(i,j) = sum_sprec(i,j) * ro_snow / ro_water
     enddo
  enddo

! c Close the Tabler surface output file if this is the end of this
! c   assimulation loop.
  if (iter.eq.max_iter) close(52)

  return
end subroutine SNOWTRAN_CODE

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine waterbal_snowtran(w_balance,prec,Qcs,&
      &  runoff,d_canopy_int,swe_depth,glacier_melt,iter,&
      &  wbal_qsubl,wbal_salt,wbal_susp,wbal_subgrid,nx,ny,&
      &  swe_depth_old,swesublim,canopy_unload,canopy_int,&
      &  canopy_int_old)

   use snowmodel_inc
   implicit none

      !include 'snowmodel.inc'

   integer iter,nx,ny,i,j

   real w_balance(nx,ny),prec(nx,ny),&
        &  Qcs(nx,ny),runoff(nx,ny),&
        &  d_canopy_int(nx,ny),swe_depth(nx,ny),&
        &  glacier_melt(nx,ny),wbal_qsubl(nx,ny),&
        &  wbal_salt(nx,ny),swe_depth_old(nx,ny),&
        &  swesublim(nx,ny),wbal_susp(nx,ny),&
        &  wbal_subgrid(nx,ny),canopy_unload(nx,ny),&
        &  canopy_int_old(nx,ny),canopy_int(nx,ny)

! c Note that the following balances should hold.  These aren't quite
! c   right, but it is a place to start.
! c   Canopy Balance (forest):
! c     canopy = sprec - unload + Qcs ==> unload = sprec - canopy + Qcs
! c
! c   Snowpack Balance (forest):
! c     swe_d = unload + rain - runoff ==>
! c       canopy + swe_d = sprec + rain + Qcs - runoff
! c     prec = sprec + rain
! c     sum_rain  = sum_sprec - sum_prec
! c
! c   Snowpack Balance (non-forest):
! c     swe_d = sprec + rain - runoff + subl + salt + susp + subgrid +
! c       glaciermelt
! c
! c   Everywhere:
! c     w_balance = sum_prec + sum_Qcs - sum_runoff + sum_subl +
! c       sum_trans - canopy_int - swe_depth + sum_glacmelt
! c
! c   The related variables that would need to be brought in are:
! c      d_canopy_int,sum_d_canopy_int,sum_unload

! c The subroutine WATERBAL_SNOWTRAN is used if the model simulation
   ! c   includes SnowTran-3D.
   do j=1,ny
      do i=1,nx
         w_balance(i,j) = swe_depth_old(i,j) - swe_depth(i,j) + &
              &      prec(i,j) - runoff(i,j) + glacier_melt(i,j) + &
              &      wbal_qsubl(i,j) + wbal_salt(i,j) + wbal_susp(i,j) + &
              &      wbal_subgrid(i,j) - swesublim(i,j) + canopy_int_old(i,j) - &
              &      canopy_int(i,j) + Qcs(i,j)

         if (abs(w_balance(i,j)).gt.1.0e-5) &
              &      print*,'water imbalance at iter =',iter,' ',w_balance(i,j)

      enddo
   enddo

   return
 end subroutine waterbal_snowtran

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine noblowsnow(nx,ny,Qsalt_max,Qsalt_maxu,&
      &  Qsalt_maxv,Qsalt,Qsalt_u,Qsalt_v,dh_salt,dh_salt_u,&
      &  dh_salt_v,conc_salt,Qsusp,Qsusp_u,Qsusp_v,dh_susp,&
      &  dh_susp_u,dh_susp_v,Qsubl,dh_subgrid)

   use caf_module, only: max_l_ny, l_ny
   use snowmodel_inc
   implicit none

      !include 'snowmodel.inc'

   integer nx,ny,i,j
   real Qsalt_max(nx,max_l_ny)
   real Qsalt_maxu(nx,max_l_ny),Qsalt_maxv(nx,max_l_ny)
   real Qsalt(nx,max_l_ny)
   real Qsalt_u(nx,max_l_ny),Qsalt_v(nx,max_l_ny)
   real dh_salt(nx,max_l_ny)
   real dh_salt_u(nx,max_l_ny),dh_salt_v(nx,max_l_ny)
   real conc_salt(nx,max_l_ny)
   real Qsusp(nx,max_l_ny)
   real Qsusp_u(nx,max_l_ny),Qsusp_v(nx,max_l_ny)
   real dh_susp(nx,max_l_ny)
   real dh_susp_u(nx,max_l_ny),dh_susp_v(nx,max_l_ny)
   real dh_subgrid(nx,max_l_ny)
   real Qsubl(nx,max_l_ny)

!   integer new_j

   conc_salt = 0.0
   Qsusp = 0.0
   Qsubl = 0.0
   dh_salt = 0.0
   dh_salt_u = 0.0
   dh_salt_v = 0.0

   do i=1,nx
      do j=1,l_ny
!         new_j = global_map(j)
         Qsalt_max(i,j) = 0.0
         Qsalt_maxu(i,j) = 0.0
         Qsalt_maxv(i,j) = 0.0
         Qsalt(i,j) = 0.0
         Qsalt_u(i,j) = 0.0
         Qsalt_v(i,j) = 0.0
         Qsusp_u(i,j) = 0.0
         Qsusp_v(i,j) = 0.0
         dh_susp(i,j) = 0.0
         dh_susp_u(i,j) = 0.0
         dh_susp_v(i,j) = 0.0
         dh_subgrid(i,j) = 0.0
      enddo
   enddo

   return
 end subroutine noblowsnow

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine bs_sublimation_only(nx,ny,dt,bs_flag,ro_water,&
      &  snow_depth,swe_depth,ro_snow_grid,soft_snow_d,Qsubl,&
      &  vegsnowd_xy,Qsubl_depth)

   use snowmodel_inc
   use caf_module, only: max_l_ny, l_ny, global_map
   implicit none

!      include 'snowmodel.inc'

   integer i,j,nx,ny

   real dt,bs_flag,ro_water,snowdmin,hard_snow_d

   real snow_depth(nx,max_l_ny)
   real swe_depth(nx,max_l_ny)
   real ro_snow_grid(nx,max_l_ny)

   real soft_snow_d(nx,max_l_ny)
   real Qsubl(nx,max_l_ny)
   real Qsubl_depth(nx,max_l_ny)

   real vegsnowd_xy(nx,max_l_ny)

! c Adjust the snow and swe depths to account for blowing-snow
! c   sublimation.

   if (bs_flag.eq.1.0) then

! c SUBLIMATION
! c Make adjustments for the case where there is no snow available
! c   on the ground (or captured within the vegetation) to be
! c   eroded.  Since Qsubl is blowing snow sublimation, don't let
! c   this sublimation reach down into the hard snow layer or into
! c   the vegsnowd, because that snow is not available to be blown
! c   around.

      do i=1,nx
         do j=1,l_ny
            hard_snow_d = snow_depth(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)

!c Convert Qsubl from sublimation flux to sublimated snow depth.
!c   Qsubl is negative here.
            Qsubl(i,j) = Qsubl(i,j) * dt / ro_snow_grid(i,j)

            if (snow_depth(i,j).gt.snowdmin) then
               if (snow_depth(i,j)+Qsubl(i,j).le.snowdmin) then
                  Qsubl(i,j) = snowdmin - snow_depth(i,j)
               endif
            else
               Qsubl(i,j) = 0.0
            endif
         enddo
      enddo

!c Account for decreases in snow depth due to sublimation.
      do i=1,nx
         do j=1,l_ny
            snow_depth(i,j) = snow_depth(i,j) + Qsubl(i,j)
            soft_snow_d(i,j) = soft_snow_d(i,j) + Qsubl(i,j)
         enddo
      enddo

   else

      do i=1,nx
         do j=1,l_ny
            Qsubl(i,j) = 0.0
         enddo
      enddo

   endif

! c Convert any snow-depth adjustments that occurred above to
! c   swe using the spatially-distributed snow density from the
! c   snowpack model.  Also convert the Qsubl_depth values from
! c   snow depth to swe depth units.
   do i=1,nx
      do j=1,l_ny
         swe_depth(i,j) = snow_depth(i,j) * ro_snow_grid(i,j) / &
              &      ro_water
         Qsubl_depth(i,j) = Qsubl(i,j) * ro_snow_grid(i,j) / &
              &      ro_water
      enddo
   enddo

   return
 end subroutine bs_sublimation_only

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine accum(snow_d,nx,ny,ro_snow,dt,ro_water,&
      &  deltax,deltay,vegsnowd_xy,Tabler_2_flag,&
      &  index_ue,index_uw,index_vn,index_vs,undef,&
      &  Qsalt_u,Qsalt_v,Qsusp_u,Qsusp_v,Qsubl,dh_salt,&
      &  dh_salt_u,dh_salt_v,dh_susp,dh_susp_u,dh_susp_v,&
      &  wbal_qsubl,wbal_salt,wbal_susp,bs_flag,&
      &  soft_snow_d,topo,topo_land,topoflag,subgrid_flag,&
      &  tabler_nn,tabler_ss,tabler_ee,tabler_ww,&
      &  tabler_ne,tabler_se,tabler_sw,tabler_nw,&
      &  uwind_grid,vwind_grid,wbal_subgrid,sum_qsubl,&
      &  sum_trans,swe_depth,snow_depth,ro_snow_grid,&
      &  dh_subgrid,tabler_dir,iter,slope_adjust,&
      &  curve_len_scale)

   use caf_module, only: me, max_l_ny, l_ny, single_gather, single_scatter,&
                        & global_map
   use snowmodel_inc
   use snowmodel_vars, only:  l_vwind_grid,&
                            & l_soft_snow_d_dep, l_snow_d_dep,&
                            & l_Qsalt_v, l_Qsusp_v, &
                            & l_topo_tmp, l_snow_d_tabler,&
                            & l_dh_dep, l_snow_d_tmp
   implicit none

!      include 'snowmodel.inc'

   integer i,j,nx,ny,iter,k,loops_snowd_smoother,II,JJ,nnx,nny,&
        &  irotate_flag, new_j

   real ro_snow,dt,deltax,deltay,bs_flag,topoflag,ro_water
   real snowdmin,hard_snow_d,subgrid_flag,tabler_dir
   real slope_adjust,xmult,curve_len_scale,Tabler_2_flag

   real snow_d(nx,max_l_ny)
   real snow_depth(nx,max_l_ny)
   real swe_depth(nx,max_l_ny)
   real ro_snow_grid(nx,max_l_ny)

   real tabler_nn(nx,max_l_ny)
   real tabler_ss(nx,max_l_ny)
   real tabler_ee(nx,max_l_ny)
   real tabler_ww(nx,max_l_ny)
   real tabler_ne(nx,max_l_ny)
   real tabler_se(nx,max_l_ny)
   real tabler_sw(nx,max_l_ny)
   real tabler_nw(nx,max_l_ny)
   real uwind_grid(nx,max_l_ny)
   real vwind_grid(nx,ny)

   real soft_snow_d(nx,max_l_ny)
   real Qsubl(nx,max_l_ny)
   real topo(nx,ny)
   real topo_land(nx,ny)

   real dh_salt(nx,max_l_ny)
   real dh_salt_u(nx,max_l_ny)
   real dh_salt_v(nx,max_l_ny)

   real dh_susp(nx,max_l_ny)
   real dh_susp_u(nx,max_l_ny)
   real dh_susp_v(nx,max_l_ny)

   real dh_subgrid(nx,max_l_ny)

   real Qsalt_u(nx,max_l_ny)
   real Qsalt_v(nx,max_l_ny)

   real Qsusp_u(nx,max_l_ny)
   real Qsusp_v(nx,max_l_ny)

   real wbal_qsubl(nx,max_l_ny)
   real wbal_salt(nx,max_l_ny)
   real wbal_susp(nx,max_l_ny)
   real wbal_subgrid(nx,max_l_ny)
   real sum_qsubl(nx,max_l_ny)
   real sum_trans(nx,max_l_ny)

   real vegsnowd_xy(nx,max_l_ny)

   integer index_ue(ny_max,2*nx_max+1)
   integer index_uw(ny_max,2*nx_max+1)
   integer index_vn(nx_max,2*ny_max+1)
   integer index_vs(nx_max,2*ny_max+1)

   real extra,space,fill,undef
   real pi,rad2deg,bsflux_dir,bs_flux_u,bs_flux_v,subgrid_dir

!c Define the required constants.
   pi = 2.0 * acos(0.0)
   rad2deg = 180.0 / pi

!c COMPUTE THE NEW SNOW DEPTH.

! c PRECIPITATION
! c Account for the addition due to snow precipitation.
! c This is now updated at the beginning of the program (day).

!c Sum the precipitation in terms of snow depth.
   if (bs_flag.eq.1.0) then

!c SALTATION
      call getnewdepth(nx,ny,deltax,deltay,Qsalt_u,&
           &    Qsalt_v,dh_salt_u,dh_salt_v,&
           &    ro_snow,dt,vegsnowd_xy,snow_d,&
           &    soft_snow_d,l_snow_d_dep,l_soft_snow_d_dep,subgrid_flag,&
           &    uwind_grid, vwind_grid)
!      call single_gather(Qsalt_v, l_Qsalt_v)

!c This saves the deposition part of the snow-depth change.
      if (subgrid_flag.eq.1.0) then
         do i=1,nx
            do j=1,l_ny
               l_dh_dep(i,j) = l_snow_d_dep(i,j)
            enddo
         enddo
      endif

      do i=1,nx
         do j=1,l_ny
            dh_salt(i,j) = dh_salt_u(i,j) + dh_salt_v(i,j)
         enddo
      enddo

!c SUSPENSION
      call getnewdepth(nx,ny,deltax,deltay,Qsusp_u,&
           &    Qsusp_v,dh_susp_u,dh_susp_v,&
           &    ro_snow,dt,vegsnowd_xy,snow_d,&
           &    soft_snow_d,l_snow_d_dep,l_soft_snow_d_dep,subgrid_flag,&
           &    uwind_grid, vwind_grid)
!      call single_gather(Qsusp_v, l_Qsusp_v)

!c This adds the suspension part the deposition to the
!c   saltation part of the snow-depth change.
      if (subgrid_flag.eq.1.0) then
         do i=1,nx
            do j=1,l_ny
               l_dh_dep(i,j) = l_dh_dep(i,j) + l_snow_d_dep(i,j)
            enddo
         enddo
      endif

      do i=1,nx
         do j=1,l_ny
            dh_susp(i,j) = dh_susp_u(i,j) + dh_susp_v(i,j)
         enddo
      enddo

!c Save a copy of the snow distribution to be used to calculate the
!c   snow distribution changes resulting from the subgrid
!c   redistribution.
      do i=1,nx
         do j=1,l_ny
            l_snow_d_tmp(i,j) = snow_d(i,j)
         enddo
      enddo

!c Run the subgrid parameterization to account for unrealistic
!c   snow accumulation spikes.
      if (subgrid_flag.eq.1.0) then

! c Do the Tabler corrections while assuming the constant snow density
! c   defined in snowmodel.par (ro_snow).  This is done to avoid all
! c   of the messy tracking of different density snows that is being
! c   blown around by the wind.  This simplification could be relaxed
! c   at some point, but it is not an easy accounting, and I don't
! c   think it is justified at this point, given everything else that
! c   we don't know.  snow_d is coming in from SnowPack where the
! c   ro_snow adjustment has already been made.  When all of the Tabler
! c   calculations are finished, then we will convert back from the
! c   SnowTran constant density convention to the SnowPack ro_snow_grid
! c   spatially distributed snow density.  Note that here, that coming
! c   from SnowPack, if the snow depth equals zero, the snow density
! c   is undefined.  This can cause problems for the case where Tabler
! c   has moved snow into a previously snow-free grid cell; that grid
! c   cell will also need a snow density assigned to it.  Right now I
! c   am assigning that value to be ro_snow from snowmodel.par.
         do i=1,nx
            do j=1,l_ny
               l_snow_d_tabler(i,j) = snow_d(i,j)
            enddo
         enddo

! c The following subgrid_1 subroutine was saved because it provides
! c   an example of how to relax the assumption of a single wind
! c   direction over the entire domain at a given time step.  This
! c   is the old subgrid routine that is no longer used or supported.
! c         if (subgrid_flag.eq.1.0) then
! c           call subgrid_1(nx,ny,snow_d_tabler,
! c    &        index_ue,index_uw,index_vn,index_vs,
! c    &        tabler_nn,tabler_ss,tabler_ee,tabler_ww,
! c    &        tabler_ne,tabler_se,tabler_sw,tabler_nw,uwind_grid,
! c    &        vwind_grid,tabler_dir)
! c         endif

! c Calculate the integrated snow-transport direction.

! c Initialize the summing arrays.
         bs_flux_u = 0.0
         bs_flux_v = 0.0

! c Sum the fluxes in the u-v directions.  Add the u-v sign to the
! c   fluxes so there is a direction associated with them.
         do i=1,nx
            do j=1,l_ny
               new_j = global_map(j)
               bs_flux_u = bs_flux_u + &
                    &          sign(Qsalt_u(i,j),uwind_grid(i,j)) + &
                    &          sign(Qsusp_u(i,j),uwind_grid(i,j))
               bs_flux_v = bs_flux_v + &
                    &          sign(Qsalt_v(i,j),vwind_grid(i,new_j)) + &
                    &          sign(Qsusp_v(i,j),vwind_grid(i,new_j))
            enddo
         enddo
         call co_sum(bs_flux_u)
         call co_sum(bs_flux_v)

! c Calculate the resulting direction.  Some compilers do not
! c   allow both u and v to be 0.0 in the atan2 computation.
         if (abs(bs_flux_u).lt.1e-10) bs_flux_u = 1e-10
         bsflux_dir = rad2deg * atan2(bs_flux_u,bs_flux_v)
         if (bsflux_dir.ge.180.0) then
            bsflux_dir = bsflux_dir - 180.0
         else
            bsflux_dir = bsflux_dir + 180.0
         endif

         print *,bsflux_dir,bs_flux_u,bs_flux_v
! c         print *,bsflux_dir,bs_flux_u,bs_flux_v
! c         print *,bsflux_dir,bs_flux_u,bs_flux_v
! c         print *,bsflux_dir,bs_flux_u,bs_flux_v

! c Decide whether tabler_dir or bsflux_dir is going to be used to do
! c   the Tabler-surface snow redistributions.
         if (tabler_dir.lt.0.0) then
            subgrid_dir = bsflux_dir
!c           print *, 'using bsflux_dir'
         else
            subgrid_dir = tabler_dir
            !c           print *, 'using tabler_dir'
         endif

! c Smooth the snow distribution to eliminate any sharp wind speed
! c   variations computed at the next time step.  Define the number
! c   of times this is done to be a function of the curvature length
! c   scale and the grid increment.  Unwanted peaks in the snow
! c   distributions can be further eliminated by adjusting the
! c   multiplication factor, xmult.  Also see "loops_windwt_smoother"
! c   in micromet_code.f.
! c         xmult = 0.15
         xmult = 0.25
!c         xmult = 0.5
         loops_snowd_smoother = nint(xmult * curve_len_scale / &
              &      (0.5 * (deltax + deltay)))

! c         print *
! c         print *, 'loops_snowd_smoother ',loops_snowd_smoother
! c         print *

!c Don't do this smoothing if the domain is arbitrarily small.
         if (nx.gt.100 .and. ny.gt.100) then
            do k=1,loops_snowd_smoother
               call smoother9(nx,l_ny,l_snow_d_tabler)
            enddo
         endif

! c Create Tabler surfaces from the snow on the ground at this
! c   time step, with just erosion taken into account (that's
! c   where we are at this point).  Add that current snow depth
! c   to topo_land.  Then use this to create the Tabler surface that
! c   will be used for this time step.  Also define this depth to be
! c   dependent on the SnowPack spatially distributed snow density,
! c   not the constant density used in SnowTran.
         do i=1,nx
            do j=1,l_ny
               new_j = global_map(j)
               l_topo_tmp(i,j) = l_snow_d_tabler(i,j) + topo_land(i,new_j)
            enddo
         enddo

! c The following loop does four things:

! c (1) Extract the Tabler surface for the direction of interest at
! c     this time step.

! c (2) The Tabler surfaces that were just generated have had topo_tmp
! c     subtracted off of them, giving just the drift profiles with
! c     things like zero drift depth on ridges and windwards slopes.
! c     So, add the snow depth, prior to any wind redistribution, to
! c     these Tabler surfaces.  This will be the maximum snow depth
! c     allowed as part of the wind redistribution.

! c (3) Set the snow-free areas equal to the snow-holding depth.

! c (4) Save the calculated Tabler surface at this time step.  You can
! c     comment this out if you don't want to write them out.

! c nn.
! c Consider N winds.
         if (subgrid_dir.gt.337.5 .and. subgrid_dir.le.360.0 .or. &
              &      subgrid_dir.ge.0.0 .and. subgrid_dir.le.22.5) then
! c           print *,'in nn'

!c Extract the Tabler surface.
            irotate_flag = 1
            call tabler_n(nx,l_ny,l_topo_tmp,tabler_nn,deltay, &
                 &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
            do i=1,nx
               do j=1,l_ny
                  tabler_nn(i,j) = l_snow_d_tabler(i,j) + tabler_nn(i,j)
                  tabler_nn(i,j) = max(tabler_nn(i,j),vegsnowd_xy(i,j))
               enddo
            enddo

!c Save the Tabler surface.
            if (Tabler_2_flag.eq.1.0) &
                 &        write(52,rec=iter) ((tabler_nn(i,j),i=1,nx),j=1,ny)

! c ne.
! c Consider NE winds.
         elseif (subgrid_dir.gt.22.5 .and. subgrid_dir.le.67.5) then
!c           print *,'in ne'

!c Extract the Tabler surface.
            irotate_flag = 2
            call tabler_e(nx,l_ny,l_topo_tmp,tabler_ne,1.41*deltax, &
                 &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
            do i=1,nx
               do j=1,l_ny
                  tabler_ne(i,j) = l_snow_d_tabler(i,j) + tabler_ne(i,j)
                  tabler_ne(i,j) = max(tabler_ne(i,j),vegsnowd_xy(i,j))
               enddo
            enddo

!c Save the Tabler surface.
            if (Tabler_2_flag.eq.1.0) &
                 &        write(52,rec=iter) ((tabler_ne(i,j),i=1,nx),j=1,ny)

!c ee.
!c Consider E winds.
         elseif (subgrid_dir.gt.67.5 .and. subgrid_dir.le.112.5) then
!c           print *,'in ee'

!c Extract the Tabler surface.
            irotate_flag = 1
            call tabler_e(nx,l_ny,l_topo_tmp,tabler_ee,deltax, &
                 &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
            do i=1,nx
               do j=1,l_ny
                  tabler_ee(i,j) = l_snow_d_tabler(i,j) + tabler_ee(i,j)
                  tabler_ee(i,j) = max(tabler_ee(i,j),vegsnowd_xy(i,j))
               enddo
            enddo

!c Save the Tabler surface.
            if (Tabler_2_flag.eq.1.0) &
                 &        write(52,rec=iter) ((tabler_ee(i,j),i=1,nx),j=1,ny)

!c se.
!c Consider SE winds.
         elseif(subgrid_dir.gt.112.5 .and. subgrid_dir.le.157.5)then
!c           print *,'in se'

!c Extract the Tabler surface.
            irotate_flag = 2
            call tabler_s(nx,l_ny,l_topo_tmp,tabler_se,1.41*deltay,&
                 &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
            do i=1,nx
               do j=1,l_ny
                tabler_se(i,j) = l_snow_d_tabler(i,j) + tabler_se(i,j)
                tabler_se(i,j) = max(tabler_se(i,j),vegsnowd_xy(i,j))
             enddo
          enddo

!c Save the Tabler surface.
          if (Tabler_2_flag.eq.1.0) &
               &        write(52,rec=iter) ((tabler_se(i,j),i=1,nx),j=1,ny)

!c ss.
!c Consider S winds.
       elseif(subgrid_dir.gt.157.5 .and. subgrid_dir.le.202.5)then
!c           print *,'in ss'

!c Extract the Tabler surface.
          irotate_flag = 1
          call tabler_s(nx,l_ny,l_topo_tmp,tabler_ss,deltay,&
               &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
          do i=1,nx
             do j=1,l_ny
                tabler_ss(i,j) = l_snow_d_tabler(i,j) + tabler_ss(i,j)
                tabler_ss(i,j) = max(tabler_ss(i,j),vegsnowd_xy(i,j))
             enddo
          enddo

!c Save the Tabler surface.
          if (Tabler_2_flag.eq.1.0) &
               &        write(52,rec=iter) ((tabler_ss(i,j),i=1,nx),j=1,ny)

!c sw.
!c Consider SW winds.
       elseif(subgrid_dir.gt.202.5 .and. subgrid_dir.le.247.5)then
!c           print *,'in sw'

!c Extract the Tabler surface.
          irotate_flag = 2
          call tabler_w(nx,l_ny,l_topo_tmp,tabler_sw,1.41*deltax,&
               &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
          do i=1,nx
             do j=1,l_ny
                tabler_sw(i,j) = l_snow_d_tabler(i,j) + tabler_sw(i,j)
                tabler_sw(i,j) = max(tabler_sw(i,j),vegsnowd_xy(i,j))
             enddo
          enddo

!c Save the Tabler surface.
          if (Tabler_2_flag.eq.1.0) &
               &        write(52,rec=iter) ((tabler_sw(i,j),i=1,nx),j=1,ny)

!c ww.
!c Consider W winds.
       elseif(subgrid_dir.gt.247.5 .and. subgrid_dir.le.292.5)then
!c           print *,'in ww'

!c Extract the Tabler surface.
          irotate_flag = 1
          call tabler_w(nx,l_ny,l_topo_tmp,tabler_ww,deltax,&
               &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
          do i=1,nx
             do j=1,l_ny
                tabler_ww(i,j) = l_snow_d_tabler(i,j) + tabler_ww(i,j)
                tabler_ww(i,j) = max(tabler_ww(i,j),vegsnowd_xy(i,j))
             enddo
          enddo

!c Save the Tabler surface.
          if (Tabler_2_flag.eq.1.0) &
               &        write(52,rec=iter) ((tabler_ww(i,j),i=1,nx),j=1,ny)

!c nw.
!c Consider NW winds.
       elseif(subgrid_dir.gt.292.5 .and. subgrid_dir.le.337.5)then
!c           print *,'in nw'

!c Extract the Tabler surface.
          irotate_flag = 2
          call tabler_n(nx,l_ny,l_topo_tmp,tabler_nw,1.41*deltay,&
               &        irotate_flag,slope_adjust)

!c Add the snow depth back on, and clip to the snow-holding depth.
          do i=1,nx
             do j=1,l_ny
                tabler_nw(i,j) = l_snow_d_tabler(i,j) + tabler_nw(i,j)
                tabler_nw(i,j) = max(tabler_nw(i,j),vegsnowd_xy(i,j))
             enddo
          enddo

!c Save the Tabler surface.
          if (Tabler_2_flag.eq.1.0) &
               &        write(52,rec=iter) ((tabler_nw(i,j),i=1,nx),j=1,ny)

       else
          print *,'subgrid_dir not found'
          stop

       endif

! c Now add the +dh snow back on, and do the sweep that does not
! c   allow any snow to exist above the Tabler surface.
       do i=1,nx
          do j=1,l_ny
             l_snow_d_tabler(i,j) = l_snow_d_tabler(i,j) + l_dh_dep(i,j)
          enddo
       enddo

! c Sweep across the domain in the direction defined by tabler_dir
! c   or bsflux_dir (see above; now defined by subgrid_dir).  This is
! c   done in the same direction across the entire domain.  If using
! c   tabler_dir, then it is the same direction for all time steps.
! c   If using bsflux_dir, then it varies for each time step depending
! c   on the dominant transport direction.  If you ever want to relax
! c   this "same direction over the entire domain" assumption, see
! c   the subgrid_1 subroutine (search subgrid_1 in this document; it
! c   is all commented out).

!c nn.
!c Consider N winds.
       if (subgrid_dir.gt.337.5 .and. subgrid_dir.le.360.0 .or. &
            &      subgrid_dir.ge.0.0 .and. subgrid_dir.le.22.5) then
!c           print *,'in nn'


!c Do the sweep in the direction of interest.
          do i=1,nx
             extra = 0.0
             do j=l_ny,1,-1
                if (l_snow_d_tabler(i,j).ge.tabler_nn(i,j)) then
                   extra = extra + l_snow_d_tabler(i,j) - tabler_nn(i,j)
                   l_snow_d_tabler(i,j) = tabler_nn(i,j)
                else
                   space = tabler_nn(i,j) - l_snow_d_tabler(i,j)
                   fill = min(extra,space)
                   l_snow_d_tabler(i,j) = l_snow_d_tabler(i,j) + fill
                   extra = extra - fill
                   if (extra.lt.0.0) print *,'extra < 0.0',extra
                endif
             enddo
          enddo

!c ne.
!c Consider NE winds.
       elseif (subgrid_dir.gt.22.5 .and. subgrid_dir.le.67.5) then
!c           print *,'in ne'

!c Do the sweep in the direction of interest.
          nny = nx + l_ny - 1
          do j=1,nny
             extra = 0.0
             do i=nx,1,-1
                JJ = j + i - nx
                if (JJ.ge.1 .and. JJ.le.l_ny) then
                   if (l_snow_d_tabler(i,JJ).ge.tabler_ne(i,JJ)) then
                      extra = extra+l_snow_d_tabler(i,JJ)-tabler_ne(i,JJ)
                      l_snow_d_tabler(i,JJ) = tabler_ne(i,JJ)
                   else
                      space = tabler_ne(i,JJ) - l_snow_d_tabler(i,JJ)
                      fill = min(extra,space)
                      l_snow_d_tabler(i,JJ) = l_snow_d_tabler(i,JJ) + fill
                      extra = extra - fill
                      if (extra.lt.0.0) print *,'extra < 0.0',extra
                   endif
                endif
             enddo
          enddo

!c ee.
!c Consider E winds.
       elseif (subgrid_dir.gt.67.5 .and. subgrid_dir.le.112.5) then
!c           print *,'in ee'

!c Do the sweep in the direction of interest.
          do j=1,l_ny
             extra = 0.0
             do i=nx,1,-1
                if (l_snow_d_tabler(i,j).ge.tabler_ee(i,j)) then
                   extra = extra + l_snow_d_tabler(i,j) - tabler_ee(i,j)
                   l_snow_d_tabler(i,j) = tabler_ee(i,j)
                else
                   space = tabler_ee(i,j) - l_snow_d_tabler(i,j)
                   fill = min(extra,space)
                   l_snow_d_tabler(i,j) = l_snow_d_tabler(i,j) + fill
                   extra = extra - fill
                   if (extra.lt.0.0) print *,'extra < 0.0',extra
                endif
             enddo
          enddo

!c se.
!c Consider SE winds.
       elseif(subgrid_dir.gt.112.5 .and. subgrid_dir.le.157.5)then
!c           print *,'in se'

!c Do the sweep in the direction of interest.
          nnx = nx + l_ny - 1
          do i=1,nnx
             extra = 0.0
             do j=1,l_ny
                II = i - j + 1
                if (II.ge.1 .and. II.le.nx) then
                   if (l_snow_d_tabler(II,j).ge.tabler_se(II,j)) then
                      extra = extra+l_snow_d_tabler(II,j)-tabler_se(II,j)
                      l_snow_d_tabler(II,j) = tabler_se(II,j)
                   else
                      space = tabler_se(II,j) - l_snow_d_tabler(II,j)
                      fill = min(extra,space)
                      l_snow_d_tabler(II,j) = l_snow_d_tabler(II,j) + fill
                      extra = extra - fill
                      if (extra.lt.0.0) print *,'extra < 0.0',extra
                   endif
                endif
             enddo
          enddo

!c ss.
!c Consider S winds.
       elseif(subgrid_dir.gt.157.5 .and. subgrid_dir.le.202.5)then
!c           print *,'in ss'

!c Do the sweep in the direction of interest.
          do i=1,nx
             extra = 0.0
             do j=1,l_ny
                if (l_snow_d_tabler(i,j).ge.tabler_ss(i,j)) then
                   extra = extra + l_snow_d_tabler(i,j) - tabler_ss(i,j)
                   l_snow_d_tabler(i,j) = tabler_ss(i,j)
                else
                   space = tabler_ss(i,j) - l_snow_d_tabler(i,j)
                   fill = min(extra,space)
                   l_snow_d_tabler(i,j) = l_snow_d_tabler(i,j) + fill
                   extra = extra - fill
                   if (extra.lt.0.0) print *,'extra < 0.0',extra
                endif
             enddo
            enddo

!c sw.
!c Consider SW winds.
         elseif(subgrid_dir.gt.202.5 .and. subgrid_dir.le.247.5)then
!c           print *,'in sw'

!c Do the sweep in the direction of interest.
            nny = nx + l_ny - 1
            do j=1,nny
               extra = 0.0
               do i=1,nx
                  JJ = j + i - nx
                  if (JJ.ge.1 .and. JJ.le.l_ny) then
                     if (l_snow_d_tabler(i,JJ).ge.tabler_sw(i,JJ)) then
                        extra = extra+l_snow_d_tabler(i,JJ)-tabler_sw(i,JJ)
                        l_snow_d_tabler(i,JJ) = tabler_sw(i,JJ)
                     else
                        space = tabler_sw(i,JJ) - l_snow_d_tabler(i,JJ)
                        fill = min(extra,space)
                        l_snow_d_tabler(i,JJ) = l_snow_d_tabler(i,JJ) + fill
                        extra = extra - fill
                        if (extra.lt.0.0) print *,'extra < 0.0',extra
                     endif
                  endif
               enddo
            enddo

!c ww.
!c Consider W winds.
         elseif(subgrid_dir.gt.247.5 .and. subgrid_dir.le.292.5)then
!c           print *,'in ww'

!c Do the sweep in the direction of interest.
            do j=1,l_ny
               extra = 0.0
               do i=1,nx
                  if (l_snow_d_tabler(i,j).ge.tabler_ww(i,j)) then
                     extra = extra + l_snow_d_tabler(i,j) - tabler_ww(i,j)
                     l_snow_d_tabler(i,j) = tabler_ww(i,j)
                  else
                     space = tabler_ww(i,j) - l_snow_d_tabler(i,j)
                     fill = min(extra,space)
                     l_snow_d_tabler(i,j) = l_snow_d_tabler(i,j) + fill
                     extra = extra - fill
                     if (extra.lt.0.0) print *,'extra < 0.0',extra
                  endif
               enddo
            enddo

!c nw.
!c Consider NW winds.
         elseif(subgrid_dir.gt.292.5 .and. subgrid_dir.le.337.5)then
!c           print *,'in nw'

!c Do the sweep in the direction of interest.
            nnx = nx + l_ny - 1
            do i=1,nnx
               extra = 0.0
               do j=l_ny,1,-1
                  II = i - j + 1
                  if (II.ge.1 .and. II.le.nx) then
                     if (l_snow_d_tabler(II,j).ge.tabler_nw(II,j)) then
                        extra = extra+l_snow_d_tabler(II,j)-tabler_nw(II,j)
                        l_snow_d_tabler(II,j) = tabler_nw(II,j)
                     else
                        space = tabler_nw(II,j) - l_snow_d_tabler(II,j)
                        fill = min(extra,space)
                        l_snow_d_tabler(II,j) = l_snow_d_tabler(II,j) + fill
                        extra = extra - fill
                        if (extra.lt.0.0) print *,'extra < 0.0',extra
                     endif
                  endif
               enddo
            enddo

         else
            print *,'subgrid_dir not found'
            stop

         endif

!c Update the snow depths using the Tabler adjustments.
         do i=1,nx
            do j=1,l_ny
               snow_d(i,j) = l_snow_d_tabler(i,j)
            enddo
         enddo

      endif
!c End subgrid_flag = 1.0.

! c Calculate the snow depth resulting from the subgrid
! c   redistribution.
      do i=1,nx
         do j=1,l_ny
            dh_subgrid(i,j) = snow_d(i,j) - l_snow_d_tmp(i,j)
         enddo
      enddo

! c SUBLIMATION
! c Make adjustments for the case where there is no snow available
! c   on the ground (or captured within the vegetation) to be
! c   eroded.  Since Qsubl is blowing snow sublimation, don't let
! c   this sublimation reach down into the hard snow layer or into
! c   the vegsnowd, because that snow is not availble to be blown
! c   around.  I don't think it matters much whether this is done
! c   before or after the Tabler subgrid redistribution.
      do i=1,nx
         do j=1,l_ny
            hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)

!c Convert Qsubl flux to sublimated snow depth.
            Qsubl(i,j) = Qsubl(i,j) * dt / ro_snow

            if (snow_d(i,j).gt.snowdmin) then
               if (snow_d(i,j)+Qsubl(i,j).le.snowdmin) then
                  Qsubl(i,j) = snowdmin - snow_d(i,j)
               endif
            else
               Qsubl(i,j) = 0.0
            endif
         enddo
      enddo

! c Account for decreases in snow depth due to sublimation.  Note
! c   that the subgrid_flag = 1.0 routines have not modified the
! c   soft snow depth in any way.  Changes must be made if this is
! c   not an appropriate assumption.
      do i=1,nx
         do j=1,l_ny
            snow_d(i,j) = snow_d(i,j) + Qsubl(i,j)
            soft_snow_d(i,j) = soft_snow_d(i,j) + Qsubl(i,j)
         enddo
      enddo

! c Find any undefined ro_snow_grid(i,j) values, that correspond to
! c   non-zero snow_d(i,j) values, and fill them in with ro_snow.
! c   This no longer seems to be required.  Any zero snow_d(i,j)
! c   values should have ro_snow in the ro_snow_grid(i,j) positions
! c   from the initial conditions that were provided.
! c       do i=1,nx
! c         do j=1,ny
! c           the checking code would go here.
! c         enddo
! c       enddo

! c Update the surface topography resulting from the snow setting
      ! c   on the land.
      if (topoflag.eq.1.0) then
         do i=1,nx
            do j=1,l_ny
               new_j = global_map(j)
! c This gives the constant density version.
! c             topo(i,j) = topo_land(i,j) + snow_d(i,j)
! c This gives the spatially distributed density version.
               l_topo_tmp(i,j) = topo_land(i,new_j) + snow_d(i,j) * &
                    &          ro_snow_grid(i,j) / ro_snow
            enddo
         enddo
         call single_gather(topo, l_topo_tmp)
      elseif (topoflag.eq.0.0) then
         do i=1,nx
            do j=1,l_ny
               new_j = global_map(j)
               l_topo_tmp(i,j) = topo_land(i,new_j)
            enddo
         enddo
         call single_gather(topo, l_topo_tmp)
      endif

   else
! c Because a Tabler surface was not generated or used (because there
! c   was no blowing snow), save this time step as undef values.
      if (Tabler_2_flag.eq.1.0) &
           &    write(52,rec=iter) ((undef,i=1,nx),j=1,ny)

! c This ends the "if blowing snow flag = 1" case.
   endif


! c MOISTURE BALANCE
   ! c Save enough information to do a moisture balance.

   do i=1,nx
      do j=1,l_ny

! c Save the sublimation in terms of swe depth.
         wbal_qsubl(i,j) = Qsubl(i,j) * ro_snow / ro_water

! c Save the saltation in terms of swe depth.
         wbal_salt(i,j) = dh_salt(i,j) * ro_snow / ro_water

!c Save the suspension in terms of swe depth.
         wbal_susp(i,j) = dh_susp(i,j) * ro_snow / ro_water

!c Save the subgrid redistribution in terms of swe depth.
         wbal_subgrid(i,j) = dh_subgrid(i,j) * ro_snow / ro_water

!c Fill summing arrays of the sublimation and transport quantities.
         sum_qsubl(i,j) = sum_qsubl(i,j) + wbal_qsubl(i,j)
         sum_trans(i,j) = sum_trans(i,j) + wbal_salt(i,j) + &
              &      wbal_susp(i,j) + wbal_subgrid(i,j)

      enddo
   enddo

   do i=1,nx
      do j=1,l_ny
! c Convert any snow-depth adjustments that occurred in SnowTran-3D
! c   to swe (using the SnowTran-3D constant snow density) so
! c   that it can be used in SNOWPACK that accounts for the
! c   time-evolution of snow density.
         swe_depth(i,j) = snow_d(i,j) * ro_snow / ro_water

! c Calculate the snow depth using the spatially-distributed snow density
! c   from the snowpack model.
         snow_depth(i,j) = swe_depth(i,j) * &
              &      ro_water / ro_snow_grid(i,j)

      enddo
   enddo

   return
 end subroutine accum

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine smoother9(nx,ny,snow)

   use snowmodel_inc
   implicit none

!      include 'snowmodel.inc'

   integer i,j,nx,ny
   real snow(nx,ny)
   real snow_tmp(nx,ny)

! c Performs a 9-point smoothing operation.

! c The result at each grid point is a weighted average of the grid
! c   point and the surrounding 8 points.  The center point receives
! c   a weight of 1.0, the points at each side and above and below
! c   receive a weight of 0.5, and corner points receive a weight of
! c   0.3.  All points are multiplied by their weights and summed,
! c   then divided by the total weight.

! c Do the interior.
   do i=2,nx-1
      do j=2,ny-1
         snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + &
              &      snow(i,j+1) + snow(i-1,j) + snow(i+1,j)) + 0.3 * &
              &      (snow(i-1,j-1) + snow(i+1,j+1) + snow(i-1,j+1) + &
              &      snow(i+1,j-1))) / 4.2
      enddo
   enddo

!c Do the sides.
   j = 1
   do i=2,nx-1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j) + &
           &    snow(i+1,j)) + 0.3 * (snow(i+1,j+1) + snow(i-1,j+1))) / 3.1
   enddo

   j = ny
   do i=2,nx-1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j) + &
           &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i-1,j-1))) / 3.1
   enddo

   i = 1
   do j=2,ny-1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) + &
           &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i+1,j+1))) / 3.1
   enddo

   i = nx
   do j=2,ny-1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) + &
           &    snow(i-1,j)) + 0.3 * (snow(i-1,j-1) + snow(i-1,j+1))) / 3.1
   enddo

!c Do the corners.
   i = 1
   j = 1
   snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i+1,j)) + &
        &  0.3 * snow(i+1,j+1)) / 2.3

   i = nx
   j = 1
   snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j)) + &
        &  0.3 * snow(i-1,j+1)) / 2.3

   i = 1
   j = ny
   snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i+1,j)) + &
        &  0.3 * snow(i+1,j-1)) / 2.3

   i = nx
   j = ny
   snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j)) + &
        &  0.3 * snow(i-1,j-1)) / 2.3

!c Return the smoothed array.
   do i=1,nx
      do j=1,ny
         snow(i,j) = snow_tmp(i,j)
      enddo
   enddo

   return
 end subroutine smoother9

 ! Called from Micromet
 subroutine smoother9_parallel(nx,ny,snow)
   use caf_module
   use snowmodel_inc
   implicit none

!      include 'snowmodel.inc'

   integer i,j,nx,ny
   real snow(nx,max_l_ny)
   real snow_tmp(nx,max_l_ny)
   real reciprocal

   snow_tmp = 0.0

! c Performs a 9-point smoothing operation.

! c The result at each grid point is a weighted average of the grid
! c   point and the surrounding 8 points.  The center point receives
! c   a weight of 1.0, the points at each side and above and below
! c   receive a weight of 0.5, and corner points receive a weight of
! c   0.3.  All points are multiplied by their weights and summed,
! c   then divided by the total weight.

   !Data is now in left_rcv and right_rcv
   call halo_exchange(snow)

   if(me > 0) then
      j = 1
      do i=2,nx-1
         !left
         snow_tmp(i,j) = (snow(i,j) + 0.5 * (left_rcv(i) + &
              &      snow(i,j+1) + snow(i-1,j) + snow(i+1,j)) + 0.3 * &
              &      (left_rcv(i-1) + snow(i+1,j+1) + snow(i-1,j+1) + &
              &      left_rcv(i+1))) /4.2
      end do
   end if

   if(me < np-1) then
      j = l_ny
      do i=2,nx-1
         !right
         snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + &
              &  right_rcv(i) + snow(i-1,j) + snow(i+1,j)) + 0.3 * &
              &  (snow(i-1,j-1) + right_rcv(i+1) + right_rcv(i-1) + &
              &   snow(i+1,j-1))) /4.2
      enddo
   end if

! c Do the interior.
   do i=2,nx-1
      do j=2,l_ny-1
         snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + &
              &      snow(i,j+1) + snow(i-1,j) + snow(i+1,j)) + 0.3 * &
              &      (snow(i-1,j-1) + snow(i+1,j+1) + snow(i-1,j+1) + &
              &      snow(i+1,j-1))) / 4.2
      enddo
   enddo

   !c Do the sides.
   if(me == 0) then
      j = 1
      do i=2,nx-1
         snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j) + &
              &    snow(i+1,j)) + 0.3 * (snow(i+1,j+1) + snow(i-1,j+1))) /3.1
      enddo
   end if

   if(me == np-1) then
      j = l_ny
      do i=2,nx-1
         snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j) + &
              &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i-1,j-1))) / 3.1
      enddo
   end if

   i = 1
   !interior
   do j=2,l_ny-1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) + &
           &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i+1,j+1))) / 3.1
   enddo

   i = nx
   !interior
   do j=2,l_ny-1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) + &
           &    snow(i-1,j)) + 0.3 * (snow(i-1,j-1) + snow(i-1,j+1))) / 3.1
   enddo

   if(me > 0) then

      !left
      j = 1
      i = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (left_rcv(i) + snow(i,j+1) + &
           &    snow(i+1,j)) + 0.3 * (left_rcv(i+1) + snow(i+1,j+1))) /3.1

      !left
      j = 1
      i = nx
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (left_rcv(i) + snow(i,j+1) + &
           &    snow(i-1,j)) + 0.3 * (left_rcv(i-1) + snow(i-1,j+1))) /3.1
   end if

   if(me < np-1) then
      !right
      j = l_ny
      i = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + right_rcv(i) + &
           &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + right_rcv(i+1))) /3.1

      !right
      j = l_ny
      i = nx
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + right_rcv(i) + &
           &    snow(i-1,j)) + 0.3 * (snow(i-1,j-1) + right_rcv(i-1))) /3.1
   end if

   !c Do the corners.
   if(me == 0) then
      i = 1
      j = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i+1,j)) + &
           &  0.3 * snow(i+1,j+1)) / 2.3

      i = nx
      j = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j)) + &
           &  0.3 * snow(i-1,j+1)) / 2.3
   end if

   if(me == np-1) then
      i = 1
      j = l_ny
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i+1,j)) + &
           &  0.3 * snow(i+1,j-1)) / 2.3

      i = nx
      j = l_ny
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j)) + &
           &  0.3 * snow(i-1,j-1)) / 2.3
   end if
!c Return the smoothed array.
   do i=1,nx
      do j=1,l_ny
         snow(i,j) = snow_tmp(i,j)
      enddo
   enddo

   return
 end subroutine smoother9_parallel

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine suspension(Utau,vonKarman,nx,ny,conc_salt,&
      &  Qsalt,Qsusp,z_0,h_star,dz,ztop,pi,&
      &  fall_vel,Ur_const,Up_const,Utau_t,Qsubl,ht_rhobs,&
      &  tair_grid,rh_grid,Qsusp_u,Qsusp_v,uwind_grid,&
      &  vwind_grid)

   use caf_module, only: l_ny, max_l_ny, global_map, single_gather
   use snowmodel_vars, only: l_Qsusp
   use snowmodel_inc
   implicit none

   !      include 'snowmodel.inc'

   integer i,j,nx,ny,nzsteps,iz

   real vonKarman,dz,ztop,fall_vel,Ur_const,Up_const
   real ht_rhobs,V_susp,V_salt,pi
   real U_p,Utau_fallvel,U_r,phistar_Cr,product,conc,z

   real Utau(nx,max_l_ny)
   real Utau_t(nx,max_l_ny)
   real uwind_grid(nx,max_l_ny)
   real vwind_grid(nx,max_l_ny)
   real z_0(nx,max_l_ny)
   real h_star(nx,max_l_ny)
   real conc_salt(nx,max_l_ny)
   real Qsalt(nx,max_l_ny)
   real Qsusp(nx,max_l_ny)
   real Qsusp_u(nx,max_l_ny)
   real Qsusp_v(nx,max_l_ny)
   real Qsubl(nx,max_l_ny)

   real tair_grid(nx,max_l_ny)
   real rh_grid(nx,max_l_ny)

!   real tmpQsusp(nx,ny)

   integer :: new_j

!c Compute the mass concentration of suspended snow according to
!c   Kind (1992).
   do j=1,l_ny
      do i=1,nx
         if (Qsalt(i,j).gt.0.0) then
            Utau_fallvel = Utau(i,j) / fall_vel
            if (h_star(i,j).eq.z_0(i,j)) h_star(i,j) = 2.0 * z_0(i,j)
            U_r = Utau(i,j)/vonKarman * log(h_star(i,j)/z_0(i,j))
            phistar_Cr = Utau(i,j)/U_r * Ur_const
            product = phistar_Cr * Utau_fallvel
            U_p = Up_const * Utau_t(i,j)

!c Compute the concentration in the saltation layer (kg/m**3).
            conc_salt(i,j) = Qsalt(i,j) / (h_star(i,j) * U_p)

            nzsteps = int((ztop - h_star(i,j)) / dz)

            Qsusp(i,j) = 0.0
            Qsubl(i,j) = 0.0

            do iz=1,nzsteps
               z = h_star(i,j) + 0.5 * dz + real(iz - 1) * dz

!c Compute the concentration of the suspended snow at height z.
               conc = conc_salt(i,j) * ((product + 1.0) * &
                    &        (z/h_star(i,j))**((-fall_vel)/(vonKarman*Utau(i,j))) - &
                    &        product)
               conc = max(conc,0.0)

!c Only do The integration if the concentration is non-zero.
               if (conc.gt.0.0) then

!c Compute the sublimation due to suspension.
                  call getsublim(z,rh_grid(i,j),tair_grid(i,j),Utau(i,j),&
                       &          z_0(i,j),V_susp,V_salt,Utau_t(i,j),ht_rhobs,1.0,pi)

!c Perform the quadrature (summation), without the constants.
                  if (z.eq.z_0(i,j)) z = 1.2 * z_0(i,j)
                  Qsusp(i,j) = Qsusp(i,j) + conc * log(z/z_0(i,j)) * dz
                  Qsubl(i,j) = Qsubl(i,j) + conc * V_susp * dz

               endif

            enddo

!c Finish the quadratures.
!c Include the constants for Qsusp.
            Qsusp(i,j) = Utau(i,j) / vonKarman * Qsusp(i,j)

!c Include the sublimation contribution due to saltation.
            z = h_star(i,j) / 2.0
            call getsublim(z,rh_grid(i,j),tair_grid(i,j),Utau(i,j),&
                 &    z_0(i,j),V_susp,V_salt,Utau_t(i,j),ht_rhobs,0.0,pi)

            Qsubl(i,j) = Qsubl(i,j) + &
                 &    V_salt * conc_salt(i,j) * h_star(i,j)

         else
            conc_salt(i,j) = 0.0
            Qsusp(i,j) = 0.0
            Qsubl(i,j) = 0.0
         endif

      enddo
   enddo

   l_Qsusp = Qsusp
!   call single_gather(tmpQsusp, l_Qsusp)

! c Separate the east-west and the north-south suspended transport
! c   components; the vector sum should equal Qsusp.
   do i=1,nx
      do j=1,l_ny
         Qsusp_u(i,j) = l_Qsusp(i,j) * abs(uwind_grid(i,j)) / &
              &      sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
         Qsusp_v(i,j) = l_Qsusp(i,j) * abs(vwind_grid(i,j)) / &
              &      sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
      enddo
   enddo

   return
 end subroutine suspension

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine saltation(Qsalt,deltax,fetch,Utau,Utau_t,nx,ny,&
      &  ro_air,gravity,vegsnowd_xy,snow_d,&
      &  Qsalt_max,Qsalt_maxu,Qsalt_maxv,deltay,Qsalt_u,Qsalt_v,&
      &  uwind_grid,vwind_grid,xmu,soft_snow_d,bc_flag)

   use caf_module, only: l_ny, max_l_ny, global_map, single_gather, single_scatter, &
        & partition_ny, me, np
   use snowmodel_inc
   use snowmodel_vars, only: l_Qsalt_v, l_vwind_grid,&
        & l_Utau
   implicit none

!      include 'snowmodel.inc'

   integer i,j,nx,ny,new_j
   integer k,istart,iend,jstart,jend

   real deltax,deltay,fetch,ro_air,gravity,dUtau,xmu,&
        &  blowby,bc_flag

   real Qsalt_max(nx,max_l_ny)
   real Qsalt_maxu(nx,max_l_ny)
   real Qsalt_maxv(nx,max_l_ny)
   real Qsalt(nx,max_l_ny)
   real Qsalt_u(nx,max_l_ny)
   real Qsalt_v(nx,max_l_ny)
   real Utau(nx,ny)
   real Utau_t(nx,max_l_ny)
   real uwind_grid(nx,max_l_ny)
   real vwind_grid(nx,ny)
   real snow_d(nx,max_l_ny)
   real soft_snow_d(nx,max_l_ny)
   real vegsnowd_xy(nx,max_l_ny)

   real scale_EW,scale_NS
   real sign1,sign2
   integer npairs, iter, salt_iter
   integer, allocatable :: index_start(:), index_end(:)
   real, allocatable :: tmp_buff(:)


! c The blowby parameter is implemented to account for the erosion
! c   of the tops of deep snow accumulations.  It corrects a
! c   deficiency in the du*/dx* < 0 formulation.  It is a number that
! c   should range from 0 to 1.0, and represents the fraction of the
! c   upwind saltation flux that is transfered farther downwind into
! c   the next grid cell.  So, the bigger the number, the less
! c   peaked the drift accumulation profile is.  blowby = 0.0 is the
! c   original model.  I am now using the Tabler surfaces to do the
! c   same kind of thing, so here I hard-code the parameter as in the
! c   original model.
       blowby = 0.0

! c Compute the maximum possible saltation flux, assuming that
! c   an abundance of snow is available at the surface.
   call single_scatter(Utau,l_Utau)
   call single_scatter(vwind_grid,l_vwind_grid)
   sync all

   do i=1,nx
      do j=1,l_ny
!c For a given wind speed, find Qsalt_max.
         Qsalt_max(i,j) = 0.68 * ro_air / gravity * &
              &      Utau_t(i,j) / l_Utau(i,j) * (l_Utau(i,j)**2 - Utau_t(i,j)**2)
         Qsalt_max(i,j) = max(Qsalt_max(i,j),0.0)

! c Now weight the max saltation flux for the u and v wind
! c   components, where the vector sum should equal Qsalt_max.
         Qsalt_maxu(i,j) = Qsalt_max(i,j) * abs(uwind_grid(i,j)) / &
              &      sqrt(uwind_grid(i,j)**2 + l_vwind_grid(i,j)**2)
         Qsalt_maxv(i,j) = Qsalt_max(i,j) * abs(l_vwind_grid(i,j)) / &
              &      sqrt(uwind_grid(i,j)**2 + l_vwind_grid(i,j)**2)
      enddo
   enddo

! c Define an upwind boundary condition.  If bc_flag = 1.0 then it is
! c   assumed that the inflow saltation flux has reached steady state.
! c   If bc_flag = 0.0 then the saltation flux is assumed to be zero.
! c   The boundary condition is implemented by initializing the arrays
! c   to Qsalt_max, and since upwind boundaries are not called in
! c   the Qsalt computation, they stay in effect for the future
! c   accumulation/erosion computation.
   if (bc_flag.eq.0.0) then
      do i=1,nx
         do j=1,l_ny
!c Zero incoming flux at the boundaries.
            Qsalt_u(i,j) = 0.0
            Qsalt_v(i,j) = 0.0
         enddo
      enddo
   elseif (bc_flag.eq.1.0) then
      do i=1,nx
         do j=1,l_ny
!c Steady-state (maximum) incoming flux at the boundaries.
            Qsalt_u(i,j) = Qsalt_maxu(i,j)
            Qsalt_v(i,j) = Qsalt_maxv(i,j)
         enddo
      enddo
   endif

!   call single_gather(Qsalt_maxv,l_Qsalt_maxv)

! c Define the scaling coefficients for Eqn. 9 in L&S 1998. Don't
! c   let them be greater than 1.0 or you will make more snow than
! c   there was before.
   scale_EW =  xmu * deltax / fetch
   scale_EW = min(1.0,scale_EW)
   scale_NS =  xmu * deltay / fetch
   scale_NS = min(1.0,scale_NS)

   allocate(index_start(nx),index_end(nx))
   index_start(:) = 0
   index_end(:) = 0

   ! !c Consider WESTERLY winds.
   do j=1,l_ny
      call get_direction_west(uwind_grid, index_start, index_end, npairs, nx,l_ny,j)
      do k=1,npairs
         istart = index_start(k)+1
         iend = index_end(k)
         do i=istart,iend
            dUtau = l_Utau(i,j) - l_Utau(i-1,j)
            if (dUtau.ge.0.0) then
               Qsalt_u(i,j) = Qsalt_u(i-1,j) + scale_EW * &
                    &          (Qsalt_maxu(i,j) - Qsalt_u(i-1,j))
            else
!c             Qsalt_u(i,j) = min(Qsalt_u(i-1,j),Qsalt_maxu(i,j))

               if (Qsalt_u(i-1,j).lt.Qsalt_maxu(i,j)) then
                  Qsalt_u(i,j) = Qsalt_u(i-1,j)
               else
                  Qsalt_u(i,j) = &
                       &            max(blowby*Qsalt_u(i-1,j),Qsalt_maxu(i,j))
               endif

            endif
         enddo !i
      enddo !k
   enddo

   index_start = 0
   index_end = 0

! !c Consider EASTERLY winds.
   do j=1,l_ny
      call get_direction_east(uwind_grid, index_start, index_end, npairs, nx, l_ny, j)
      do k=1,npairs
         iend = index_start(k)
         istart = index_end(k)-1
         do i=istart,iend,-1
            dUtau = l_Utau(i,j) - l_Utau(i+1,j)
            if (dUtau.ge.0.0) then
               Qsalt_u(i,j) = Qsalt_u(i+1,j) + scale_EW * &
                    &          (Qsalt_maxu(i,j) - Qsalt_u(i+1,j))
            else
!c             Qsalt_u(i,j) = min(Qsalt_u(i+1,j),Qsalt_maxu(i,j))

               if (Qsalt_u(i+1,j).lt.Qsalt_maxu(i,j)) then
                  Qsalt_u(i,j) = Qsalt_u(i+1,j)
               else
                  Qsalt_u(i,j) = &
                       &            max(blowby*Qsalt_u(i+1,j),Qsalt_maxu(i,j))
               endif

            endif
         enddo
      enddo
   enddo


   deallocate(index_start,index_end)
   allocate(index_start(ny),index_end(ny))

   index_start = 0
   index_end = 0

   allocate(tmp_buff(nx))


!   do iter = 1,np
   salt_iter = ceiling(10000 / (deltay * max_l_ny))
   salt_iter = max(2,salt_iter)
   do iter = 1, min(np,salt_iter)
      tmp_buff = 0.0
      if (iter.gt.1) then
        sync all
        if(me /= 0) tmp_buff(:) = l_Qsalt_v(:,partition_ny(me-1)) [me] !Remember each processor is me-1
      endif
!c Consider SOUTHERLY winds.
      do i=1,nx
         call get_direction_south(vwind_grid, index_start, index_end, npairs, nx, ny, i)
         do k=1,npairs
            jstart = index_start(k)+1
            jend = index_end(k)
            if(jstart == 1 .and. jend >= jstart) then
               ! Use the buffer
               new_j = global_map(jstart)
               dUtau = Utau(i,new_j) - Utau(i,global_map(jstart-1))
               if (dUtau.ge.0.0) then
                  Qsalt_v(i,jstart) = tmp_buff(i) + scale_NS * &
                       &          (Qsalt_maxv(i,jstart) - tmp_buff(i))
               else
                  if (tmp_buff(i) .lt. Qsalt_maxv(i,jstart)) then
                     Qsalt_v(i,jstart) = tmp_buff(i)
                  else
                     Qsalt_v(i,jstart) = &
                          &            max(blowby*tmp_buff(i),Qsalt_maxv(i,jstart))
                  endif
               endif
               jstart = jstart + 1
            end if

            do j=jstart,jend
               new_j = global_map(j)
               dUtau = Utau(i,new_j) - Utau(i,new_j-1)
               if (dUtau.ge.0.0) then
                  Qsalt_v(i,j) = Qsalt_v(i,j-1) + scale_NS * &
                       &          (Qsalt_maxv(i,j) - Qsalt_v(i,j-1))
               else
                  !c             Qsalt_v(i,j) = min(Qsalt_v(i,j-1),Qsalt_maxv(i,j))

                  if (Qsalt_v(i,j-1) .lt. Qsalt_maxv(i,j)) then
                     Qsalt_v(i,j) = Qsalt_v(i,j-1)
                  else
                     Qsalt_v(i,j) = &
                          &            max(blowby*Qsalt_v(i,j-1),Qsalt_maxv(i,j))
                  endif
               endif
               !if(i == 4 .and. j==22) write(*,*) me, Qsalt_v(i,j), Qsalt_maxv(i,new_j),dUtau,Qsalt_v(i,j-1)
            enddo
         enddo
      enddo
   end do

   index_start = 0
   index_end = 0
   tmp_buff = 0.0

!   sync all

!   do iter = 1,np
   do iter = 1,min(np,salt_iter)
      if (iter.gt.1) then
        sync all
        if(me /= np-1) tmp_buff(:) = l_Qsalt_v(:,1) [me+2]!Remember each processor is me-1
      endif
      !c Consider NORTHERLY winds.
      do i=1,nx
         call get_direction_north(vwind_grid, index_start, index_end, npairs, nx, ny, i)
         do k=1, npairs
            jend = index_start(k)
            jstart = index_end(k)-1
            if(jstart == l_ny .and. jend <= jstart) then
               new_j = global_map(jstart)
               dUtau = Utau(i,new_j) - Utau(i,global_map(jstart+1))
               if (dUtau.ge.0.0) then
                  Qsalt_v(i,jstart) = tmp_buff(i) + scale_NS * &
                       &          (Qsalt_maxv(i,jstart) - tmp_buff(i))
               else
                  if (tmp_buff(i).lt.Qsalt_maxv(i,jstart)) then
                     Qsalt_v(i,jstart) = tmp_buff(i)
                  else
                     Qsalt_v(i,jstart) = &
                          &            max(blowby*tmp_buff(i),Qsalt_maxv(i,jstart))
                  endif
               endif
               jstart = jstart - 1
            end if

            do j=jstart,jend,-1
               new_j = global_map(j)
               dUtau = Utau(i,new_j) - Utau(i,global_map(j+1))
               if (dUtau.ge.0.0) then
                  Qsalt_v(i,j) = Qsalt_v(i,j+1) + scale_NS * &
                       &          (Qsalt_maxv(i,j) - Qsalt_v(i,j+1))
               else
                  !c             Qsalt_v(i,j) = min(Qsalt_v(i,j+1),Qsalt_maxv(i,j))

                  if (Qsalt_v(i,j+1).lt.Qsalt_maxv(i,j)) then
                     Qsalt_v(i,j) = Qsalt_v(i,j+1)
                  else
                     Qsalt_v(i,j) = &
                          &            max(blowby*Qsalt_v(i,j+1),Qsalt_maxv(i,j))
                  endif

               endif
            enddo
         enddo
      enddo
   enddo

!   sync all


   deallocate(tmp_buff)

!c Combine the u and v components to yield the total saltation flux
   !c   at each grid cell.
   do j=1,l_ny
      do i=1,nx
         Qsalt(i,j) = Qsalt_u(i,j) + Qsalt_v(i,j)
      enddo
   enddo

! c Adjust Qsalt to account for the availablity of snow for transport;
! c   taking into consideration whether there is snow on the ground,
   ! c   the holding depth of the vegetation, etc..
   do j=1,l_ny
      do i=1,nx
         if (snow_d(i,j).le.vegsnowd_xy(i,j)) Qsalt(i,j) = 0.0
         if (soft_snow_d(i,j).le.0.0) Qsalt(i,j) = 0.0
      enddo
   enddo

   return
 end subroutine saltation

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine solveUtau(Utau,ht_windobs,windspd_grid,C_z,vonKarman, &
      &  gravity,z_0,h_star,h_const,vegsnowd_xy,snow_d, &
      &  snow_z0,veg_z0,bs_flag,nx,ny,Utau_t,soft_snow_d)

   use caf_module, only: max_l_ny, l_ny, global_map
   use snowmodel_inc
   implicit none

!      include 'snowmodel.inc'

   integer i,j,nx,ny,new_j

   real bs_flag,guess,sfrac,vonKarman,ht_windobs,C_z,gravity
   real h_const,snow_z0,Utautmp,windtmp,wind_max
   real threshold,threshold_flag,z_0_tmp

   real Utau(nx,max_l_ny)
   real Utau_t(nx,max_l_ny)
   real windspd_grid(nx,max_l_ny)
   real z_0(nx,max_l_ny)
   real h_star(nx,max_l_ny)

   real snow_d(nx,max_l_ny)
   real soft_snow_d(nx,max_l_ny)
   real veg_z0(nx,max_l_ny)
   real vegsnowd_xy(nx,max_l_ny)

! c Initially set the blowing snow flag to no blowing snow
! c   (bs_flag = 0.0).  Then, if snow is found to blow in any
! c   domain grid cell, set the flag to on (bs_flag = 1.0).
   bs_flag = 0.0
!   sync all
!c Build the Utau array.
   guess = 0.1
   do j=1,l_ny
      do i=1,nx

! c Determine whether snow is saltating (this influences how Utau
! c   and z_0 are computed).
         if (snow_d(i,j).le.vegsnowd_xy(i,j)) then

!c Saltation will not occur.
            sfrac = snow_d(i,j) / max(vegsnowd_xy(i,j),veg_z0(i,j))
            z_0(i,j) = sfrac * snow_z0 + (1.0 - sfrac) * veg_z0(i,j)
            z_0_tmp = min(0.25*ht_windobs,z_0(i,j))
            Utau(i,j) = windspd_grid(i,j) * &
                 &      vonKarman / log(ht_windobs/z_0_tmp)
            h_star(i,j) = z_0(i,j) * h_const / C_z
         elseif (soft_snow_d(i,j).le.0.0) then
!c Saltation will not occur.
            z_0(i,j) = snow_z0
            Utau(i,j) = windspd_grid(i,j) * &
                 &      vonKarman / log(ht_windobs/z_0(i,j))
            h_star(i,j) = z_0(i,j)
         else
! c Saltation may occur.  Test for that possibility by assuming that
! c   saltation is present, solving for Utau and z_0, and comparing
! c   whether Utau exceeds Utau_t.  If it does not, set z_0 to that
! c   of snow and recompute Utau.

! c To help insure that the iteration converges, set the minimum
! c   wind speed to be 1.0 m/s, and the maximum wind speed to be
! c   30 m/s at 10-m height.
            windtmp = max(1.0,windspd_grid(i,j))
            wind_max = 30.0 * log(ht_windobs/snow_z0)/log(10.0/snow_z0)
            windtmp = min(windtmp,wind_max)

! c For u* over 0.6, use the relation z0 = 0.00734 u* - 0.0022,
! c   instead of Equation (5) in Liston and Sturm (1998).  Note that
! c   for windspeeds greater than about 35 m/s this will have to be
! c   modified for the solution algorithm to converge (because the
! c   roughness length will start to be higher than the obs height!).
            threshold = 0.6/vonKarman * log(ht_windobs/0.0022)
            if (windtmp.le.threshold) then
               threshold_flag = 1.0
            else
               threshold_flag = 2.0
            endif

            call solve1(Utautmp,guess,ht_windobs,windtmp,C_z,vonKarman,&
                 &      gravity,threshold_flag)

            if (Utautmp.gt.Utau_t(i,j)) then

!c We have saltation.
               Utau(i,j) = Utautmp
               z_0(i,j) = C_z * Utau(i,j)**2 / (2.0 * gravity)
               h_star(i,j) = h_const * Utau(i,j)**2 / (2.0 * gravity)
               bs_flag = 1.0
            else

! c We do not have saltation, but the vegetation is covered by snow.
! c   Because we have determined that we do not have saltation, make
! c   sure Utau does not exceed Utau_t.
               z_0(i,j) = snow_z0
               Utau(i,j) = windspd_grid(i,j) * &
                    &        vonKarman / log(ht_windobs/z_0(i,j))
               Utau(i,j) = min(Utau(i,j),Utau_t(i,j))
               h_star(i,j) = z_0(i,j) * h_const / C_z

            endif
         endif

      enddo
   enddo

   call co_max(bs_flag)

   return
 end subroutine solveUtau

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine solve1(xnew,guess,z,windtmp,C_z,vonKarman,&
      &  gravity,threshold_flag)

   use snowmodel_inc
   implicit none

   integer i,maxiter

   real xnew,guess,z,windtmp,C_z,vonKarman,tol,old,gravity
   real fprime,funct,threshold_flag

   tol = 1.0e-3
   maxiter = 20
   old = guess

   if (threshold_flag.eq.1.0) then

      do i=1,maxiter
         fprime = - 1.0 + 2.0 / old * windtmp * vonKarman * &
              &      (log(z) - log(C_z/(2.0*gravity)) - 2.0*log(old))**(-2)
         funct = - old + windtmp * vonKarman * &
              &      (log(z) - log(C_z/(2.0*gravity)) - 2.0*log(old))**(-1)
         xnew = old - funct/fprime
         if (abs(xnew - old).lt.tol) return
         old = xnew
      end do

   elseif (threshold_flag.eq.2.0) then

      old = 0.6
      do i=1,maxiter
         fprime = - 1.0 + windtmp * vonKarman * &
              &      0.00734 / (0.00734 * old - 0.0022) * &
              &      (log(z) - log(0.00734 * old - 0.0022))**(-2)
         funct = - old + windtmp * vonKarman * &
              &      (log(z) - log(0.00734 * old - 0.0022))**(-1)
         xnew = old - funct/fprime
         if (abs(xnew - old).lt.tol) return
         old = xnew
      end do

   endif

   print *,'max iteration exceeded when solving for Utau, Utau=',old

   return
 end subroutine solve1

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine getsublim(z,rh,tair,Utau,&
      &  z_0,V_susp,V_salt,Utau_t,ht_rhobs,flag,pi)

   implicit none

   real pi,ro_ice,xM,R,R_dryair,vonKarman,visc_air,h_s
   real xlamdaT,D,ro_sat,rh_offset,sigma
   real alfa,rbar_r,xmbar,rbar,u_z,x_r,wbar,flag
   real z,rh,tair,Utau,z_0,V_susp,V_salt,Utau_t,ht_rhobs
   real V_r,xN_r,xNu,xSh,tmp1,tmp2,top,bottom,V_rsalt

   ro_ice = 917.0
   xM = 18.01
   R = 8313.
   R_dryair = 287.
   vonKarman = 0.4
   visc_air = 13.e-6
   h_s = 2.838e6

!c     xlamdaT = 0.00063 * tair + 0.0673
   xlamdaT = 0.024
   D = 2.06e-5 * (tair/273.0)**(1.75)
!c     ro_sat = 0.622 * 10.0**(11.40 - 2353./tair) / (R_dryair * tair)
   ro_sat = 0.622 / (R_dryair * tair) * &
        &  610.78 * exp(21.875 * (tair - 273.15) / (tair - 7.66))

! c Assume that the rh varies according to a modification to
! c   Pomeroy's humidity variation with height equation.
   rh_offset = 1.0 - 0.027 * log(ht_rhobs)
   sigma = (0.01 * rh - 1.0) * (rh_offset + 0.027 * log(z))
   sigma = min(0.0,sigma)
   sigma = max(-1.0,sigma)

   alfa = 4.08 + 12.6 * z
   rbar_r = 4.6e-5 * z**(-0.258)
   xmbar = 4.0/3.0 * pi * ro_ice * rbar_r**3 * &
        &  (1.0 + 3.0/alfa + 2.0/alfa**2)
   rbar = ((3.0 * xmbar) / (4.0 * pi * ro_ice))**(0.33)
   u_z = Utau/vonKarman * log(z/z_0)
   x_r = 0.005 * u_z**(1.36)
   wbar = 1.1e7 * rbar**(1.8)

   if (flag.eq.1.0) then

! c Compute the sublimation loss rate coefficient for the suspension
! c   layer.
      V_r = wbar + 3.0 * x_r * cos(pi/4.0)
      xN_r = 2.0 * rbar * V_r / visc_air
      xNu = 1.79 + 0.606 * xN_r**(0.5)
      xSh = xNu
      tmp1 = (h_s * xM)/(R * tair) - 1.0
      tmp2 = xlamdaT * tair * xNu
      top = 2.0 * pi * rbar * sigma
      bottom = h_s/tmp2 * tmp1 + 1.0/(D * ro_sat * xSh)
      V_susp = (top/bottom)/xmbar
      V_salt = 0.0

   elseif (flag.eq.0.0) then

! c Compute the sublimation loss rate coefficient for the saltation
! c   layer.
      V_rsalt = 0.68 * Utau + 2.3 * Utau_t
      xN_r = 2.0 * rbar * V_rsalt / visc_air
      xNu = 1.79 + 0.606 * xN_r**(0.5)
      xSh = xNu
      tmp1 = (h_s * xM)/(R * tair) - 1.0
      tmp2 = xlamdaT * tair * xNu
      top = 2.0 * pi * rbar * sigma
      bottom = h_s/tmp2 * tmp1 + 1.0/(D * ro_sat * xSh)
      V_salt = (top/bottom)/xmbar
      V_susp = 0.0

   endif

   return
 end subroutine getsublim

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine get_direction_west(uwind_grid, index_start, index_end, npairs, nx,ny,new_j)
   use caf_module, only: l_ny, max_l_ny, global_map
   use snowmodel_inc

   real, intent(in) ::  uwind_grid(nx,ny)
   integer, intent(inout) :: index_start(nx), index_end(nx)
   integer, intent(out) :: npairs
   integer, intent(in) :: nx,ny,new_j
   real sign1,sign2
   integer :: i

   npairs = 0
! Sweep looking for WESTERLY winds, looking for positive numbers.

   if (uwind_grid(1,new_j).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.gt.0.0) then
      npairs = 1
      index_start(npairs) = 1
      ! index_uw(j,2) =  1
   else
      npairs = 0
   endif
   do i=2,nx

      if (uwind_grid(i-1,new_j).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif
      if (uwind_grid(i,new_j).le.0.0) then
         sign2 = -1.0
      else
         sign2 = 1.0
      endif

      if (sign2.ne.sign1) then
         !c We have a sign change.
         if (sign2.gt.0.0) then
            ! c We have gone from negative to positive, indicating the start
            ! c   of a new positive group.
            npairs = npairs + 1
            index_start(npairs) = i
         else
            !c We have gone from positive to negative, indicating the end of
            !c   the group.
            index_end(npairs) = i - 1
         endif
      endif
   enddo

   if (uwind_grid(nx,new_j).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.gt.0.0) then
      index_end(npairs) = nx
   endif

 end subroutine get_direction_west

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine get_direction_east(uwind_grid, index_start, index_end, npairs, nx,ny,new_j)
   use caf_module, only: l_ny, max_l_ny, global_map
   use snowmodel_inc

   real, intent(in) ::  uwind_grid(nx,ny)
   integer, intent(inout) :: index_start(nx), index_end(nx)
   integer, intent(out) :: npairs
   integer, intent(in) :: nx,ny,new_j
   real sign1,sign2
   integer :: i

   if (uwind_grid(1,new_j).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.lt.0.0) then
      npairs = 1
      index_start(npairs) = 1
   else
      npairs = 0
   endif
   do i=2,nx

      if (uwind_grid(i-1,new_j).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif
      if (uwind_grid(i,new_j).le.0.0) then
         sign2 = -1.0
      else
         sign2 = 1.0
      endif

      if (sign2.ne.sign1) then
         !c We have a sign change.
         if (sign2.lt.0.0) then
            ! c We have gone from positive to negative, indicating the start
            ! c   of a new negative group.
            npairs = npairs + 1
            index_start(npairs) = i
         else
            ! c We have gone from negative to positive, indicating the end of
            ! c   the group.
            index_end(npairs) = i - 1
         endif
      endif
   enddo

   if (uwind_grid(nx,new_j).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.lt.0.0) then
      index_end(npairs) = nx
   endif

 end subroutine get_direction_east

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine get_direction_south(vwind_grid, index_start, index_end, npairs, nx, ny, i)
   use caf_module, only: l_ny, max_l_ny, global_map, me
   use snowmodel_inc

   real, intent(in) ::  vwind_grid(nx,ny)
   integer, intent(inout) :: index_start(ny), index_end(ny)
   integer, intent(out) :: npairs
   integer, intent(in) :: nx,ny,i
   real sign1,sign2
   integer :: j, new_j, start_j

   npairs = 0

   if(me == 0) then
      start_j = 1
   else
      start_j = 0
   end if

   if (vwind_grid(i,global_map(start_j)).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.gt.0.0) then
      npairs = 1
      if(me == 0) then
         index_start(npairs) = 1 !start_j
      else
         index_start(npairs) = 0
      end if
   else
      npairs = 0
   endif

   if(me == 0) then
      start_j = 2
   else
      start_j = 1
   end if

   do j=start_j,l_ny

      new_j = global_map(j)

      if (vwind_grid(i,new_j-1).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif
      if (vwind_grid(i,new_j).le.0.0) then
         sign2 = -1.0
      else
         sign2 = 1.0
      endif

      if (sign2.ne.sign1) then
         !c We have a sign change.
         if (sign2.gt.0.0) then
            ! c We have gone from negative to positive, indicating the start
            ! c   of a new positive group.
            npairs = npairs + 1
            index_start(npairs) = j
         else
            ! c We have gone from positive to negative, indicating the end of
            ! c   the group.
            index_end(npairs) = j - 1
         endif
      endif
   enddo

   if (vwind_grid(i,global_map(l_ny)).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.gt.0.0) then
      index_end(npairs) = l_ny
   endif

 end subroutine get_direction_south

 ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine get_direction_north(vwind_grid, index_start, index_end, npairs, nx, ny, i)
   use caf_module, only: l_ny, max_l_ny, global_map, me, np
   use snowmodel_inc

   real, intent(in) ::  vwind_grid(nx,ny)
   integer, intent(inout) :: index_start(ny), index_end(ny)
   integer, intent(out) :: npairs
   integer, intent(in) :: nx,ny,i
   real sign1,sign2
   integer :: j, new_j, start_j

   npairs = 0

   if(me == 0) then
      start_j = 1
   else
      start_j = 0
   end if

   if (vwind_grid(i,global_map(start_j)).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.lt.0.0) then
      npairs = 1
      index_start(npairs) = 1
   else
      npairs = 0
   endif

   start_j = start_j + 1

   do j=start_j,l_ny

      new_j = global_map(j)

      if (vwind_grid(i,new_j-1).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif
      if (vwind_grid(i,new_j).le.0.0) then
         sign2 = -1.0
      else
         sign2 = 1.0
      endif

      if (sign2.ne.sign1) then
         !c We have a sign change.
         if (sign2.lt.0.0) then
            ! c We have gone from positive to negative, indicating the start
            ! c   of a new negative group.
            npairs = npairs + 1
            index_start(npairs) = j
         else
            ! c We have gone from negative to positive, indicating the end of
            ! c   the group.
            index_end(npairs) = j - 1
         endif
      endif
   enddo

   if (vwind_grid(i,global_map(l_ny)).le.0.0) then
      sign1 = -1.0
   else
      sign1 = 1.0
   endif

   if (sign1.lt.0.0) then
      index_end(npairs) = l_ny
      if(me /= np-1) then
         if (vwind_grid(i,global_map(l_ny)+1).le.0.0) then
            index_end(npairs) = l_ny+1 !interval goes to right neighbour
         end if
      end if
   endif

 end subroutine get_direction_north

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine getdirection(nx,ny,uwind_grid,vwind_grid,index_ue,&
      &  index_uw,index_vn,index_vs)

   use snowmodel_inc
   implicit none

      ! include 'snowmodel.inc'

   integer i,j,nx,ny,npairs

   real sign1,sign2

   integer index_ue(ny_max,2*nx_max+1)
   integer index_uw(ny_max,2*nx_max+1)
   integer index_vn(nx_max,2*ny_max+1)
   integer index_vs(nx_max,2*ny_max+1)

   real uwind_grid(nx,ny)
   real vwind_grid(nx,ny)

! c Index whether the winds are blowing east or west.  The first
! c   column of the index array is the number of pairs of begining
! c   and ending array index of blocks of wind running in the same
! c   direction.

! c Sweep looking for WESTERLY winds, looking for positive numbers.
   do j=1,ny

      if (uwind_grid(1,j).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.gt.0.0) then
         npairs = 1
         index_uw(j,2) =  1
      else
         npairs = 0
      endif
      do i=2,nx

         if (uwind_grid(i-1,j).le.0.0) then
            sign1 = -1.0
         else
            sign1 = 1.0
         endif
         if (uwind_grid(i,j).le.0.0) then
            sign2 = -1.0
         else
            sign2 = 1.0
         endif

         if (sign2.ne.sign1) then
!c We have a sign change.
            if (sign2.gt.0.0) then
! c We have gone from negative to positive, indicating the start
! c   of a new positive group.
               npairs = npairs + 1
               index_uw(j,npairs*2) = i
            else
!c We have gone from positive to negative, indicating the end of
!c   the group.
               index_uw(j,npairs*2+1) = i - 1
            endif
         endif
      enddo

      if (uwind_grid(nx,j).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.gt.0.0) then
         index_uw(j,npairs*2+1) = nx
      endif
      index_uw(j,1) = npairs
   enddo

! c     do j=1,ny
! c       print 30, (index_uw(j,k),k=1,index_uw(j,1)*2+1)
! c     enddo
! c     print *
! c     print *

! c Sweep looking for EASTERLY winds, looking for negative numbers.
   do j=1,ny

      if (uwind_grid(1,j).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.lt.0.0) then
         npairs = 1
         index_ue(j,2) = 1
      else
         npairs = 0
      endif
      do i=2,nx

         if (uwind_grid(i-1,j).le.0.0) then
            sign1 = -1.0
         else
            sign1 = 1.0
         endif
         if (uwind_grid(i,j).le.0.0) then
            sign2 = -1.0
         else
            sign2 = 1.0
         endif

         if (sign2.ne.sign1) then
!c We have a sign change.
            if (sign2.lt.0.0) then
! c We have gone from positive to negative, indicating the start
! c   of a new negative group.
               npairs = npairs + 1
               index_ue(j,npairs*2) = i
            else
! c We have gone from negative to positive, indicating the end of
! c   the group.
               index_ue(j,npairs*2+1) = i - 1
            endif
         endif
      enddo

      if (uwind_grid(nx,j).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.lt.0.0) then
         index_ue(j,npairs*2+1) = nx
      endif
      index_ue(j,1) = npairs
   enddo

! c     do j=1,ny
! c       print 30, (index_ue(j,k),k=1,index_ue(j,1)*2+1)
! c     enddo
! c     print *
! c     print *

! c Sweep looking for SOUTHERLY winds, looking for positive numbers.
   do i=1,nx

      if (vwind_grid(i,1).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.gt.0.0) then
         npairs = 1
         index_vs(i,2) = 1
      else
         npairs = 0
      endif
      do j=2,ny

         if (vwind_grid(i,j-1).le.0.0) then
            sign1 = -1.0
         else
            sign1 = 1.0
         endif
         if (vwind_grid(i,j).le.0.0) then
            sign2 = -1.0
         else
            sign2 = 1.0
         endif

         if (sign2.ne.sign1) then
!c We have a sign change.
            if (sign2.gt.0.0) then
! c We have gone from negative to positive, indicating the start
! c   of a new positive group.
               npairs = npairs + 1
               index_vs(i,npairs*2) = j
            else
! c We have gone from positive to negative, indicating the end of
! c   the group.
               index_vs(i,npairs*2+1) = j - 1
            endif
         endif
      enddo

      if (vwind_grid(i,ny).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.gt.0.0) then
         index_vs(i,npairs*2+1) = ny
      endif
      index_vs(i,1) = npairs
   enddo

! c     do i=1,nx
! c       print 30, (index_vs(i,k),k=1,index_vs(i,1)*2+1)
! c     enddo
! c     print *
! c     print *

! c Sweep looking for NORTHERLY winds, looking for negative numbers.
   do i=1,nx

      if (vwind_grid(i,1).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.lt.0.0) then
         npairs = 1
         index_vn(i,2) = 1
      else
         npairs = 0
      endif
      do j=2,ny

         if (vwind_grid(i,j-1).le.0.0) then
            sign1 = -1.0
         else
            sign1 = 1.0
         endif
         if (vwind_grid(i,j).le.0.0) then
            sign2 = -1.0
         else
            sign2 = 1.0
         endif

         if (sign2.ne.sign1) then
!c We have a sign change.
            if (sign2.lt.0.0) then
! c We have gone from positive to negative, indicating the start
! c   of a new negative group.
               npairs = npairs + 1
               index_vn(i,npairs*2) = j
            else
! c We have gone from negative to positive, indicating the end of
! c   the group.
               index_vn(i,npairs*2+1) = j - 1
            endif
         endif
      enddo

      if (vwind_grid(i,ny).le.0.0) then
         sign1 = -1.0
      else
         sign1 = 1.0
      endif

      if (sign1.lt.0.0) then
         index_vn(i,npairs*2+1) = ny
      endif
      index_vn(i,1) = npairs
   enddo

! c     do i=1,nx
! c       print 30, (index_vn(i,k),k=1,index_vn(i,1)*2+1)
! c     enddo
! c     print *
! c     print *
! c 30  format(20i4)

   return
 end subroutine getdirection

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine getnewdepth(nx,ny,deltax,deltay,Qsalt_u,&
      &  Qsalt_v,dh_salt_u,dh_salt_v,&
      &  ro_snow,dt,vegsnowd_xy,snow_d,&
      &  soft_snow_d,snow_d_dep,soft_snow_d_dep,subgrid_flag,&
      &  uwind_grid,vwind_grid)

! c If subgrid_flag = 1.0, this routine separates the snow
! c   accumulation and erosion contributions to the snowpack
! c   evolution so they can be dealt with separately in the
   ! c   subgrid_flag = 1.0 approach and routines.

   use caf_module, only: l_ny, max_l_ny, single_gather, single_scatter, &
        & partition_ny, me, np
   use snowmodel_inc
   use snowmodel_vars, only: l_Qsalt_v

   implicit none

   !include 'snowmodel.inc'

   integer i,j,nx,ny
   integer k,istart,iend,jstart,jend

   real deltax,deltay,ro_snow,dt,dQsalt,snowdmin
   real hard_snow_d,weight_u,weight_v,eps


   real uwind_grid(nx,max_l_ny)
   real vwind_grid(nx,ny)

   real Qsalt_u(nx,max_l_ny)
   real Qsalt_v(nx,max_l_ny)

   real snow_d(nx,max_l_ny)
   real soft_snow_d(nx,max_l_ny)
   real dh_salt_u(nx,max_l_ny)
   real dh_salt_v(nx,max_l_ny)
   real vegsnowd_xy(nx,max_l_ny)
   real snow_d_dep(nx,max_l_ny)
   real soft_snow_d_dep(nx,max_l_ny)

   real subgrid_flag

   integer npairs, iter, salt_iter
   integer, allocatable :: index_start(:), index_end(:)
   real, allocatable :: tmp_buff(:)


! c Define an upwind boundary condition for saltation (here I have
! c   assumed that the transport is in equilibrium).
   do i=1,nx
      do j=1,l_ny
         dh_salt_u(i,j) = 0.0
         dh_salt_v(i,j) = 0.0
      enddo
   enddo

   allocate(index_start(nx),index_end(nx))
   index_start(:) = 0
   index_end(:) = 0


!c Consider WESTERLY winds.
   do j=1,l_ny
!      new_j = global_map(j)
      call get_direction_west(uwind_grid, index_start, index_end, npairs, nx,l_ny,j)
      do k=1, npairs
         istart = index_start(k)+1
         iend = index_end(k)
         do i=istart,iend
            dQsalt = Qsalt_u(i,j) - Qsalt_u(i-1,j)
            dh_salt_u(i,j) = (- dt) / ro_snow * dQsalt / deltax

            ! c Make adjustments for the case where there is no snow available
            ! c   on the ground (or captured within the vegetation) to be
            ! c   eroded.
            hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
            if (snow_d(i,j).gt.snowdmin) then
               if (snow_d(i,j) + dh_salt_u(i,j).le.snowdmin) then
                  dh_salt_u(i,j) = snowdmin - snow_d(i,j)
                  Qsalt_u(i,j) = Qsalt_u(i-1,j) - dh_salt_u(i,j) * &
                       &            ro_snow * deltax / dt
               endif
            else
               Qsalt_u(i,j) = 0.0
               dh_salt_u(i,j) = 0.0
            endif
         enddo
      enddo
   enddo

   index_start = 0
   index_end = 0

!c Consider EASTERLY winds.
   do j=1,l_ny
!      new_j = global_map(j)
      call get_direction_east(uwind_grid, index_start, index_end, npairs, nx,l_ny,j)
      do k=1,npairs
         iend = index_start(k)
         istart = index_end(k)-1
         do i=istart,iend,-1
            dQsalt = Qsalt_u(i,j) - Qsalt_u(i+1,j)
            dh_salt_u(i,j) = (- dt) / ro_snow * dQsalt / deltax

            ! c Make adjustments for the case where there is no snow available
            ! c   on the ground (or captured within the vegetation) to be
            ! c   eroded.
            hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
            snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
            if (snow_d(i,j).gt.snowdmin) then
               if (snow_d(i,j) + dh_salt_u(i,j).le.snowdmin) then
                  dh_salt_u(i,j) = snowdmin - snow_d(i,j)
                  Qsalt_u(i,j) = Qsalt_u(i+1,j) - dh_salt_u(i,j) * &
                       &            ro_snow * deltax / dt
               endif
            else
               Qsalt_u(i,j) = 0.0
               dh_salt_u(i,j) = 0.0
            endif
         enddo
      enddo
   enddo

   deallocate(index_start,index_end)
   allocate(index_start(ny),index_end(ny))

   index_start = 0
   index_end = 0


!   call single_scatter(Qsalt_v,l_Qsalt_v)

   allocate(tmp_buff(nx))

!   do iter = 1, np
   salt_iter = ceiling(10000 / (deltay * max_l_ny))
   salt_iter = max(2,salt_iter)
   do iter = 1,min(np,salt_iter)
      sync all
      tmp_buff = 0.0
      l_Qsalt_v = Qsalt_v
      if(me /= 0) tmp_buff(:) = l_Qsalt_v(:,partition_ny(me-1)) [me]!Remember each processor is me-1
      !c Consider SOUTHERLY winds.
      do i=1,nx
         call get_direction_south(vwind_grid, index_start, index_end, npairs, nx, ny, i)
         do k=1,npairs
            jstart = index_start(k)+1
            jend = index_end(k)
            if(jstart == 1 .and. jend >= jstart) then
               !use buffer
!               new_j = global_map(jstart)
               dQsalt = Qsalt_v(i,jstart) - tmp_buff(i)!Qsalt_v(i,j-1)
               dh_salt_v(i,jstart) = (- dt) / ro_snow * dQsalt / deltay

               ! c Make adjustments for the case where there is no snow available
               ! c   on the ground (or captured within the vegetation) to be
               ! c   eroded.
               hard_snow_d = snow_d(i,jstart) - soft_snow_d(i,jstart)
               snowdmin = max(vegsnowd_xy(i,jstart),hard_snow_d)
               if (snow_d(i,jstart).gt.snowdmin) then
                  if (snow_d(i,jstart) + dh_salt_v(i,jstart).le.snowdmin) then
                     dh_salt_v(i,jstart) = snowdmin - snow_d(i,jstart)
                     Qsalt_v(i,jstart) = tmp_buff(i) - dh_salt_v(i,jstart) * &
                          &            ro_snow * deltay / dt
                  endif
               else
                  Qsalt_v(i,jstart) = 0.0
                  dh_salt_v(i,jstart) = 0.0
               endif
               jstart = jstart + 1
            end if

            do j=jstart,jend
!               new_j = global_map(j)
               dQsalt = Qsalt_v(i,j) - Qsalt_v(i,j-1)
               dh_salt_v(i,j) = (- dt) / ro_snow * dQsalt / deltay

               ! c Make adjustments for the case where there is no snow available
               ! c   on the ground (or captured within the vegetation) to be
               ! c   eroded.
               hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
               snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
               if (snow_d(i,j).gt.snowdmin) then
                  if (snow_d(i,j) + dh_salt_v(i,j).le.snowdmin) then
                     dh_salt_v(i,j) = snowdmin - snow_d(i,j)
                     Qsalt_v(i,j) = Qsalt_v(i,j-1) - dh_salt_v(i,j) * &
                          &            ro_snow * deltay / dt
                  endif
               else
                  Qsalt_v(i,j) = 0.0
                  dh_salt_v(i,j) = 0.0
               endif
            enddo
         enddo
      enddo
   enddo

   index_start = 0
   index_end = 0
   tmp_buff = 0.0

!   sync all

!   do iter = 1, np
   do iter = 1,min(np,salt_iter)
      sync all
      l_Qsalt_v = Qsalt_v
      if(me /= np-1) tmp_buff(:) = l_Qsalt_v(:,1) [me+2]!Remember each processor is me-1
      !c Consider NORTHERLY winds.
      do i=1,nx
         call get_direction_north(vwind_grid, index_start, index_end, npairs, nx, ny, i)
         do k=1,npairs
            jend = index_start(k)
            jstart = index_end(k)-1
            if(jstart == l_ny .and. jend <= jstart) then
!               new_j = global_map(jstart)
               dQsalt = Qsalt_v(i,jstart) - tmp_buff(i)!Qsalt_v(i,j+1)
               dh_salt_v(i,jstart) = (- dt) / ro_snow * dQsalt / deltay
               ! c Make adjustments for the case where there is no snow available
               ! c   on the ground (or captured within the vegetation) to be
               ! c   eroded.
               hard_snow_d = snow_d(i,jstart) - soft_snow_d(i,jstart)
               snowdmin = max(vegsnowd_xy(i,jstart),hard_snow_d)
               if (snow_d(i,jstart).gt.snowdmin) then
                  if (snow_d(i,jstart) + dh_salt_v(i,jstart).le.snowdmin) then
                     dh_salt_v(i,jstart) = snowdmin - snow_d(i,jstart)
                     Qsalt_v(i,jstart) = tmp_buff(i) - dh_salt_v(i,jstart) * &
                          &            ro_snow * deltay / dt
                  endif
               else
                  Qsalt_v(i,jstart) = 0.0
                  dh_salt_v(i,jstart) = 0.0
               endif
               jstart = jstart - 1
            end if

            do j=jstart,jend,-1
!               new_j = global_map(j)
               dQsalt = Qsalt_v(i,j) - Qsalt_v(i,j+1)
               dh_salt_v(i,j) = (- dt) / ro_snow * dQsalt / deltay
               ! c Make adjustments for the case where there is no snow available
               ! c   on the ground (or captured within the vegetation) to be
               ! c   eroded.
               hard_snow_d = snow_d(i,j) - soft_snow_d(i,j)
               snowdmin = max(vegsnowd_xy(i,j),hard_snow_d)
               if (snow_d(i,j).gt.snowdmin) then
                  if (snow_d(i,j) + dh_salt_v(i,j).le.snowdmin) then
                     dh_salt_v(i,j) = snowdmin - snow_d(i,j)
                     Qsalt_v(i,j) = Qsalt_v(i,j+1) - dh_salt_v(i,j) * &
                          &            ro_snow * deltay / dt
                  endif
               else
                  Qsalt_v(i,j) = 0.0
                  dh_salt_v(i,j) = 0.0
               endif
            enddo
         enddo
      enddo
   end do

   sync all

!   call single_gather(Qsalt_v,l_Qsalt_v)

   deallocate(tmp_buff)

! c Update the snow depth changes due to saltation transport from the
! c   the east and west, and north and south.  Also correct dh_salt_u
! c   and dh_salt_v to account for the minimum snow depth.
   eps = 1e-6
   do i=1,nx
      do j=1,l_ny
         weight_u = abs(dh_salt_u(i,j)) / &
              &      (abs(dh_salt_u(i,j)) + abs(dh_salt_v(i,j)) + eps)

         weight_v = abs(dh_salt_v(i,j)) / &
              &      (abs(dh_salt_u(i,j)) + abs(dh_salt_v(i,j)) + eps)

         dh_salt_u(i,j) = weight_u * dh_salt_u(i,j)
         dh_salt_v(i,j) = weight_v * dh_salt_v(i,j)

         if (subgrid_flag.eq.0.0) then
            snow_d(i,j) = snow_d(i,j) + dh_salt_u(i,j) + dh_salt_v(i,j)

            soft_snow_d(i,j) = soft_snow_d(i,j) + dh_salt_u(i,j) + &
                 &        dh_salt_v(i,j)

         elseif (subgrid_flag.eq.1.0) then

!c Just do the erosion.
            snow_d(i,j) = snow_d(i,j) + min(0.0,dh_salt_u(i,j)) + &
                 &        min(0.0,dh_salt_v(i,j))

            soft_snow_d(i,j) = soft_snow_d(i,j) + &
                 &        min(0.0,dh_salt_u(i,j))+ min(0.0,dh_salt_v(i,j))

!c And save an array of what was deposited during this time step.
            snow_d_dep(i,j) = max(0.0,dh_salt_u(i,j)) + &
                 &        max(0.0,dh_salt_v(i,j))

            soft_snow_d_dep(i,j) =  &
                 &        max(0.0,dh_salt_u(i,j))+ max(0.0,dh_salt_v(i,j))
         endif

      enddo
   enddo

   return
 end subroutine getnewdepth

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! c     subroutine subgrid_1(nx,ny,snow_d,
! c    &  index_ue,index_uw,index_vn,index_vs,
! c    &  tabler_nn,tabler_ss,tabler_ee,tabler_ww,
! c    &  tabler_ne,tabler_se,tabler_sw,tabler_nw,uwind_grid,
! c    &  vwind_grid,tabler_dir)

! c This subroutine forces SnowTran-3D's snow accumluation profiles
! c   to be bounded by the equilibrium topographic drift catchment
! c   profiles observed and modeled by Tabler (1975).

! c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
! c   topographic catchments.  Proceedings of the 43rd Annual Western
! c   Snow Conference, San Diego, California, 87-97.

! c     implicit none

! c     include 'snowmodel.inc'

! c     integer i,j,nx,ny
! c     integer k,istart,iend,jstart,jend

! c     real snow_d_extra,snow_sfc,tabler,tabler_dir

! c     integer index_ue(ny_max,2*nx_max+1)
! c     integer index_uw(ny_max,2*nx_max+1)
! c     integer index_vn(nx_max,2*ny_max+1)
! c     integer index_vs(nx_max,2*ny_max+1)

! c     real snow_d(nx_max,ny_max)
! c     real snow_d1(nx_max,ny_max)
! c     real snow_d2(nx_max,ny_max)
! c     real tabler_nn(nx_max,ny_max)
! c     real tabler_ss(nx_max,ny_max)
! c     real tabler_ee(nx_max,ny_max)
! c     real tabler_ww(nx_max,ny_max)
! c     real tabler_ne(nx_max,ny_max)
! c     real tabler_se(nx_max,ny_max)
! c     real tabler_sw(nx_max,ny_max)
! c     real tabler_nw(nx_max,ny_max)
! c     real uwind_grid(nx_max,ny_max)
! c     real vwind_grid(nx_max,ny_max)

! c     real weight_u(nx_max,ny_max)
! c     real weight_v(nx_max,ny_max)

! c This is just a summary of all of the possibilities.
! cc          if(winddir(i,j).gt.337.5.or.winddir(i,j).le.22.5)then
! cc            tabler = tabler_nn(i,j)
! cc          elseif(winddir(i,j).gt.22.5.and.winddir(i,j).le.67.5)then
! cc            tabler = tabler_ne(i,j)
! cc          elseif(winddir(i,j).gt.67.5.and.winddir(i,j).le.112.5)then
! cc            tabler = tabler_ee(i,j)
! cc          elseif(winddir(i,j).gt.112.5.and.winddir(i,j).le.157.5)then
! cc            tabler = tabler_se(i,j)
! cc          elseif(winddir(i,j).gt.157.5.and.winddir(i,j).le.202.5)then
! cc            tabler = tabler_ss(i,j)
! cc          elseif(winddir(i,j).gt.202.5.and.winddir(i,j).le.247.5)then
! cc            tabler = tabler_sw(i,j)
! cc          elseif(winddir(i,j).gt.247.5.and.winddir(i,j).le.292.5)then
! cc            tabler = tabler_ww(i,j)
! cc          elseif(winddir(i,j).gt.292.5.and.winddir(i,j).le.337.5)then
! cc            tabler = tabler_nw(i,j)
! cc          endif

! c Create a copy of the incoming snow depth distribution.  Also define
! c   the u and v weighting functions.
! c     do j=1,ny
! c       do i=1,nx
! c         snow_d1(i,j) = snow_d(i,j)
! c         snow_d2(i,j) = snow_d(i,j)

! c         weight_u(i,j) = abs(uwind_grid(i,j)) /
! c    &          sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
! c         weight_v(i,j) = abs(vwind_grid(i,j)) /
! c    &          sqrt(uwind_grid(i,j)**2 + vwind_grid(i,j)**2)
! c       enddo
! c     enddo

! c Consider WESTERLY winds.
! c     do j=1,ny
! c       do k=1,index_uw(j,1)
! c         istart = index_uw(j,k*2)+1
! c         iend = index_uw(j,k*2+1)
! c         do i=istart,iend

! c           if(tabler_dir.gt.337.5.and.tabler_dir.le.360.0.or.
! c    &        tabler_dir.ge.0.0.and.tabler_dir.le.22.5)then
! c             tabler = tabler_nn(i,j)
! c           elseif(tabler_dir.gt.157.5.and.tabler_dir.le.202.5)then
! c             tabler = tabler_ss(i,j)
! c           elseif(tabler_dir.gt.202.5.and.tabler_dir.le.247.5)then
! c             tabler = tabler_sw(i,j)
! c           elseif(tabler_dir.gt.247.5.and.tabler_dir.le.292.5)then
! c             tabler = tabler_ww(i,j)
! c           elseif(tabler_dir.gt.292.5.and.tabler_dir.le.337.5)then
! c             tabler = tabler_nw(i,j)
! c           endif

! c           snow_sfc = tabler

! c           if (snow_d1(i,j).gt.snow_sfc) then
! c             snow_d_extra = (snow_d1(i,j) - snow_sfc) * weight_u(i,j)
! c             snow_d1(i,j) = snow_d1(i,j) - snow_d_extra
! c             if (i.lt.nx) then
! c               snow_d1(i+1,j) = snow_d1(i+1,j) + snow_d_extra
! c             else
! c               snow_d1(i,j) = snow_d1(i,j)
! c             endif
! c           endif

! c         enddo
! c       enddo
! c     enddo

! c Consider EASTERLY winds.
! c     do j=1,ny
! c       do k=1,index_ue(j,1)
! c         iend = index_ue(j,k*2)
! c         istart = index_ue(j,k*2+1)-1
! c         do i=istart,iend,-1

! c           if(tabler_dir.gt.337.5.and.tabler_dir.le.360.0.or.
! c    &        tabler_dir.ge.0.0.and.tabler_dir.le.22.5)then
! c             tabler = tabler_nn(i,j)
! c           elseif(tabler_dir.gt.22.5.and.tabler_dir.le.67.5)then
! c             tabler = tabler_ne(i,j)
! c           elseif(tabler_dir.gt.67.5.and.tabler_dir.le.112.5)then
! c             tabler = tabler_ee(i,j)
! c           elseif(tabler_dir.gt.112.5.and.tabler_dir.le.157.5)then
! c             tabler = tabler_se(i,j)
! c           elseif(tabler_dir.gt.157.5.and.tabler_dir.le.202.5)then
! c             tabler = tabler_ss(i,j)
! c           endif

! c           snow_sfc = tabler

! c           if (snow_d1(i,j).gt.snow_sfc) then
! c             snow_d_extra = (snow_d1(i,j) - snow_sfc) * weight_u(i,j)
! c             snow_d1(i,j) = snow_d1(i,j) - snow_d_extra
! c             if (i.gt.1) then
! c               snow_d1(i-1,j) = snow_d1(i-1,j) + snow_d_extra
! c             else
! c               snow_d1(i,j) = snow_d1(i,j)
! c             endif
! c           endif
! c         enddo
! c       enddo
! c     enddo

! c Consider SOUTHERLY winds.
! c     do i=1,nx
! c       do k=1,index_vs(i,1)
! c         jstart = index_vs(i,k*2)+1
! c         jend = index_vs(i,k*2+1)
! c         do j=jstart,jend

! c           if(tabler_dir.gt.67.5.and.tabler_dir.le.112.5)then
! c             tabler = tabler_ee(i,j)
! c           elseif(tabler_dir.gt.112.5.and.tabler_dir.le.157.5)then
! c             tabler = tabler_se(i,j)
! c           elseif(tabler_dir.gt.157.5.and.tabler_dir.le.202.5)then
! c             tabler = tabler_ss(i,j)
! c           elseif(tabler_dir.gt.202.5.and.tabler_dir.le.247.5)then
! c             tabler = tabler_sw(i,j)
! c           elseif(tabler_dir.gt.247.5.and.tabler_dir.le.292.5)then
! c             tabler = tabler_ww(i,j)
! c           endif

! c           snow_sfc = tabler

! c           if (snow_d2(i,j).gt.snow_sfc) then
! c             snow_d_extra = (snow_d2(i,j) - snow_sfc) * weight_v(i,j)
! c             snow_d2(i,j) = snow_d2(i,j) - snow_d_extra
! c             if (j.lt.ny) then
! c               snow_d2(i,j+1) = snow_d2(i,j+1) + snow_d_extra
! c             else
! c               snow_d2(i,j) = snow_d2(i,j)
! c             endif
! c           endif
! c         enddo
! c       enddo
! c     enddo

! c Consider NORTHERLY winds.
! c     do i=1,nx
! c       do k=1,index_vn(i,1)
! c         jend = index_vn(i,k*2)
! c         jstart = index_vn(i,k*2+1)-1
! c         do j=jstart,jend,-1

! c           if(tabler_dir.gt.337.5.and.tabler_dir.le.360.0.or.
! c    &        tabler_dir.ge.0.0.and.tabler_dir.le.22.5)then
! c             tabler = tabler_nn(i,j)
! c           elseif(tabler_dir.gt.22.5.and.tabler_dir.le.67.5)then
! c             tabler = tabler_ne(i,j)
! c           elseif(tabler_dir.gt.67.5.and.tabler_dir.le.112.5)then
! c             tabler = tabler_ee(i,j)
! c           elseif(tabler_dir.gt.247.5.and.tabler_dir.le.292.5)then
! c             tabler = tabler_ww(i,j)
! c           elseif(tabler_dir.gt.292.5.and.tabler_dir.le.337.5)then
! c             tabler = tabler_nw(i,j)
! c           endif

! c           snow_sfc = tabler

! c           if (snow_d2(i,j).gt.snow_sfc) then
! c             snow_d_extra = (snow_d2(i,j) - snow_sfc) * weight_v(i,j)
! c             snow_d2(i,j) = snow_d2(i,j) - snow_d_extra
! c             if (j.gt.1) then
! c               snow_d2(i,j-1) = snow_d2(i,j-1) + snow_d_extra
! c             else
! c               snow_d2(i,j) = snow_d2(i,j)
! c             endif
! c           endif

! c         enddo
! c       enddo
! c     enddo

! c Update the snow depths resulting from these redistributions.
! c     do j=1,ny
! c       do i=1,nx
! c         snow_d(i,j) = snow_d1(i,j) * weight_u(i,j) +
! c    &      snow_d2(i,j) * weight_v(i,j)
! c       enddo
! c     enddo

! c Clean up the boundaries.  Make the boundary values equal to
! c   the values just inside the boundaries.
! c     do i=2,nx-1
! c       snow_d(i,1) = snow_d(i,2)
! c       snow_d(i,ny) = snow_d(i,ny-1)
! c     enddo
! c     do j=1,ny
! c       snow_d(1,j) = snow_d(2,j)
! c       snow_d(nx,j) = snow_d(nx-1,j)
! c     enddo

! c     return
! c     end

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine tabler_3d(nx,ny,topo_land,deltax,deltay,&
      &  tabler_ww,tabler_ee,tabler_ss,&
      &  tabler_nn,tabler_ne,tabler_se,&
      &  tabler_sw,tabler_nw,slope_adjust)

! c This subroutine uses Tabler (1975) to define equilibrium profiles
! c   for the topographic drift catchments, for the case of winds
! c   from the EAST, WEST, NORTH, and SOUTH, and anywhere inbetween.

   use snowmodel_inc
   implicit none

!      include 'snowmodel.inc'

   integer nx,ny,irotate_flag

   real deltax,deltay,slope_adjust

   real topo_land(nx,ny)

   real tabler_nn(nx,ny)
   real tabler_ss(nx,ny)
   real tabler_ee(nx,ny)
   real tabler_ww(nx,ny)

   real tabler_ne(nx,ny)
   real tabler_se(nx,ny)
   real tabler_sw(nx,ny)
   real tabler_nw(nx,ny)

! c Here we generate maximum snow accumulation surfaces for n, ne, e,
! c   se, s, sw, w, and nw winds.  I call these "tabler surfaces".
! c
! c They are valid for the wind direction ranges: N=337.5-22.5,
! c   NE=22.5-67.5, E=67.5-112.5, SE=112.5-157.5, S=157.5-202.5,
! c   SW=202.5-247.5, W=247.5-292.5, and NW=292.5-337.5.
! c
! c These Tabler Surfaces define a "potential" snow surface that
! c   represents the maximum possible snow-accumulation depth from winds
! c   coming from these directions.
! c
! c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
! c   topographic catchments.  Proceedings of the 43rd Annual Western
! c   Snow Conference, San Diego, California, 87-97.


! c Consider N winds.
   irotate_flag = 1
   call tabler_n(nx,ny,topo_land,tabler_nn,deltay,&
        &  irotate_flag,slope_adjust)

!c Consider NE winds.
   irotate_flag = 2
   call tabler_e(nx,ny,topo_land,tabler_ne,1.41*deltax,&
        &  irotate_flag,slope_adjust)

!c Consider E winds.
   irotate_flag = 1
   call tabler_e(nx,ny,topo_land,tabler_ee,deltax, &
        &  irotate_flag,slope_adjust)

!c Consider SE winds.
   irotate_flag = 2
   call tabler_s(nx,ny,topo_land,tabler_se,1.41*deltay, &
        &  irotate_flag,slope_adjust)

!c Consider S winds.
   irotate_flag = 1
   call tabler_s(nx,ny,topo_land,tabler_ss,deltay,&
        &  irotate_flag,slope_adjust)

!c Consider SW winds.
   irotate_flag = 2
   call tabler_w(nx,ny,topo_land,tabler_sw,1.41*deltax, &
        &  irotate_flag,slope_adjust)

!c Consider W winds.
   irotate_flag = 1
   call tabler_w(nx,ny,topo_land,tabler_ww,deltax, &
        &  irotate_flag,slope_adjust)

!c Consider NW winds.
   irotate_flag = 2
   call tabler_n(nx,ny,topo_land,tabler_nw,1.41*deltay, &
        &  irotate_flag,slope_adjust)


   return
 end subroutine tabler_3d

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine tabler_w(nx,ny,topo_land,tabler_ww,deltax, &
      &  irotate_flag,slope_adjust)

! c This subroutine uses Tabler (1975) to define equilibrium profiles
! c   for the topographic drift catchments, for the case of winds
! c   from the WEST and SOUTHWEST.

! c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
! c   topographic catchments.  Proceedings of the 43rd Annual Western
! c   Snow Conference, San Diego, California, 87-97.
   use snowmodel_inc
   implicit none

   !   include 'snowmodel.inc'

   integer nx,ny,i,j,istart,iend,JJ,irotate_flag,nny,nnx,&
        &  ii,iii,maxlines
   real deltax,y,slope_adjust,xmax_slope,dx,test,&
        &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

   real topo_land(nx,ny)
   real tabler_ww(nx,ny)
   real tabler(nx_max)
   real topo_line(nx_max)
   real drift_start_topo(nx_max)
   integer, parameter :: nx_max_tabler = nx_max*1000
   real topo_1m(nx_max_tabler)
   real tabler_1m(nx_max_tabler)
   real drift_start_topo_1m(nx_max_tabler)

! c This program:
! c   1) takes the coarse model topography and generates 1.0-m grid
! c        increment topo lines;
! c   2) uses those to generate the Tabler surface profiles; and
! c   3) extracts the profiles at the coarse model grid cells.

!c Fill the snow array with topo data.
   do j=1,ny
      do i=1,nx
         tabler_ww(i,j) = topo_land(i,j)
      enddo
   enddo

! c Define the length of the j-looping, depending on whether there
! c   is any rotation done to get 45 deg winds.
   if (irotate_flag.eq.2) then
      nny = nx+ny-1
   else
      nny = ny
   endif

!c Required parameters.
   if (irotate_flag.eq.2) then
      dx = 1.41
      test = amod(deltax/1.41,dx/1.41)
   else
      dx = 1.0
      test = amod(deltax,dx)
   endif

   xmax_slope = -0.20

! c This Tabler program has not been made general enough to deal
! c   with deltax and deltay values that are not evenly divisible
! c   by 1.0 (or 1.41 for diagonal profiles).
   if (abs(test).gt.1.0e10-5) then
      print *,'To generate the Tabler surfaces, deltax and deltay'
      print *,'  must be evenly divisible by 1.0 or 1.41.'
      print *,'  deltax = ',deltax
      stop
   endif

! c Define the number of 1-m grid cells in each model grid cell, and
! c   calculate how many of these are in the entire nx domain.
   nnx = nint(deltax/dx)
   maxlines = (nx - 1) * nnx + 1

! c Define the starting and ending points of the line we are going to
! c   work with.  This is done to make sure we are never looking
! c   outside the data for numbers to work with.
   istart = 1
   if (irotate_flag.eq.2) then
      iend = maxlines - (32+11+10+11)
   else
      iend = maxlines - (45+15+15+15)
   endif

!c Extract the line we are going to work with.
   do j=1,nny

      if (irotate_flag.eq.2) then
         do i=1,nx
            JJ = j + i - nx
            drift_start_topo(i) = 0.0
            if (JJ.le.0) then
               tabler(i) = tabler_ww(1-j+nx,1)
               topo_line(i) = tabler_ww(1-j+nx,1)
            elseif (JJ.gt.ny) then
               tabler(i) = tabler(i-1)
               topo_line(i) = tabler(i-1)
            else
               tabler(i) = tabler_ww(i,JJ)
               topo_line(i) = tabler_ww(i,JJ)
            endif
         enddo
      else
         do i=1,nx
            tabler(i) = tabler_ww(i,j)
            topo_line(i) = topo_land(i,j)
            drift_start_topo(i) = 0.0
         enddo
      endif

! c To build the 1.0 m line, use linear interpolation between the
! c   model topo data.  Include the end point.
      do i=1,nx-1
         do ii=1,nnx
            iii = (i - 1) * nnx + ii
            x1 = 0.0
            x = real(ii - 1) * dx
            y2 = topo_line(i+1)
            y1 = topo_line(i)
            topo_1m(iii) = y1 + ((y2 - y1)/deltax) * (x - x1)
         enddo
      enddo
      topo_1m((nx - 1) * nnx + 1) = topo_line(nx)

! c Use this topo array to be the starting point for generating the
! c   Tabler surfaces.
      do i=1,maxlines
         tabler_1m(i) = topo_1m(i)
         drift_start_topo_1m(i) = 0.0
      enddo

! c Run the Tabler model.
      do i=istart,iend
         if (irotate_flag.eq.2) then
            t1 = tabler_1m(i)
            t2 = tabler_1m(i+31)
            t3 = tabler_1m(i+31+11)
            t4 = tabler_1m(i+31+21)
            t5 = tabler_1m(i+31+32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i+32) = max(topo_1m(i+32),&
                 &        tabler_1m(i+31) + y * slope_adjust * dx)
         else
            t1 = tabler_1m(i)
            t2 = tabler_1m(i+44)
            t3 = tabler_1m(i+44+15)
            t4 = tabler_1m(i+44+30)
            t5 = tabler_1m(i+44+45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i+45) = max(topo_1m(i+45),&
                 &        tabler_1m(i+44) + y * slope_adjust * dx)
         endif
      enddo

!c Extract the profile at the model grid points.
      do i=1,nx
         ii = (i - 1) * nnx + 1
         tabler(i) = tabler_1m(ii)
         drift_start_topo(i) = drift_start_topo_1m(ii)
      enddo

!c Use the 1-D arrays to fill in the 2-D tabler-surface array.
      do i=1,nx
         if (irotate_flag.eq.2) then
            JJ = j + i - nx
            if (JJ.ge.1 .and. JJ.le.ny) then
               tabler_ww(i,JJ) = tabler(i) + drift_start_topo(i)
            endif
         else
            tabler_ww(i,j) = tabler(i) + drift_start_topo(i)
         endif
      enddo

   enddo

!c Convert snow_traps back to actual snow depths instead of
!c   depth plus topography.
   do j=1,ny
      do i=1,nx
         tabler_ww(i,j) = tabler_ww(i,j) - topo_land(i,j)
      enddo
   enddo

   return
 end subroutine tabler_w

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine tabler_e(nx,ny,topo_land,tabler_ee,deltax,&
      &  irotate_flag,slope_adjust)

! c This subroutine uses Tabler (1975) to define equilibrium profiles
! c   for the topographic drift catchments, for the case of winds
! c   from the EAST and NORTHEAST.

! c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
! c   topographic catchments.  Proceedings of the 43rd Annual Western
! c   Snow Conference, San Diego, California, 87-97.

   use snowmodel_inc
   implicit none

   !      include 'snowmodel.inc'

   integer nx,ny,i,j,istart,iend,JJ,irotate_flag,nny,nnx,&
        &  ii,iii,maxlines
   real deltax,y,slope_adjust,xmax_slope,dx,test,&
        &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

   real topo_land(nx,ny)
   real tabler_ee(nx,ny)
   real tabler(nx_max)
   real topo_line(nx_max)
   real drift_start_topo(nx_max)

   integer, parameter :: nx_max_tabler = nx_max*1000
   real topo_1m(nx_max_tabler)
   real tabler_1m(nx_max_tabler)
   real drift_start_topo_1m(nx_max_tabler)

! c This program:
! c   1) takes the coarse model topography and generates 1.0-m grid
! c        increment topo lines;
! c   2) uses those to generate the Tabler surface profiles; and
! c   3) extracts the profiles at the coarse model grid cells.

!c Fill the snow array with topo data.
   do j=1,ny
      do i=1,nx
         tabler_ee(i,j) = topo_land(i,j)
      enddo
   enddo

! c Define the length of the j-looping, depending on whether there
! c   is any rotation done to get 45 deg winds.
   if (irotate_flag.eq.2) then
      nny = nx+ny-1
   else
      nny = ny
   endif

!c Required parameters.
   if (irotate_flag.eq.2) then
      dx = 1.41
      test = amod(deltax/1.41,dx/1.41)
   else
      dx = 1.0
      test = amod(deltax,dx)
   endif

   xmax_slope = -0.20

! c This Tabler program has not been made general enough to deal
! c   with deltax and deltay values that are not evenly divisible
! c   by 1.0 (or 1.41 for diagonal profiles).
   if (abs(test).gt.1.0e10-5) then
      print *,'To generate the Tabler surfaces, deltax and deltay'
      print *,'  must be evenly divisible by 1.0 or 1.41.'
      print *,'  deltax = ',deltax
      stop
   endif

!c Define the number of 1-m grid cells in each model grid cell, and
!c   calculate how many of these are in the entire nx domain.
   nnx = nint(deltax/dx)
   maxlines = (nx - 1) * nnx + 1

! c Define the starting and ending points of the line we are going to
! c   work with.  This is done to make sure we are never looking
! c   outside the data for numbers to work with.
   istart = maxlines
   if (irotate_flag.eq.2) then
      iend = 1 + (32+11+10+11)
   else
      iend = 1 + (45+15+15+15)
   endif

!c Extract the line we are going to work with.
   do j=1,nny

      if (irotate_flag.eq.2) then
         do i=1,nx
            JJ = j + i - nx
            drift_start_topo(i) = 0.0
            if (JJ.le.0) then
               tabler(i) = tabler_ee(1-j+nx,1)
               topo_line(i) = tabler_ee(1-j+nx,1)
            elseif (JJ.gt.ny) then
               tabler(i) = tabler(i-1)
               topo_line(i) = tabler(i-1)
            else
               tabler(i) = tabler_ee(i,JJ)
               topo_line(i) = tabler_ee(i,JJ)
            endif
         enddo
      else
         do i=1,nx
            tabler(i) = tabler_ee(i,j)
            topo_line(i) = topo_land(i,j)
            drift_start_topo(i) = 0.0
         enddo
      endif

! c To build the 1.0 m line, use linear interpolation between the
! c   model topo data.  Include the end point.
      do i=1,nx-1
         do ii=1,nnx
            iii = (i - 1) * nnx + ii
            x1 = 0.0
            x = real(ii - 1) * dx
            y2 = topo_line(i+1)
            y1 = topo_line(i)
            topo_1m(iii) = y1 + ((y2 - y1)/deltax) * (x - x1)
         enddo
      enddo
      topo_1m((nx - 1) * nnx + 1) = topo_line(nx)

! c Use this topo array to be the starting point for generating the
! c   Tabler surfaces.
      do i=1,maxlines
         tabler_1m(i) = topo_1m(i)
         drift_start_topo_1m(i) = 0.0
      enddo

!c Run the Tabler model.
      do i=istart,iend,-1
         if (irotate_flag.eq.2) then
            t1 = tabler_1m(i)
            t2 = tabler_1m(i-31)
            t3 = tabler_1m(i-31-11)
            t4 = tabler_1m(i-31-21)
            t5 = tabler_1m(i-31-32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i-32) = max(topo_1m(i-32),&
                 &        tabler_1m(i-31) + y * slope_adjust * dx)
         else
            t1 = tabler_1m(i)
            t2 = tabler_1m(i-44)
            t3 = tabler_1m(i-44-15)
            t4 = tabler_1m(i-44-30)
            t5 = tabler_1m(i-44-45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(i-45) = max(topo_1m(i-45),&
                 &        tabler_1m(i-44) + y * slope_adjust * dx)
         endif
      enddo

!c Extract the profile at the model grid points.
      do i=1,nx
         ii = (i - 1) * nnx + 1
         tabler(i) = tabler_1m(ii)
         drift_start_topo(i) = drift_start_topo_1m(ii)
      enddo

!c Use the 1-D arrays to fill in the 2-D tabler-surface array.
      do i=1,nx
         if (irotate_flag.eq.2) then
            JJ = j + i - nx
            if (JJ.ge.1 .and. JJ.le.ny) then
               tabler_ee(i,JJ) = tabler(i) + drift_start_topo(i)
            endif
         else
            tabler_ee(i,j) = tabler(i) + drift_start_topo(i)
         endif
      enddo

   enddo

! c Convert snow_traps back to actual snow depths instead of
! c   depth plus topography.
   do j=1,ny
      do i=1,nx
         tabler_ee(i,j) = tabler_ee(i,j) - topo_land(i,j)
      enddo
   enddo

   return
 end subroutine tabler_e

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine tabler_s(nx,ny,topo_land,tabler_ss,deltay,&
      &  irotate_flag,slope_adjust)

! c This subroutine uses Tabler (1975) to define equilibrium profiles
! c   for the topographic drift catchments, for the case of winds
! c   from the SOUTH and SOUTHEAST.

! c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
! c   topographic catchments.  Proceedings of the 43rd Annual Western
! c   Snow Conference, San Diego, California, 87-97.

   use snowmodel_inc
   implicit none

!      include 'snowmodel.inc'

   integer nx,ny,i,j,jstart,jend,II,irotate_flag,nny,nnx,&
        &  jj,jjj,maxlines
   real deltay,y,slope_adjust,xmax_slope,dy,test,&
        &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

   real topo_land(nx,ny)
   real tabler_ss(nx,ny)
   real tabler(ny_max)
   real topo_line(ny_max)
   real drift_start_topo(ny_max)

   integer,parameter :: ny_max_tabler=ny_max*1000
   real topo_1m(ny_max_tabler)
   real tabler_1m(ny_max_tabler)
   real drift_start_topo_1m(ny_max_tabler)

! c This program:
! c   1) takes the coarse model topography and generates 1.0-m grid
! c        increment topo lines;
! c   2) uses those to generate the Tabler surface profiles; and
! c   3) extracts the profiles at the coarse model grid cells.

!c Fill the snow array with topo data.
   do j=1,ny
      do i=1,nx
         tabler_ss(i,j) = topo_land(i,j)
      enddo
   enddo

! c Define the length of the i-looping, depending on whether there
! c   is any rotation done to get 45 deg winds.
   if (irotate_flag.eq.2) then
      nnx = nx+ny-1
   else
      nnx = nx
   endif

!c Required parameters.
   if (irotate_flag.eq.2) then
      dy = 1.41
      test = amod(deltay/1.41,dy/1.41)
   else
      dy = 1.0
      test = amod(deltay,dy)
   endif

   xmax_slope = -0.20

!c This Tabler program has not been made general enough to deal
!c   with deltax and deltay values that are not evenly divisible
!c   by 1.0 (or 1.41 for diagonal profiles).
   if (abs(test).gt.1.0e10-5) then
      print *,'To generate the Tabler surfaces, deltax and deltay'
      print *,'  must be evenly divisible by 1.0 or 1.41.'
      print *,'  deltay = ',deltay
      stop
   endif

! c Define the number of 1-m grid cells in each model grid cell, and
! c   calculate how many of these are in the entire nx domain.
   nny = nint(deltay/dy)
   maxlines = (ny - 1) * nny + 1

! c Define the starting and ending points of the line we are going to
! c   work with.  This is done to make sure we are never looking
! c   outside the data for numbers to work with.
   jstart = 1
   if (irotate_flag.eq.2) then
      jend = maxlines - (32+11+10+11)
   else
      jend = maxlines - (45+15+15+15)
   endif

! c Extract the line we are going to work with.
   do i=1,nnx
      if (irotate_flag.eq.2) then
         do j=1,ny
            II = i - j + 1
            drift_start_topo(j) = 0.0
            if (II.le.0) then
               tabler(j) = tabler(j-1)
               topo_line(j) = tabler(j-1)
            elseif (II.gt.nx) then
               tabler(j) = tabler_ss(nx,i-nx+1)
               topo_line(j) = tabler_ss(nx,i-nx+1)
            else
               tabler(j) = tabler_ss(II,j)
               topo_line(j) = tabler_ss(II,j)
            endif
         enddo
      else
         do j=1,ny
            tabler(j) = tabler_ss(i,j)
            topo_line(j) = topo_land(i,j)
            drift_start_topo(j) = 0.0
         enddo
      endif

! c To build the 1.0 m line, use linear interpolation between the
! c   model topo data.  Include the end point.
      do j=1,ny-1
         do jj=1,nny
            jjj = (j - 1) * nny + jj
            x1 = 0.0
            x = real(jj - 1) * dy
            y2 = topo_line(j+1)
            y1 = topo_line(j)
            topo_1m(jjj) = y1 + ((y2 - y1)/deltay) * (x - x1)
         enddo
      enddo
      topo_1m((ny - 1) * nny + 1) = topo_line(ny)

!c Use this topo array to be the starting point for generating the
!c   Tabler surfaces.
      do j=1,maxlines
         tabler_1m(j) = topo_1m(j)
         drift_start_topo_1m(j) = 0.0
      enddo

!c Run the Tabler model.
      do j=jstart,jend
         if (irotate_flag.eq.2) then
            t1 = tabler_1m(j)
            t2 = tabler_1m(j+31)
            t3 = tabler_1m(j+31+11)
            t4 = tabler_1m(j+31+21)
            t5 = tabler_1m(j+31+32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j+32) = max(topo_1m(j+32),&
                 &        tabler_1m(j+31) + y * slope_adjust * dy)
         else
            t1 = tabler_1m(j)
            t2 = tabler_1m(j+44)
            t3 = tabler_1m(j+44+15)
            t4 = tabler_1m(j+44+30)
            t5 = tabler_1m(j+44+45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j+45) = max(topo_1m(j+45),&
                 &        tabler_1m(j+44) + y * slope_adjust * dy)
         endif
      enddo

!c Extract the profile at the model grid points.
      do j=1,ny
         jj = (j - 1) * nny + 1
         tabler(j) = tabler_1m(jj)
         drift_start_topo(j) = drift_start_topo_1m(jj)
      enddo

!c Use the 1-D arrays to fill in the 2-D tabler-surface array.
      do j=1,ny
         if (irotate_flag.eq.2) then
            II = i - j + 1
            if (II.ge.1 .and. II.le.nx) then
               tabler_ss(II,j) = tabler(j) + drift_start_topo(j)
            endif
         else
            tabler_ss(i,j) = tabler(j) + drift_start_topo(j)
         endif
      enddo

   enddo

!c Convert snow_traps back to actual snow depths instead of
!c   depth plus topography.
   do j=1,ny
      do i=1,nx
         tabler_ss(i,j) = tabler_ss(i,j) - topo_land(i,j)
      enddo
   enddo

   return
 end subroutine tabler_s

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine tabler_n(nx,ny,topo_land,tabler_nn,deltay,&
      &  irotate_flag,slope_adjust)

! c This subroutine uses Tabler (1975) to define equilibrium profiles
! c   for the topographic drift catchments, for the case of winds
! c   from the NORTH and NORTHWEST.

! c Tabler, R. D., 1975: Predicting profiles of snowdrifts in
! c   topographic catchments.  Proceedings of the 43rd Annual Western
! c   Snow Conference, San Diego, California, 87-97.

   use snowmodel_inc
   implicit none

!      include 'snowmodel.inc'

   integer nx,ny,i,j,jstart,jend,II,irotate_flag,nny,nnx,&
     &  jj,jjj,maxlines
   real deltay,y,slope_adjust,xmax_slope,dy,test,&
        &  x,x1,x2,x3,x4,y1,y2,t1,t2,t3,t4,t5

   real topo_land(nx,ny)
   real tabler_nn(nx,ny)
   real tabler(ny_max)
   real topo_line(ny_max)
   real drift_start_topo(ny_max)

   integer, parameter :: ny_max_tabler = ny_max*1000
   real topo_1m(ny_max_tabler)
   real tabler_1m(ny_max_tabler)
   real drift_start_topo_1m(ny_max_tabler)

! c This program:
! c   1) takes the coarse model topography and generates 1.0-m grid
! c        increment topo lines;
! c   2) uses those to generate the Tabler surface profiles; and
! c   3) extracts the profiles at the coarse model grid cells.

! c Fill the snow array with topo data.
   do j=1,ny
      do i=1,nx
         tabler_nn(i,j) = topo_land(i,j)
      enddo
   enddo

! c Define the length of the i-looping, depending on whether there
! c   is any rotation done to get 45 deg winds.
   if (irotate_flag.eq.2) then
      nnx = nx+ny-1
   else
      nnx = nx
   endif

!c Required parameters.
   if (irotate_flag.eq.2) then
      dy = 1.41
      test = amod(deltay/1.41,dy/1.41)
   else
      dy = 1.0
      test = amod(deltay,dy)
   endif

   xmax_slope = -0.20

! c This Tabler program has not been made general enough to deal
! c   with deltax and deltay values that are not evenly divisible
! c   by 1.0 (or 1.41 for diagonal profiles).
   if (abs(test).gt.1.0e10-5) then
      print *,'To generate the Tabler surfaces, deltax and deltay'
      print *,'  must be evenly divisible by 1.0 or 1.41.'
      print *,'  deltay = ',deltay
      stop
   endif

! c Define the number of 1-m grid cells in each model grid cell, and
! c   calculate how many of these are in the entire nx domain.
   nny = nint(deltay/dy)
   maxlines = (ny - 1) * nny + 1

! c Define the starting and ending points of the line we are going to
! c   work with.  This is done to make sure we are never looking
! c   outside the data for numbers to work with.
   jstart = maxlines
   if (irotate_flag.eq.2) then
      jend = 1 + (32+11+10+11)
   else
      jend = 1 + (45+15+15+15)
   endif

!c Extract the line we are going to work with.
   do i=1,nnx
      if (irotate_flag.eq.2) then
         do j=1,ny
            II = i - j + 1
            drift_start_topo(j) = 0.0
            if (II.le.0) then
               tabler(j) = tabler(j-1)
               topo_line(j) = tabler(j-1)
            elseif (II.gt.nx) then
               tabler(j) = tabler_nn(nx,i-nx+1)
               topo_line(j) = tabler_nn(nx,i-nx+1)
            else
               tabler(j) = tabler_nn(II,j)
               topo_line(j) = tabler_nn(II,j)
            endif
         enddo
      else
         do j=1,ny
            tabler(j) = tabler_nn(i,j)
            topo_line(j) = topo_land(i,j)
            drift_start_topo(j) = 0.0
         enddo
      endif

! c To build the 1.0 m line, use linear interpolation between the
! c   model topo data.  Include the end point.
      do j=1,ny-1
         do jj=1,nny
            jjj = (j - 1) * nny + jj
            x1 = 0.0
            x = real(jj - 1) * dy
            y2 = topo_line(j+1)
            y1 = topo_line(j)
            topo_1m(jjj) = y1 + ((y2 - y1)/deltay) * (x - x1)
         enddo
      enddo
      topo_1m((ny - 1) * nny + 1) = topo_line(ny)

!c Use this topo array to be the starting point for generating the
!c   Tabler surfaces.
      do j=1,maxlines
         tabler_1m(j) = topo_1m(j)
         drift_start_topo_1m(j) = 0.0
      enddo

!c Run the Tabler model.
      do j=jstart,jend,-1
         if (irotate_flag.eq.2) then
            t1 = tabler_1m(j)
            t2 = tabler_1m(j-31)
            t3 = tabler_1m(j-31-11)
            t4 = tabler_1m(j-31-21)
            t5 = tabler_1m(j-31-32)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j-32) = max(topo_1m(j-32),&
                 &        tabler_1m(j-31) + y * slope_adjust * dy)
         else
            t1 = tabler_1m(j)
            t2 = tabler_1m(j-44)
            t3 = tabler_1m(j-44-15)
            t4 = tabler_1m(j-44-30)
            t5 = tabler_1m(j-44-45)

            x1 = (t2 - t1) / 45.0
            x2 = max((t3 - t2) / 15.0,xmax_slope)
            x3 = max((t4 - t3) / 15.0,xmax_slope)
            x4 = max((t5 - t4) / 15.0,xmax_slope)

            y = 0.25*x1 + 0.55*x2 + 0.15*x3 + 0.05*x4

            tabler_1m(j-45) = max(topo_1m(j-45),&
                 &        tabler_1m(j-44) + y * slope_adjust * dy)
         endif
      enddo

!c Extract the profile at the model grid points.
      do j=1,ny
         jj = (j - 1) * nny + 1
         tabler(j) = tabler_1m(jj)
         drift_start_topo(j) = drift_start_topo_1m(jj)
      enddo

!c Use the 1-D arrays to fill in the 2-D tabler-surface array.
      do j=1,ny
         if (irotate_flag.eq.2) then
            II = i - j + 1
            if (II.ge.1 .and. II.le.nx) then
               tabler_nn(II,j) = tabler(j) + drift_start_topo(j)
            endif
         else
            tabler_nn(i,j) = tabler(j) + drift_start_topo(j)
         endif
      enddo

   enddo

!c Convert snow_traps back to actual snow depths instead of
!c   depth plus topography.
   do j=1,ny
      do i=1,nx
         tabler_nn(i,j) = tabler_nn(i,j) - topo_land(i,j)
      enddo
   enddo

   return
 end subroutine tabler_n

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine surface_snow_1(Tair,windspd_2m,prec,ro_soft_snow,&
      &    Utau_t,ro_soft_snow_old,dt,ro_nsnow)

   implicit none

   real C,alfa,ro_min,ro_max,prec,Tair,ro_nsnow,dt,Tf,&
        &  Utau_t,ro_soft_snow_old,ro_soft_snow,windspd_2m

!c Define the density rate coefficients.
   C = 0.10

!c Define alfa.
   alfa = 0.2

! c Define the minimum and maximum snow density that should be
! c   simulated.
   ro_min = 50.0
   ro_max = 450.0

!c Freezing temperature.
   Tf = 273.15

! c Calculate the new snow density.  First calculate this under the
! c   assumption of no wind using the standard SnowModel formulation,
! c   then calculate the offset for wind speeds > 5 m/s.
   if (prec.gt.0.0) then
      ro_soft_snow = ro_nsnow
      ro_soft_snow = min(ro_soft_snow,ro_max)
      ro_soft_snow = max(ro_soft_snow,ro_min)
      if (ro_soft_snow.le.300.0) then
         Utau_t = 0.10 * exp(0.003 * ro_soft_snow)
      else
         Utau_t = 0.005 * exp(0.013 * ro_soft_snow)
      endif
      ro_soft_snow_old = ro_soft_snow
   else
      call surface_snow_2(ro_soft_snow_old,ro_soft_snow,Utau_t, &
           &    dt,Tair,windspd_2m,C,ro_max,ro_min,alfa,Tf)
      ro_soft_snow_old = ro_soft_snow
   endif

   return
 end subroutine surface_snow_1

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine surface_snow_2(ro_soft_snow_old,ro_soft_snow,Utau_t,&
      &  dt,Tair,windspd_2m,C,ro_max,ro_min,alfa,Tf)

   implicit none

   real Tf,Tsnow,Tair,ro_soft_snow,dt,A1,A2,B,U,C,&
           &  windspd_2m,Utau_t,ro_soft_snow_old,ro_max,ro_min,alfa

   A1 = 0.0013
   A2 = 0.021
   B = 0.08

! c Evolve the near-surface snow density under the influence of
! c   temperature and snow-transporting wind speeds.

! c Assume that the near-surface snow temperature equals the air
! c   temperature, but is not above the melting point.
!       Tsnow = min(Tf,Tair)

! c Update the snow density of the soft snow layer.  Eliminate the
! c   wind speed influence for speeds below 5 m/s, but account for it
! c   if speeds are >= 5 m/s.
   if (windspd_2m.ge.5.0) then
      U = 5.0 + 15.0 * (1.0 - exp(-(alfa*(windspd_2m - 5.0))))
   else
      U = 1.0
   endif

   ro_soft_snow = ro_soft_snow_old + dt * &
        &  (C * A1 * U * ro_soft_snow_old * &
        &  exp((- B)*(Tf-Tsnow)) * exp((- A2)*ro_soft_snow_old))

!c Bound the calculated density.
   ro_soft_snow = min(ro_max,ro_soft_snow)
   ro_soft_snow = max(ro_min,ro_soft_snow)

!c Calculate the snow threshold friction velocity.
   if (ro_soft_snow.le.300.0) then
      Utau_t = 0.10 * exp(0.003 * ro_soft_snow)
   else
      Utau_t = 0.005 * exp(0.013 * ro_soft_snow)
   endif

   return
 end subroutine surface_snow_2

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
