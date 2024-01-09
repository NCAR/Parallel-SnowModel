To process (re)analysis met forcing for your SnowModel run,
you have to work through each of the 7 steps in this directory
(you can skip steps 5 and 6 if you don't need them).  Specifically,
you have to:

1_topo_lonlat
  Don't need to do anything. These are inputs associated with the
  specific (re)analysis.

2_define_points
  Run "1_points_from_ll_to_proj.script".

  Edit "2_pts_sm_domain.f" between "BEGIN USER INPUT" and "END
  USER INPUT".

  Compile and run "2_pts_sm_domain.f".

3_figs
  Run "met_points.gs" in GrADS to plot the met points relative to
  your simulation domain, to see if you like the results.

4_maxiter_offset
  Edit "start_end_dates_maxiter_ioffset.f" between "BEGIN USER
  INPUT" and "END USER INPUT" to define your temporal domain.

  Compile and run "start_end_dates_maxiter_ioffset.f".

5_mk_gdat_optional
  Edit "mk_SM_gdat_met_file.f" between "BEGIN USER INPUT" and
  "END USER INPUT".

  Compile and run "mk_SM_gdat_met_file.f".

6_fix_drizzle_optional
  Edit "fix_drizzle.f" to define: nlayers_target, iimo_start,
  iidy_start, iimo_end, iidy_end.

  Compile and run "fix_drizzle.f".

7_mk_mm
  Edit "mk_micromet_XXXXX.f" between "BEGIN USER INPUT" and
  "END USER INPUT".

  Compile and run "mk_micromet_XXXXX.f".

