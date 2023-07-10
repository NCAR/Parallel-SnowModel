c snowmodel.inc

      integer nx_max,ny_max,nstns_max,nvegtypes,max_time_steps,
     &  max_obs_dates,nz_max,n_print_vars

c THESE FIRST THREE PARAMETER VALUES OFTEN NEED TO BE CHANGED
c   FOR BIG AND LONG MODEL SIMULATIONS.

c nx_max and ny_max define the maximum nx and ny grid points
c   that you can have in your model run.  If you want to run a
c   larger domain, then you have to change these numbers and
c   recompile the code.
      parameter (nx_max=1000,ny_max=1000)

c max_time_steps defines the maximum number of time steps that
c   will be used in the current compliled version of the code.
c   If you want to run a longer time domain, then you have to
c   change this number and recompile the code.
      parameter (max_time_steps=8784)

c nstns_max is the maximum number of stations that can be used
c   in the data assimilation routines.
      parameter (nstns_max=10000)

c max_obs_dates is used in the data assimilation routines.  It's
c   value must be greater than the number of observation dates
c   in the entire simulation + (plus) the number of years in the
c   simulation.  For example, for a 6-year simulation with two
c   observation dates in each year, you would set max_obs_dates
c   to be at least = 18.
      parameter (max_obs_dates=12)

c If you are running the multi-layer snow model (even with a single
c   layer) nz_max must be at least one greater than max_layers in
c   snowmodel.par.  This is because the model will build a new layer
c   with the new snowfall and then it is merged with the layer below
c   if you only want a single snow layer.  If you are running
c   SnowModel's original single layer model, nz_max can be 1 (but if
c   nz_max=2 it will avoid a warning message if you are compiling
c   the code with gfortran).
c     parameter (nz_max=25)
      parameter (nz_max=2)

c This is the number of print variables that are controlled by
c   the variable list in snowmodel.par.
      parameter (n_print_vars=30)

c nvegtypes is the number of vegetation types used in the model
c   run.  If you change this then you have made some big changes
c   in the codes' vegetation description.
      parameter (nvegtypes=30)
