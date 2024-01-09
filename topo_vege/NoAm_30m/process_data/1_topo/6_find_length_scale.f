c 5_find_length_scale.f

c This program figures out which scaling factor to use.  This
c   allows the code to use pre-gridded, coarser resolution data
c   for coarser resolution SnowModel simulations.  For example,
c   if you want a 10-km SnowModel grid, there is no reason to
c   start your re-gridding with a 30-m dem dataset.  This is
c   done through the iscale parameter.

      implicit none

      integer iscale
      real deltax,deltay,deltaxy

      character*2 cscale
      character*9 scale

c Figure out what iscale should be.  iscale defines whether you
c   want coarser-scale outputs.  For example, iscale = 4 will
c   give you 4 arc-sec data outputs.

c Read the SnowModel grid increment.
      open (21,file='outputs/1_grid_increment.dat')
      read (21,*) deltax,deltay
      deltaxy = min(deltax,deltay)

c Here I am defining iscale to provide at least 5 lat-lon grid
c   cells for ever SnowModel grid cell that is ultimately
c   created.  Note that iscale must divide evenly into 3600.
      if (deltaxy.le.150.0) then
        iscale = 1
      elseif (deltaxy.le.300.0) then
        iscale = 2
      elseif (deltaxy.le.450.0) then
        iscale = 3
      elseif (deltaxy.le.600.0) then
        iscale = 4
      elseif (deltaxy.le.750.0) then
        iscale = 5
      elseif (deltaxy.le.900.0) then
        iscale = 6
      elseif (deltaxy.le.1200.0) then
        iscale = 8
      elseif (deltaxy.le.1500.0) then
        iscale = 10
      elseif (deltaxy.le.1800.0) then
        iscale = 12
      elseif (deltaxy.le.2250.0) then
        iscale = 15
      elseif (deltaxy.le.3000.0) then
        iscale = 20
      elseif (deltaxy.le.4500.0) then
        iscale = 30
      elseif (deltaxy.gt.4500.0) then
        iscale = 30
      endif

c Turn iscale into a two-character string.
      write(cscale,'(i2.2)') iscale

c Define the file path associated with this scaling factor.
      scale = 'scale_'//cscale//'/'

c Save the data.
      open (31,file='outputs/6_grid_scale.dat')
      write (31,91) scale

   91 format (a9)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

