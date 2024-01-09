c pts_sm_domain.f

c This program tosses all points outside the SnowModel
c   simulation domain plus a small border of met forcing
c   points around the SnowModel simulation domain that
c   will be used during the interpolation to near-border
c   SnowModel points.

c This example uses ERA5 reanalysis.  If you are using
c   some other (re)analysis dataset, you will have to change
c   the nnx and nny parameter definitions below.

c ERA5.
      parameter (nnx=1440,nny=321)

      real topo(nnx,nny)
      real x(nnx,nny)
      real y(nnx,nny)

c BEGIN USER INPUT.
c BEGIN USER INPUT.
c BEGIN USER INPUT.

c Define the path where the topo-vege grid information is located.
      character path1*(*)
      parameter (path1=
     &  '../../../topo_vege/NoAm_30m/')

      character path2*(*)
      parameter (path2= '../1_topo_lonlat/')

c Point to the coordinate and topography data files.
      open (21,file='met_points_proj.dat')
      open (22,file=path2//'era5_topo.txt')

c Define whether you want to extract a single point that is
c   closest to the center of your simulation domain (single
c   = 1.0), or multiple points (single = 0.0) in and around
c   your SnowModel simulation domain.
      single = 0.0
c     single = 1.0

c Define whether you want to thin the points as you go north
c   (thin = 1.0), or not (thin = 0.0).  This is critical for
c   met data that are on a lat-lon grid at high latitudes, to
c   make sure the 'met station' grid points are evenly spaced
c   in x and y.
c     thin = 0.0
      thin = 1.0

c Define whether you want to eliminate ocean 'met stations'
c   (ocean_mask = 1.0), or not (ocean_mask = 0.0).
c NOTE: this does not work yet, but it could be used to provide
c   some guidance of what might be done if you are interested
c   in removing ocean points to save space in the MicroMet
c   file.
      ocean_mask = 0.0

c END USER INPUT.
c END USER INPUT.
c END USER INPUT.

c Provide the grid information. This includes the grid
c   increments in degrees.
      delta_lon_deg = 0.25
      delta_lat_deg = 0.25
      xlat_start = 10.0

c Find the approximate latitude of the center of the simulation
c   domain.
      if (thin.eq.0.0) then
        open (23,file=path1//'SM_dxdy_cornerll_proj_INPUTS.dat')
        read (23,*)
        read (23,*)
        read (23,*) dummy1,yll_lat
        read (23,*) dummy1,dummy2
        read (23,*) dummy1,yul_lat
        read (23,*) dummy1,dummy2
        center_lat = (yll_lat + yul_lat) / 2.0

c Adjust the longitude increment for the latitude of the domain
c   center.
        pi = 2.0 * acos(0.0)
        deg2rad = pi / 180.0
        print *
        write (*,96) delta_lon_deg
        delta_lon_deg = delta_lon_deg * cos(center_lat*deg2rad)
        write (*,97) delta_lon_deg
      endif

  96  format ('the original delta_lon_deg =',f8.4)
  97  format ('at this latitude, delta_lon_deg =',f8.4)

c Read in the SnowModel grid information from the topo-vege
c   processing file.
      open (31,file=path1//'SM_domain_config_info_OUTPUTS.dat')

      read (31,101) nx
      read (31,101) ny
      read (31,102) deltax
      read (31,102) deltay
      read (31,103) xmn
      read (31,103) ymn

  101 format (9x,i8)
  102 format (9x,f8.1)
  103 format (9x,f12.2)

      print *,nx
      print *,ny
      print *,deltax
      print *,deltay
      print *,xmn
      print *,ymn

c Name the output files that will contain the met analyses
c   points that you want to use in the SnowModel simulations.
c   The second file is just a duplicate that can/will be used
c   to plot the points in a figure.  The third file lists how
c   many reanalysis grid points you are going to be using in
c   your SnowModel simulation.
      open (81,file='ij_xy_topo.dat')
      open (82,file='../3_figs/pts.dat')
      open (83,file='npts.dat')
      open (84,file='../3_figs/SM_info.dat')

c Save some information required for plotting the results.
      write(84,*) nx
      write(84,*) ny
      write(84,*) deltax
      write(84,*) deltay
      write(84,*) xmn
      write(84,*) ymn

c Read in the elevation and projected coordinate data.
      do j=1,nny
        do i=1,nnx
          read (21,*) x(i,j),y(i,j)
          read (22,*) topo(i,j)
        enddo
      enddo

c Do the multi-points calculations.
      if (single.eq.0.0) then 

c Define the width of the met 'station' boundary outside the
c   SnowModel simulation domain, in analysis model lat-lon grid
c   increments.
        edge_scale = 1.1

c Distance in a degree of latitude, in meters.
        dist_per_deg = 111133.

c Undefined value.
        undef = -9999.0

c Define the border around the SnowModel simulation domain where
c   we want to include analysis met forcing grid points.
        extra_e = edge_scale * delta_lon_deg * dist_per_deg
        extra_w = edge_scale * delta_lon_deg * dist_per_deg
        extra_n = edge_scale * delta_lat_deg * dist_per_deg
        extra_s = edge_scale * delta_lat_deg * dist_per_deg

        xmin = xmn - extra_w
        ymin = ymn - extra_s
        xmax = xmn + nx * deltax + extra_e
        ymax = ymn + ny * deltay + extra_n

c Toss all of the met forcing points that are outside the
c   simulation domain, plus a small border around the
c   simulation domain boundaries.  This is done by setting
c   any analysis-grid topo value that we are not interested
c   in to undef.
        do j=1,nny
          do i=1,nnx
            if (x(i,j).lt.xmin) then
              topo(i,j) = undef
            endif
            if (x(i,j).gt.xmax) then
              topo(i,j) = undef
            endif
            if (y(i,j).lt.ymin) then
              topo(i,j) = undef
            endif
            if (y(i,j).gt.ymax) then
              topo(i,j) = undef
            endif
          enddo
        enddo

c If required, thin the analysis grid points as you go north in
c   the domain.  Note that it is critical for the x and y spacing
c   be approximately equal, or the MicroMet interpolations will
c   not be as accurate or representative as they could and should
c   be (for example, the 5 nearest met forcing grid points will
c   all be very close to each other in the E-W direction and will
c   not include any points to the N or S).
        if (thin.eq.1.0) then
          call thin_north (nnx,nny,topo,delta_lat_deg,undef,
     &      dist_per_deg,delta_lon_deg,xlat_start)
        endif

c Toss all ocean points that are within some distance of the
c   coast.
        if (ocean_mask.eq.1.0) then
          print *,'ocean_mask = 1.0 has not been implemented yet'
          stop
          call mask_ocean_points (nnx,nny,topo,x,y,undef,
     &      nx,ny,deltax,deltay,xmn,ymn)
        endif

c Save the valid grid point coordinates a text file.  Include
c   the elevation data.
        icount = 0
        do j=1,nny
          do i=1,nnx
            if( topo(i,j).ne.undef) then
              icount = icount + 1
              write(81,98) icount,i,j,x(i,j),y(i,j),topo(i,j)
              write(82,98) icount,i,j,x(i,j),y(i,j),topo(i,j)
            endif
          enddo
        enddo

c Save npts (the number of analyses grid points you are going
c   to use in the SnowModel simulation).
        write(83,99) icount
        print *,icount
        print *

      elseif (single.eq.1.0) then

        icount = 1

c Find the center coords of the SnowModel simulation domain.
        xmin = xmn
        ymin = ymn
        xmax = xmn + nx * deltax
        ymax = ymn + ny * deltay
        x1 = xmin + (xmax - xmin) / 2.0
        y1 = ymin + (ymax - ymin) / 2.0

c Find the minimum distance between the center of the SnowModel
c   simulation domain and the available (re)analyes grid cells.
        dist_old = 1.0e10
        do j=1,nny
          do i=1,nnx
            x2 = x(i,j)
            y2 = y(i,j)
            dist = sqrt((x2-x1)**2 + (y2-y1)**2)
            if (dist.lt.dist_old) then
              dist_old = dist
              ii = i
              jj = j
            endif
          enddo
        enddo

c Save the grid point.
        write(81,98) icount,ii,jj,x(ii,jj),y(ii,jj),topo(ii,jj)
        write(82,98) icount,ii,jj,x(ii,jj),y(ii,jj),topo(ii,jj)

c Save npts (the number of analyses grid points you are going
c   to use in the SnowModel simulation).
        write(83,99) icount

      endif

  98  format(3i5,2f16.1,f8.1)
  99  format(i10)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine thin_north (nnx,nny,topo,delta_lat_deg,undef,
     &  dist_per_deg,delta_lon_deg,xlat_start)

      implicit none

      integer nnx,nny,i,j
      real topo(nnx,nny),dx_meters(nny)
      integer ithin(nny)
      real dy_meters,delta_lat_deg,dist_per_deg,undef,xlat,
     &  delta_lon_deg,pi,deg2rad,xlat_start

      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0

c Calculate the latitudinal distance (m).
      dy_meters = delta_lat_deg * dist_per_deg

c Calculate how frequently the points need to be tossed.
      do j=1,nny
        xlat = xlat_start + (real(j) - 1.0) * delta_lat_deg
        dx_meters(j) = delta_lon_deg * dist_per_deg * cos(xlat*deg2rad)
        ithin(j) = nint(dy_meters / dx_meters(j))
        ithin(j) = max(1,ithin(j))
        if (xlat.eq.90.0) ithin(j) = ithin(j-1)
c       print *,j,dx_meters(j),dy_meters,ithin(j),xlat
      enddo

c Thin the points by setting the topo to undef.
      do j=1,nny
        do i=1,nnx
          if (mod(i,ithin(j)).ne.0) then
            topo(i,j) = undef
          endif
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mask_ocean_points (nnx,nny,topo,x,y,undef,
     &  nx,ny,deltax,deltay,xmn,ymn)

c NOTE: this has not been successfully implemented yet; but I
c   have done similar things to the code below for simulations
c   over Greenland, etc.

      implicit none

      integer nx,ny,i,j,nnx,nny,ismoother_loops,k
      real undef,deltax,deltay,xmn,ymn,topo_min

      real topo(nnx,nny),x(nnx,nny),y(nnx,nny)
      real topo1(nx,ny),vege1(nx,ny),vege2(nx,ny)

      integer ii(nnx,nny),jj(nnx,nny)

c BEGIN USER INPUT.
c BEGIN USER INPUT.
c BEGIN USER INPUT.

c This depends on grid resolution.
      ismoother_loops = 500

c Read in the veg data for the domain.
      open (51,file=
     &  '/data1/working/caribou/topo_veg/caribou_topo_veg_500m.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny)
      read (51,rec=2) ((vege1(i,j),i=1,nx),j=1,ny)

c END USER INPUT.
c END USER INPUT.
c END USER INPUT.

c Create a land-ocean elevation mask over the domain.  Save a
c   copy with the oceans set to undef.
      do j=1,ny
        do i=1,nx
          if (vege1(i,j).eq.24.0) then
            topo1(i,j) = 0.0
            vege2(i,j) = undef
          else
            topo1(i,j) = 1.0
            vege2(i,j) = 1.0
          endif
        enddo
      enddo

c Smooth this block mask so it puts non-zero topo out in the
c   oceans.
      do k=1,ismoother_loops
        print *,k
        call smoother9(nx,ny,topo1)
      enddo

c Use this information to remove analysis grid points from the
c   ocean part of the simulation domain.

c Convert the x and y locations to (ii,jj) locations on
c   the SnowModel (nx,ny) grid.
      do j=1,nny
        do i=1,nnx
          ii(i,j) = 1 + nint((x(i,j) - xmn) / deltax)
          jj(i,j) = 1 + nint((y(i,j) - ymn) / deltay)
        enddo 
      enddo 

c Now toss all of the points in the smoothed topo array that
c   are below some topo level.  You may have to plot the topo
c   result to see what this level should be.  And it does depend
c   on the value of ismoother_loops, grid resolution, etc.
      topo_min = 0.01

c For this special domain, just do it for the west and
c   north sides of the domain.
c     do j=101,nny
c       do i=1,nnx-100
      do j=1,nny
        do i=1,nnx
          if (topo(i,j).ne.undef) then
            if (topo1(ii(i,j),jj(i,j)).lt.topo_min) then
              topo(i,j) = undef
            endif
          endif
        enddo
      enddo

c Save the mask.
      open (41,file='coastal_mask.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny)
      write (41,rec=1) ((topo1(i,j),i=1,nx),j=1,ny)
      write (41,rec=2) ((vege2(i,j),i=1,nx),j=1,ny)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine smoother9(nx,ny,snow)

      implicit none

      integer i,j,nx,ny
      real snow(nx,ny)
      real snow_tmp(nx,ny)

c Performs a 9-point smoothing operation.

c The result at each grid point is a weighted average of the grid
c   point and the surrounding 8 points.  The center point receives
c   a weight of 1.0, the points at each side and above and below
c   receive a weight of 0.5, and corner points receive a weight of
c   0.3.  All points are multiplied by their weights and summed,
c   then divided by the total weight.

c Do the interior.
      do i=2,nx-1
        do j=2,ny-1
          snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) +
     &      snow(i,j+1) + snow(i-1,j) + snow(i+1,j)) + 0.3 *
     &      (snow(i-1,j-1) + snow(i+1,j+1) + snow(i-1,j+1) +
     &      snow(i+1,j-1))) / 4.2
        enddo
      enddo

c Do the sides.
      j = 1
      do i=2,nx-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j) +
     &    snow(i+1,j)) + 0.3 * (snow(i+1,j+1) + snow(i-1,j+1))) / 3.1
      enddo

      j = ny
      do i=2,nx-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j) +
     &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i-1,j-1))) / 3.1
      enddo

      i = 1
      do j=2,ny-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) +
     &    snow(i+1,j)) + 0.3 * (snow(i+1,j-1) + snow(i+1,j+1))) / 3.1
      enddo

      i = nx
      do j=2,ny-1
        snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i,j+1) +
     &    snow(i-1,j)) + 0.3 * (snow(i-1,j-1) + snow(i-1,j+1))) / 3.1
      enddo

c Do the corners.
      i = 1
      j = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i+1,j)) +
     &  0.3 * snow(i+1,j+1)) / 2.3

      i = nx
      j = 1
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j+1) + snow(i-1,j)) +
     &  0.3 * snow(i-1,j+1)) / 2.3

      i = 1
      j = ny
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i+1,j)) +
     &  0.3 * snow(i+1,j-1)) / 2.3

      i = nx
      j = ny
      snow_tmp(i,j) = (snow(i,j) + 0.5 * (snow(i,j-1) + snow(i-1,j)) +
     &  0.3 * snow(i-1,j-1)) / 2.3

c Return the smoothed array.
      do i=1,nx
        do j=1,ny
          snow(i,j) = snow_tmp(i,j)
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

