c 1_find_subdomain_ll.f

c This program finds the lat-lon coordinates of a North American
c   subdomain that contains the SnowModel simulation domain
c   within it, and that is made up of 1-degree by 1-degree data
c   blocks.

c This must be compiled with gfortran, pgf90, or pgf95 (pgf77
c   won't work because of the floor and ceiling calls).

      character*250 proj_string

c Read in the approximate corner lat-lon coordinates of the
c   SnowModel simulation domain.
      open (21,file='../../SM_dxdy_cornerll_proj_INPUTS.dat')

c Read past the first two lines of the INPUT file.
      read (21,*)
      read (21,*)

c Read the grid increment.
      read (21,*) deltax,deltay
      read (21,*)

c Read in the SnowModel simulation domain corner coords.
      read (21,*) xll_lon,yll_lat
      read (21,*) xur_lon,yur_lat
      read (21,*)

c Read in the proj.4 string.
      read (21,96) proj_string

c Define the max extent of these coords.
      x_min = xll_lon
      x_max = xur_lon

      y_min = yll_lat
      y_max = yur_lat

c     print *,x_min,x_max,y_min,y_max

c Find the integer lat-lon coords that this domain is contained
c   within.
      lon_min = floor(x_min) - 1
      lon_max = ceiling(x_max) + 1

      lat_min = floor(y_min) - 1
      lat_max = ceiling(y_max) + 1

c     print *,lon_min,lon_max,lat_min,lat_max

      open (31,file='outputs/1_subdomain_ll_extent.dat')
      open (32,file='outputs/1_subdomain_ll_corners.dat')

      write (31,98) lon_min,lat_min,lon_max,lat_max

      write (32,99) lon_min,lat_min
      write (32,99) lon_max,lat_min
      write (32,99) lon_min,lat_max
      write (32,99) lon_max,lat_max

c Having these as individual files with nothing else in then
c   makes using the data in scripts easier in the following
c   steps.
      open (41,file='outputs/1_grid_increment.dat')
      write (41,97) deltax,deltay

      open (51,file='outputs/1_proj_string.dat')
      write (51,96) proj_string

      open (61,file='outputs/1_ll_sm_corners.dat')
      write (61,95) xll_lon,yll_lat
      write (61,95) xur_lon,yll_lat
      write (61,95) xll_lon,yur_lat
      write (61,95) xur_lon,yur_lat

   95 format (2f12.6)
   96 format (a250)
   97 format (2f8.1)
   98 format (4i8)
   99 format (2i8)

      end

