c B_find_subdomain_ll.f

c This program finds the lat-lon coordinates of a North American
c   subdomain that contains the SnowModel simulation domain
c   within it, and that is made up of 1-degree by 1-degree data
c   blocks.

c This must be compiled with gfortran, pgf90, or pgf95 (pgf77
c   won't work because of the floor and ceiling calls).

c Read in the approximate corner lat-lon coordinates of the
c   SnowModel simulation domain.
      open (21,file='outputs/4_subdomain_ll_corners_tmp.dat')

c Read in the SnowModel simulation domain corner coords.
      read (21,*) xll_lon,yll_lat
      read (21,*) xlr_lon,ylr_lat
      read (21,*) xul_lon,yul_lat
      read (21,*) xur_lon,yur_lat

c Find the max extent of these coords.
      x_min = min(xll_lon,xul_lon)
      x_max = max(xlr_lon,xur_lon)

      y_min = min(yll_lat,ylr_lat)
      y_max = max(yul_lat,yur_lat)

c     print *,x_min,x_max,y_min,y_max

c Find the integer lat-lon coords that this domain is contained
c   within.
      lon_min = floor(x_min) - 1
      lon_max = ceiling(x_max) + 1

      lat_min = floor(y_min) - 1
      lat_max = ceiling(y_max) + 1

c     print *,lon_min,lon_max,lat_min,lat_max

      open (31,file='outputs/5_subdomain_ll_extent_final.dat')
      open (32,file='outputs/5_subdomain_ll_corners_final.dat')

      write (31,98) lon_min,lat_min,lon_max,lat_max

      write (32,99) lon_min,lat_min
      write (32,99) lon_max,lat_min
      write (32,99) lon_min,lat_max
      write (32,99) lon_max,lat_max

   98 format (4i8)
   99 format (2i8)

      end

