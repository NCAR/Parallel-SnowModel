c 6_mk_ll_tile_list.f

c This program identifies the 1-degree by 1-degree blocks of
c   topo data that cover the SnowModel simulation domain you
c   are interested in.

c Because this was set up to work for my North American topo
c   dataset, it assumes the integer degrees span -180 to 0
c   for W to E, and 0 to 90 for S to N.  If you are working
c   in the Eastern or Southern Hemispheres, you may have to
c   make the code a little more general (see the "mx_start =
c   abs(lon_right) + 1, etc." calculations below).

c These lat-lon coordinate bounds are read in from the file:
c   /outputs/subdomain_ll_extent_final.dat

      implicit none

      integer lon_left,lon_right,lat_bottom,lat_top,
     &  mx_start,mx_end,my_start,my_end,mx,my,mmx,mmy

      character*200 fname
      character*2 lat
      character*3 lon

c Read in the lat-lon coordinate bounds for this subdomain.
      open (21,file='outputs/5_subdomain_ll_extent_final.dat')
      read (21,*) lon_left,lat_bottom,lon_right,lat_top

      print *
      print *,'lon_left =  ',lon_left
      print *,'lon_right = ',lon_right
      print *,'lat_bottom =',lat_bottom
      print *,'lat_top =   ',lat_top

      mx_start = abs(lon_right) + 1
      mx_end = abs(lon_left)
      my_start = lat_bottom + 1
      my_end = lat_top

c     print *, mx_start,mx_end,my_start,my_end

c mx = number of 1 degree blocks e-w.
c my = number of 1 degree blocks n-s.
      mx = mx_end - mx_start + 1
      my = my_end - my_start + 1

c Calculate the array dimensions for the coarse grid.
      print *
      print *,'mx =, my =',mx,my

c Open the output file.
      open (31,file='outputs/7_ll_file_list.dat')

c Loop through each 1 degree block.
      do mmy=my_start,my_end
        do mmx=mx_start,mx_end

c Find the appropriate file name.
          write(lon,'(i3.3)') mmx
          write(lat,'(i2.2)') mmy

          fname ='n'//lat//'w'//lon

c         write (*,92) -mmx,mmy,fname

          write (31,91) fname

        enddo
      enddo

      close (31)

  91  format (a7)
c 92  format ('working on lon, lat (upper left corner):',2i6,4x,a7)

      end

