c 3_mk_lat_lon.f

c This is the SnowModel grid.
      parameter (nx=310,ny=200)

      real glon(nx,ny)
      real glat(nx,ny)

c Open the output files.
      open (unit=51,file='grid_lon.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny,
     &  status='replace')

      open (unit=52,file='grid_lat.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny,
     &  status='replace')

c Read in the lat-lon data for the SnowModel grid.
      open (31,file='ll_coords.txt')
      do j=1,ny
        do i=1,nx
          read (31,*) glon(i,j),glat(i,j)
        enddo
      enddo

      write (51,rec=1) ((glon(i,j),i=1,nx),j=1,ny)
      write (52,rec=1) ((glat(i,j),i=1,nx),j=1,ny)

      end

