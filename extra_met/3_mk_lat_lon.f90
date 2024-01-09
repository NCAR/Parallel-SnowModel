PROGRAM mk_lat_lon_3
  implicit none

! c This is the SnowModel grid.
      integer, parameter :: nx=6916,ny=9107

      real glon(nx,ny)
      real glat(nx,ny)
      integer :: i,j

! c Open the output files.
      open (unit=51,file='grid_lon.gdat',&
     &  form='unformatted',access='direct',recl=4*nx,&
     &  status='replace')

      open (unit=52,file='grid_lat.gdat',&
     &  form='unformatted',access='direct',recl=4*nx,&
     &  status='replace')


! c Read in the lat-lon data for the SnowModel grid.
      open (31,file='ll_coords.txt')
      do j=1,ny
        do i=1,nx
          read (31,*) glon(i,j),glat(i,j)
        enddo
      enddo

      print *, 'Finished reading'

      do j=1,ny
        write (51,rec=j) (glon(i,j),i=1,nx)
      enddo

      print *,'Finished writing glon'

      do j=1,ny
        write (52,rec=j) (glat(i,j),i=1,nx)
      enddo

      print *, 'Finished writing glat'

      close(51)
      close(52)
      close(31)

end program
