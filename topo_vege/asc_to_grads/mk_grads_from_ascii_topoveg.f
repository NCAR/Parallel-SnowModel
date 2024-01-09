c mk_grads_from_ascii_topoveg.f

c This program reads ARC/INFO ascii gridded vegetation and
c   topography data files and outputs a grads binary direct
c   access data file for use in SnowModel.

      parameter (iheader=6)
      parameter (nx=31,ny=31)

      real topo(nx,ny)
      real vege(nx,ny)

      open (21,file='SnowModel_test_dem.asc')
      open (22,file='SnowModel_test_veg.asc')

      do k=1,iheader
        read (21,*)
        read (22,*)
      enddo

c Read the data in as real numbers, and do the yrev.
      do j=ny,1,-1
        read (21,*) (topo(i,j),i=1,nx)
        read (22,*) (vege(i,j),i=1,nx)
      enddo

      open (31,file='SnowModel_test_topo_veg_200m.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny)

      write (31,rec=1) ((topo(i,j),i=1,nx),j=1,ny)
      write (31,rec=2) ((vege(i,j),i=1,nx),j=1,ny)

      end

