c 1_merge_topo_vege.f

c This program merges the topo and vege .asc files into one .gdat
c   file for use in SnowModel.

      parameter (nx_max=60000,ny_max=35000)

      real topo(nx_max,ny_max)
      real veg1(nx_max,ny_max)
      real veg2(nx_max,ny_max)

      parameter (ntypes=25)
      integer num_in_class(ntypes)
      integer i_length,trailing_blanks

      character*39 fname_veg
      character*250 fname_in

c This is going to be used to compare against the actual input
c   file that you are using, to make sure the conversions from
c   this file to the SnowModel landcover classes are correct.
c   If you are using a different input file, you will have to
c   search and change the two occurances of '39'.  39 is the
c   length of "fname_veg".
      data fname_veg /'NA_NALCMS_2015_LC_30m_LAEA_mmu5pix_.tif'/

c First check to see if you are working with the NALCMS dataset.
      open (81,file='../../input_files_info/vege_filename.dat')
      read (81,98) fname_in
      i_length = 250 - trailing_blanks(fname_in)
      i_start = i_length - 39 + 1

c     print *,fname_in(i_start:)
c     print *,fname_veg

      if (fname_in(i_start:).ne.fname_veg) then
        print *,'*** The program has stopped because the'
        print *,'vege class number conversions are not/may'
        print *,'not be correct for this dataset.  See:'
        print *,'/3_merge_topo_vege/merge_topo_vege.f. ***'
        print *
        stop
      endif

   98 format (a250)

      undef = -9999.0

c Read in nx and ny.
      open (21,file=
     &  '../1_topo/outputs/3_nxny_dxdy_xmnymn.dat')
      read (21,99) nx
      read (21,99) ny

  99  format (9x,i8)

c Read in the .flt topo data.
      open (31,file='../1_topo/outputs/10_topo_proj_cropped.flt',
     &  form='unformatted',access='direct',recl=4*nx)

c Read the data in as real numbers, and do the yrev.
      do j=1,ny
        jj = ny + 1 - j
        read (31,rec=j) (topo(i,jj),i=1,nx)
      enddo

c For Alaska, the original NED topography dataset has an
c   approximately 150-m wavelength processing feature in x
c   and y.  SnowModel simulations (see the Dall sheep project:
c   dallsheep_snowmodel_0 run) suggest 10 iterations of the
c   smoother is required to correct the problem.  It looks to
c   me like this kind of smoothing is required for SnowModel
c   simulations at resolutions below ~150-m grid increment.
c     do k=1,10
c       call smoother9(nx,ny,topo)
c     enddo

c Check for any negative topography values.  If there are any
c   negative topo values, clip them to zero and write them to
c   a file.
      open (91,file='negative_topo_values.dat')

      do j=1,ny
        do i=1,nx
          if (topo(i,j).eq.undef) then
            topo(i,j) = 0.0
          endif
          if (topo(i,j).lt.0.0) then
c           print *,'neg topo values found',i,j,topo(i,j)
            write (91,*) i,j,topo(i,j)
            topo(i,j) = 0.0
          endif
        enddo
      enddo

c Read in the .flt vege data.
      open (41,file='../2_vege/outputs/1_vege_proj_cropped.flt',
     &  form='unformatted',access='direct',recl=4*nx)

c Read the data in as real numbers, and do the yrev.
      do j=1,ny
        jj = ny + 1 - j
        read (41,rec=j) (veg1(i,jj),i=1,nx)
      enddo

c Convert from the NALCMS classes to SnowModel classes.  Check
c   for any undefined land cover values.  If you find any, assume
c   they are ocean and write a list of those to a file so you can
c   look at them if you want.
c     open (92,file='undef_vege_values.dat')

      do j=1,ny
        do i=1,nx
          veg2(i,j) = undef

c Ocean = 0.
          if (veg1(i,j).eq.0.0) veg2(i,j) = 24.0
          if (veg1(i,j).eq.1.0) veg2(i,j) = 1.0
          if (veg1(i,j).eq.2.0) veg2(i,j) = 1.0
          if (veg1(i,j).eq.3.0) veg2(i,j) = 2.0
          if (veg1(i,j).eq.4.0) veg2(i,j) = 2.0
          if (veg1(i,j).eq.5.0) veg2(i,j) = 2.0
          if (veg1(i,j).eq.6.0) veg2(i,j) = 3.0
          if (veg1(i,j).eq.7.0) veg2(i,j) = 8.0
          if (veg1(i,j).eq.8.0) veg2(i,j) = 6.0
          if (veg1(i,j).eq.9.0) veg2(i,j) = 12.0
          if (veg1(i,j).eq.10.0) veg2(i,j) = 12.0
          if (veg1(i,j).eq.11.0) veg2(i,j) = 10.0
          if (veg1(i,j).eq.12.0) veg2(i,j) = 14.0
          if (veg1(i,j).eq.13.0) veg2(i,j) = 18.0
          if (veg1(i,j).eq.14.0) veg2(i,j) = 17.0
          if (veg1(i,j).eq.15.0) veg2(i,j) = 22.0
          if (veg1(i,j).eq.16.0) veg2(i,j) = 18.0
          if (veg1(i,j).eq.17.0) veg2(i,j) = 21.0
          if (veg1(i,j).eq.18.0) veg2(i,j) = 19.0
          if (veg1(i,j).eq.19.0) veg2(i,j) = 20.0

          if (veg2(i,j).eq.undef) then
c           print *,'veg undef found at ',i,j,veg1(i,j)
c           print *,'  setting to ocean '
c           write (92,*) i,j,veg1(i,j)
            veg2(i,j) = 24.0
c           stop
          endif
        enddo
      enddo

c Count the number of grid cells in each number-class that is used.
      do k=1,ntypes
         num_in_class(k) = 0
      enddo

      do j=1,ny
        do i=1,nx
          idat = int(veg2(i,j))
          num_in_class(idat) = num_in_class(idat) + 1
        enddo
      enddo

      open (51,file='classes.txt')
      do k=1,ntypes
        write (51,*) k,num_in_class(k)
      enddo

c Save the data in a GrADS file.
      open (71,file='../../topo_vege.gdat',
     &  form='unformatted',access='direct',recl=4*nx,
     &  status='replace')

      do j=1,ny
        write (71,rec=j) (topo(i,j),i=1,nx)
      enddo

      do j=1,ny
        write (71,rec=j+ny) (veg2(i,j),i=1,nx)
      enddo

c Create a grads .ctl file for this topo_vege.gdat file.
      call mk_topo_vege_ctl (nx,ny)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mk_topo_vege_ctl (nx,ny)

      implicit none

      integer nx,ny

      open (61,file='../../topo_vege.ctl')

      write (61,51)
      write (61,52)
      write (61,53)
      write (61,54) nx
      write (61,55) ny
      write (61,56)
      write (61,57)
      write (61,58)
      write (61,59)
      write (61,60)
      write (61,61)

      close (61)

   51 format ('DSET ^topo_vege.gdat')
   52 format ('TITLE topography and landcover for SnowModel')
   53 format ('UNDEF -9999.0')
   54 format ('XDEF ',i5,' LINEAR 1.0 1.0')
   55 format ('YDEF ',i5,' LINEAR 1.0 1.0')
   56 format ('ZDEF 1 LINEAR 1 1')
   57 format ('TDEF 1 LINEAR 1jan9999 1yr')
   58 format ('VARS 2')
   59 format ('topo 0 0 topography (m)')
   60 format ('vege 0 0 land cover (SnowModel classes)')
   61 format ('ENDVARS')

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

      integer function trailing_blanks(input_string)

      implicit none

      integer k
      character*(*) input_string

c Count the number of blanks following the last non-blank
c   character.
      trailing_blanks = 0
      do k=len(input_string),1,-1
        if (input_string(k:k).eq.' ') then
          trailing_blanks = trailing_blanks + 1
        else
          return
        endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

