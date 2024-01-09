c A_write_coords_info.f

c This program is used if you have provided the exact coords
c   of the SnowModel domain in the SnowModel projection of
c   interest.  It just writes out the information required
c   so the reprojection script (#10) can run in the same way
c   that it is run if you provided approximate lat-lon
c   boundaries of the simulation domain you are interested
c   in.

      implicit none

      integer nx,ny

      real deltax,deltay,x_min,y_min,x_max,y_max,xmn,ymn

      character*250 proj_string

c Read in the SnowModel simulation domain corner coords.
      open (21,file='../../SM_dxdy_cornerll_proj_INPUTS.dat')

c Read past the first two lines of the INPUT file.
      read (21,*)
      read (21,*)

c Read the grid increment.
      read (21,*) deltax,deltay
      read (21,*)

c Read in the SnowModel simulation domain corner coords.  Note
c   that (x_min,y_min) is the lower-left corner, and (x_max,y_max)
c   is the upper-right corner.
      read (21,*) x_min,y_min
      read (21,*) x_max,y_max
      read (21,*)

c Read in the proj.4 string.
      read (21,96) proj_string

c Note that I assume here that your coords are evenly divisible
c   by deltax and deltay, and the correct nx and ny are
c   produced.
      nx = nint((x_max - x_min) / deltax)
      ny = nint((y_max - y_min) / deltay)

c     print *,nx,ny,x_min,x_max,y_min,y_max

c This is the center of the lower left grid cell in the SnowModel
c   simulation domain.  This goes in the .par file.
      xmn = x_min + 0.5 * deltax
      ymn = y_min + 0.5 * deltay

c SAVE ALL OF THIS INFORMATION IN SIMPLE TEXT FILES.

c Save this information as individual files with nothing else
c   in them.  This makes using the data in scripts easier in the
c   steps that follow.

c Create a file that contains the SnowModel domain information.
c   Save this in two places.
      open (31,file='outputs/3_nxny_dxdy_xmnymn.dat')
      write (31,101) nx
      write (31,102) ny
      write (31,103) deltax
      write (31,104) deltay
      write (31,105) xmn
      write (31,106) ymn
      
      open (32,file='../../SM_domain_config_info_OUTPUTS.dat')
      write (32,101) nx
      write (32,102) ny
      write (32,103) deltax
      write (32,104) deltay
      write (32,105) xmn
      write (32,106) ymn

  101 format ('nx =     ',i8)
  102 format ('ny =     ',i8)
  103 format ('deltax = ',f8.1)
  104 format ('deltay = ',f8.1)
  105 format ('xmn    = ',f12.2)
  106 format ('ymn    = ',f12.2)

      open (41,file='outputs/1_grid_increment.dat')
      write (41,97) deltax,deltay
   97 format (2f8.1)

      open (51,file='outputs/1_proj_string.dat')
      write (51,96) proj_string
   96 format (a250)

      open (61,file='outputs/3_cropping_coords.dat')
      write (61,99) x_min,y_min,x_max,y_max
   99 format (4f12.2)

c Write the corner coordinates of the SnowModel simulation
c   domain, in the projection of interest.
      open (71,file='outputs/2_corners_proj.dat')
      write (71,95) x_min,y_min
      write (71,95) x_max,y_min
      write (71,95) x_min,y_max
      write (71,95) x_max,y_max

c This is used to convert back to lat-lon to define the required
c   ll topo subdomain to extract.
      open (52,file='outputs/3_cropping_coords_tmp.dat')
      write (52,95) x_min,y_min
      write (52,95) x_max,y_min
      write (52,95) x_min,y_max
      write (52,95) x_max,y_max
   95 format (2f12.2)

      end

