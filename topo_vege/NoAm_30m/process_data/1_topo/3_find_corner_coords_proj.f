c 3_find_corner_coords_proj.f

c This program finds the corner coordinates of the SnowModel
c   simulation domain, based on the lat-lon coordinate
c   information you provided in /1_topo/1_define_corners/
c   in ll_corners.dat.

c This requires deltax and deltay values.
      open (41,file='outputs/1_grid_increment.dat')
      read (41,*) deltax,deltay

c Read in the approximate corner coordinates of the SnowModel
c   simulation domain, in the projection of interest.
      open (21,file='outputs/2_corners_proj.dat')

      read (21,*) x_ll,y_ll
      read (21,*) x_lr,y_lr
      read (21,*) x_ul,y_ul
      read (21,*) x_ur,y_ur

c Find the max extent of these coords.
      x_min = min(x_ll,x_ul)
      x_max = max(x_lr,x_ur)

      y_min = min(y_ll,y_lr)
      y_max = max(y_ul,y_ur)

c     print *,x_min,x_max,y_min,y_max

c Find the corner coords that are integer-divisable by the
c   grid increment.  If they are, then the difference should
c   be too (giving reasonable nx and ny values).
c     xx_min = deltax * real(nint(x_min / deltax))
c     xx_max = deltax * real(nint(x_max / deltax))
c     yy_min = deltay * real(nint(y_min / deltay))
c     yy_max = deltay * real(nint(y_max / deltay))

c This version always makes things a little bigger (it always
c   rounds "up", never "down"; it makes the domain bigger,
c   not smaller).  Using this also means you have to compile
c   with the gfortran compiler.
      xx_min = deltax * real(floor(x_min / deltax))
      xx_max = deltax * real(ceiling(x_max / deltax))
      yy_min = deltay * real(floor(y_min / deltay))
      yy_max = deltay * real(ceiling(y_max / deltay))

      nx = nint((xx_max - xx_min) / deltax)
      ny = nint((yy_max - yy_min) / deltay)

c     print *,nx,ny,xx_min,xx_max,yy_min,yy_max

c This is the center of the lower left grid cell in the SnowModel
c   simulation domain.  This goes in the .par file.
      xmn = xx_min + 0.5 * deltax
      ymn = yy_min + 0.5 * deltay

c Convert the min-max values to values for the edges of the
c   grid cells, like in GIS conventions (this is what is
c   assumed in the gdal calcuations the next steps).
c     xxx_min = xx_min - 0.5 * deltax
c     xxx_max = xx_max + 0.5 * deltax
c     yyy_min = yy_min - 0.5 * deltay
c     yyy_max = yy_max + 0.5 * deltay
      xxx_min = xx_min
      xxx_max = xx_max
      yyy_min = yy_min
      yyy_max = yy_max

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

      open (51,file='outputs/3_cropping_coords.dat')
      write (51,99) xxx_min,yyy_min,xxx_max,yyy_max

c This is used to convert back to lat-lon to define the required
c   ll topo subdomain to extract.
      open (52,file='outputs/3_cropping_coords_tmp.dat')
      write (52,98) xxx_min,yyy_min
      write (52,98) xxx_max,yyy_min
      write (52,98) xxx_min,yyy_max
      write (52,98) xxx_max,yyy_max

   98 format (2f12.2)
   99 format (4f12.2)

c Output the coordinates of the center of the northern boundary.
      open (53,file='outputs/3_northern_center_coords.dat')
      write (53,98) (xxx_max-xxx_min)/2.0+xxx_min,yyy_max

      end

