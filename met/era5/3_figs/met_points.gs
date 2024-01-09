 'reinit'
 'set display color white'
 'c'

 'set mproj scaled'
 'set grads off'
 'set mpdraw off'

 'open ../../../topo_vege/NoAm_30m/topo_vege.ctl'

**********************************

* define the SnowModel domain
  fname=SM_info.dat

  ret = read(fname)
  rec = sublin(ret,2)
  nx = subwrd(rec,1)

  ret = read(fname)
  rec = sublin(ret,2)
  ny = subwrd(rec,1)

  ret = read(fname)
  rec = sublin(ret,2)
  dx = subwrd(rec,1)

  ret = read(fname)
  rec = sublin(ret,2)
  dy = subwrd(rec,1)

  ret = read(fname)
  rec = sublin(ret,2)
  xmn = subwrd(rec,1)

  ret = read(fname)
  rec = sublin(ret,2)
  ymn = subwrd(rec,1)

* say nx
* say ny
* say dx
* say dy
* say xmn
* say ymn

  fname = pts.dat
  set_parea(nx,ny,fname,dx,dy,xmn,ymn)

 'set mpdraw off'
 'set mproj scaled'
 'set gxout grfill'
 'set grads off'

 'set xlab off'
 'set ylab off'

 'd topo'

  points(fname,dx,dy,xmn,ymn)

************************************
************************************

 'gprint fig_files/met_points'

************************************
************************************
************************************
************************************

  function points(fname,dx,dy,xmn,ymn)

 'set line 1 1 6'
while (1)
  ret = read(fname)
  rc = sublin(ret,1)
  if (rc>0) 
    if (rc!=2) 
      say 'File I/O Error'
      return
    endif
    break
  endif
  rec = sublin(ret,2)
  xstn = subwrd(rec,4)
  ystn = subwrd(rec,5)

* convert these 'meters' coords to i,j to be compatible with
*   the i,j .ctl file
  ii=1+(xstn-xmn)/dx
  jj=1+(ystn-ymn)/dy

* say ii
* say jj

  'q ll2xy 'ii' 'jj
  x1=subwrd(result,1)
  y1=subwrd(result,2)

  'set line 1'
  'draw mark 3 'x1' 'y1' 0.10'
endwhile

  rc=close(fname)

  return

************************************
************************************

     function set_parea(nx,ny,fname,dx,dy,xmn,ymn)

* Adjust parea to make the aspect ratio correct.
*   This should work regardles of whether you
*   are plotting in landscape or portrait.
*    'set mproj scaled'
*     set_parea(nx,ny)

* find the min and max extent of the met points
     min_max_points(fname,dx,dy,xmn,ymn)

*    say nx' 'ny
*    say _dx_pts' '_dy_pts

     scale_x = nx / _dx_pts
     scale_y = ny / _dy_pts
     if (scale_x < scale_y)
       scale = scale_x
     else
       scale = scale_y
     endif
*    say scale_x' 'scale_y' 'scale

     min_boundary = 1.0
    'query gxinfo'
     rec2 = sublin(result,2)
     gxlolim = 0.0 + min_boundary
     gxhilim = subwrd(rec2,4) - min_boundary
     gylolim = 0.0 + min_boundary
     gyhilim = subwrd(rec2,6) - min_boundary

     gxmax = scale * (gxhilim - gxlolim)
     gymax = scale * (gyhilim - gylolim)
     gxcent = (gxhilim + gxlolim) / 2
     gycent = (gyhilim + gylolim) / 2
     gxlolim = gxcent - gxmax / 2
     gxhilim = gxcent + gxmax / 2
     gylolim = gycent - gymax / 2
     gyhilim = gycent + gymax / 2

     eps = 0.0001
     ratio = nx / ny + eps

     if (ratio>1)
* x is the biggest.
       xlo = gxlolim
       xhi = gxhilim
       ylo = gycent - gxmax / ratio / 2
       yhi = gycent + gxmax / ratio / 2
       if (ylo<gylolim)
         xlo = gxcent - gymax * ratio / 2
         xhi = gxcent + gymax * ratio / 2
         ylo = gylolim
         yhi = gyhilim
       endif
     else
* y is the biggest.
       xlo = gxcent - gymax * ratio / 2
       xhi = gxcent + gymax * ratio / 2
       ylo = gylolim
       yhi = gyhilim
     endif

    'set parea 'xlo' 'xhi' 'ylo' 'yhi

     return

************************************
************************************

  function min_max_points(fname,dx,dy,xmn,ymn)

  xmin_old = 10000
  xmax_old = 0
  ymin_old = 10000
  ymax_old = 0

* say

while (1)
  ret = read(fname)
  rc = sublin(ret,1)
  if (rc>0) 
    if (rc!=2) 
      say 'File I/O Error'
      return
    endif
    break
  endif
  rec = sublin(ret,2)
  xstn = subwrd(rec,4)
  ystn = subwrd(rec,5)

* convert these 'meters' coords to i,j to be compatible with
*   the i,j .ctl file
* ii=1+(xstn-xmn)/dx
* jj=1+(ystn-ymn)/dy
  x=1+(xstn-xmn)/dx
  y=1+(ystn-ymn)/dy
* say x' 'y

  if (x < xmin_old) ; xmin = x ; xmin_old = xmin ; endif
  if (x > xmax_old) ; xmax = x ; xmax_old = xmax ; endif
  if (y < ymin_old) ; ymin = y ; ymin_old = ymin ; endif
  if (y > ymax_old) ; ymax = y ; ymax_old = ymax ; endif

endwhile

* say
* say xmin' 'xmax
* say ymin' 'ymax
* say

  _dx_pts = xmax - xmin
  _dy_pts = ymax - ymin

  rc=close(fname)

  return

************************************
************************************
