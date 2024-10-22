c mk_3hrly.f

c Read in the hourly data and write out 3-hrly data.  This
c   is just done to create a workable 3-hrly met file.  It
c   is not really the way I normally do things.

c     implicit none

      integer nlines,nstns,k,n,j

      parameter (nlines=21288)
      parameter (nstns=1)

      integer iyr(nlines,nstns)
      integer imo(nlines,nstns)
      integer idy(nlines,nstns)
      real xhr(nlines,nstns)

      real ta(nlines,nstns)
      real rh(nlines,nstns)
      real ws(nlines,nstns)
      real wd(nlines,nstns)
      real pp(nlines,nstns)

      real northing(nlines,nstns)
      real easting(nlines,nstns)
      real elev(nlines,nstns)
      integer id(nlines,nstns)

      real undef

      undef = -9999.0

      pi = 2.0 * acos(0.0)
      degtrad = pi / 180.0
      radtdeg = 180.0 / pi

      nhrs_in_3hrs = 3

c Input files.
      open (41,file='met_station_408_hrly.dat')

c Read the data.
      do j=1,nstns
        do n=1,7656
          read (40+j,*) iyr(n,j),imo(n,j),idy(n,j),xhr(n,j),id(n,j),
     &      easting(n,j),northing(n,j),elev(n,j),
     &      ta(n,j),rh(n,j),ws(n,j),wd(n,j),pp(n,j)
        enddo
        close (40+j)
      enddo

c Output files.
      open (41,file='met_station_408_3hrly.dat')

c Create the daily data.
      do j=1,nstns
        do k=1,319*8

c Initialize the summing arrays.
          temp = 0.0
          relh = 0.0
          wspd = 0.0
          prec = 0.0
          uwnd = 0.0
          vwnd = 0.0

          miss_ta = 0
          miss_rh = 0
          miss_ws = 0
          miss_wd = 0
          miss_pp = 0

c Summing loop.
          do i=1,nhrs_in_3hrs
            kk = (k-1) * nhrs_in_3hrs + i

            if (ta(kk,j).eq.undef) miss_ta = 1
            if (rh(kk,j).eq.undef) miss_rh = 1
            if (ws(kk,j).eq.undef) miss_ws = 1
            if (wd(kk,j).eq.undef) miss_wd = 1
            if (pp(kk,j).eq.undef) miss_pp = 1

            temp = temp + ta(kk,j)
            relh = relh + rh(kk,j)
            wspd = wspd + ws(kk,j)
            prec = prec + pp(kk,j)

            uwnd = uwnd - sin(wd(kk,j)*degtrad)
            vwnd = vwnd - cos(wd(kk,j)*degtrad)
          enddo

c Calculate the daily averages.
          wspd = wspd / real(nhrs_in_3hrs)
          temp = temp / real(nhrs_in_3hrs)
          relh = relh / real(nhrs_in_3hrs)

c Leave the precipitation as a sum of the hourly values.
c         prec = prec

c Deal with the wind direction.  Some compilers don't allow both
c   variables to be zero in the atan2 calculation.
          if (abs(uwnd).lt.1e-10) uwnd = 1e-10
          wdir = radtdeg * atan2(uwnd,vwnd)
          if (wdir.ge.180.0) then
            wdir = wdir - 180.0
          else
            wdir = wdir + 180.0
          endif

          if (miss_ta.eq.1) temp = undef
          if (miss_rh.eq.1) relh = undef
          if (miss_ws.eq.1) wspd = undef
          if (miss_wd.eq.1) wdir = undef
          if (miss_pp.eq.1) prec = undef

c Deal with the time stamp.
c         xhr1 = xhr(kk-1,j)
c         idy1 = idy(kk-1,j)
c         imo1 = imo(kk-1,j)
c         iyr1 = iyr(kk-1,j)
          xhr1 = xhr(kk,j)
          idy1 = idy(kk,j)
          imo1 = imo(kk,j)
          iyr1 = iyr(kk,j)

c Increase the precipitation a little.
          prec = 1.4 * prec

c Clean up the air temperatures so we have a well defined
c   accumulation season.
          if (imo1.eq.10) temp = temp+5.0
          if (imo1.eq.10) temp = max(temp,5.0)
          if (imo1.eq.11) temp = temp-5.0
          if (imo1.eq.12) temp = temp-5.0
          if (imo1.eq.1) temp = temp-5.0
          if (imo1.eq.2) temp = temp-5.0
          if (imo1.eq.3) temp = temp-5.0
          if (imo1.eq.4) temp = temp-5.0

c Reduce the wind speed a little to make the data
c   assimilation more exact.
          wspd = 0.5 * wspd

c Save the data in the MicroMet required format.
          write (40+j,80) iyr1,imo1,idy1,xhr1,
     &      id(kk,j),easting(kk,j),northing(kk,j),elev(kk,j),
     &      temp,relh,wspd,wdir,prec
        enddo
        close (40+j)
      enddo

  80  format(i5,i3,i3,f6.2,i6,2f12.1,f8.1,5f9.2)

      end

