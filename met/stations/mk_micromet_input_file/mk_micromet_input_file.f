c mk_micromet_input_file.f

c This program reads the individual met station data files, merges
c   them into one file, and puts the station count before each
c   hours worth of data.

c It assumes that you have data increments of hourly to daily, and
c   defines enough space for 1000 stations with one year worth of
c   hourly data.  maxlines must be changed in three places if this
c   is exceeded.

      parameter (maxlines=1000000)
      parameter (maxntimes=10000)
      parameter (nfiles=2)
      parameter (nvars=13)

      integer iyr(maxlines)
      integer imo(maxlines)
      integer idy(maxlines)
      integer istn(maxlines)

      real xhr(maxlines)
      real x(maxlines)
      real y(maxlines)
      real elev(maxlines)
      real tair(maxlines)
      real rh(maxlines)
      real spd(maxlines)
      real dir(maxlines)
      real prec(maxlines)

      real var_1(maxlines,nvars)

      integer istncount(maxntimes)

      character*80 fnamein(nfiles),fnameout

      data fnamein/'met_station_101.dat',
     &             'met_station_408.dat'/

      data fnameout/'SnowModel_test_met.dat'/

c Open the output file.
      open (11,file=fnameout)

c READ ALL OF THE MET FILES INTO A SINGLE TIME-VARIABLE ARRAY.

c   Note, the data are assumed to be in the following order:
c    iyr(k),imo(k),idy(k),xhr(k),istn(k),x(k),y(k),elev(k),
c    tair(k),rh(k),spd(k),dir(k),prec(k)
      irecs = 0
      do nfile=1,nfiles
        open (21,file=fnamein(nfile))

        do k=irecs+1,maxlines
          read (21,*,end=99) (var_1(k,i),i=1,nvars)
          irecs = irecs + 1
        enddo
  99    continue

        close (21)
      enddo

c SORT THIS TIME-VARIABLE ARRAY ON THE YEAR, MONTH, DAY, HOUR, AND
c   STATION ID.
      call sort_met(irecs,var_1)

c Extract the sorted variables into their named arrays.
      do k=1,irecs
        iyr(k) = nint(var_1(k,1))
        imo(k) = nint(var_1(k,2))
        idy(k) = nint(var_1(k,3))
        xhr(k) = var_1(k,4)
        istn(k) = nint(var_1(k,5))
        x(k) = var_1(k,6)
        y(k) = var_1(k,7)
        elev(k) = var_1(k,8)
        tair(k) = var_1(k,9)
        rh(k) = var_1(k,10)
        spd(k) = var_1(k,11)
        dir(k) = var_1(k,12)
        prec(k) = var_1(k,13)
      enddo

c TAKE THE SORTED MET-DATA ARRAY, AND WRITE A STATION COUNT BEFORE
c   EACH TIME SLICE.

c This allows the program to processes the last time slice.
      xhr(irecs+1) = xhr(irecs) + 1.0

c Loop through all of the xhr's and build an array of the number of
c   stations for each hour.
      ntimes = 0
      xhr_test = xhr(1)
      idy_test = idy(1)
      icount = 0
      do k=1,irecs+1
        if (xhr(k).eq.xhr_test .and. idy(k).eq.idy_test) then
          icount = icount + 1
        else
          ntimes = ntimes + 1
          istncount(ntimes) = icount
          xhr_test = xhr(k)
          idy_test = idy(k)
          icount = 1
        endif

      enddo

c Loop through all of the data, writing a station count before each
c   hour's data write.
      k = 1
      do kk=1,ntimes

        write (11,980) istncount(kk)

        do kkk=1,istncount(kk)
          write (11,990) iyr(k),imo(k),idy(k),xhr(k),istn(k),
     &      x(k),y(k),elev(k),tair(k),rh(k),spd(k),dir(k),prec(k)
          k = k + 1
        enddo
      enddo

  980 format(i5)
  990 format(i5,i3,i3,f6.2,i6,2f10.1,f8.1,5f9.2)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sort_met(irecs,var_1)

      implicit none

      integer irecs,maxlines,nvars,k,kkk,i,icount,icount_1,
     &  icount_2,ii,isegment

      parameter (maxlines=1000000)
      parameter (nvars=13)

      real var_1(maxlines,nvars)
      real var_2(maxlines,nvars)
      real sort_var(maxlines)

      integer istart(maxlines)
      integer iend(maxlines)
      integer index(maxlines)
      integer irecs_in_segment(maxlines)

c Sort on the year.
c Extract the sort variable and define the index array.
      do k=1,irecs
        sort_var(k) = var_1(k,1)
        index(k) = k
      enddo

c Get the new order (index).
      call isort(sort_var,index,irecs)

c Re-order all of the variables.
      do k=1,irecs
        do i=1,nvars
          var_2(k,i) = var_1(index(k),i)
        enddo
      enddo

c Re-set the variable array.
      do k=1,irecs
        do i=1,nvars
          var_1(k,i) = var_2(k,i)
        enddo
      enddo

c Initialize the start and end indexing arrays.
      do k=1,irecs
        istart(k) = 0
        iend(k) = 0
      enddo

c Sort on the month, day, hour, and station.
      do ii=2,5

c First loop through and find the start and end of the common
c   date segments.  Identify the new segments, while keeping the
c   old segments.
        istart(1) = 1
        iend(irecs) = irecs
        do k=2,irecs
          if (var_1(k-1,ii-1).ne.var_1(k,ii-1)) then
            istart(k) = k
            iend(k-1) = k - 1
          endif
        enddo

c Extract the number of segments and the segment lengths.
        isegment = 0
        icount = 0
        do k=1,irecs
          if (iend(k).ne.0) then
            isegment = isegment + 1
            icount = icount + 1
            irecs_in_segment(isegment) = icount
            icount = 0
          else
            icount = icount + 1
          endif
        enddo

c Extract the sort variable and define the index array.
        icount_1 = 0
        icount_2 = 0
        do kkk=1,isegment
          do k=1,irecs_in_segment(kkk)
            icount_1 = icount_1 + 1
            sort_var(k) = var_1(icount_1,ii)
            index(k) = icount_1
          enddo

c Get the new order (index).
          call isort(sort_var,index,irecs_in_segment(kkk))

c Re-order all of the variables.
          do k=1,irecs_in_segment(kkk)
            icount_2 = icount_2 + 1
            do i=1,nvars
              var_2(icount_2,i) = var_1(index(k),i)
            enddo
          enddo
        enddo

c Re-set the variable array.
        do k=1,irecs
          do i=1,nvars
            var_1(k,i) = var_2(k,i)
          enddo
        enddo

      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine isort(x,iy,n)

      implicit none

c Use insertion sort to sort an array of values, and keep track of
c   the original variable positions.  Force the algorithm to sort
c   the array in ascending order instead of decreasing order.

c X - array of values to be sorted.
c IY - after the sort, IY(J) contains the original position of the
c      value X(J) in the unsorted X array.
c N - number of values in array X to be sorted.

      integer n,maxlines
      parameter (maxlines=1000000)
      real x(maxlines)
      integer iy(maxlines)

      real temp
      integer i,j,k,itemp

c Force this to sort in ascending order.
      do i=1,n
        x(i) = - x(i)
      enddo

      do i=2,n
        if (x(i).gt.x(i-1)) then
          do j=i-2,1,-1
            if(x(i).lt.x(j)) go to 70
          enddo
          j=0
  70      temp = x(i)
          itemp = iy(i)
          do 90 k=i,j+2,-1
            iy(k) = iy(k-1)
  90        x(k) = x(k-1)
          x(j+1) = temp
          iy(j+1) = itemp
        endif
      enddo

c Adjust back for the ascending order sort.
      do i=1,n
        x(i) = - x(i)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

