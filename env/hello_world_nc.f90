program hello_world_nc  
! This is a basic hello_world program for coarray fortran and netcdf  
  use netcdf
  implicit none 
  
  integer :: me,np,i
  me = this_image() - 1
  np = num_images()
  
  do i=0,np-1
    if (me.eq.i) then 
      print*,'NetCDF worked! I am process number:',me,'out of',np
    endif 
  enddo 
end program