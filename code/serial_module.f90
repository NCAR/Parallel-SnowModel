!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! serial_module.f90
module serial_module
  implicit none

  public

  integer :: me, np, global_ny, global_nx
  integer :: rem, l_ny, max_l_ny, max_l_nx,buff_index
  integer :: neww_lny, new_jj,smooth_index
  integer,parameter :: k15 = selected_int_kind(15)
  integer(kind=k15) :: byte_index_r1,byte_index_r2
  integer :: smooth_index_s
  integer :: salt_iter
  integer, allocatable :: partition_ny(:)
  integer, allocatable :: prefix_sum(:)
  integer, allocatable :: smooth_partition(:)

! Used for topo_data wind field generation in micromet
  real, allocatable :: l_topo_buff_n(:,:), l_topo_buff_s(:,:)
  real, allocatable :: topo_buffer(:,:)

! Used for smoother_array_2d parallelization
  real, allocatable :: l_smoother_buff(:,:)

! Used for saltation subroutine in SnowTran
  real, allocatable :: l_qsalt_buff(:), l_utau_buff(:)

! Used to create wind_buffer in SnowTran
  real, allocatable :: l_vwind_north(:), l_vwind_south(:)
  real, allocatable :: vwind_buffer(:,:)

! Used for halo exchange
  real, allocatable :: left_rcv(:), right_rcv(:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partitioning()
    use snowmodel_vars
    use git_version, only: git_branch,git_hash
    implicit none
    
    integer :: inc 
    real :: deltaxy
    
  90 format ('(i',i2,')')

    me = 0
    np = 1

! Print the number of processors and github information to output file   
    if (me.eq.0) then 
      print*, 'Number of processors =',np
      print*, 'Github branch name = ',git_branch
      print*, 'Github commit hash = ',git_hash
    endif

    allocate(partition_ny(0:np-1))
    allocate(prefix_sum(0:np-1))
    partition_ny = 0
    prefix_sum = 0

! Wind interpolation smoothing index and windowsize designation
    allocate(smooth_partition(0:np-1))


    rem = mod(ny,np)
    l_ny = (ny-rem)/np
    max_l_ny = l_ny

    global_ny = ny ! Global ny allocation variable
    global_nx = nx ! Global nx allocation variable

    deltaxy = 0.5 * (deltax + deltay)

    inc = max(1,nint(curve_len_scale/deltaxy))


    if(rem.gt.0) max_l_ny = max_l_ny + 1
    if(me.lt.rem) l_ny = l_ny + 1

    if (np.gt.1) then
      if (me.eq.0) then
        buff_index = l_ny + 1 !Buffer index for vwind_grid in SnowTran

        neww_lny = l_ny + inc !Buffer y-dimension for topo in micromet
        new_jj = 0 !Starting index for topo buffer in micromet
      elseif(me.eq.np-1) then
        buff_index = l_ny + 1 !Buffer index for vwind_grid in SnowTran

        neww_lny = l_ny + inc !Buffer y-dimension for topo in micromet
        new_jj = inc !Starting index for topo buffer in micromet
      else
        buff_index = l_ny + 2 !Buffer index for vwind_grid in SnowTran

        neww_lny = l_ny + (2*inc) !Buffer y-dimension for topo in micromet
        new_jj = inc !Starting index for topo buffer in micromet
      endif
    else
      buff_index = l_ny !Buffer index for vwind_grid in SnowTran

      neww_lny = l_ny !Buffer y-dimension for topo in micromet
      new_jj = 0 !Starting index for topo buffer in micromet
    endif
    
    call partition()
    call local_ny(l_ny,ny)
    call allocation(inc)

  end subroutine partitioning

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine processor_ID()
    use snowmodel_vars
    implicit none

    me = 0
    np = 1

  end subroutine processor_ID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partition()
    use snowmodel_vars
    implicit none
    integer :: ii,idx

    do idx = 0, np-1
       partition_ny(idx) = (ny-rem)/np
       if(idx.lt.rem) then
          partition_ny(idx) = partition_ny(idx) + 1
       end if
    end do

    do ii = 1, np-1
       prefix_sum(ii) = partition_ny(ii-1) + prefix_sum(ii-1)
    end do
    
    byte_index_r1 = (4*prefix_sum(me)*nx)+1 ! When binary files have 1 record
    byte_index_r2 = (nx*global_ny*4)+(4*prefix_sum(me)*nx)+1 ! When binary files have 2 records

  end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocation(inc)
    use snowmodel_vars
    implicit none
    
    integer, intent(in) :: inc ! Indicates whether local array should be created 
                               !  based on print_parallel parameter
    real,dimension(nx,max_l_ny) :: zero_real
    real,dimension(nx,max_l_ny,nz_max) :: zero_real_three_d

    max_l_nx = nx

    zero_real = 0.0
    zero_real_three_d = 0.0

! Halo exchange buffer allocation
    allocate(left_rcv(nx),right_rcv(nx))
    left_rcv = 0.0
    right_rcv = 0.0

! SnowTran buffer allocations
    allocate(l_qsalt_buff(nx))
    allocate(l_utau_buff(nx))
    allocate(l_vwind_north(nx))
    allocate(l_vwind_south(nx))
    l_utau_buff = 0.0
    l_qsalt_buff = 0.0
    l_vwind_north = 0.0
    l_vwind_south = 0.0

! Topo buffers
    allocate(l_topo_buff_n(nx,inc))
    allocate(l_topo_buff_s(nx,inc))
    allocate(vwind_buffer(nx,buff_index))
    allocate(topo_buffer(nx,neww_lny))
    l_topo_buff_n = 0.0
    l_topo_buff_s = 0.0
    vwind_buffer = 0.0
    topo_buffer = 0.0
    
! Wind smoothing buffers 
    allocate(l_smoother_buff(nx,max_l_ny))
    l_smoother_buff = 0.0
    
    call saltation_iteration()

  end subroutine allocation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine local_buffer_1D(local_array,local_buffer,north_buff_co,south_buff_co,nx,new_lny)

    real, intent(in) ::  local_array(:,:)
    real, intent(inout) :: south_buff_co(:)
    real, intent(inout) :: north_buff_co(:)
    integer, intent(in) :: nx, new_lny
    real, intent(out) :: local_buffer(:,:)


    if (np.gt.1) then
      continue
    else
      local_buffer(1:nx,1:l_ny) = local_array(1:nx,1:l_ny)
    endif

  end subroutine local_buffer_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine topo_buffer_1D(local_array, local_buffer,nx,inc,&
                        &   local_array_s, local_array_n)

    real, intent(in) ::  local_array(:,:)
    integer, intent(in) :: nx, inc
    real, intent(inout) :: local_array_n(:,:)
    real, intent(inout) :: local_array_s(:,:)
    real, intent(out) :: local_buffer(:,:)

    if (np.gt.1) then
      continue
    else
      local_buffer(1:nx,1:l_ny) = local_array(1:nx,1:l_ny)
    endif

  end subroutine topo_buffer_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine buffer_north(local_array, local_buffer)

    real, intent(in) ::  local_array(:,:)
    real, intent(out) :: local_buffer(:)


  end subroutine buffer_north

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine buffer_south(local_array, local_buffer,index_end)

    real, intent(in) ::  local_array(:,:)
    integer, intent(in) :: index_end
    real, intent(out) :: local_buffer(:)


  end subroutine buffer_south

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function global_map(j) result(new_j)
    implicit none
    integer, intent(in) :: j
    integer :: new_j

    new_j = j + prefix_sum(me)

  end function global_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine synchronization()
    implicit none


  end subroutine synchronization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine local_ny(ny_local, ny_global)
    implicit none
    integer, intent(out) :: ny_global
    integer, intent(in) :: ny_local


  end subroutine local_ny

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine maximum(num)
    implicit none
    real, intent(inout) :: num


  end subroutine maximum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine summation(num)
    implicit none
    real, intent(inout) :: num


  end subroutine summation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The halo exchange is done using "puts" instead of "gets"
  subroutine halo_exchange(array_2d)
    implicit none

    real, intent(in) :: array_2d(max_l_nx,max_l_ny)
    integer :: nx


  end subroutine halo_exchange

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partition_smoother(nx,ny,forcing_dx,deltax,smooth_index,windowsize)
    implicit none
    integer,intent(in) :: nx,ny
    real, intent(in) :: forcing_dx,deltax
    integer, intent(inout) :: smooth_index,windowsize
    integer :: resid,i


! north loop
    smooth_partition(:) = partition_ny(:)
    windowsize = int(forcing_dx / deltax / 2.0)
    resid = windowsize
! create index counter
    smooth_index = 0
    smooth_index_s = 1
    do i = 0,np-1
      smooth_index = smooth_index + smooth_partition(i)
      if (i.lt.me) then
        smooth_index_s = smooth_index_s + smooth_partition(i)
      endif
    enddo
    
  ! allocate local buffer
    ! allocate(smoother_buff(nx,smooth_index))
    ! allocate(l_smoother_buff(nx,max_l_ny))



  end subroutine partition_smoother

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine distribute_buffer(local_array,local_buffer,nx,ny)
    real, intent(in) ::  local_array(:,:)
    integer, intent(in) :: nx,ny
    real, intent(inout) :: local_buffer(:,:)

    call smoother_corray(local_array,nx,ny)
    call window_distribute(nx,ny,local_buffer)

  end subroutine distribute_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine smoother_corray(local_array,nx,ny)
    real, intent(in) ::  local_array(:,:)
    integer, intent(in) :: nx,ny

    l_smoother_buff(:,:) = local_array(:,:)
  end subroutine smoother_corray

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine window_distribute(nx,ny,local_buffer)
    implicit none
    integer,intent(in) :: nx,ny
    real,intent(inout) :: local_buffer(:,:)
    integer :: i,start_index,end_index,local_start
    
    local_buffer(:,:) = l_smoother_buff(:,:)



  end subroutine window_distribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine localize(buffer_array,local_array,nx,ny)

    real, intent(inout) ::  local_array(:,:),buffer_array(:,:)
    integer, intent(in) :: nx,ny

    local_array(:,:) = 0.0
    local_array(:,:) = buffer_array(:,:)

    buffer_array(:,:) = 0.0

  end subroutine localize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine saltation_iteration()
  use snowmodel_vars
  
  salt_iter = 1

end subroutine saltation_iteration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module serial_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
