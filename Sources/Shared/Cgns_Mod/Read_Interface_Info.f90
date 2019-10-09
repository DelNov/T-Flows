!==============================================================================!
  subroutine Cgns_Mod_Read_Interface_Info(base, block, interface)
!------------------------------------------------------------------------------!
!   Reads general interface info                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer             :: base, block, interface
!-----------------------------------[Locals]-----------------------------------!
  integer             :: base_id           ! base index number
  integer             :: block_id          ! block index number
  integer             :: int_id            ! interface index number
  character(len=80)   :: int_name          ! name of the interface
  character(len=80)   :: int_name_short    ! shortened name of the interface
  integer             :: int_type          ! interface type
  integer             :: int_location      ! interface location
  integer             :: int_ptset_type    ! interface node/cell placement
  integer             :: int_n_nodes       ! interface nodes or cells count
  character(len=80)   :: donor_name        ! name of the donor
  integer             :: donor_block_type  ! donor block type
  integer             :: donor_data_type   ! donor data type
  integer             :: donor_ptset_type  ! donor node/cell placement
  integer             :: donor_n_nodes     ! donor nodes or cells count
  integer,allocatable :: point_list(:)     ! interface nodes or cells list
  integer             :: error, i
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  int_id   = interface

  ! Get interface info
  call Cg_Conn_Info_F(file_id,           & !(in )
                      base_id,           & !(in )
                      block_id,          & !(in )
                      int_id,            & !(in )
                      int_name,          & !(out)
                      int_location,      & !(out)
                      int_type,          & !(out)
                      int_ptset_type,    & !(out)
                      int_n_nodes,       & !(out)
                      donor_name,        & !(out)
                      donor_block_type,  & !(out)
                      donor_ptset_type,  & !(out)
                      donor_data_type,   & !(out)
                      donor_n_nodes,     & !(out)
                      error)               !(out)

  if (error .ne. 0) then
    print *, '# Failed to read interface info'
    call Cg_Error_Exit_F()
    stop
  end if

  !'For a int_ptset_type of PointRange, int_n_nodes is always two.
  ! For a int_ptset_type of PointList, int_n_nodes is the number of
  ! points in the PointList'

  if (int_n_nodes .eq. 2) then
    print *, '# Interface is set as PointRange -> Unsupported'
    stop
    ! To do: int_n_nodes .eq. 2
  end if

  ! Cut substring from last 'dom-'
  int_name_short = int_name(index(int_name,'dom-', &
    back = .true.):)

  ! allocate space for interface specified as PointList or PointRange
  allocate( point_list(int_n_nodes) ); point_list = 0

  ! Read interface data
  call Cg_Conn_Read_Short_F(file_id,        & !(in )
                            base_id,        & !(in )
                            block_id,       & !(in )
                            int_id,         & !(in )
                            point_list(:),  & !(out)
                            error)            !(out)

  if (error .ne. 0) then
    print *, '# Failed to read interface PointList'
    call Cg_Error_Exit_F()
    stop
  end if

  cgns_base(base) % block(block) % interface(int_id) % n_nodes = &
    int_n_nodes

  ! List of faces, which belongs to this interface
  allocate(cgns_base(base) % block(block) % interface(int_id) % &
    point_list(1:int_n_nodes))

  ! Array which points at section number
  allocate( cgns_base(base) % block(block) % interface(int_id) % belongs_to_sect(1:int_n_nodes) )

  cgns_base(base) % block(block) % interface(int_id) % belongs_to_sect(:) = 0

  !allocate(cgns_base(base) % block(block) % interface(int_id) % &
  !  assigned(1:int_n_nodes))
  !cgns_base(base) % block(block) % interface(int_id) % assigned(:) = .false.

  ! Fetch received parameters
  cgns_base(base) % block(block) % interface(int_id) % &
    point_list(:) = point_list(:)

  cgns_base(base_id) % block(block_id) % interface(int_id) % name = &
    trim(int_name_short)

  if(verbose) then
    print '(a)',     ' #======================================================='
    print '(a, a32)',' #       Interface name: ', trim(int_name)
    print '(a, a32)',' #       Shortened to:   ', trim(int_name_short)
    print '(a, a32)',' #       Donor name:     ', trim(donor_name)
    print '(a, i31)', ' #       Interface nodes: ', &
             cgns_base(base) % block(block) % interface(int_id) % n_nodes
    print '(a,a30)', ' #       Interface Extent: ',  'PointList'
    print '(a)',     ' #       Point list (sample): '
    print '(a,a12,6i7)', ' # ', ' ', (point_list(i), &
      i = 1, min(6,int_n_nodes))
    print '(a)',     ' #-------------------------------------------------------'
  end if

  ! Free memory
  deallocate(point_list)

  end subroutine
