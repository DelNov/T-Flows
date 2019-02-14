!==============================================================================!
  subroutine Cgns_Mod_Read_Bnd_Conds_Info(base, block, bc)
!------------------------------------------------------------------------------!
!   Reads boundary condition info.                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block, bc
!-----------------------------------[Locals]-----------------------------------!
  integer             :: base_id        ! base index number
  integer             :: block_id       ! block index number
  integer             :: bc_id          ! block index number
  character(len=80)   :: bc_name        ! name of the boundary condition
  integer             :: bc_type        ! boundary condition type
  integer             :: bc_ptset_type  ! boundary node/cell placement
  integer             :: bc_n_nodes     ! boundary nodes or cells number
  integer,allocatable :: point_list(:)  ! boundary nodes or cells list
  integer             :: NormalIndex(3)
  integer             :: NormalListFlag
  integer             :: NormalDataType
  integer,allocatable :: NormalList(:)  ! normal vector pointing into the block
  integer             :: bc_n_datasets
  integer             :: error
  integer             :: one = 1         ! go figure :-(
  integer             :: i, color, first_point, last_point
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  bc_id    = bc

  ! Position yourself at ZoneBC
  call Cg_Goto_F(file_id,     &  ! cgns file index number
                 base,        &  ! base index number
                 error,       &  ! error status
                 'Zone_t',    &  ! node of block type
                 block,       &  ! block index number
                 'ZoneBC_t',  &  ! search for node "ZoneBC_t"
                 one,         &  ! ???
                 'end')          ! indicates end of call
  if (error.ne.0) then
    print "(a)", " # Failed to navigate to ZoneBC node"
    call Cg_Error_Exit_F()
  endif

  ! Get boundary condition info
  call Cg_Boco_Info_F(file_id,         & !(in )
                      base_id,         & !(in )
                      block_id,        & !(in )
                      bc_id,           & !(in )
                      bc_name,         & !(out)
                      bc_type,         & !(out)
                      bc_ptset_type,   & !(out)
                      bc_n_nodes,      & !(out)
                      NormalIndex,     & !(out)
                      NormalListFlag,  & !(out)
                      NormalDataType,  & !(out)
                      bc_n_datasets,   & !(out)
                      error)             !(out)
  if (error .ne. 0) then
    print "(a)", " # Failed to read boundary conditions info"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % bnd_cond(bc) % name = trim(bc_name)
  ! "For a ptset_type of PointRange, npnts is always two"
  ! "For a ptset_type of PointList, npnts is the number of points or elements
  !  in the list"

  ! Fill up the boundary condition names
  color = 0
  do i = 1, cnt_bnd_conds
    if(bnd_cond_names(i) .eq. trim(bc_name)) then
      color = i
      goto 1
    end if
  end do
  cnt_bnd_conds = cnt_bnd_conds + 1
  bnd_cond_names(cnt_bnd_conds) = trim(bc_name)
  color = cnt_bnd_conds
1 continue
  cgns_base(base) % block(block) % bnd_cond(bc) % color = color

  ! Allocate space for b.c. specified as PointList or PointRange
  allocate( point_list(bc_n_nodes) ); point_list(:) = 0
  allocate( NormalList(bc_n_nodes) ); NormalList(:) = 0

  ! Read boundary condition data and normals
  call Cg_Boco_Read_F(file_id,        & ! (in )
                      base_id,        & ! (in )
                      block_id,       & ! (in )
                      bc_id,          & ! (in )
                      point_list(:),  & ! (out)
                      NormalList(:),  & ! (out)
                      error)            ! (out)

  if (error .ne. 0) then
    print "(a)", " # Failed to read boundary conditions PointList"
    call Cg_Error_Exit_F()
  endif

  !-------------------------------------------------------!
  !   Copy b.c. point_list to Cgns_Block_Type structure   !
  !-------------------------------------------------------!
  if (trim(PointSetTypeName(bc_ptset_type)) .eq. 'PointList') then

    cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes = bc_n_nodes

    ! Boundary points (faces from ZoneBC)
    allocate( cgns_base(base) % block(block) % bnd_cond(bc) % &
      point_list( bc_n_nodes ) )
    cgns_base(base) % block(block) % bnd_cond(bc) % point_list(:) = &
      point_list(:)

    ! Array which points at section number
    allocate( cgns_base(base) % block(block) % bnd_cond(bc) % &
      belongs_to_sect( bc_n_nodes ) )
    cgns_base(base) % block(block) % bnd_cond(bc) % belongs_to_sect(:) = 0

    !allocate( cgns_base(base) % block(block) % bnd_cond(bc) % &
    !  assigned( bc_n_nodes ) )
    !cgns_base(base) % block(block) % bnd_cond(bc) % assigned(:) = .false.
  elseif (trim(PointSetTypeName(bc_ptset_type)) .eq. 'PointRange' .or. &
          trim(PointSetTypeName(bc_ptset_type)) .eq. 'ElementRange' ) then

    first_point = min(point_list(1), point_list(2))
    last_point  = max(point_list(1), point_list(2))

    cgns_base(base) % block(block) % bnd_cond(bc) % &
      n_nodes = last_point - first_point + 1

    ! Boundary points (faces from ZoneBC)
    allocate( cgns_base(base) % block(block) % bnd_cond(bc) % &
      point_list( last_point - first_point + 1 ) )

    do i = first_point, last_point
      cgns_base(base) % block(block) % bnd_cond(bc) % point_list( &
        i - first_point + 1) = i
    end do

    ! Array which points at section number
    allocate( cgns_base(base) % block(block) % bnd_cond(bc) % &
      belongs_to_sect( last_point - first_point + 1 ) )
    cgns_base(base) % block(block) % bnd_cond(bc) % belongs_to_sect(:) = 0

    !allocate( cgns_base(base) % block(block) % bnd_cond(bc) % assigned( &
    !    last_point - first_point + 1 ) )
    !cgns_base(base) % block(block) % bnd_cond(bc) % assigned(:) = .false.

  end if

  if(verbose) then
    print "(a)", " #======================================================="
    print "(a,a23)", " #       Boundary condition name: ", &
             trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
    print "(a,i22)", ' #-------------------------------------------------------'
    print "(a,i22)", " #       Boundary condition index: ", bc
    print "(a,i22)", " #       Boundary condition color: ", color
    print "(a,i22)", " #       Boundary condition nodes: ", &
             cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
    print "(a,a21)", " #       Boundary condition Extent: ",   &
             trim(PointSetTypeName(bc_ptset_type))
    print "(a)",     ' #       Point list (sample): '
    print "(a,a12,6i7)", " # ", " ", (point_list(i), &
      i = 1, min(6,bc_n_nodes))

  end if

  ! Free memory
  deallocate(point_list)
  deallocate(NormalList)

  end subroutine