!==============================================================================!
  subroutine Cgns_Mod_Read_Bnd_Conds_Info(base, block, bc)
!------------------------------------------------------------------------------!
!   Reads boundary condition info.                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block, bc
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id         ! base index number
  integer           :: block_id        ! block index number
  integer           :: bc_id           ! block index number
  character(len=80) :: bc_name         ! name of the boundary condition
  integer           :: bc_type         ! boundary condition type
  integer           :: bc_ptset_type   ! boundary node/cell placement
  integer           :: bc_n_nodes      ! boundary nodes or cells
  integer           :: NormalIndex(3)
  integer           :: NormalListFlag
  integer           :: bc_data_type
  integer           :: bc_n_datasets
  integer           :: error
  integer           :: one = 1         ! go figure :-(
  integer           :: i, color
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
    print *,"# Failed to navigate to ZoneBC node"
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
                      bc_data_type,    & !(out)
                      bc_n_datasets,   & !(out)
                      error)             !(out)
  if (error .ne. 0) then
    print *,"# Failed to read boundary conditions info"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % bnd_cond(bc) % name    = trim(bc_name)
  !"For a ptset_type of PointRange, npnts is always two"
  cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes = bc_n_nodes

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

  if(verbose) then
    print *, '#       ----------------------------------------'
    print *, "#       Boundary condition name:   ",   &
             trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
    print *, '#       ----------------------------------------'
    print *, "#       Boundary condition index:  ", bc
    print *, "#       Boundary condition color:  ", color
    print *, "#       Boundary condition nodes:  ",   &
             cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
    print *, "#       Boundary condition Extent: ",   &
             PointSetTypeName(bc_ptset_type)

  end if

  end subroutine