!==============================================================================!
  subroutine Save_Surf(surf, time_step)
!------------------------------------------------------------------------------!
!   Writes surface vertices in CGNS file format (for VisIt and Paraview)       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  integer                 :: time_step
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: name_out
  integer       :: base
  integer       :: block
  integer       :: solution
  integer       :: field
  integer       :: coord
  integer       :: sect
  integer       :: i
  real          :: n_array(surf % n_verts)
  real          :: c_array(surf % n_elems)
!==============================================================================!

  if(surf % n_verts < 1) return

  call Cpu_Timer_Mod_Start('Save_Cgns_Surf')

  if(this_proc < 2) then

    call File_Mod_Set_Name(name_out,               &
                           time_step = time_step,  &
                           appendix  = '-surf',    &
                           extension = '.cgns') ! -> name_out

    ! Mostly repeats Save_Cgns_Cells

    file_mode = CG_MODE_WRITE
    call Cgns_Mod_Open_File(name_out, file_mode) ! -> file_id

    call Cgns_Mod_Initialize_Counters

    !-----------------!
    !                 !
    !   Bases block   !
    !                 !
    !-----------------!
    n_bases = 1
    allocate(cgns_base(n_bases))

    base = 1
    cgns_base(base) % name = 'Base 1'
    cgns_base(base) % cell_dim = 2 ! 2 -> surfaces cells, 3 -> volume cells
    cgns_base(base) % phys_dim = 3 ! 3 points to define a vector

    call Cgns_Mod_Write_Base_Info(base)

    !-----------------!
    !                 !
    !   Zones block   !
    !                 !
    !-----------------!

    cgns_base(base) % n_blocks = 1
    allocate(cgns_base(base) % block(cgns_base(base) % n_blocks))

    block = 1
    cgns_base(base) % block(block) % name = 'Zone 1'
    cgns_base(base) % block(block) % mesh_info(1) = surf % n_verts
    cgns_base(base) % block(block) % mesh_info(2) = surf % n_elems
    cgns_base(base) % block(block) % mesh_info(3) = 0

    call Cgns_Mod_Write_Block_Info(base, block)

    !-----------------------!
    !                       !
    !   Coordinates block   !
    !                       !
    !-----------------------!

    coord = 1
    cgns_base(base) % block(block) % coord_name(coord) = 'CoordinateX'

    coord = 2
    cgns_base(base) % block(block) % coord_name(coord) = 'CoordinateY'

    coord = 3
    cgns_base(base) % block(block) % coord_name(coord) = 'CoordinateZ'

    ! Actually write grid coordinates in DB
    do coord = 1, cgns_base(base) % phys_dim
      call Write_Surf_Coordinate_Array(base, block, coord, surf)
    end do

    ! Mostly repeats Save_Cgns_Cells

    !-----------------------------!
    !                             !
    !   Cells connections block   !
    !                             !
    !-----------------------------!
    cnt_cells = 0

    cgns_base(base) % block(block) % n_sects = 1
    allocate(cgns_base(base) % block(block) % section( &
      cgns_base(base) % block(block) % n_sects))

    sect = 1
    cgns_base(base) % block(block) % section(sect) % name = 'Elem_Triangles'
    cgns_base(base) % block(block) % section(sect) % cell_type = TRI_3
    cgns_base(base) % block(block) % section(sect) % first_cell = 1
    cgns_base(base) % block(block) % section(sect) % last_cell = surf % n_elems

    ! Mostly repeats Write_Section_Connections
    call Write_Surf_Section_Connections(base, block, sect, surf)

    ! Mostly repeats Save_Results

    !--------------------!
    !                    !
    !   Solution block   !
    !                    !
    !--------------------!

    cgns_base(base) % block(block) % n_solutions = 2

    allocate(cgns_base(base) % block(block) % solution( &
      cgns_base(base) % block(block) % n_solutions))

    ! First solution block is for node(vertex) centered data
    solution = 1
    cnt_cells = 0

    cgns_base(base) % block(block) % solution(solution) % name = &
      'PointData'
    cgns_base(base) % block(block) % solution(solution) % sol_type = Vertex
    call Cgns_Mod_Write_Solution_Info(base, block, solution)

    ! Mostly repeats Write_Field

    do i = 1, surf % n_verts
      n_array(i) = real(i)
    end do
    call Write_Surf_Field_In_Nodes(base, block, solution, field, surf, &
                                   n_array, 'Index')

    n_array = real(surf % vert(1:surf % n_verts) % nne)
    call Write_Surf_Field_In_Nodes(base, block, solution, field, surf, &
                                   n_array, 'NeighboursN')

    ! First solution block is for cell centered data
    solution = 2
    cnt_cells = 0

    cgns_base(base) % block(block) % solution(solution) % name = &
      'CellData'
    cgns_base(base) % block(block) % solution(solution) % sol_type = CellCenter
    call Cgns_Mod_Write_Solution_Info(base, block, solution)

    c_array = real(surf % elem(1:surf % n_elems) % nne)
    call Write_Field_In_Cells(base, block, solution, field, surf, &
                              c_array, 'NeighboursC')

    c_array = surf % elem(1:surf % n_elems) % nx
    call Write_Field_In_Cells(base, block, solution, field, surf, &
                              c_array, 'SurfaceNormalsX')

    c_array = surf % elem(1:surf % n_elems) % ny
    call Write_Field_In_Cells(base, block, solution, field, surf, &
                              c_array, 'SurfaceNormalsY')

    c_array = surf % elem(1:surf % n_elems) % nz
    call Write_Field_In_Cells(base, block, solution, field, surf, &
                              c_array, 'SurfaceNormalsZ')

    ! Close DB
    call Cgns_Mod_Close_File

  end if !(this_proc < 2) then

  if (this_proc < 2) print *, '# Saved surface to ', trim(name_out)

  deallocate(cgns_base)

  call Cpu_Timer_Mod_Stop('Save_Cgns_Surf')

  contains

!==============================================================================!
  subroutine Write_Surf_Coordinate_Array(base, block, coord, surf)
!------------------------------------------------------------------------------!
!   Writes surf grid coordinates (RealDouble)                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                 :: base, block, coord
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  integer       :: base_id         ! base index number
  integer       :: block_id        ! block index number
  integer       :: coord_id        ! coord index number
  character(SL) :: coord_name
  integer       :: i               ! lower range index
  integer       :: j               ! upper range index
  integer       :: error           ! error status
  real          :: coordinates(surf % n_verts)
!==============================================================================!

  ! Set input parameters
  base_id    = base
  block_id   = block
  coord_id   = coord
  coord_name = cgns_base(base_id) % block(block_id) % coord_name(coord_id)

  i = 1
  j = cgns_base(base_id) % block(block_id) % mesh_info(1)

  ! Fetch received parameters
  select case (coord)
    case (1)
      coordinates = surf % vert(1:surf % n_verts) % x_n
    case (2)
      coordinates = surf % vert(1:surf % n_verts) % y_n
    case (3)
      coordinates = surf % vert(1:surf % n_verts) % z_n
  end select

  ! Write grid coordinates
  call Cg_Coord_Write_F(file_id,      & !(in )
                        base_id,      & !(in )
                        block_id,     & !(in )
                        RealDouble,   & !(in )
                        coord_name,   & !(in )
                        coordinates,  & !(in )
                        coord_id,     & !(out)
                        error)          !(out)

  if (error .ne. 0) then
    print *, '#         Failed to write: ', trim(coord_name)
    call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose) then
    print *, '#         Coord array: ', coord_name
  end if
  if(verbose.and.coord_id.eq.1) then
    print *, '#         Number of nodes: ', j - i + 1
    print *, '#         First node:', i
    print *, '#         Last node: ', j
  end if

  end subroutine

!==============================================================================!
  subroutine Write_Surf_Section_Connections(base, block, sect, surf)
!------------------------------------------------------------------------------!
!   Writes surf elements connection for sect                                   !
!------------------------------------------------------------------------------!
!   Each node in zone/block must have unique id                                !
!   https://cgns.github.io/CGNS_docs_current/midlevel/grid.html                !
!                                                                              !
!   Connections:  | TRI_3 | QUAD_4 |                                           !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                 :: base, block, sect
  type(Surf_Type), target :: surf
!-----------------------------------[Locals]-----------------------------------!
  integer       :: base_id       ! base index number
  integer       :: block_id      ! block index number
  integer       :: sect_id       ! element section index
  character(SL) :: sect_name     ! name of the Elements_t node
  integer       :: cell_type     ! types of elements in the section
  integer       :: first_cell    ! index of first element
  integer       :: last_cell     ! index of last element
  integer       :: n_bnd         ! index of last boundary element
  integer       :: error
  integer       :: n_nodes, e, cnt, i, j
  integer       :: cell_n(4, surf % n_elems)
!==============================================================================!

  ! Set input parameters
  base_id   = base
  block_id  = block
  sect_id   = sect

  sect_name = trim(cgns_base(base_id)%block(block_id)%section(sect_id)%name)
  cell_type = cgns_base(base_id)%block(block_id)%section(sect_id)%cell_type
  i         = cgns_base(base_id)%block(block_id)%section(sect_id)%first_cell
  j         = cgns_base(base_id)%block(block_id)%section(sect_id)%last_cell

  ! cells of cell_type
  cnt = j - i + 1
  if ( cnt .ne. 0 ) then

    ! first and last cells have to be shifted according to previous sections
    first_cell = 1 + cnt_cells
    last_cell  = first_cell - 1 + cnt

    n_bnd = 0 ! unsorted boundary elements

    ! Allocate memory
    if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) n_nodes = 3
    if ( ElementTypeName(cell_type) .eq. 'QUAD_4') n_nodes = 4

    ! Convert T-FlowS -> CGNS [same as VTK]
    i = 1
    do e = 1, surf % n_elems
      !if (surf % surf_n_nodes(c).eq.3 .and. n_nodes.eq.3 ) then ! tri
      cell_n(1, i) = surf % elem(e) % i
      cell_n(2, i) = surf % elem(e) % j
      cell_n(3, i) = surf % elem(e) % k
      i = i + 1
    end do

    ! Write element data
    call Cg_Section_Write_F(file_id,                   & !(in )
                            base_id,                   & !(in )
                            block_id,                  & !(in )
                            sect_name,                 & !(in )
                            cell_type,                 & !(in )
                            first_cell,                & !(in )
                            last_cell,                 & !(in )
                            n_bnd,                     & !(in )
                            cell_n(1:n_nodes, 1:cnt),  & !(in )
                            sect_id,                   & !(out)
                            error)                       !(out)

    if (error .ne. 0) then
      print*, ' #         Failed to write ', trim(sect_name), ' connections'
      call Cg_Error_Exit_F()
    endif

    ! Print some info
    if(verbose ) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:    ', sect_id
      print *, '#         Cell section type:   ', ElementTypeName(cell_type)
      print *, '#         Number of cells: ', cnt
      print *, '#         First cell:', first_cell
      print *, '#         Last cell: ', last_cell
    end if

    cnt_cells = cnt_cells + cnt
  end if

  end subroutine

!==============================================================================!
  subroutine Write_Surf_Field_In_Nodes(base, block, solution, field, surf, &
    input_array, input_name)
!------------------------------------------------------------------------------!
!   Writes field to solution node and sets its field                           !
!   Solution type: Vertex                                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                 :: base, block, solution, field
  type(Surf_Type), target :: surf
  real                    :: input_array(surf % n_verts)
  character(len=*)        :: input_name
!-----------------------------------[Locals]-----------------------------------!
  integer                 :: base_id     ! base index number
  integer                 :: block_id    ! block index number
  integer                 :: solution_id ! solution index
  integer                 :: field_id    ! field index
  character(SL)           :: field_name  ! name of the FlowSolution_t node
  integer                 :: cnt            ! cells of sect_id
  real                    :: field_array(surf % n_verts) ! field array
  integer                 :: i, j, k, v
  integer                 :: error
!==============================================================================!

  ! Set input parameters
  base_id     = base
  block_id    = block
  solution_id = solution

  field_name = trim(input_name)

  ! Find first and last cells of sect_id
  ! cells of sect_id
  cnt = surf % n_verts

  if (cnt.ne.0) then

    i = 1 + cnt_cells
    j = i + cnt - 1

    ! Copy input array to field_array
    k = 1
    do v = 1, surf % n_verts
      field_array(k) = input_array(v)
      k = k + 1
    end do

    !--------------------------------------------!
    !   Writing cells assocciated with sect_id   !
    !--------------------------------------------!

    ! Add field to FlowSolution_t node for sect_id
    call Cg_Field_Partial_Write_F(file_id,      & !(in )
                                  base_id,      & !(in )
                                  block_id,     & !(in )
                                  solution_id,  & !(in )
                                  RealDouble,   & !(in )
                                  field_name,   & !(in )
                                  i,            & !(in )
                                  j,            & !(in )
                                  field_array,  & !(in )
                                  field_id,     & !(out)
                                  error)          !(out)

    if (error .ne. 0) then
      print *, '# Failed to write field: ', trim(field_name)
      call Cg_Error_Exit_F()
    endif

    cnt_cells = cnt_cells + cnt

    ! Print some info
    if(verbose ) then
      print *, '#           ---------------------------------'
      print *, '#           Field name: ',   field_name
      print *, '#           Field idx:    ', field_id
      print *, '#           ---------------------------------'
      print *, '#           First cell:', i
      print *, '#           Last cell: ', j
    end if

  end if

  cnt_cells = 0

  end subroutine

!==============================================================================!
  subroutine Write_Field_In_Cells(base, block, solution, field, surf, &
    input_array, input_name)
!------------------------------------------------------------------------------!
!   Writes field to solution node and sets its field                           !
!   Solution type: Face centered                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                 :: base, block, solution, field
  type(Surf_Type), target :: surf
  real                    :: input_array(surf % n_elems)
  character(len=*)        :: input_name
!-----------------------------------[Locals]-----------------------------------!
  integer       :: base_id        ! base index number
  integer       :: block_id       ! block index number
  integer       :: solution_id    ! solution index
  integer       :: field_id       ! field index
  character(SL) :: field_name     ! name of the FlowSolution_t node
  integer       :: cnt            ! cells of sect_id
  real          :: field_array(surf % n_elems) ! field array
  integer       :: i, j, k, e
  integer       :: error
!==============================================================================!

  ! Set input parameters
  base_id     = base
  block_id    = block
  solution_id = solution
  field_id    = field

  field_name = trim(input_name)

  !---------------------!
  !   Mapping 1:n_elems !
  !---------------------!

  ! Find first and last cells of sect_id
  ! cells of sect_id
  cnt = surf % n_elems

  if (cnt.ne.0) then

    i = 1 + cnt_cells
    j = i + cnt - 1

    ! Copy input array to field_array
    k = 1
    do e = 1, surf % n_elems
      field_array(k) = input_array(e)
      k = k + 1
    end do

    !--------------------------------------------!
    !   Writing cells assocciated with sect_id   !
    !--------------------------------------------!

    ! Add field to FlowSolution_t node for sect_id
    call Cg_Field_Partial_Write_F(file_id,      & !(in )
                                  base_id,      & !(in )
                                  block_id,     & !(in )
                                  solution_id,  & !(in )
                                  RealDouble,   & !(in )
                                  field_name,   & !(in )
                                  i,            & !(in )
                                  j,            & !(in )
                                  field_array,  & !(in )
                                  field_id,     & !(out)
                                  error)          !(out)

    if (error .ne. 0) then
      print *, '# Failed to write field ', trim(field_name)
      call Cg_Error_Exit_F()
    endif

    cnt_cells = cnt_cells + cnt

    ! Print some info
    if(verbose ) then
      print *, '#           ---------------------------------'
      print *, '#           Field name: ',   field_name
      print *, '#           Field idx:    ', field_id
      print *, '#           ---------------------------------'
      print *, '#           First cell:', i
      print *, '#           Last cell: ', j
    end if

  end if

  cnt_cells = 0

  end subroutine

  end subroutine
