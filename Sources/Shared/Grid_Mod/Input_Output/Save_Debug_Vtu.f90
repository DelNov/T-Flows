!==============================================================================!
  subroutine Save_Debug_Vtu(Grid, append,                                 &
                                  inside_cell, inside_name,               &
                                  scalar_cell, scalar_node, scalar_name,  &
                                  vector_cell, vector_node, vector_name,  &
                                  tensor_cell, tensor_node, tensor_name,  &
                                  tensor_comp,                            &
                                  plot_inside)
!------------------------------------------------------------------------------!
!>  The Save_Debug_Vtu subroutine is designed for writing .vtu files with
!>  additional debugging capabilities, especially useful from Process, but it
!>  can certainly be used from any other sub-program.  It can handle both
!>  cell-based and node-based data, offering a more detailed visualization of
!>  various scalar, vector, and tensor variables within the domain.
!>  Despite similarities with the Save_Vtu_Cells subroutine used in
!>  pre-processing, Save_Debug_Vtu has its distinct usage scenario during the
!>  processing stages. The decision to merge them would depend on how much the
!>  duplicated code affects maintenance and whether the combined complexity
!>  remains manageable.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid    !! grid under consideration
  character(*)     :: append  !! appendix to file name

  !-------------!
  !   Scalars   !
  !-------------!
  real, optional :: inside_cell(1:Grid % n_cells)
    !! scalar variables defined only on inside cell, which
    !! is only the right hand side vector "b" at this point
  character(*), optional :: inside_name
    !! name assigned to inside_cell

  real, optional :: scalar_cell(-Grid % n_bnd_cells:Grid % n_cells)
    !! scalar variables defined on both boundary and inside cells,
    !! like most variables in T-Flows
  real, optional :: scalar_node(1:Grid % n_nodes)
    !! scalar variables defined on nodes
  character(*), optional :: scalar_name
    !! name assigned to either scalar_cell or scalar_node, implying that
    !! scalar_cell and scalar_node can't be processed in the same call

  !-------------!
  !   Vectors   !
  !-------------!
  real, optional :: vector_cell(-Grid % n_bnd_cells:Grid % n_cells, 3)
    !! vector variables defined on both boundary and inside cells,
    !! like velocities and gradients in T-Flows, for example
  real, optional :: vector_node(1:Grid % n_nodes, 3)
    !! vector variables defined on nodes
  character(*), optional :: vector_name
    !! name assigned to either vector_cell or vector_node, implying that
    !! vector_cell and vector_node can't be processed in the same call

  !----------------------------------------------!
  !   Tensors                                    !
  !                                              !
  !   - symmetric tensors are stored in order:   !
  !     xx, yy, zz, xy, yz, xz:                  !
  !       | xx xy xz |   | 1 4 6 |               !
  !       |    yy yz | = |   2 5 |               !
  !       |       zz |   |     3 |               !
  !                                              !
  !   - non-symmetric tensors are stored in      !
  !     this order, explicitly defined:          !
  !       | xx xy xz |   | 1 2 3 |               !
  !       | yx yy yz | = | 4 5 6 |               !
  !       | zx zy zz |   | 7 8 9 |               !
  !----------------------------------------------!
  real, optional :: tensor_cell(-Grid % n_bnd_cells:Grid % n_cells, *)
    !! tensor variables defined on both boundary and inside cells
  real, optional :: tensor_node(1:Grid % n_nodes, *)
    !! tensor variables defined on nodes
  character(*), optional :: tensor_name
    !! name assigned to either tensor_cell or tensor_node, implying that
    !! tensor_cell and tensor_node can't be processed in the same call
  integer, optional :: tensor_comp
    !! number of tensor components, can be six for symmetric
    !! tensors or nine for non-symmetric tensors

  logical, optional :: plot_inside
    !! plot inside (true) or on booundary (false)
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)       :: data_size
  integer           :: c, n, s, i_fac, data_offset, cell_offset, fu
  integer           :: n_conns, n_polyg
  integer           :: cs, ce, nc, s1, s2
  logical           :: inside
  real              :: dist1, dist2
  real, allocatable :: values_cell(:)
  character(SL)     :: values_name, name_out, str1, str2
!==============================================================================!

  call Profiler % Start("Save_Debug_Vtu")

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  ! If inside cells (such as the b array) or other scalar cells are present
  ! allocate memory for values_cell (which subsequently take either of them)
  if(present(inside_cell) .or. present(scalar_cell)) then
    allocate(values_cell(-Grid % n_bnd_cells:Grid % n_cells));
    values_cell(:) = 0.0
  end if

  if(present(inside_cell)) then
    values_cell(1:Grid % n_cells) = inside_cell(1:Grid % n_cells)
    values_name = inside_name
  else if(present(scalar_cell)) then
    values_cell(-Grid % n_bnd_cells:Grid % n_cells) =  &
    scalar_cell(-Grid % n_bnd_cells:Grid % n_cells)
    values_name = scalar_name
  else
    ! An error trap here would be lovely
  end if

  ! Fetch the name for node-based variables
  if(present(scalar_node)) then
    values_name = scalar_name
  end if

  ! Set initial value for inside (which, if .true., means plotting inside cells)
  inside = .true.

  ! Take the value of argument, if provided
  if(present(plot_inside)) then
    inside = plot_inside
  end if

  ! Set start and ending cell depending if you save inside or boundary cells
  if(inside) then
    cs = 1
    ce = Grid % n_cells
  else
    cs = -Grid % n_bnd_cells
    ce = -1
  end if
  nc = ce - cs + 1

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do c = cs, ce
    n_conns = n_conns + abs(Grid % cells_n_nodes(c))
  end do

  ! Count face data for polyhedral cells, you will need it later
  n_polyg = 0
  do c = cs, ce
    if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
      n_polyg = n_polyg + 1                  ! add one for number of polyfaces
      do i_fac = 1, Grid % cells_n_faces(c)  ! add all faces and their nodes
        s = Grid % cells_f(i_fac, c)
        n = Grid % faces_n_nodes(s)
        n_polyg = n_polyg + 1 + n
      end do
    end if
  end do

  !----------------------!
  !                      !
  !   Create .vtu file   !
  !                      !
  !----------------------!
  call File % Set_Name(name_out,                                &
                       appendix  = '-'//trim(append),           &
                       processor = (/This_Proc(), N_Procs()/),  &
                       extension = '.vtu')
  call File % Open_For_Writing_Binary(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu) IN_0 // '<?xml version="1.0"?>'             // LF
  write(fu) IN_0 // '<VTKFile type="UnstructuredGrid"'  //  &
                    ' version="0.1"'                    //  &
                    ' byte_order="LittleEndian">'       // LF
  write(fu) IN_1 // '<UnstructuredGrid>'                // LF
  write(str1, '(i0.0)') Grid % n_nodes
  write(str2, '(i0.0)') nc
  write(fu) IN_2 // '<Piece NumberOfPoints="' // trim(str1) // '"' //  &
                    ' NumberOfCells="' // trim(str2) // '">'       // LF
  data_offset = 0

  !-----------!
  !           !
  !   Nodes   !
  !           !
  !-----------!
  write(str1, '(i1)') data_offset
  write(fu) IN_3 // '<Points>'                       // LF
  write(fu) IN_4 // '<DataArray type='//floatp       //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + Grid % n_nodes * RP * 3  ! prepare for next

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! First write all cells' nodes (a.k.a. connectivity)
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp         //  &
                    ' Name="connectivity"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP  ! prepare for next

  ! Cells' offsets
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp         //  &
                    ' Name="offsets"'                //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + nc * IP  ! prepare for next

  ! Cells' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp         //  &
                    ' Name="types"'                  //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + nc * IP  ! prepare for next

  ! For polyhedral grids, save faces and face offsets
  if(Grid % polyhedral) then

    ! Write polyhedral cells' faces
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="faces"'                  //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + n_polyg * IP  ! prepare for next

    ! Write polyhedral cells' faces offsets
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="faceoffsets"'            //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + nc * IP  ! prepare for next

  end if  ! is Grid polyhedral

  !----------------------!
  !   The end of cells   !
  !----------------------!
  write(fu) IN_3 // '</Cells>' // LF

  !----------------!
  !   Point data   !
  !----------------!
  if(     present(scalar_node)  &
     .or. present(vector_node)  &
     .or. present(tensor_node)) then
    write(fu) IN_3 // '<PointData>' // LF
  end if

  ! Additional node-based scalar array
  if(present(scalar_node)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//floatp             //  &
                      ' Name="'// trim(values_name) // '"'   //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + Grid % n_nodes * RP  ! prepare for next
  end if

  ! Additional node-based vector array
  if(present(vector_node)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//floatp             //  &
                      ' Name="'// trim(vector_name) // '"'   //  &
                      ' NumberOfComponents="3"'              //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + Grid % n_nodes * RP * 3  ! prepare for next
  end if

  ! Additional node-based tensor array
  if(present(tensor_node)) then
    write(str1, '(i0.0)') tensor_comp
    write(str2, '(i0.0)') data_offset
    if(tensor_comp .eq. 6) then
      write(fu) IN_4 // '<DataArray type='//floatp                   //  &
                        ' Name="'// trim(tensor_name) // '"'         //  &
                        ' NumberOfComponents="' // trim(str1) // '"' //  &
                        ' format="appended"'                         //  &
                        ' offset="' // trim(str2) // '">'            // LF
    else
      write(fu) IN_4 // '<DataArray type='//floatp                   //  &
                        ' Name="'// trim(tensor_name) // '"'         //  &
                        ' NumberOfComponents="' // trim(str1) // '"' //  &
                        ' ComponentName0="XX"'                       //  &
                        ' ComponentName1="XY"'                       //  &
                        ' ComponentName2="XZ"'                       //  &
                        ' ComponentName3="YX"'                       //  &
                        ' ComponentName4="YY"'                       //  &
                        ' ComponentName5="YZ"'                       //  &
                        ' ComponentName6="ZX"'                       //  &
                        ' ComponentName7="ZY"'                       //  &
                        ' ComponentName8="ZZ"'                       //  &
                        ' format="appended"'                         //  &
                        ' offset="' // trim(str2) // '">'            // LF
    end if
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + Grid % n_nodes * RP * tensor_comp
  end if

  if(     present(scalar_node)  &
     .or. present(vector_node)  &
     .or. present(tensor_node)) then
    write(fu) IN_3 // '</PointData>' // LF
  end if

  !---------------!
  !               !
  !   Cell data   !
  !               !
  !---------------!
  write(fu) IN_3 // '<CellData>' // LF

  ! Processor i.d.
  if(allocated(Grid % Comm % cell_proc)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="Processor"'              //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + nc * IP  ! prepare for next
  end if

  ! Additional cell scalar
  if(present(inside_cell) .or.  &
     present(scalar_cell)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//floatp             //  &
                      ' Name="'// trim(values_name) // '"'   //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + nc * RP  ! prepare for next
  end if

  ! Additional cell vector
  if(present(vector_cell)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//floatp             //  &
                      ' Name="'// trim(vector_name) // '"'   //  &
                      ' NumberOfComponents="3"'              //  &
                      ' format="appended"'                   //  &
                      ' offset="' // trim(str1)       //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + nc * RP * 3  ! prepare for next
  end if

  ! Additional cell tensor
  if(present(tensor_cell)) then
    write(str1, '(i0.0)') tensor_comp
    write(str2, '(i0.0)') data_offset
    if(tensor_comp .eq. 6) then
      write(fu) IN_4 // '<DataArray type='//floatp                   //  &
                        ' Name="'// trim(tensor_name) // '"'         //  &
                        ' NumberOfComponents="' // trim(str1) // '"' //  &
                        ' format="appended"'                         //  &
                        ' offset="' // trim(str2) // '">'            // LF
    else
      write(fu) IN_4 // '<DataArray type='//floatp                   //  &
                        ' Name="'// trim(tensor_name) // '"'         //  &
                        ' NumberOfComponents="' // trim(str1) // '"' //  &
                        ' ComponentName0="XX"'                       //  &
                        ' ComponentName1="XY"'                       //  &
                        ' ComponentName2="XZ"'                       //  &
                        ' ComponentName3="YX"'                       //  &
                        ' ComponentName4="YY"'                       //  &
                        ' ComponentName5="YZ"'                       //  &
                        ' ComponentName6="ZX"'                       //  &
                        ' ComponentName7="ZY"'                       //  &
                        ' ComponentName8="ZZ"'                       //  &
                        ' format="appended"'                         //  &
                        ' offset="' // trim(str2) // '">'            // LF
    end if
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + nc * RP * tensor_comp  ! prepare for next
  end if

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu) IN_3 // '</CellData>'         // LF
  write(fu) IN_2 // '</Piece>'            // LF
  write(fu) IN_1 // '</UnstructuredGrid>' // LF

  !-------------------!
  !                   !
  !   Appended data   !
  !                   !
  !-------------------!
  write(fu) IN_0 // '<AppendedData encoding="raw">' // LF
  write(fu) '_'

  !-----------!
  !   Nodes   !
  !-----------!
  data_size = int(Grid % n_nodes * RP * 3, SP)
  write(fu) data_size
  do n = 1, Grid % n_nodes
    write(fu) Grid % xn(n), Grid % yn(n), Grid % zn(n)
  end do

  !-----------!
  !   Cells   !
  !-----------!

  ! Cells' nodes
  data_size = int(n_conns * IP, SP)
  write(fu) data_size

  do c = cs, ce

    !---------------------------!
    !   Plotting inside cells   !
    !---------------------------!
    if(inside) then

      ! Tetrahedral, pyramid, wedge and hexahedral cells
      if( any( Grid % cells_n_nodes(c) .eq. (/4,5,6,8/)  ) ) then
        write(fu) (Grid % cells_n(1:Grid % cells_n_nodes(c), c))-1

      ! Polyhedral cells
      else if(Grid % cells_n_nodes(c) .lt. 0) then
        write(fu) (Grid % cells_n(1:-Grid % cells_n_nodes(c), c))-1

      end if

    !-----------------------------!
    !   Plotting boundary cells   !
    !-----------------------------!
    else

      ! All cell types
      write(fu) (Grid % cells_n(1:Grid % cells_n_nodes(c), c))-1

    end if  ! if inside
  end do

  ! Cells' offsets
  data_size = int(nc * IP, SP)
  write(fu) data_size
  cell_offset = 0
  do c = cs, ce
    cell_offset = cell_offset + abs(Grid % cells_n_nodes(c))
    write(fu) cell_offset
  end do

  ! Cells' types
  data_size = int(nc * IP, SP)
  write(fu) data_size
  if(inside) then
    do c = cs, ce
      if(Grid % cells_n_nodes(c) .eq. 4) write(fu) VTK_TETRA
      if(Grid % cells_n_nodes(c) .eq. 8) write(fu) VTK_HEXAHEDRON
      if(Grid % cells_n_nodes(c) .eq. 6) write(fu) VTK_WEDGE
      if(Grid % cells_n_nodes(c) .eq. 5) write(fu) VTK_PYRAMID
      if(Grid % cells_n_nodes(c) .lt. 0) write(fu) VTK_POLYHEDRON
    end do
  else
    do c = cs, ce
      if(Grid % cells_n_nodes(c) .eq. 3) then
        write(fu) VTK_TRIANGLE
      else if(Grid % cells_n_nodes(c) .eq. 4) then
        write(fu) VTK_QUAD
      else
        write(fu) VTK_POLYGON
      end if
    end do
  end if

  ! For polyhedral grids, save faces and face offsets
  if(Grid % polyhedral) then

    ! Write polyhedral cells' faces
    data_size = int(n_polyg * IP, SP)
    write(fu) data_size
    do c = cs, ce
      if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
        write(fu) Grid % cells_n_faces(c)      ! write number of its polyfaces
        do i_fac = 1, Grid % cells_n_faces(c)  ! and all polyfaces
          s = Grid % cells_f(i_fac, c)
          if(Grid % faces_s(s) .ne. 0) then    ! face has a shadow, if it ...
            s1 = s                             ! ... is closer, plot that!
            s2 = Grid % faces_s(s)
            dist1 = Math % Distance(                              &
                    Grid % xc(c),  Grid % yc(c),  Grid % zc(c),   &
                    Grid % xf(s1), Grid % yf(s1), Grid % zf(s1))
            dist2 = Math % Distance(                              &
                    Grid % xc(c),  Grid % yc(c),  Grid % zc(c),   &
                    Grid % xf(s2), Grid % yf(s2), Grid % zf(s2))
            if(dist1 < dist2) s = s1
            if(dist2 < dist1) s = s2
          end if
          n = Grid % faces_n_nodes(s)
          write(fu) n, Grid % faces_n(1:n, s)-1
        end do
      end if
    end do

    ! Write polyhedral cells' faces offsets
    data_size = int(Grid % n_cells * IP, SP)
    write(fu) data_size
    cell_offset = 0
    do c = 1, Grid % n_cells
      if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
        cell_offset = cell_offset + 1          ! to store number of polyfaces
        do i_fac = 1, Grid % cells_n_faces(c)  ! to store polyfaces
          s = Grid % cells_f(i_fac, c)
          n = Grid % faces_n_nodes(s)
          cell_offset = cell_offset + 1 + n    ! number of nodes and nodes
        end do
        write(fu) cell_offset                  ! write the current offset
      else
        write(fu) -1             ! not a polyhedron, offsets are not needed
      end if
    end do

  end if  ! if Grid % polyhedral

  !----------------!
  !   Point data   !
  !----------------!
  if(present(scalar_node)) then
    data_size = int(Grid % n_nodes * RP, SP)
    write(fu) data_size
    do n = 1, Grid % n_nodes
      write(fu) scalar_node(n)
    end do
  end if

  if(present(vector_node)) then
    data_size = int(Grid % n_nodes * RP * 3, SP)
    write(fu) data_size
    do n = 1, Grid % n_nodes
      write(fu) vector_node(n, 1), vector_node(n, 2), vector_node(n, 3)
    end do
  end if

  if(present(tensor_node)) then
    data_size = int(Grid % n_nodes * RP * tensor_comp, SP)
    write(fu) data_size

    ! Symmetric tensor
    if(tensor_comp .eq. 6) then
      do n = 1, Grid % n_nodes
        write(fu) tensor_node(n, 1), tensor_node(n, 2), tensor_node(n, 3),  &
                  tensor_node(n, 4), tensor_node(n, 5), tensor_node(n, 6)
      end do

    ! Non-symmetric tensor
    else if(tensor_comp .eq. 9) then
      do n = 1, Grid % n_nodes
        write(fu) tensor_node(n, 1), tensor_node(n, 2), tensor_node(n, 3),  &
                  tensor_node(n, 4), tensor_node(n, 5), tensor_node(n, 6),  &
                  tensor_node(n, 7), tensor_node(n, 8), tensor_node(n, 9)
      end do
    end if
  end if

  !---------------!
  !   Cell data   !
  !---------------!

  ! Processor i.d.
  if(allocated(Grid % Comm % cell_proc)) then
    data_size = int(nc * IP, SP)
    write(fu) data_size
    do c = cs, ce
      write(fu) Grid % Comm % cell_proc(c)
    end do
  end if

  ! Additional cell data
  if(present(scalar_cell) .or.   &
     present(inside_cell)) then
    data_size = int(nc * RP, SP)
    write(fu) data_size
    do c = cs, ce
      write(fu) values_cell(c)
    end do
  end if

  if(present(vector_cell)) then
    data_size = int(nc * RP * 3, SP)
    write(fu) data_size
    do c = cs, ce
      write(fu) vector_cell(c, 1), vector_cell(c, 2), vector_cell(c, 3)
    end do
  end if

  if(present(tensor_cell)) then
    data_size = int(nc * RP * tensor_comp, SP)
    write(fu) data_size

    ! Symmetric tensor
    if(tensor_comp .eq. 6) then
      do c = cs, ce
        write(fu) tensor_cell(c, 1), tensor_cell(c, 2), tensor_cell(c, 3),  &
                  tensor_cell(c, 4), tensor_cell(c, 5), tensor_cell(c, 6)
      end do

    ! Non-symmetric tensor
    else if(tensor_comp .eq. 9) then
      do c = cs, ce
        write(fu) tensor_cell(c, 1), tensor_cell(c, 2), tensor_cell(c, 3),  &
                  tensor_cell(c, 4), tensor_cell(c, 5), tensor_cell(c, 6),  &
                  tensor_cell(c, 7), tensor_cell(c, 8), tensor_cell(c, 9)
      end do
    end if
  end if

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  !-----------------------!
  !                       !
  !   Create .pvtu file   !
  !                       !
  !-----------------------!
  if(.not. allocated(Grid % Comm % cell_proc)) return

  ! Create it only from subdomain 1, when decomposed
  if(maxval(Grid % Comm % cell_proc(:)) > 1 .and. This_Proc() .eq. 1) then

    call File% Set_Name(name_out, appendix='-'//trim(append),  &
                        extension='.pvtu')
    call File % Open_For_Writing_Ascii(name_out, fu)

    ! Header
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(fu,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'

    ! This section must be present
    write(fu,'(a,a)') IN_2, '<PPoints>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp  // &
                            ' NumberOfComponents="3"/>'
    write(fu,'(a,a)') IN_2, '</PPoints>'

    ! Data section is not mandatory, but very useful
    write(fu,'(a,a)') IN_2, '<PCellData>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//intp//' Name="Processor"/>'
    if(present(scalar_cell) .or.  &
       present(inside_cell)) then
      write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp    //  &
                              ' Name="'// trim(values_name)  // '"/>'
    end if
    if(present(scalar_node)) then
      write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp    //  &
                              ' Name="'// trim(scalar_name)  // '"/>'
    end if
    if(present(vector_cell)) then
      write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp    //  &
                              ' NumberOfComponents="3"'      //  &
                              ' Name="'// trim(vector_name)  // '"/>'
    end if
    if(present(tensor_cell)) then
      write(str1, '(i0.0)') tensor_comp
      if(tensor_comp .eq. 6) then
        write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp                //  &
                                ' NumberOfComponents="'//trim(str1) // '"' //  &
                                ' Name="'// trim(tensor_name) // '"/>'
      else
        write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp                //  &
                                ' NumberOfComponents="'//trim(str1) // '"' //  &
                                ' Name="'// trim(tensor_name) // '"'       //  &
                                ' ComponentName0="XX"'                     //  &
                                ' ComponentName1="XY"'                     //  &
                                ' ComponentName2="XZ"'                     //  &
                                ' ComponentName3="YX"'                     //  &
                                ' ComponentName4="YY"'                     //  &
                                ' ComponentName5="YZ"'                     //  &
                                ' ComponentName6="ZX"'                     //  &
                                ' ComponentName7="ZY"'                     //  &
                                ' ComponentName8="ZZ"' // '/>'             // LF
      end if
    end if
    write(fu,'(a,a)') IN_2, '</PCellData>'

    ! Write out the names of all the pieces
    do n = 1, N_Procs()
      call File % Set_Name(name_out, appendix='-'//trim(append),  &
                           processor=(/n, N_Procs()/), extension='.vtu')
      write(fu, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
    end do

    ! Footer
    write(fu, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(fu, '(a,a)') IN_0, '</VTKFile>'

    close(fu)

  end if

  call Profiler % Stop("Save_Debug_Vtu")

  end subroutine
