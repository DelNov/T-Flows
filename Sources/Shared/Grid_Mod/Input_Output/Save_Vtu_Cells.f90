!==============================================================================!
  subroutine Save_Vtu_Cells(Grid, sub, n_nodes_sub, n_cells_sub)
!------------------------------------------------------------------------------!
!>  Writes cells in vtu file format, used only from pre-processing stages,
!>  meaing from Generate, Convert and Divide.  Depending on the arguments sent
!>  to the function, it can write a whole domain (Generate and Convert), or
!>  an individual a sub-domain (from Divide).
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine begins with some initial setup, such as setting precision !
!     and counting connections and polyhedral cells.                           !
!   * Memory allocation for local buffers is performed next.                   !
!   * It opens a .vtu file for writing and starts constructing the file        !
!     content, which includes header information, node data, cell data, and    !
!     cell-related information like cell types, offsets, etc.                  !
!   * For polyhedral grids, additional steps are included to handle faces and  !
!     their offsets.                                                           !
!   * Cell data, such as processor ID, number of nodes, wall distace, cell     !
!     volume, cell inertia, and cell porosity, is appended.                    !
!   * The subroutine concludes by closing the VTU file and, if necessary,      !
!     creating a .pvtu file for multi-processor runs.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid      !! grid being processes
  integer, intent(in) :: sub(1:2)  !! subdomain and total number of subdomains
  integer, intent(in) :: n_nodes_sub  !! number of nodes in a (sub)domain
  integer, intent(in) :: n_cells_sub  !! number of cells in a (sub)domain
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)          :: data_size
  integer              :: c, s, i_fac, data_offset, cell_offset, fu, s1, s2
  integer              :: n_conns, n_polyg, i, j, n, n1, n2
  real                 :: dist1, dist2
  character(SL)        :: name_out
  character(DL*2)      :: str1, str2
  real,    allocatable :: r_buffer(:)
  integer, allocatable :: i_buffer(:)
!==============================================================================!

  ! A very rudimentary test to check if regions are set
  Assert(Grid % region % f_cell(Grid % n_regions) .eq.              1)
  Assert(Grid % region % l_cell(Grid % n_regions) .eq. Grid % n_cells)

  call Profiler % Start('Save_Vtu_Cells')

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) n_conns = n_conns + abs(Grid % cells_n_nodes(c))
  end do

  ! Count face data for polyhedral cells, you will need it later
  n_polyg = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
      if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
        n_polyg = n_polyg + 1                  ! add one for number of polyfaces
        do i_fac = 1, Grid % cells_n_faces(c)  ! add all faces and their nodes
          s = Grid % cells_f(i_fac, c)
          n = Grid % faces_n_nodes(s)
          n_polyg = n_polyg + 1 + n
        end do
      end if
    end if
  end do

  ! Allocate memory for local buffers
  allocate(r_buffer(max(Grid % n_nodes * 3, Grid % n_cells * 6)))
  n1 = 0
  do c = Cells_In_Domain()
    n1 = n1 + max(abs(Grid % cells_n_nodes(c)), Grid % cells_n_faces(c))
  end do

  ! Still counting
  n2 = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
      if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
        n2 = n2 + 1
        do i_fac = 1, Grid % cells_n_faces(c)  ! and all polyfaces
          s = Grid % cells_f(i_fac, c)
          n = Grid % faces_n_nodes(s)
          n2 = n2 + 1                          ! to store number of nodes
          n2 = n2 + n                          ! to store each node
        end do
      end if
    end if
  end do

  allocate(i_buffer(max(Grid % n_nodes, max(n1,n2))))

  !------------------------!
  !   Open the .vtu file   !
  !------------------------!
  call File % Set_Name(name_out, processor=sub, extension='.vtu')
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
  write(fu) IN_1 // '<UnstructuredGrid>' // LF
  write(str1, '(i0.0)') n_nodes_sub
  write(str2, '(i0.0)') n_cells_sub
  write(fu) IN_2 // '<Piece NumberOfPoints="' // trim(str1) // '"' //  &
                    ' NumberOfCells="' // trim(str2) // '">'       // LF
  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  write(str1, '(i1)') data_offset
  write(fu) IN_3 // '<Points>'                       // LF
  write(fu) IN_4 // '<DataArray type='//floatp       //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + n_nodes_sub * RP * 3  ! prepare for next

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! Cells' nodes
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
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  ! Cells' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp         //  &
                    ' Name="types"'                  //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

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
    data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  end if  ! is Grid polyhedral

  !----------------------!
  !   The end of cells   !
  !----------------------!
  write(fu) IN_3 // '</Cells>' // LF

  !---------------!
  !               !
  !   Cell data   !
  !               !
  !---------------!
  write(fu) IN_3 // '<CellData>' // LF

  ! Cell processor
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp            //  &
                    ' Name="Grid Processor [1]"'        //  &
                    ' format="appended"'                //  &
                    ' offset="' // trim(str1) //'">'    // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  ! Number of nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp            //  &
                    ' Name="Grid Number Of Nodes [1]"'  //  &
                    ' format="appended"'                //  &
                    ' offset="' // trim(str1) //'">'    // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

  ! Wall distance
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//floatp          //  &
                    ' Name="Grid Wall Distance [m]"'    //  &
                    ' format="appended"'                //  &
                    ' offset="' // trim(str1) //'">'    // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * RP  ! prepare for next

  ! Cell volume
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//floatp          //  &
                    ' Name="Grid Cell Volume [m^3]"'    //  &
                    ' format="appended"'                //  &
                    ' offset="' // trim(str1) //'">'    // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * RP  ! prepare for next

  ! Cell centers
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//floatp          //  &
                    ' Name="Grid Cell Centers [m]"'     //  &
                    ' NumberOfComponents="3"'           //  &
                    ' format="appended"'                //  &
                    ' offset="' // trim(str1) //'">'    // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * RP * 3  ! prepare for next

  ! Cell inertia
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//floatp          //  &
                    ' Name="Grid Cell Inertia [m^2]"'   //  &
                    ' NumberOfComponents="6"'           //  &
                    ' format="appended"'                //  &
                    ' offset="' // trim(str1) //'">'    // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * RP * 6  ! prepare for next

  ! Cell porosity region
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp                 //  &
                    ' Name="Grid Cell Porosity Region [1]"'  //  &
                    ' format="appended"'                     //  &
                    ' offset="' // trim(str1) //'">'         // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_cells_sub * IP  ! prepare for next

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
  data_size = int(n_nodes_sub * RP * 3, SP)
  write(fu) data_size
  i = 0
  do n = 1, Grid % n_nodes
    if(Grid % new_n(n) .ne. 0) then
      i=i+1;  r_buffer(i) = Grid % xn(n)
      i=i+1;  r_buffer(i) = Grid % yn(n)
      i=i+1;  r_buffer(i) = Grid % zn(n)
    end if
  end do
  write(fu) r_buffer(1:i)

  !-----------!
  !   Cells   !
  !-----------!

  ! Cells' nodes
  data_size = int(n_conns * IP, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then

      ! Tetrahedral, pyramid, wedge and hexahedral cells
      if( any( Grid % cells_n_nodes(c) .eq. (/4,5,6,8/)  ) ) then
        do j = 1, Grid % cells_n_nodes(c)
          i=i+1;  i_buffer(i) = Grid % new_n(Grid % cells_n(j, c)) - 1
        end do

      ! Polyhedral cells
      else if(Grid % cells_n_nodes(c) < 0) then
        do j = 1, -Grid % cells_n_nodes(c)
          i=i+1;  i_buffer(i) = Grid % new_n(Grid % cells_n(j, c)) - 1
        end do

      else
        print *, '# Unsupported cell type with ',  &
                    Grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        stop
      end if
    end if
  end do
  write(fu) i_buffer(1:i)

  ! Cells' offsets
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  i = 0
  cell_offset = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      cell_offset = cell_offset + abs(Grid % cells_n_nodes(c))
      i=i+1;  i_buffer(i) = cell_offset
    end if
  end do
  write(fu) i_buffer(1:i)

  ! Cells' types
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1
      if(Grid % cells_n_nodes(c) .eq. 4) i_buffer(i) = VTK_TETRA
      if(Grid % cells_n_nodes(c) .eq. 8) i_buffer(i) = VTK_HEXAHEDRON
      if(Grid % cells_n_nodes(c) .eq. 6) i_buffer(i) = VTK_WEDGE
      if(Grid % cells_n_nodes(c) .eq. 5) i_buffer(i) = VTK_PYRAMID
      if(Grid % cells_n_nodes(c) .lt. 0) i_buffer(i) = VTK_POLYHEDRON
    end if
  end do
  write(fu) i_buffer(1:i)

  ! For polyhedral grids, save faces and face offsets
  if(Grid % polyhedral) then

    ! Write polyhedral cells' faces
    data_size = int(n_polyg * IP, SP)
    write(fu) data_size
    i = 0
    do c = Cells_In_Domain()
      if(Grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
        if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
          i=i+1
          i_buffer(i) = Grid % cells_n_faces(c)  ! write number of its polyfaces
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
            i=i+1;  i_buffer(i) = n
            do j = 1, n
              i=i+1;  i_buffer(i) = Grid % new_n(Grid % faces_n(j, s))-1
            end do
          end do
        end if
      end if
    end do
    write(fu) i_buffer(1:i)

    ! Write polyhedral cells' faces offsets
    data_size = int(Grid % n_cells * IP, SP)
    write(fu) data_size
    i = 0
    cell_offset = 0
    do c = Cells_In_Domain()
      if(Grid % new_c(c) .ne. 0) then            ! cell is in this subdomain
        if(Grid % cells_n_nodes(c) .lt. 0) then  ! found a polyhedron
          cell_offset = cell_offset + 1          ! to store number of polyfaces
          do i_fac = 1, Grid % cells_n_faces(c)  ! to store polyfaces
            s = Grid % cells_f(i_fac, c)
            n = Grid % faces_n_nodes(s)
            cell_offset = cell_offset + 1 + n    ! number of nodes and nodes
          end do
          i=i+1;  i_buffer(i) = cell_offset      ! write the current offset
        else
          i=i+1;  i_buffer(i) = -1  ! not a polyhedron, offsets are not needed
        end if
      end if
    end do
    write(fu) i_buffer(1:i)

  end if  ! if Grid % polyhedral

  !---------------!
  !   Cell data   !
  !---------------!

  ! Cell processor
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1;  i_buffer(i) = Grid % Comm % cell_proc(c)
    end if
  end do
  write(fu) i_buffer(1:i)

  ! Number of nodes
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1;  i_buffer(i) = abs(Grid % cells_n_nodes(c))
    end if
  end do
  write(fu) i_buffer(1:i)

  ! Wall distance
  data_size = int(n_cells_sub * RP, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1;  r_buffer(i) = Grid % wall_dist(c)
    end if
  end do
  write(fu) r_buffer(1:i)

  ! Cell volume
  data_size = int(n_cells_sub * RP, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1;  r_buffer(i) = Grid % vol(c)
    end if
  end do
  write(fu) r_buffer(1:i)

  ! Cell centers
  data_size = int(n_cells_sub * RP * 3, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1;  r_buffer(i) = Grid % xc(c)
      i=i+1;  r_buffer(i) = Grid % yc(c)
      i=i+1;  r_buffer(i) = Grid % zc(c)
    end if
  end do
  write(fu) r_buffer(1:i)

  ! Cell inertia
  data_size = int(n_cells_sub * RP * 6, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1;  r_buffer(i) = Grid % ixx(c)
      i=i+1;  r_buffer(i) = Grid % iyy(c)
      i=i+1;  r_buffer(i) = Grid % izz(c)
      i=i+1;  r_buffer(i) = Grid % ixy(c)
      i=i+1;  r_buffer(i) = Grid % iyz(c)
      i=i+1;  r_buffer(i) = Grid % ixz(c)
    end if
  end do
  write(fu) r_buffer(1:i)

  ! Cell porosity
  data_size = int(n_cells_sub * IP, SP)
  write(fu) data_size
  i = 0
  do c = Cells_In_Domain()
    if(Grid % new_c(c) .ne. 0) then
      i=i+1;  i_buffer(i) = Grid % por(c)
    end if
  end do
  write(fu) i_buffer(1:i)

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  !-----------------------!
  !                       !
  !   Create .pvtu file   !
  !                       !
  !-----------------------!

  ! Create it only from subdomain 1, when decomposed
  Assert(maxval(Grid % Comm % cell_proc(:)) .eq. sub(2))
  if(sub(2) > 1 .and. sub(1) .eq. 1) then

    call File % Set_Name(name_out, extension='.pvtu')
    call File % Open_For_Writing_Ascii(name_out, fu)

    ! Header
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(fu,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'

    ! This section must be present
    write(fu,'(a,a)') IN_2, '<PPoints>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp  //  &
                           ' NumberOfComponents="3"/>'
    write(fu,'(a,a)') IN_2, '</PPoints>'

    write(fu,'(a,a)') IN_2, '<PCells>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//intp//' Name="connectivity"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//intp//' Name="offsets"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//intp//' Name="types"/>'

    if(Grid % polyhedral) then
      write(fu,'(a,a)') IN_3, '<PDataArray type='//intp//' Name="faces"/>'
      write(fu,'(a,a)') IN_3, '<PDataArray type='//intp//' Name="faceoffsets"/>'
    end if

    write(fu,'(a,a)') IN_2, '</PCells>'

    ! Data section is not mandatory, but very useful
    write(fu,'(a,a)') IN_2, '<PCellData>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//intp        //  &
                            ' Name="Grid Processor [1]"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//intp        //  &
                            ' Name="Grid Number Of Nodes [1]"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp      //  &
                            ' Name="Grid Wall Distance [m]"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp      //  &
                            ' Name="Grid Cell Volume [m^3]"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp      //  &
                            ' Name="Grid Cell Centers [m]"'  //  &
                            ' NumberOfComponents="3"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//floatp      //  &
                            ' Name="Grid Cell Inertia [m^2]"'//  &
                            ' NumberOfComponents="6"/>'
    write(fu,'(a,a)') IN_3, '<PDataArray type='//intp        //  &
                            ' Name="Grid Cell Porosity Region [1]"/>'
    write(fu,'(a,a)') IN_2, '</PCellData>'

    ! Write out the names of all the pieces
    do n = 1, sub(2)
      call File % Set_Name(name_out, processor=(/n, sub(2)/), extension='.vtu')
      write(fu, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
    end do

    ! Footer
    write(fu, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(fu, '(a,a)') IN_0, '</VTKFile>'

    close(fu)

  end if

  call Profiler % Stop('Save_Vtu_Cells')

  end subroutine
