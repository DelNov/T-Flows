!==============================================================================!
  subroutine Save_Vtu_Faces(grid)
!------------------------------------------------------------------------------!
!   Writes boundary condition .faces.vtu file.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)   :: data_size
  integer       :: c2, n, s, cell_offset, data_offset, n_conns, fu
  character(SL) :: name_out, str1, str2
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: IP = DP  ! int. precision is double precision
  integer, parameter :: RP = DP  ! real precision is double precision
!==============================================================================!

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do s = 1, grid % n_faces
    n_conns = n_conns + grid % faces_n_nodes(s)
  end do

  !------------------------!
  !   Open the .vtu file   !
  !------------------------!
  call File_Mod_Set_Name(name_out, extension='.faces.vtu')
  call File_Mod_Open_File_For_Writing_Binary(name_out, fu)

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
  write(str1, '(i0.0)') grid % n_nodes
  write(str2, '(i0.0)') grid % n_faces
  write(fu) IN_2 // '<Piece NumberOfPoints="' // trim(str1) // '"'  //  &
                    ' NumberOfCells="'        // trim(str2) // '">' // LF
  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  write(str1, '(i1)') data_offset
  write(fu) IN_3 // '<Points>'                       // LF
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  write(fu) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + grid % n_nodes * RP * 3  ! prepare for next

  !-----------!
  !   Faces   !
  !-----------!
  write(fu) IN_3 // '<Cells>' // LF

  ! Faces' nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="connectivity"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP             ! prepare for next

  ! Faces' offsets
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="offsets"'                //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + grid % n_faces * IP      ! prepare for next

  ! Faces' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="types"'                  //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + grid % n_faces * IP      ! prepare for next

  !----------------------!
  !   The end of faces   !
  !----------------------!
  write(fu) IN_3 // '</Cells>' // LF

  !---------------!
  !   Face data   !
  !---------------!
  write(fu) IN_3 // '<CellData Scalars="scalars" vectors="velocity">' // LF

  ! Boundary conditions
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="BoundaryConditions"'     //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + grid % n_faces * IP      ! prepare for next

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
  data_size = grid % n_nodes * RP * 3
  write(fu) data_size
  do n = 1, grid % n_nodes
    write(fu) grid % xn(n), grid % yn(n), grid % zn(n)
  end do

  !-----------!
  !   Faces   !
  !-----------!

  ! Faces' nodes
  data_size = n_conns * IP
  write(fu) data_size
  do s = 1, grid % n_faces
    n = grid % faces_n_nodes(s)
    write(fu) grid % faces_n(1:n,s)-1
  end do

  ! Faces' offsets
  data_size = grid % n_faces * IP
  write(fu) data_size
  cell_offset = 0
  do s = 1, grid % n_faces
    cell_offset = cell_offset + grid % faces_n_nodes(s)
    write(fu) cell_offset
  end do

  ! Faces' types
  data_size = grid % n_faces * IP
  write(fu) data_size
  do s = 1, grid % n_faces
    if(grid % faces_n_nodes(s) .eq. 4) then
      write(fu) VTK_QUAD
    else if(grid % faces_n_nodes(s) .eq. 3) then
      write(fu) VTK_TRIANGLE
    else
      write(fu) VTK_POLYGON
    end if
  end do

  ! Boundary conditions
  data_size = grid % n_faces * IP
  write(fu) data_size
  do s = 1, grid % n_faces
    c2 = grid % faces_c(2,s)

    if(c2 < 0) write(fu) grid % bnd_cond % color(c2)
    if(c2 > 0) write(fu) 0
  end do

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  end subroutine
