!==============================================================================!
  subroutine Save_Vtu_Faces(Grid, plot_shadows, phi_f)
!------------------------------------------------------------------------------!
!   Writes boundary condition .faces.vtu or shadow .shadow.vtu file.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)  :: Grid
  logical, optional :: plot_shadows  ! plot shadow faces
  real,    optional :: phi_f(1:Grid % n_faces)
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)   :: data_size
  integer       :: c2, n, s, s_f, s_l, cell_offset, data_offset, n_conns, fu
  character(SL) :: name_out, ext, str1, str2
  real          :: mag
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: IP = DP  ! int. precision is double precision
  integer, parameter :: RP = DP  ! real precision is double precision
!==============================================================================!

  ! Starting and ending counters; file extension
  s_f = 1
  s_l = Grid % n_faces
  ext = '.faces.vtu'

  ! Fix counters and file extension if you are plotting shadows
  ! (Keep in mind that minval and maxval return HUGE if no mask is matched)
  if(present(plot_shadows)) then
    if( plot_shadows ) then
      if( any(Grid % faces_s(1:Grid % n_faces) .ne. 0) ) then
        s_f = minval(Grid % faces_s(1:Grid % n_faces),  &
                mask=Grid % faces_s(1:Grid % n_faces) .ne. 0)
        s_l = maxval(Grid % faces_s(1:Grid % n_faces),  &
                mask=Grid % faces_s(1:Grid % n_faces) .ne. 0)
        ext = '.shadows.vtu'
      else
        print *, '# NOTE: No shadow faces in this domain, nothing to plot!'
        return
      end if
    end if
  end if

  ! Count connections in this subdomain, you will need it later
  n_conns = 0
  do s = s_f, s_l
    n_conns = n_conns + Grid % faces_n_nodes(s)
  end do

  !------------------------!
  !   Open the .vtu file   !
  !------------------------!
  call File % Set_Name(name_out, processor=this_proc, extension=trim(ext))
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
  write(str1, '(i0.0)') Grid % n_nodes
  write(str2, '(i0.0)') (s_l-s_f+1)
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
  data_offset = data_offset + SP + Grid % n_nodes * RP * 3  ! prepare for next

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
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Faces' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="types"'                  //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

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
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Optional face variable
  if(present(phi_f)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                      ' Name="FaceVariable"'           //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + (s_l-s_f+1) * RP      ! prepare for next
  end if

  ! Number of nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="GridNumberOfNodes"'      //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Surface vectors
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridSurfaceVectors"'     //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

  ! Surface normals
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridSurfaceNormals"'     //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

  ! Connection vectors
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type="Float64"'      //  &
                    ' Name="GridConnectionVectors"'  //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

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
  !   Faces   !
  !-----------!

  ! Faces' nodes
  data_size = int(n_conns * IP, SP)
  write(fu) data_size
  do s = s_f, s_l
    n = Grid % faces_n_nodes(s)
    write(fu) Grid % faces_n(1:n,s)-1
  end do

  ! Faces' offsets
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  cell_offset = 0
  do s = s_f, s_l
    cell_offset = cell_offset + Grid % faces_n_nodes(s)
    write(fu) cell_offset
  end do

  ! Faces' types
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  do s = s_f, s_l
    if(Grid % faces_n_nodes(s) .eq. 4) then
      write(fu) VTK_QUAD
    else if(Grid % faces_n_nodes(s) .eq. 3) then
      write(fu) VTK_TRIANGLE
    else
      write(fu) VTK_POLYGON
    end if
  end do

  ! Boundary conditions
  ! (Check c1 and c2 for shadow faces, seems to be something messed up)
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  do s = s_f, s_l
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      write(fu) Grid % bnd_cond % color(c2)
    else
      write(fu) 0
    end if
  end do

  ! Optional face variable
  ! (Check c1 and c2 for shadow faces, seems to be something messed up)
  if(present(phi_f)) then
    data_size = int((s_l-s_f+1) * IP, SP)
    write(fu) data_size
    do s = s_f, s_l
      write(fu) phi_f(s)
    end do
  end if

  ! Number of nodes
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  do s = s_f, s_l
    write(fu) Grid % faces_n_nodes(s)
  end do

  ! Surface vectors
  data_size = int((s_l-s_f+1) * RP * 3, SP)
  write(fu) data_size
  do s = s_f, s_l
    write(fu) Grid % sx(s), Grid % sy(s), Grid % sz(s)
  end do

  ! Surface normals
  data_size = int((s_l-s_f+1) * RP * 3, SP)
  write(fu) data_size
  do s = s_f, s_l
    mag = sqrt(Grid % sx(s)**2 + Grid % sy(s)**2 + Grid % sz(s)**2)
    write(fu) Grid % sx(s) / mag, Grid % sy(s) / mag, Grid % sz(s) / mag
  end do

  ! Connection vectors
  data_size = int((s_l-s_f+1) * RP * 3, SP)
  write(fu) data_size
  do s = s_f, s_l
    write(fu) Grid % dx(s), Grid % dy(s), Grid % dz(s)
  end do

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  end subroutine
