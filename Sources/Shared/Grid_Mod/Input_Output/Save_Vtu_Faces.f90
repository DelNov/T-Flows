!==============================================================================!
  subroutine Save_Vtu_Faces(Grid, sub, plot_inside, plot_shadows,  &
                                  volume_flux)
!------------------------------------------------------------------------------!
!>  Writes boundary condition .faces.vtu or shadow .shadow.vtu file.
!>  It is called primarily from Generate and Convert, and the most usefult use
!>  of the .vtu files created by this function is to check the prescribed
!>  boundary conditions, as well as periodic boundary conditions.  The optoinal
!>  argument (volume_flux) can be used to plot face-based volumetric flux
!>  (unit: [m^3/s]) and face-normal velocity (unit: [m/s])
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine handles both regular faces (boundary conditions) and      !
!     shadow faces (if plot_shadows is true).                                  !
!   * The subroutine opens a .vtu file for writing and constructs the file     !
!     content, including headers, node data, face data, and additional         !
!     face-related information.                                                !
!   * The face data includes the face's nodes, offsets, types, and other       !
!     optional variables like real and integer face variables.                 !
!   * If shadow faces are to be plotted, the subroutine adjusts its            !
!     processing accordingly.                                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid          !! grid being processes
  integer, intent(in) :: sub(1:2)      !! 1=subdomaintotal number of subdomains
  logical,   optional :: plot_inside   !! plot faces inside the domain
  logical,   optional :: plot_shadows  !! plot shadow faces
  real,      optional :: volume_flux(1:Grid % n_faces)  !! volumetric flux
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)          :: data_size
  integer              :: c2, n, s, s_f, s_l, cell_offset, data_offset, n_conns
  integer              :: fu, i, j
  character(SL)        :: name_out, ext, str1, str2
  real                 :: mag
  real,    allocatable :: r_buffer(:)
  integer, allocatable :: i_buffer(:)
!==============================================================================!

  call Profiler % Start('Save_Vtu_Faces')

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  ! Starting and ending counters; file extension,
  ! this is the default behavior: plot all faces
  s_f = 1
  s_l = Grid % n_faces
  ext = '.faces.vtu'

  ! Allocate memory for local buffers
  allocate(r_buffer(max(Grid % n_nodes, Grid % n_faces) * 3))
  n = 0
  do s = s_f, s_l
    n = n + Grid % faces_n_nodes(s)
  end do
  allocate(i_buffer(max(grid % n_nodes * 3, grid % n_faces * 3, n)))

  ! Fix counters and file extension if you are plotting the boundary only
  if(present(plot_inside)) then
    if( .not. plot_inside ) then
      s_f = 1
      s_l = Grid % n_bnd_cells
      ext = '.bnd.faces.vtu'
    end if
  end if

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
  call File % Set_Name(name_out,         &
                       processor = sub,  &
                       extension = trim(ext))
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
  write(fu) IN_4 // '<DataArray type='//floatp       //  &
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
  write(fu) IN_4 // '<DataArray type='//intp         //  &
                    ' Name="connectivity"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP             ! prepare for next

  ! Faces' offsets
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp         //  &
                    ' Name="offsets"'                //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Faces' types
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp         //  &
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
  write(fu) IN_3 // '<CellData>' // LF

  ! Boundary conditions
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp           //  &
                    ' Name="Grid Boundary Conditions"' //  &
                    ' format="appended"'               //  &
                    ' offset="' // trim(str1) //'">'   // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Optional volumetric flux
  if(present(volume_flux)) then
    write(str1, '(i0.0)') data_offset
    write(fu) IN_4 // '<DataArray type='//floatp        //  &
                      ' Name="Volumetric Flux [m^3/s]"' //  &
                      ' format="appended"'              //  &
                      ' offset="' // trim(str1) //'">'  // LF
    write(fu) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + (s_l-s_f+1) * RP      ! prepare for next
  end if

  ! Number of nodes
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//intp         //  &
                    ' Name="Grid Number Of Nodes"'   //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * IP      ! prepare for next

  ! Surface vectors
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//floatp       //  &
                    ' Name="Grid Surface Vectors"'   //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

  ! Surface normals
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//floatp       //  &
                    ' Name="Grid Surface Normals"'   //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(fu) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + (s_l-s_f+1) * RP * 3  ! prepare for next

  ! Connection vectors
  write(str1, '(i0.0)') data_offset
  write(fu) IN_4 // '<DataArray type='//floatp         //  &
                    ' Name="Grid Connection Vectors"'  //  &
                    ' NumberOfComponents="3"'          //  &
                    ' format="appended"'               //  &
                    ' offset="' // trim(str1) //'">'   // LF
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
  i = 0
  do n = 1, Grid % n_nodes
    i=i+1;  r_buffer(i) = Grid % xn(n)
    i=i+1;  r_buffer(i) = Grid % yn(n)
    i=i+1;  r_buffer(i) = Grid % zn(n)
  end do
  write(fu) r_buffer(1:i)

  !-----------!
  !   Faces   !
  !-----------!

  ! Faces' nodes
  data_size = int(n_conns * IP, SP)
  write(fu) data_size
  i = 0
  do s = s_f, s_l
    n = Grid % faces_n_nodes(s)
    do j = 1, n
      i=i+1;  i_buffer(i) = Grid % faces_n(j,s)-1
    end do
  end do
  write(fu) i_buffer(1:i)

  ! Faces' offsets
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  i = 0
  cell_offset = 0
  do s = s_f, s_l
    cell_offset = cell_offset + Grid % faces_n_nodes(s)
    i=i+1;  i_buffer(i) = cell_offset
  end do
  write(fu) i_buffer(1:i)

  ! Faces' types
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  i = 0
  do s = s_f, s_l
    if(Grid % faces_n_nodes(s) .eq. 4) then
      i=i+1;  i_buffer(i) = VTK_QUAD
    else if(Grid % faces_n_nodes(s) .eq. 3) then
      i=i+1;  i_buffer(i) = VTK_TRIANGLE
    else
      i=i+1;  i_buffer(i) = VTK_POLYGON
    end if
  end do
  write(fu) i_buffer(1:i)

  ! Boundary conditions
  ! (Check c1 and c2 for shadow faces, seems to be something messed up)
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  i = 0
  do s = s_f, s_l
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      i=i+1;  i_buffer(i) = Grid % region % at_cell(c2)
    else
      i=i+1;  i_buffer(i) = 0
    end if
  end do
  write(fu) i_buffer(1:i)

  call Profiler % Start('Save_Vtu_Faces (optional real - optimize!)')

  ! Optional real face variable
  ! (Check c1 and c2 for shadow faces, seems to be something messed up)
  if(present(volume_flux)) then
    data_size = int((s_l-s_f+1) * RP, SP)
    write(fu) data_size
    do s = s_f, s_l
      write(fu) volume_flux(s)
    end do
  end if

  call Profiler % Stop('Save_Vtu_Faces (optional real - optimize!)')

  ! Number of nodes
  data_size = int((s_l-s_f+1) * IP, SP)
  write(fu) data_size
  i = 0
  do s = s_f, s_l
    i=i+1;  i_buffer(i) = Grid % faces_n_nodes(s)
  end do
  write(fu) i_buffer(1:i)

  ! Surface vectors
  data_size = int((s_l-s_f+1) * RP * 3, SP)
  write(fu) data_size
  i = 0
  do s = s_f, s_l
    i=i+1;  r_buffer(i) = Grid % sx(s)
    i=i+1;  r_buffer(i) = Grid % sy(s)
    i=i+1;  r_buffer(i) = Grid % sz(s)
  end do
  write(fu) r_buffer(1:i)

  ! Surface normals
  data_size = int((s_l-s_f+1) * RP * 3, SP)
  write(fu) data_size
  i = 0
  do s = s_f, s_l
    mag = sqrt(Grid % sx(s)**2 + Grid % sy(s)**2 + Grid % sz(s)**2)
    i=i+1;  r_buffer(i) = Grid % sx(s) / (mag+TINY)
    i=i+1;  r_buffer(i) = Grid % sy(s) / (mag+TINY)
    i=i+1;  r_buffer(i) = Grid % sz(s) / (mag+TINY)
  end do
  write(fu) r_buffer(1:i)

  ! Connection vectors
  data_size = int((s_l-s_f+1) * RP * 3, SP)
  write(fu) data_size
  i = 0
  do s = s_f, s_l
    i=i+1;  r_buffer(i) = Grid % dx(s)
    i=i+1;  r_buffer(i) = Grid % dy(s)
    i=i+1;  r_buffer(i) = Grid % dz(s)
  end do
  write(fu) r_buffer(1:i)

  write(fu) LF // IN_0 // '</AppendedData>' // LF
  write(fu) IN_0 // '</VTKFile>' // LF

  close(fu)

  call Profiler % Stop('Save_Vtu_Faces')

  end subroutine
