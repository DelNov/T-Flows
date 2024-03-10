!==============================================================================!
  subroutine Save_Debug_Front_Vtu(Front, time_step)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to create visual representations of a front
!>  object for debugging and analysis purposes. This subroutine generates the
!>  .vtu files only for sequential runs, but for parallel runs also creates
!>  a .pvtu file.  In addition to writing vertex coordinates (which is clearly
!>  a must), the subroutine also writes indices, curvatures, and neighboring
!>  elements and curvatures for detailed examination.
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Front_Type), target :: Front      !! parent class
  integer                   :: time_step  !! current time step
!----------------------------------[Locals]------------------------------------!
  type(Vert_Type), pointer :: Vert
  type(Grid_Type), pointer :: Grid
  integer                  :: v, e, n, s, n_int, offset, fu
  integer,     allocatable :: n_elems(:)
  character(SL)            :: name_out
!==============================================================================!

  ! Take alias
  Grid => Front % pnt_grid

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  !-----------------------------!
  !                             !
  !   Create .front.pvtu file   !
  !                             !
  !-----------------------------!
  if(Parallel_Run()) then

    ! Work out how many elements are in each processor
    allocate(n_elems(N_Procs()))
    n_elems(:)         = 0
    n_elems(This_Proc()) = Front % n_elems
    call Global % Sum_Int_Array(N_Procs(), n_elems)

    if(First_Proc()) then
      call File % Set_Name(name_out,             &
                           time_step=time_step,  &
                           appendix='-front',    &
                           extension='.pvtu',    &
                           domain=Grid % rank)
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

      ! Write out the names of all the pieces
      do n = 1, N_Procs()
        if(n_elems(n) > 0) then
          call File % Set_Name(name_out,                    &
                               time_step=time_step,         &
                               appendix='-front',           &
                               processor=(/n, N_Procs()/),  &
                               extension='.vtu')
          write(fu, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out), '"/>'
        end if
      end do

      ! Footer
      write(fu, '(a,a)') IN_1, '</PUnstructuredGrid>'
      write(fu, '(a,a)') IN_0, '</VTKFile>'

      close(fu)
    end if
  end if

  if(Front % n_verts < 1) return

  !----------------------------!
  !                            !
  !   Create .front.vtu file   !
  !                            !
  !----------------------------!

  call File % Set_Name(name_out,                                &
                       processor = (/This_Proc(), N_Procs()/),  &
                       time_step = time_step,                   &
                       appendix  = '-front',                    &
                       extension = '.vtu',                      &
                       domain=Grid % rank)
  call File % Open_For_Writing_Ascii(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//&
                          'byte_order="LittleEndian">'
  write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'

  write(fu,'(a,a,i0.0,a,i0.0,a)')   &
              IN_2, '<Piece NumberOfPoints="', Front % n_verts,  &
                         '" NumberOfCells ="', Front % n_elems, '">'

  !------------------------!
  !                        !
  !   Vertex coordinates   !
  !                        !
  !------------------------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                          ' NumberOfComponents'       //  &
                          '="3" format="ascii">'
  do v = 1, Front % n_verts
    Vert => Front % Vert(v)
    write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                IN_5, Vert % x_n, Vert % y_n, Vert % z_n
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Points>'

  !----------------!
  !                !
  !   Point data   !
  !                !
  !----------------!
  write(fu,'(a,a)') IN_3, '<PointData>'

  !--------------------!
  !   Particle i.d.s   !
  !--------------------!
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="Index" ' //  &
                          'format="ascii">'
  do v = 1, Front % n_verts
    write(fu,'(a,i9)') IN_5, v
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !--------------------------!
  !   Number of neighbours   !
  !--------------------------!
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="Neighbours" '     //  &
                          'format="ascii">'
  do v = 1, Front % n_verts
    write(fu,'(a,i9)') IN_5, Front % Vert(v) % nne
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !-----------------------------!
  !   Curvatures at the nodes   !
  !-----------------------------!
  write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                          ' Name="NodeCurv" '         //  &
                         ' format="ascii">'
  do v = 1, Front % n_verts
    Vert => Front % Vert(v)
    write(fu,'(a,1pe16.6e4)') IN_5, Vert % curv
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  write(fu,'(a,a)') IN_3, '</PointData>'

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Cells>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="connectivity"'    //  &
                          ' format="ascii">'
  ! Cell topology
  do e = 1, Front % n_elems
    write(fu,'(a,99i9)') IN_5, Front % Elem(e) % v(1:Front % Elem(e) % nv)-1
  end do

  ! Cell offsets
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="offsets"' //  &
                          ' format="ascii">'
  offset = 0
  do e = 1, Front % n_elems
    offset = offset + Front % Elem(e) % nv
    write(fu,'(a,i9)') IN_5, offset
  end do

  ! Cell types
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="types"' //  &
                          ' format="ascii">'
  do e = 1, Front % n_elems
    write(fu,'(a,i9)') IN_5, VTK_POLYGON
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !               !
  !   Cell data   !
  !               !
  !---------------!

  ! Beginning of cell data
  write(fu,'(a,a)') IN_3, '<CellData>'

  !-------------------------------------!
  !   Number of neighbouring elements   !
  !-------------------------------------!
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="Neighbours"'      //  &
                          ' format="ascii">'
  do e = 1, Front % n_elems
    write(fu,'(a,i9)') IN_5, Front % Elem(e) % nne
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !---------------------!
  !   Surface normals   !
  !---------------------!
  write(fu,'(4a)') IN_4,                                                   &
                 '<DataArray type='//floatp//' Name="ElementNormals" ' //  &
                 ' NumberOfComponents="3" format="ascii">'
  do e = 1, Front % n_elems
    write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
              IN_5, Front % Elem(e) % nx,           &
                    Front % Elem(e) % ny,           &
                    Front % Elem(e) % nz
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !-------------------!
  !   Element areas   !
  !-------------------!
  write(fu,'(4a)') IN_4,                                                &
                 '<DataArray type='//floatp//' Name="ElementArea" ' //  &
                 ' format="ascii">'
  do e = 1, Front % n_elems
    write(fu,'(a,1pe16.6e4)') IN_5, Front % Elem(e) % area
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !-------------------------!
  !   Element coordinates   !
  !-------------------------!
  write(fu,'(4a)') IN_4,                           &
                 '<DataArray type='//floatp    //  &
                 ' Name="ElementCoordinates" ' //  &
                 ' NumberOfComponents="3" format="ascii">'
  do e = 1, Front % n_elems
    write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
              IN_5, Front % Elem(e) % xe,           &
                    Front % Elem(e) % ye,           &
                    Front % Elem(e) % ze
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  !------------------------!
  !   Surface curvatures   !
  !------------------------!
  write(fu,'(4a)') IN_4,                                                &
                 '<DataArray type='//floatp//' Name="ElementCurv" ' //  &
                 ' format="ascii">'
  do e = 1, Front % n_elems
    write(fu,'(a,1pe16.6e4)') IN_5, Front % Elem(e) % curv
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'

  ! End of cell data
  write(fu,'(a,a)') IN_3, '</CellData>'

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_2, '</Piece>'
  write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(fu,'(a,a)') IN_0, '</VTKFile>'
  close(fu)

  !------------------------------------!
  !                                    !
  !   Create .intersections.vtu file   !
  !                                    !
  !------------------------------------!

  ! Count number of intersections in a very simple way
  n_int = 0
  do s = 1, Grid % n_faces
    if(.not. Math % Approx_Real(Front % xs(s), 0.0) .and.  &
       .not. Math % Approx_Real(Front % ys(s), 0.0) .and.  &
       .not. Math % Approx_Real(Front % zs(s), 0.0)) then
      n_int = n_int + 1
    end if
  end do
  print *, '# Number of intersecton = ', n_int

  call File % Set_Name(name_out,                    &
                       time_step=time_step,         &
                       appendix ='-intersections',  &
                       extension='.vtu',            &
                       domain=Grid % rank)
  call File % Open_For_Writing_Ascii(name_out, fu)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//&
                          'byte_order="LittleEndian">'
  write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'

  write(fu,'(a,a,i0.0,a)')   &
              IN_2, '<Piece NumberOfPoints="', n_int,  &
                         '" NumberOfCells="0">'

  !-------------------------------!
  !                               !
  !   Intersections coordinates   !
  !                               !
  !-------------------------------!
  write(fu,'(a,a)') IN_3, '<Points>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                          ' NumberOfComponents="3"'   //  &
                          ' format="ascii">'
  do s = 1, Grid % n_faces
    if(Front % intersects_face(s)) then
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')       &
                  IN_5, Front % xs(s), Front % ys(s), Front % zs(s)
    end if
  end do
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Points>'

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------!
  write(fu,'(a,a)') IN_3, '<Cells>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="connectivity"'    //  &
                         ' format="ascii">'
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="offsets"' //  &
                         ' format="ascii">'
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="types"' //  &
                         ' format="ascii">'
  write(fu,'(a,a)') IN_4, '</DataArray>'
  write(fu,'(a,a)') IN_3, '</Cells>'

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  write(fu,'(a,a)') IN_2, '</Piece>'
  write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(fu,'(a,a)') IN_0, '</VTKFile>'
  close(fu)

  end subroutine
