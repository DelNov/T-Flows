!==============================================================================!
  subroutine Save_Vtu_Front(Results, Front, time_step, domain)
!------------------------------------------------------------------------------!
!   Writes surface vertices in VTU file format (for VisIt and Paraview)        !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type)      :: Results
  type(Front_Type), target :: Front
  integer                  :: time_step
  integer,        optional :: domain
!----------------------------------[Locals]------------------------------------!
  type(Vert_Type), pointer :: Vert
  integer                  :: v, e     ! vertex and element counters
  integer                  :: offset, n, f8, f9
  character(SL)            :: name_out_8, name_out_9
!-----------------------------[Local parameters]-------------------------------!
  character(len= 0)  :: IN_0 = ''         ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  !----------------------------!
  !                            !
  !   Create .front.vtu file   !
  !                            !
  !----------------------------!

  call File % Set_Name(name_out_8,             &
                       time_step = time_step,  &
                       appendix  = '-front',   &
                       extension = '.pvtu',    &
                       domain    = domain)
  call File % Set_Name(name_out_9,             &
                       processor = this_proc,  &
                       time_step = time_step,  &
                       appendix  = '-front',   &
                       extension = '.vtu',     &
                       domain    = domain)

  if(this_proc .eq. 1) then
    call File % Open_For_Writing_Binary(name_out_8, f8)
  end if
  call File % Open_For_Writing_Ascii(name_out_9, f9)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  if(this_proc .eq. 1)  then
    write(f8) IN_0 // '<?xml version="1.0"?>'              // LF
    write(f8) IN_0 // '<VTKFile type="PUnstructuredGrid">' // LF
    write(f8) IN_1 // '<PUnstructuredGrid GhostLevel="1">' // LF
  end if

  write(f9,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(f9,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//&
                          'byte_order="LittleEndian">'
  write(f9,'(a,a)') IN_1, '<UnstructuredGrid>'

  write(f9,'(a,a,i0.1,a,i0.1,a)')   &
              IN_2, '<Piece NumberOfPoints="', Front % n_verts,  &
                         '" NumberOfCells ="', Front % n_elems, '">'

  !------------------------!
  !                        !
  !   Vertex coordinates   !
  !                        !
  !------------------------!
  if(this_proc .eq. 1)  then
    write(f8) IN_3 // '<PPoints>'                                 // LF
    write(f8) IN_4 // '<PDataArray type='//floatp//''             //  &
                      ' NumberOfComponents="3" format="ascii"/>'  // LF
    write(f8) IN_3 // '</PPoints>'                                // LF
  end if

  write(f9,'(a,a)') IN_3, '<Points>'
  write(f9,'(a,a)') IN_4, '<DataArray type='//floatp   //  &
                          ' NumberOfComponents'        //  &
                          '="3" format="ascii">'
  do v = 1, Front % n_verts
    Vert => Front % Vert(v)
    write(f9, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                IN_5, Vert % x_n, Vert % y_n, Vert % z_n
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'
  write(f9,'(a,a)') IN_3, '</Points>'

  !----------------!
  !                !
  !   Point data   !
  !                !
  !----------------!
  if(this_proc .eq. 1) then
    write(f8) IN_3 // '<PPointData>' // LF
  end if
  write(f9,'(a,a)') IN_3, '<PointData>'

  !--------------------!
  !   Particle i.d.s   !
  !--------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//intp//' Name="Index" '    //  &
                      'format="ascii"/>'                             // LF
  end if
  write(f9,'(a,a)') IN_4, '<DataArray type='//intp//' Name="Index" ' //  &
                          'format="ascii">'
  do v = 1, Front % n_verts
    write(f9,'(a,i9)') IN_5, v
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  !--------------------------!
  !   Number of neighbours   !
  !--------------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//intp//' Name="Neighbours" '  //  &
                      'format="ascii"/>'                                // LF
  end if
  write(f9,'(a,a)') IN_4, '<DataArray type='//intp//' Name="Neighbours" '  //  &
                          'format="ascii">'
  do v = 1, Front % n_verts
    write(f9,'(a,i9)') IN_5, Front % Vert(v) % nne
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  !-----------------------------!
  !   Curvatures at the nodes   !
  !-----------------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//floatp//' Name="NodeCurv" ' //  &
                      ' format="ascii"/>'                              // LF
  end if
  write(f9,'(a,a)') IN_4, '<DataArray type='//floatp//' Name="NodeCurv" ' // &
                         ' format="ascii">'
  do v = 1, Front % n_verts
    Vert => Front % Vert(v)
    write(f9,'(a,1pe16.6e4)') IN_5, Vert % curv
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  if(this_proc .eq. 1) then
    write(f8) IN_3 // '</PPointData>' // LF
  end if
  write(f9,'(a,a)') IN_3, '</PointData>'

  !-----------!
  !           !
  !   Cells   !
  !           !
  !-----------!
  if(this_proc .eq. 1) then
    write(f8) IN_3 // '<PCells>' // LF
    write(f8) IN_4 // '<PDataArray type='//intp//' Name="connectivity"'  //  &
                      ' format="ascii"/>'                                // LF
  end if
  write(f9,'(a,a)') IN_3, '<Cells>'
  write(f9,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                          ' Name="connectivity"'    //  &
                          ' format="ascii">'
  ! Cell topology
  do e = 1, Front % n_elems
    write(f9,'(a,99i9)') IN_5, Front % elem(e) % v(1:Front % elem(e) % nv)-1
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  ! Cell offsets
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//intp//' Name="offsets"'  //  &
                      ' format="ascii"/>'                           // LF
  end if
  write(f9,'(a,a)') IN_4, '<DataArray type='//intp//' Name="offsets"'  //  &
                          ' format="ascii">'
  offset = 0
  do e = 1, Front % n_elems
    offset = offset + Front % elem(e) % nv
    write(f9,'(a,i9)') IN_5, offset
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  ! Cell types
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//intp//' Name="types"'  //  &
                      ' format="ascii"/>'                         // LF
  end if
  write(f9,'(a,a)') IN_4, '<DataArray type='//intp//' Name="types"'  //  &
                          ' format="ascii">'
  do e = 1, Front % n_elems
    write(f9,'(a,i9)') IN_5, VTK_POLYGON
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  if(this_proc .eq. 1) then
    write(f8) IN_3 // '</PCells>' // LF
  end if
  write(f9,'(a,a)') IN_3, '</Cells>'

  !---------------!
  !               !
  !   Cell data   !
  !               !
  !---------------!

  ! Beginning of cell data
  if(this_proc .eq. 1) then
    write(f8) IN_3 // '<PCellData>'  // LF
  end if
  write(f9,'(a,a)') IN_3, '<CellData>'

  !-------------------------------------!
  !   Number of neighbouring elements   !
  !-------------------------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//intp//' Name="Neighbours"'  //  &
                      ' format="ascii"/>'                              // LF
  end if
  write(f9,'(a,a)') IN_4, '<DataArray type='//intp//' Name="Neighbours"'  //  &
                          ' format="ascii">'
  do e = 1, Front % n_elems
    write(f9,'(a,i9)') IN_5, Front % elem(e) % nne
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  !---------------------!
  !   Surface normals   !
  !---------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//floatp                 //  &
                      ' Name="ElementNormals" '                   //  &
                      ' NumberOfComponents="3" format="ascii"/>'  // LF
  end if
  write(f9,'(4a)') IN_4,                                                   &
                 '<DataArray type='//floatp//' Name="ElementNormals" ' //  &
                 ' NumberOfComponents="3" format="ascii">'
  do e = 1, Front % n_elems
    write(f9, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
              IN_5, Front % elem(e) % nx,           &
                    Front % elem(e) % ny,           &
                    Front % elem(e) % nz
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  !-------------------!
  !   Element areas   !
  !-------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//floatp  //  &
                      ' Name="ElementArea" '       //  &
                      ' format="ascii"/>'          // LF
  end if
  write(f9,'(4a)') IN_4,                                                 &
                 '<DataArray type='//floatp//' Name="ElementArea" '  //  &
                 ' format="ascii">'
  do e = 1, Front % n_elems
    write(f9,'(a,1pe16.6e4)') IN_5, Front % elem(e) % area
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  !-------------------------!
  !   Element coordinates   !
  !-------------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//floatp                 //  &
                      ' Name="ElementCoordinates"'                //  &
                      ' NumberOfComponents="3" format="ascii"/>'  // LF
  end if
  write(f9,'(4a)') IN_4,                            &
                 '<DataArray type='//floatp     //  &
                 ' Name="ElementCoordinates" '  //  &
                 ' NumberOfComponents="3" format="ascii">'
  do e = 1, Front % n_elems
    write(f9, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
              IN_5, Front % elem(e) % xe,           &
                    Front % elem(e) % ye,           &
                    Front % elem(e) % ze
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  !------------------------!
  !   Surface curvatures   !
  !------------------------!
  if(this_proc .eq. 1) then
    write(f8) IN_4 // '<PDataArray type='//floatp  //  &
                      ' Name="ElementCurv" '       //  &
                      ' format="ascii"/>'          // LF
  end if
  write(f9,'(4a)') IN_4,                                                &
                 '<DataArray type='//floatp//' Name="ElementCurv" ' //  &
                 ' format="ascii">'
  do e = 1, Front % n_elems
    write(f9,'(a,1pe16.6e4)') IN_5, Front % elem(e) % curv
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'

  ! End of cell data
  if(this_proc .eq. 1) then
    write(f8) IN_3 // '</PCellData>' // LF
  end if
  write(f9,'(a,a)') IN_3, '</CellData>'

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  if(this_proc .eq. 1) then
    do n = 1, n_proc
      call File % Set_Name(name_out_9,             &
                           processor = n,          &
                           time_step = time_step,  &
                           appendix  = '-front',   &
                           extension = '.vtu',     &
                           domain    = domain)
      write(f8) IN_2 // '<Piece Source="', trim(name_out_9), '"/>' // LF
    end do
    write(f8) IN_1 // '</PUnstructuredGrid>' // LF
    write(f8) IN_0 // '</VTKFile>'           // LF
    close(f8)
  end if

  write(f9,'(a,a)') IN_2, '</Piece>'
  write(f9,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(f9,'(a,a)') IN_0, '</VTKFile>'
  close(f9)

  end subroutine
