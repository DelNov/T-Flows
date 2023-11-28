!==============================================================================!
  subroutine Save_Vtu_Surf(Results, Surf)
!------------------------------------------------------------------------------!
!>  The Save_Vtu_Surf subroutine is responsible for writing surface vertices
!>  in .vtu file format, compatible with visualization tools like VisIt and
!>  ParaView.  This subroutine primarily deals with surfacing structures
!>  (most often created as interfaces between two phases in VOF), capturing
!>  their  geometric details and properties for visualization purposes.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initial setup: Ensures the necessary precision for plotting is set and   !
!     checks if the number of vertices to be processed is sufficient.          !
!   * File creation: A .surf.vtu file is created, specifically for storing     !
!     surface-related data.                                                    !
!   * XML header: Constructs the necessary XML structure for the VTU file,     !
!     including the number of vertices and elements involved.                  !
!   * Vertex data: Iteratively writes the coordinates of each vertex, along    !
!     with other vertex-specific information like the number of neighbors and  !
!     curvature.                                                               !
!   * Element data: Writes the connectivity of elements, their types, and      !
!     other attributes like neighboring element count, surface normals, and    !
!     curvatures.                                                              !
!   * Finalization: Closes the XML structure and the file, completing the      !
!     process of saving the surface data.                                      !
!------------------------------------------------------------------------------!
!   Workflow                                                                   !
!                                                                              !
!   * File setup: The subroutine begins by preparing the file name and opening !
!     the file for writing, restricted to the first processor.                 !
!   * Writing header: It then writes the XML header, specifying the VTU format !
!     and other metadata.                                                      !
!   * Processing vertices: For each vertex in the surface structure, its       !
!     coordinates and additional data (e.g., curvature) are written.           !
!   * Processing elements: The subroutine then iterates over each element in   !
!     the surface structure, writing data about its connectivity, type, and    !
!     other attributes.                                                        !
!   * Completing data write: After finishing with vertices and elements, the   !
!     subroutine writes the end tags for the XML structure and closes the file.!
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type),     intent(in) :: Results  !! parent class
  type(Surf_Type), target, intent(in) :: Surf     !! surface object to save
!----------------------------------[Locals]------------------------------------!
  type(Vert_Type), pointer :: Vert
  integer                  :: v, e     ! vertex and element counters
  integer                  :: offset, fu
  character(SL)            :: name_out
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  if(Surf % n_verts < 1) return

  !---------------------------!
  !                           !
  !   Create .surf.vtu file   !
  !                           !
  !---------------------------!

  if(First_Proc()) then

    call File % Set_Name(name_out,                      &
                         time_step = Time % Curr_Dt(),  &
                         appendix  = '-surf',           &
                         extension = '.vtu')
    call File % Open_For_Writing_Ascii(name_out, fu)

    !------------!
    !            !
    !   Header   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid"'  //  &
                            ' version="0.1"'                    //  &
                            ' byte_order="LittleEndian">'
    write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'

    write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                IN_2, '<Piece NumberOfPoints="', Surf % n_verts,  &
                           '" NumberOfCells ="', Surf % n_elems, '">'

    !------------------------!
    !                        !
    !   Vertex coordinates   !
    !                        !
    !------------------------!
    write(fu,'(a,a)') IN_3, '<Points>'
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                            ' NumberOfComponents'       //  &
                            '="3" format="ascii">'
    do v = 1, Surf % n_verts
      Vert => Surf % Vert(v)
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
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                            ' Name="Index" '          //  &
                            'format="ascii">'
    do v = 1, Surf % n_verts
      write(fu,'(a,i9)') IN_5, v
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !--------------------------!
    !   Number of neighbours   !
    !--------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                            ' Name="Neighbours" '     // &
                            'format="ascii">'
    do v = 1, Surf % n_verts
      write(fu,'(a,i9)') IN_5, Surf % Vert(v) % nne
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-----------------------------!
    !   Curvatures at the nodes   !
    !-----------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                            ' Name="NodeCurv" '         //  &
                            ' format="ascii">'
    do v = 1, Surf % n_verts
      Vert => Surf % Vert(v)
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
    do e = 1, Surf % n_elems
      write(fu,'(a,3i9)')          &
         IN_5,                     &
         Surf % Elem(e) % v(1)-1,  &
         Surf % Elem(e) % v(2)-1,  &
         Surf % Elem(e) % v(3)-1
    end do

    ! Cell offsets
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="offsets"' //  &
                            ' format="ascii">'
    offset = 0
    do e = 1, Surf % n_elems
      offset = offset + 3
      write(fu,'(a,i9)') IN_5, offset
    end do

    ! Cell types
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="types"' //  &
                            ' format="ascii">'
    do e = 1, Surf % n_elems
      write(fu,'(a,i9)') IN_5, VTK_TRIANGLE
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
    do e = 1, Surf % n_elems
      write(fu,'(a,i9)') IN_5, Surf % Elem(e) % nne
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !---------------------!
    !   Surface normals   !
    !---------------------!
    write(fu,'(4a)') IN_4,                                                   &
                   '<DataArray type='//floatp//' Name="ElementNormals" ' //  &
                   ' NumberOfComponents="3" format="ascii">'
    do e = 1, Surf % n_elems
      write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
            IN_5, Surf % Elem(e) % nx, Surf % Elem(e) % ny, Surf % Elem(e) % nz
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !------------------------!
    !   Surface curvatures   !
    !------------------------!
    write(fu,'(4a)') IN_4,                                                &
                   '<DataArray type='//floatp//' Name="ElementCurv" ' //  &
                   ' format="ascii">'
    do e = 1, Surf % n_elems
      write(fu,'(a,1pe16.6e4)') IN_5, Surf % Elem(e) % curv
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
  end if

  end subroutine
