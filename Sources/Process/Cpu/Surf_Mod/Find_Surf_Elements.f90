!==============================================================================!
  subroutine Find_Surf_Elements(Surf)
!------------------------------------------------------------------------------!
!>  Thus subroutine establishes the relationships between elements and their
!>  adjacent sides in a surface mesh. This process involves iterating through
!>  each side and determining the adjacent elements, ensuring that each element
!>  is correctly linked to its neighbors and sides, essential for the integrity
!>  and accuracy of the surface representation.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias assignment: Sets up aliases for vertices (Vert), sides (side),     !
!     elements (Elem), and their respective counts (nv, ns, ne) for simplified !
!     reference in the code.                                                   !
!   * Connectivity establishment:                                              !
!     - Iterates through all the sides in the surface mesh.                    !
!     - For each side, identifies adjacent elements (referred to as element a  !
!       and element b) and their associated vertices.                          !
!     - Checks if the vertices of the side match with the vertices of the      !
!       adjacent elements.                                                     !
!     - If a match is found, updates the element's information with the side   !
!       and the neighboring element. This includes updating the side and       !
!       neighboring element arrays (s and e) and their respective counts (ns   !
!       and nne).                                                              !
!   * Element a and Element b processing:                                      !
!     - Handles both elements on either side of a side (element a and b)       !
!       ensuring that both elements' connectivity is established.              !
!     - For each element, checks and updates connectivity for all three        !
!       possible pairs of vertices (ensuring all edges of the element are      !
!       considered).                                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: eb, ea, c, d, s
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  Vert => Surf % Vert
  side => Surf % side
  Elem => Surf % Elem

  !-------------------------------!
  !   Find elements' neighbours   !
  !-------------------------------!
  do s = 1, ns
    c  = side(s) % c
    d  = side(s) % d
    ea = side(s) % ea
    eb = side(s) % eb

    ! Element a
    if(ea > 0) then
      if(Elem(ea) % v(2) .eq. c .and. Elem(ea) % v(3) .eq. d  .or. &
         Elem(ea) % v(3) .eq. c .and. Elem(ea) % v(2) .eq. d) then
        Elem(ea) % ns = Elem(ea) % ns + 1
        Elem(ea) % s(1) = s
        if(eb .gt. 0) then
          Elem(ea) % nne = Elem(ea) % nne + 1
          Elem(ea) % e(1) = eb
        end if
      end if

      if(Elem(ea) % v(1) .eq. c .and. Elem(ea) % v(3) .eq. d  .or. &
         Elem(ea) % v(3) .eq. c .and. Elem(ea) % v(1) .eq. d) then
        Elem(ea) % ns = Elem(ea) % ns + 1
        Elem(ea) % s(2) = s
        if(eb .gt. 0) then
          Elem(ea) % nne = Elem(ea) % nne + 1
          Elem(ea) % e(2) = eb
        end if
      end if

      if(Elem(ea) % v(1) .eq. c .and. Elem(ea) % v(2) .eq. d  .or. &
         Elem(ea) % v(2) .eq. c .and. Elem(ea) % v(1) .eq. d) then
        Elem(ea) % ns = Elem(ea) % ns + 1
        Elem(ea) % s(3) = s
        if(eb .gt. 0) then
          Elem(ea) % nne = Elem(ea) % nne + 1
          Elem(ea) % e(3) = eb
        end if
      end if
    end if  ! ea > 0

    ! Element b
    if(eb > 0) then
      if(Elem(eb) % v(2) .eq. c .and. Elem(eb) % v(3) .eq. d  .or. &
         Elem(eb) % v(3) .eq. c .and. Elem(eb) % v(2) .eq. d) then
        Elem(eb) % ns = Elem(eb) % ns + 1
        Elem(eb) % s(1) = s
        if(ea .gt. 0) then
          Elem(eb) % nne = Elem(eb) % nne + 1
          Elem(eb) % e(1) = ea
        end if
      end if

      if(Elem(eb) % v(1) .eq. c .and. Elem(eb) % v(3) .eq. d  .or. &
         Elem(eb) % v(3) .eq. c .and. Elem(eb) % v(1) .eq. d) then
        Elem(eb) % ns = Elem(eb) % ns + 1
        Elem(eb) % s(2) = s
        if(ea .gt. 0) then
          Elem(eb) % nne = Elem(eb) % nne + 1
          Elem(eb) % e(2) = ea
        end if
      end if

      if(Elem(eb) % v(1) .eq. c .and. Elem(eb) % v(2) .eq. d  .or. &
         Elem(eb) % v(2) .eq. c .and. Elem(eb) % v(1) .eq. d) then
        Elem(eb) % ns = Elem(eb) % ns + 1
        Elem(eb) % s(3) = s
        if(ea .gt. 0) then
          Elem(eb) % nne = Elem(eb) % nne + 1
          Elem(eb) % e(3) = ea
        end if
      end if
    end if

  end do

  end subroutine
