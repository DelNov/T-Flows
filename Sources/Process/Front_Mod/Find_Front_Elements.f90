!==============================================================================!
  subroutine Find_Front_Elements(Front)
!------------------------------------------------------------------------------!
!>  This subroutine establishes connectivity between sides and elements of a
!>  front object.  This process ensures that each element in the front object
!>  is correctly linked to its neighboring elements and sides, crucial for
!>  accurate front representation and analysis.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias assignment: Sets up aliases for vertices (Vert), sides (side),     !
!     elements (Elem), and their respective counts (nv, ns, ne) for simplified !
!     reference in the code.                                                   !
!   * Iterating through sides:                                                 !
!     - The subroutine iterates through each side in the front object. For     !
!       each side, it identifies the adjacent elements (referred to as         !
!       elements a and b) and their associated vertices.                       !
!   * Connectivity establishment:                                              !
!     - For each side, the subroutine determines the elements on either side   !
!       of it. It then identifies the corresponding vertices of these elements.!
!     - The subroutine checks if the vertices of these elements match the      !
!       side's vertices.                                                       !
!   * Updating element information:                                            !
!     - Upon finding a match, the subroutine updates the element's information !
!       with the identified side and the neighboring elements. This includes   !
!       updating the side and neighboring element arrays (s and e) and their   !
!       respective counts (ns and nne).                                        !
!   * Processing for both sides of each element:                               !
!     - The subroutine processes both elements adjacent to a side (elements a  !
!       and b), ensuring comprehensive establishment of connectivity.          !
!     - For each element, it checks and updates connectivity for all pairs of  !
!       vertices. This ensures that all edges of each element are considered   !
!       in the connectivity establishment.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: eb, ea, c, d, s, i_ver, v1, v2
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ns   => Front % n_sides
  ne   => Front % n_elems
  Vert => Front % Vert
  side => Front % side
  Elem => Front % Elem

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

      do i_ver = 1, Elem(ea) % nv

        ! Get first and second vertex
        v1 = Elem(ea) % v(i_ver)
        if(i_ver < Elem(ea) % nv) then
          v2 = Elem(ea) % v(i_ver+1)
        else
          v2 = Elem(ea) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          Elem(ea) % ns  = Elem(ea) % ns  + 1
          Elem(ea) % s(Elem(ea) % ns)  = s
          if(eb .gt. 0) then
            Elem(ea) % nne = Elem(ea) % nne + 1
            Elem(ea) % e(Elem(ea) % nne) = eb
          end if
        end if
      end do

    end if  ! ea > 0

    ! Element b
    if(eb > 0) then

      do i_ver = 1, Elem(eb) % nv

        ! Get first and second vertex
        v1 = Elem(eb) % v(i_ver)
        if(i_ver < Elem(eb) % nv) then
          v2 = Elem(eb) % v(i_ver+1)
        else
          v2 = Elem(eb) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          Elem(eb) % ns  = Elem(eb) % ns  + 1
          Elem(eb) % s(Elem(eb) % ns)  = s
          if(ea .gt. 0) then
            Elem(eb) % nne = Elem(eb) % nne + 1
            Elem(eb) % e(Elem(eb) % nne) = ea
          end if
        end if
      end do

    end if  ! eb > 0

  end do

  end subroutine
