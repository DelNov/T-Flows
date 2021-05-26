!==============================================================================!
  subroutine Find_Surf_Elements(Surf, verbose)
!------------------------------------------------------------------------------!
!   Compresses sides' list                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  logical                  :: verbose
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
