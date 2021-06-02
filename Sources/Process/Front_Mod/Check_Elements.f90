!==============================================================================!
  subroutine Check_Elements(Front)
!------------------------------------------------------------------------------!
!   Finds connectivity for sides and elements                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
!-----------------------------------[Locals]-----------------------------------!
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: ne
  integer                  :: e, sum_ijk, sum_cd, i_ver, i_s
!==============================================================================!

  ! Take aliases
  ne   => Front % n_elems
  side => Front % side
  Elem => Front % Elem

  ! Checking
  do e = 1, ne
    sum_ijk = 0
    sum_cd  = 0

    do i_ver = 1, Elem(e) % nv
      sum_ijk = sum_ijk + Elem(e) % v(i_ver)
    end do

    do i_s = 1, Elem(e) % ns
      sum_cd = sum_cd + side(Elem(e) % s(i_s)) % c  &
                      + side(Elem(e) % s(i_s)) % d
    end do

    if( sum_cd / sum_ijk .ne. 2 ) then
      print *, '# ERROR in forming elements'' neighbours!'
      stop
    end if
  end do

  end subroutine
