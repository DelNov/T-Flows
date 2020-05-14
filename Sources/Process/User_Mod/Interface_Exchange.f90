!==============================================================================!
  subroutine User_Mod_Interface_Exchange(inter, flow, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)        :: inter(MD, MD)
  type(Field_Type),    target :: flow(MD)
  integer                     :: n_dom
!-----------------------------------[Locals]-----------------------------------!
  integer :: n1, n2
!==============================================================================!

  do d1 = 1, n_dom
    do d2 = 1, n_dom
      call Interface_Mod_Exchange(inter,         &
                                  flow(d1) % t,  &
                                  flow(d2) % t,  &
                                  d1, d2)
    end do
  end do

  end subroutine
