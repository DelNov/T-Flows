!==============================================================================!
  subroutine Source_Zeta_K_Eps_Zeta_F(grid, n_step)
!------------------------------------------------------------------------------!
!   Calculate source terms in equation for v2                                  !
!   Term which is negative is put on left hand side in diagonal of             !
!   martix of coefficient                                                      !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n_step
!----------------------------------[Locals]------------------------------------!
  integer :: c
!==============================================================================!
!   In transport equation for zeta two source terms exist which have form:     !
!                                                                              !
!   int( f22 * kin * dV  - ( eps / zeta) * dV )                                !
!                                                                              !
!   First term can appear as positive and as negative as well so depend of     !
!   sign of term , it is placed on left or right hand side.  Second, negative  !
!   source term is added to main diagonal left hand side coefficient matrix    !
!   in order to increase stability of solver                                   !
!------------------------------------------------------------------------------!

  if(turbulence_model .eq. K_EPS_ZETA_F) then

    ! Positive source term 
    ! The first option in treating the source is making computation very 
    ! sensitive to initial condition while the second one can lead to 
    ! instabilities for some cases such as flow around cylinder. That is why we
    ! choose this particular way to the add source term.
    do c = 1, grid % n_cells
      if(n_step > 500) then
        b(c) = b(c) + f22 % n(c) * grid % vol(c) * density
      else
        b(c) = b(c) + max(0.0, f22 % n(c)*grid % vol(c)) * density
        a % val(a % dia(c)) = a % val(a % dia(c))                  &
                            + max(0.0, -f22 % n(c) * grid % vol(c) &
                            / (zeta % n(c) + TINY)) * density
      end if      
      a % val(a % dia(c)) =  a % val(a % dia(c))      &
                          + grid % vol(c) * p_kin(c)  &
                          / (kin % n(c) + TINY) 
    end do

  end if

  end subroutine
