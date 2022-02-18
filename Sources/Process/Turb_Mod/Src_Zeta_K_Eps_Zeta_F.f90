!==============================================================================!
  subroutine Turb_Mod_Src_Zeta_K_Eps_Zeta_F(turb, Nat, curr_dt)
!------------------------------------------------------------------------------!
!   Calculate source terms in equation for zeta.                               !
!   Term which is negative is put on left hand side in diagonal of             !
!   martix of coefficient                                                      !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Turb_Type),   target :: turb
  type(Native_Type), target :: Nat
  integer, intent(in)       :: curr_dt
!----------------------------------[Locals]------------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: kin, eps, f22, zeta
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c
!==============================================================================!
!   In transport equation for zeta two source terms exist which have form:     !
!                                                                              !
!    /                                                                         !
!    |                     eps                                                 !
!    | ( f22 * kin * dV - ---- * dV )                                          !
!    |                    zeta                                                 !
!    /                                                                         !
!                                                                              !
!   First term can appear as positive and as negative as well so depend of     !
!   sign of term , it is placed on left or right hand side.  Second, negative  !
!   source term is added to main diagonal left hand side coefficient matrix    !
!   in order to increase stability of solver                                   !
!------------------------------------------------------------------------------!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Nat % Alias_Native(A, b)

  ! Positive source term 
  ! The first option in treating the source is making computation very 
  ! sensitive to initial condition while the second one can lead to 
  ! instabilities for some cases such as Flow around cylinder. That is why we
  ! choose this particular way to the add source term.
  do c = 1, Grid % n_cells
    if(curr_dt > 500) then
      b(c) = b(c) + f22 % n(c) * Grid % vol(c) * Flow % density(c)
    else
      b(c) = b(c) + max(0.0, f22 % n(c)*Grid % vol(c)) * Flow % density(c)
      A % val(A % dia(c)) = A % val(A % dia(c))                  &
                          + max(0.0, -f22 % n(c) * Grid % vol(c) &
                          / (zeta % n(c) + TINY)) * Flow % density(c)
    end if
    A % val(A % dia(c)) =  A % val(A % dia(c))             &
                        + Grid % vol(c) * turb % p_kin(c)  &
                        / (kin % n(c) + TINY) 
  end do

  end subroutine
