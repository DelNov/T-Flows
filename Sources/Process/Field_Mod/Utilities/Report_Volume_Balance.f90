!==============================================================================!
  subroutine Report_Volume_Balance(Flow, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Opens file for volume balance reporting.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  type(Solver_Type), target :: Sol
  integer, intent(in)       :: curr_dt, ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  real, contiguous, pointer :: b(:)
  real                      :: src
!==============================================================================!

  if(Flow % report_vol_balance) then

    ! Take aliases
    Grid => Flow % pnt_grid
    b    => Sol % Nat % b % val

    ! Compute global volume imbalance
    src = sum(b(1:Grid % n_cells - Grid % Comm % n_buff_cells))
    call Comm_Mod_Global_Sum_Real(src)

    ! Print it out
    if(First_Proc()) then
      if(Flow % p_m_coupling == PISO) then
        write(Flow % fuvbr, '(a,i6.6,a,i3.3,a,i2.2,es13.5)')  &
                            '# ts-', curr_dt,  &
                            ' oi-', ini,      &
                            ' ii-', Flow % i_corr, src
      else
        write(Flow % fuvbr, '(a,i6.6,a,i3.3,es13.5)')  &
                            '# ts-', curr_dt,  &
                            ' oi-',  ini, src
      end if
    end if

  end if

  end subroutine
