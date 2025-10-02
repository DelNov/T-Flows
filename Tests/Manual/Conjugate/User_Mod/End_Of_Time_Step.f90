!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: n_stat_t  ! start time step for Turb. stat.
  integer, intent(in)         :: n_stat_p  ! start time step for Swarm. stat.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: t
  integer                  :: c1, c2, s, reg
  real                     :: nu, area
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  call Flow % Alias_Energy(t)

  !------------------------------------!
  !   Compute average Nusselt number   !
  !------------------------------------!
  if(Grid % name(1:6) .eq. 'FLUID') then

    ! Initialize variables for computing average Nusselt number
    nu   = 0.0
    area = 0.0

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. WALL) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          area = area + Grid % s(s)
          nu   = nu + Grid % s(s)                 &
                    * abs(t % n(c2) - t % n(c1))  &
                    / Grid % d(s)
        end do  ! faces in this region
      end if    ! region is at the wall; it is a wall region
    end do      ! through regions

    !-----------------------------------------------!
    !   Integrate (summ) heated area, and heat up   !
    !-----------------------------------------------!
    call Global % Sum_Real(area)
    call Global % Sum_Real(nu)

    !-------------------------------------------------!
    !   Compute averaged Nussel number and print it   !
    !-------------------------------------------------!
    nu = nu / area

    if(First_Proc()) then
      print '(a)',        ' #==========================================='
      print '(a)',        ' # Output from user function, Nusselt number!'
      print '(a)',        ' #-------------------------------------------'
      print '(a,es12.4)', ' # Toral  area    : ', area
      print '(a,es12.4)', ' # Nusselt number : ', nu
    end if

  end if  ! domain is middle

  end subroutine
