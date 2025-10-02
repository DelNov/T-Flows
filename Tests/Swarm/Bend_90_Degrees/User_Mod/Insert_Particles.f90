!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i, j, k, c, n_parts_in_buffers
  integer                  :: n_rows, n_new, n_old
  real                     :: rx, ry, rz, r, theta, d_r, d_theta, x, z
!==============================================================================!

  ! Take alias(es)
  Grid => Flow % pnt_grid

  ! Re-set the diameters, since what you have in backup could be different
  Swarm % Particle(:) % d =  Swarm % diameter

  !----------------------------------------!
  !   Number of rows in radial direction   !
  !----------------------------------------!

  ! Number of particles = 1 + 6 + 12 + ... (n_rows - 1) * 6
  !                     = 1 + 6 * (1 + 2 + 3 + ... + n_rows)
  !                     = 1 + 3 * n_rows * (n_rows-1)
  n_rows = 15

  !-------------------------------------------!
  !                                           !
  !   1st time step of particle computation   !
  !                                           !
  !-------------------------------------------!
  n_old = Swarm % n_particles   ! old number of particles

  if(Time % Curr_Dt() .eq. 1001 .or.  &
     Time % Curr_Dt() .eq. 1501 .or.  &
     Time % Curr_Dt() .eq. 2001 .or.  &
     Time % Curr_Dt() .eq. 2501) then  ! should be after the flow is developed

    ! First Particle in the center
    k = n_old + 1
    call Swarm % Particle(k) % Insert_At(0.0, 0.0999, 0.0,    &
                                         n_parts_in_buffers,  &
                                         Flow=Flow)

    d_r = 0.009 / (n_rows - 2)

    ! Place the particles where you want them
    ! Theta loop
    do i = 1, n_rows - 1
      r = i * d_r
      d_theta = 2.0 * PI / (i * 6)
      do j = 1, i * 6
        theta = (j-1) * d_theta
        x = r * cos(theta)
        z = r * sin(theta)

        k = k + 1
        call Swarm % Particle(k) % Insert_At(x, 0.0999, z,        &
                                             n_parts_in_buffers,  &
                                             Flow=Flow)
      end do
    end do

    n_new = k  ! new number of particles

    if(n_new <= Swarm % max_particles) then
      if(First_Proc()) then
        print *, '# @User_Mod_Insert_Particles: inserted',  &
                 n_new - n_old, ' particles'
      end if
    else
      if(First_Proc()) then
        print *, '# @User_Mod_Insert_Particles: too many particles'
        call Global % End_Parallel
        stop
      end if
    end if

    !------------------------------------!
    !   Update the number of particles   !
    !------------------------------------!
    Swarm % n_particles = n_new

  end if

  end subroutine
