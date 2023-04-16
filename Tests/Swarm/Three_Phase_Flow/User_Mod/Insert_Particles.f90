!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!----------------------------------[Locals]------------------------------------!
  integer :: i, j, k, p, n_parts_in_buffers, c
  real    :: x, y, z, xo, yo, zo, xn, yn, zn
  real    :: dx, dy, dz, mx, my, mz
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NI = 25, NJ = 25, NK = 2
!==============================================================================!

  !----------------------------------------------------!
  !   Initialize particles only in the 1st time step   !
  !----------------------------------------------------!
  if(Time % Curr_Dt() .eq. 1) then

    ! Leave 10% margin
    xo = -0.005
    yo = -0.005
    xn = +0.005
    yn = +0.005
    zo =  0.005
    zn =  0.015
    dx =  (xn-xo) / real(NI-1)
    dy =  (yn-yo) / real(NJ-1)
    dz =  (zn-zo) / real(NK-1)

    ! Place particles where you want them
    p = 0
    do i = 1, NI
      do j = 1, NJ
        do k = 1, NK

          p = p + 1

          ! Placing particles (only at the 1st time step)
          x = xo + (i-1) * dx
          y = yo + (j-1) * dy
          z = zo + (k-1) * dz

          mx = 0
          my = 0
          call random_number(mx);  mx = (mx - 0.5) * dx * 0.8
          call random_number(my);  my = (my - 0.5) * dy * 0.8

          call Swarm % Particle(p) % Insert_At(x+mx, y+my, z,       &
                                               n_parts_in_buffers,  &
                                               Vof=Vof)

          ! Cell is invalid
          if(Swarm % Particle(p) % cell .eq. -1) then
            print '(a,i6,a,i6,a)', ' # PANIC: Particle ', p,            &
                                   ' from processor ',    This_Proc(),  &
                                   ' couldn''t be located!'
            print *, '# Check initial placement of particles.'
            print *, '# This error is critical, exiting!'
            call Global % End_Parallel
            stop
          end if

        end do
      end do
    end do

    !--------------------------------------------------!
    !   Update number of particles in the simulation   !
    !--------------------------------------------------!
    Swarm % n_particles = p

  end if

  end subroutine
