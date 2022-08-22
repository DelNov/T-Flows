!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n     ! time step
  real,    intent(in)      :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i, j, k, p, n_parts_in_buffers, c
  real                     :: x, y, z, xo, yo, zo, xn, yn, zn
  real                     :: dx, dy, dz, mx, my, mz
  real                     :: rx_o, ry_o, rz_o
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NI = 25, NJ = 25, NK = 2
!==============================================================================!

  Grid => Flow % pnt_grid

  !----------------------------------------------------!
  !   Initialize particles only in the 1st time step   !
  !----------------------------------------------------!
  if(n .eq. 1) then

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
          Swarm % Particle(p) % x_n = x + mx
          Swarm % Particle(p) % y_n = y + my
          Swarm % Particle(p) % z_n = z

          Swarm % Particle(p) % x_o = Swarm % Particle(p) % x_n
          Swarm % Particle(p) % y_o = Swarm % Particle(p) % y_n
          Swarm % Particle(p) % z_o = Swarm % Particle(p) % z_n

          ! Searching for the closest cell and node to place the moved particle
          call Swarm % Particle(p) % Find_Nearest_Cell(n_parts_in_buffers)
          call Swarm % Particle(p) % Find_Nearest_Node()

          ! Cell is invalid
          if(Swarm % Particle(p) % cell .eq. -1) then
            print '(a,i6,a,i6,a)', ' # PANIC: Particle ', p,          &
                                   ' from processor ',    this_proc,  &
                                   ' couldn''t be located!'
            print *, '# Check initial placement of particles.'
            print *, '# This error is critical, exiting!'
            call Comm_Mod_End
            stop
          end if

          ! Index of the closest cell for interpolation
          c = Swarm % Particle(p) % cell

          ! Vector connecting new particle position and cell center
          rx_o = Swarm % Particle(p) % x_o - Grid % xc(c)
          ry_o = Swarm % Particle(p) % y_o - Grid % yc(c)
          rz_o = Swarm % Particle(p) % z_o - Grid % zc(c)

          ! Value of smoothed vof at the old position
          Swarm % Particle(p) % smooth_o = Vof % smooth % n(c)         &
                                         + Vof % smooth % x(c) * rx_o  &
                                         + Vof % smooth % y(c) * ry_o  &
                                         + Vof % smooth % z(c) * rz_o
        end do
      end do
    end do

    !--------------------------------------------------!
    !   Update number of particles in the simulation   !
    !--------------------------------------------------!
    Swarm % n_particles = p

  end if

  end subroutine
