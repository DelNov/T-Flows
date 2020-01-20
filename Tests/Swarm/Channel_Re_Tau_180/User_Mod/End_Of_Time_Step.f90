!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, turb, swarm, n_stat_p, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type,  &
                       viscosity, density, conductivity, heat_transfer
  use Var_Mod,   only: Var_Type
  use Const_Mod, only: PI
  use Comm_Mod,  only: Comm_Mod_Global_Max_Real, this_proc
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  type(Swarm_Type), target :: swarm
  integer                  :: n, n_stat_p     ! time step
  real                     :: time            ! physical time
!----------------------------------[Locals]------------------------------------!
  integer                        :: i, j, k
  real                           :: L1, L2, L3             ! domain dimensions 
  real                           :: c1, c2, c3             ! random variables 
!  integer, parameter             :: N_P     =    16        ! number of particles
  type(Var_Type),  pointer       :: u, v, w, t
  type(Grid_Type), pointer       :: grid
  integer                        :: c, eddy, dir, npb = 0
  real                           :: lo, xo(4), yo(4),                          &
                                    rx, ry, rz,                                &
                                    zo, ro, xc, yc, zc, vc, wc, sig_x, sig_yz, &
                                    rmin, rmax, sg, lx, ly, lz, vmax
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t


  !----------------------!
  !                      !
  !   Particle-related   !
  !                      !
  !----------------------!

  ! random variables
  c1 = 1.0
  c2 = 1.0
  c3 = 1.0

  ! domain size 
  L1 = 6.28 !streamwise 
  L2 = 3.14 !spanwise
  L3 = 2.0  !wall-normal

!  print *, 'End_Of_Time_Step; cnt_d = ', swarm % cnt_d

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 340010) then     ! should be after the flow is developed

    ! Initializing both deposition and departure counters
    swarm % cnt_d = 0
    swarm % cnt_e = 0
    swarm % cnt_r = 0

    ! Browsing through all introduced particles
    do k = 1, swarm % n_particles

        ! Initalizing particle position
        !swarm % particle(k) % x_n = 0.0 
        !swarm % particle(k) % y_n = 0.0 
        !swarm % particle(k) % z_n = 0.0 

        ! Generating random locations for particle
        call random_number(c1)
        call random_number(c2)
        call random_number(c3)

        ! Initalizing particle position
        swarm % particle(k) % x_n = (L1 * c1) + swarm % particle(k) % x_n
        swarm % particle(k) % y_n = (L2 * c2) + swarm % particle(k) % y_n
        swarm % particle(k) % z_n = (L3 * c3) + swarm % particle(k) % z_n

        ! you essentially moved them a lot (from 0, 0, 0)
        swarm % particle(k) % cell = 0
        swarm % particle(k) % node = 0
        swarm % particle(k) % proc = 0
        swarm % particle(k) % buff = 0

        swarm % particle(k) % x_o = swarm % particle(k) % x_n
        swarm % particle(k) % y_o = swarm % particle(k) % y_n
        swarm % particle(k) % z_o = swarm % particle(k) % z_n

        ! Searching for the closest cell and node to place the moved particle
        call Swarm_Mod_Find_Nearest_Cell(swarm, k, npb)
        call Swarm_Mod_Find_Nearest_Node(swarm, k)

        c = swarm % particle(k) % cell

        ! Set initial particle velocities
        rx = swarm % particle(k) % x_n - grid % xc(c)
        ry = swarm % particle(k) % y_n - grid % yc(c)
        rz = swarm % particle(k) % z_n - grid % zc(c)

        ! Compute velocities at the particle position from velocity gradients
        swarm % particle(k) % u    &
           = flow % u % n(c)       &  ! u velocity at the new time step (% n)
           + flow % u % x(c) * rx  &  ! u % x is gradient du/dx
           + flow % u % y(c) * ry  &  ! u % y is gradient du/dy
           + flow % u % z(c) * rz     ! u % x is gradient du/dz

        swarm % particle(k) % v    &
           = flow % v % n(c)       &  ! v velocity at the new time step (% n)
           + flow % v % x(c) * rx  &  ! v % x is gradient dv/dx
           + flow % v % y(c) * ry  &  ! v % y is gradient dv/dy
           + flow % v % z(c) * rz     ! v % x is gradient dv/dz

        swarm % particle(k) % w    &
           = flow % w % n(c)       &  ! w velocity at the new time step (% n)
           + flow % w % x(c) * rx  &  ! w % x is gradient dw/dx
           + flow % w % y(c) * ry  &  ! w % y is gradient dw/dy
           + flow % w % z(c) * rz     ! w % x is gradient dw/dz

    end do
  end if

  !----------------------!
  !   2nd time step on   !
  !----------------------!
  if(n .gt. 340010) then     ! should be started after the flow is fully developed
    call Swarm_Mod_Advance_Particles(swarm, turb, n_stat_p, n)
  end if

  if(this_proc < 2) then
    write(*,'(a,i4,a,i4,a,i4,a,i4)')                 &
             " # particles: ", swarm % n_particles,  &
             " trapped:   ",   swarm % cnt_d,        &
             " escaped:   ",   swarm % cnt_e,        &
             " reflected: ",   swarm % cnt_r
!    print *, 'particle statistics begins at: ', n_stat_p
!    stop
  end if

  !------------------------!
  !                        !
  !   Turbulence-related   !
  !                        !
  !------------------------!

  ! If not time for disturbing the velocity field, return
  if(mod(n, 120) .ne. 0) return

  ! If too late to disturb, get out too
  if(n > 1200) return

  ! Print a message
  if(this_proc < 2) then
    print *, '# Superimposing random eddies on top of velocity field!'
  end if

  ! Minimum and maximum size of eddies
  rmin = 0.2
  rmax = 0.6

  ! Size of the computational domain
  lx = 6.28
  ly = 3.14
  lz = 2.0

  !-------------------------------!
  !   Browse through all eddies   !
  !-------------------------------!
  do eddy = 1, 60

    ! Random direction of the vortex
    call random_number(sg);
    if(sg < 0.5) then
      sg = -1.0
    else
      sg = +1.0
    end if

    ! Determine random position of a vortex
    call random_number(ro);     ro    = rmin + (rmax-rmin)*ro  ! rmin -> rmax
    call random_number(xo(1));  xo(1) = xo(1) * lx
    call random_number(yo(1));  yo(1) = yo(1) * ly
    call random_number(zo);     zo = ro + (lz - 2.0*ro) * zo

    ! Handle periodicity; that is: copy eddie in periodic directions
    xo(2:4) = xo(1)
    yo(2:4) = yo(1)
    if(xo(1) > lx/2.0) xo(3) = xo(1) - lx
    if(xo(1) < lx/2.0) xo(3) = xo(1) + lx
    if(yo(1) > ly/2.0) yo(2) = yo(1) - ly
    if(yo(1) < ly/2.0) yo(2) = yo(1) + ly
    xo(4) = xo(3)
    yo(4) = yo(2)

    ! Length of the eddy is six times the diameter
    lo = ro * 6.0

    sig_yz = ro / 2.0
    sig_x  = lo / 2.0

    ! Superimpose eddies on the velocity field
    do dir = 1, 4
      do c = 1, grid % n_cells
        xc = grid % xc(c)
        yc = grid % yc(c)
        zc = grid % zc(c)
        wc = sg * ( (yc-yo(dir))/ro )
        vc = sg * ( (zo-zc     )/ro )

        !--------------------------------------------+
        !   Gaussian distribution:                   !
        !                                - (x-a)^2   !
        !                  1           ^ ---------   !
        !   f(x) = ------------------ e  2 sigma^2   !
        !          sqrt(2 PI sigma^2)                !
        !                                            !
        !                                            !
        !          exp[-(x-a)^2 / (2 sigma^2)]       !
        !   f(x) = ---------------------------       !
        !              sqrt(2 PI sigma^2)            !
        !                                            !
        !          exp[-0.5 ((x-a) / sigma)^2]       !
        !   f(x) = ---------------------------       !
        !               sigma sqrt(2 PI)             !
        !--------------------------------------------!
        vc = vc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((yc-yo(dir))/sig_yz)**2)
        vc = vc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((zc-zo)     /sig_yz)**2)

        wc = wc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((yc-yo(dir))/sig_yz)**2)
        wc = wc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((zc-zo)     /sig_yz)**2)

        vc = vc / (sig_x *sqrt(PI+PI))*exp(-0.5*((xc-xo(dir))/sig_x)**2)
        wc = wc / (sig_x *sqrt(PI+PI))*exp(-0.5*((xc-xo(dir))/sig_x)**2)

        ! Superimposed those fluctuations on spanwise and normal velocity comp.
        v % n(c) = v % n(c) + vc
        v % o(c) = v % o(c) + vc
        w % n(c) = w % n(c) + wc
        w % o(c) = w % o(c) + wc
      end do
    end do
  end do

  vmax = 0
  do c = 1, grid % n_cells
    vmax = max(vmax, abs(v % n(c)))
    vmax = max(vmax, abs(w % n(c)))
  end do
  call Comm_Mod_Global_Max_Real(vmax)
  do c = 1, grid % n_cells
    v % n(c) = v % n(c) / vmax / 5.0
    v % o(c) = v % o(c) / vmax / 5.0
    w % n(c) = w % n(c) / vmax / 5.0
    w % o(c) = w % o(c) / vmax / 5.0
  end do

  end subroutine
