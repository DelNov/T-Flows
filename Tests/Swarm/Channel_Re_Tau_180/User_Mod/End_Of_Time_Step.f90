!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, turb, mult, swarm, n, first_dt)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type
  use Var_Mod,   only: Var_Type
  use Const_Mod, only: PI
  use Comm_Mod,  only: Comm_Mod_Global_Max_Real, this_proc
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target     :: flow
  type(Turb_Type),       target     :: turb
  type(Multiphase_Type), target     :: mult
  type(Swarm_Type),      target     :: swarm
  integer,               intent(in) :: n         ! current time step
  integer,               intent(in) :: first_dt  ! 1st time step of this sim.
!----------------------------------[Locals]------------------------------------!
  real                           :: l1, l2, l3   ! domain dimensions
  real                           :: c1, c2, c3   ! random variables
  type(Var_Type),  pointer       :: u, v, w, t
  type(Grid_Type), pointer       :: grid
  type(Particle_Type), pointer   :: part
  character(len=80)              :: result_name
  integer                        :: c, eddy, dir, npb = 0, nn
  integer                        :: ss, oo, n_b, n_bp, fu, n_bin1, n_bin2, temp
  integer                        :: i, j, k, n_stat_p, r, s, ii, mark, n_test
  integer, allocatable           :: bin_count(:)
  real, allocatable              :: rep(:), delta(:), bin(:)
  real                           :: lo, xo(4), yo(4),                          &
                                    rx, ry, rz, Re_tau, ss0, ss1, ss00, ss11,  &
                                    zo, ro, xc, yc, zc, vc, wc, sig_x, sig_yz, &
                                    rmin, rmax, sg, lx, ly, lz, vmax, max_rep, &
                                    level
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

  ! Reading starting time for swarm statistics from control file
  call Control_Mod_Starting_Time_Step_For_Swarm_Statistics &
       (n_stat_p, verbose=.true.)

  ! Reynolds number should be passed from Save_Results and number of bins should
  ! be defined in control file, also same for n_bin... (it's okey for now!) 
  Re_tau =    142  ! operating shear Reynolds number (a little bit above what
                   ! we have so we don't lose any particle in counting!) 
  n_b    =     64  ! number of bins across half of the channel 
  n_bin1 = 342237  ! time at which we should collect bins info t+=675
  n_bin2 = 343729  ! time at which we should collect bins info t+=1125
                   ! ..particle concentration should be equivalent to t+=1000).
  n_test = 350000  ! testing (for swarm statistics)

  ! Allocating some arrays for bins
  allocate(rep(swarm % n_particles)); rep = 0.0
  allocate(bin(n_b)); bin = 0.0     ! bin distance from wall
  allocate(delta(n_b - 1)); delta = 0.0 ! bin thickness
  allocate(bin_count(swarm % n_particles)); bin_count = 0

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
  l1 = 6.28 !streamwise 
  l2 = 3.14 !spanwise
  l3 = 2.0  !wall-normal

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 100001) then     ! should be after the flow is developed

    ! Initializing both deposition and departure counters
    swarm % cnt_d = 0
    swarm % cnt_e = 0
    swarm % cnt_r = 0

    ! Browsing through all introduced particles
    do k = 1, swarm % n_particles

        ! Initalizing particle position (already initialized in
        ! Swarm_Mod_Allocate)
        swarm % particle(k) % x_n = 0.0 
        swarm % particle(k) % y_n = 0.0 
        swarm % particle(k) % z_n = 0.0 

        ! Generating random locations for particle
        call random_number(c1)
        call random_number(c2)
        call random_number(c3)

        ! Initalizing particle position
        swarm % particle(k) % x_n = (l1 * c1) + swarm % particle(k) % x_n
        swarm % particle(k) % y_n = (l2 * c2) + swarm % particle(k) % y_n
        swarm % particle(k) % z_n = (l3 * c3) + swarm % particle(k) % z_n

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
  if(n .gt. 150001) then     ! should be started after the flow is fully developed
    call Swarm_Mod_Advance_Particles(swarm, n, n_stat_p)
    if(this_proc < 2) then
      write(*,'(a,i7,a,i4,a,i4,a,i4)')                 &
               " # particles: ", swarm % n_particles,  &
               " trapped:   ",   swarm % cnt_d,        &
               " escaped:   ",   swarm % cnt_e,        &
               " reflected: ",   swarm % cnt_r
    end if
  end if

  ! Moved these five lines from Main_Pro.f90
  ! The call assumes that one uses dynamic LES model.
  ! Can we really always make such an assumption?
  if(mult % model .eq. LAGRANGIAN_PARTICLES) then
    if(swarm % subgrid_scale_model .eq. BROWNIAN_FUKAGATA) then
      call Turb_Mod_Vis_T_Dynamic(turb)
    end if
  end if

!!<<<<<<<<<<< Binning Simulation latest 26th Feb 2020 >>>>>>>>>>>>>
!  ! Chebyshev polynomials to compute slice thickness 
!  if(n .eq. 340001 .or.  n .eq. n_bin1   .or. & 
!     n .eq. n_bin2 .or.  n .eq. 344975   .or. &  ! t+ = 1000, 1500
!     n .eq. 346630 .or.  n .eq. 348290   .or. &  ! t+ = 2000, 2500
!     n .eq. 349945 .or.  n .eq. 353260   .or. &  ! t+ = 3000, 4000
!     n .eq. 360000 .or.  n .eq. 370000) then
!!  nn =  n - 340000 - 1 ! since we started having particle distribution from...
!                       ! ... n = 340001! 
!!  if(mod(nn , 2500) .eq. 0 .or. n .eq. n_bin .and. n .le. 347500) then
!    if(this_proc < 2) then
!
!      level = 0.0 
!      temp  = 0 
!      ! browse though the bins 
!      do ss = 1, n_b 
!
!        ! calculate the thickness of this bin and the following one
!        bin(ss) = 0.5 * Re_tau * (1.0 - cos(PI*(ss - 1)/(n_b - 1)))
!        ss0 = bin(ss) 
!        ss1 = (ss0 * 3.0E-5) / 0.004255 
!        bin_count(ss) = 0
!
!        ! browse through the particles 
!        !do oo = 1, swarm % n_particles
!        !  part => swarm % particle(oo)
!        !  if(part % z_n .le. ss1 .and. ss .lt. n_b) then
!        !    bin_count(ss) = bin_count(ss) + 1
!        !  end if 
!        !  rep(oo) = part % re
!        !end do
!
!        ! browse through the particles 
!        do oo = 1, swarm % n_particles
!          part => swarm % particle(oo)
!          if(part % z_n .gt. level .and. part % z_n .le. ss1) then
!            bin_count(ss) = bin_count(ss) + 1
!          end if 
!          rep(oo) = part % re
!        end do
!        level = ss1
!
!        max_rep =  maxval(rep) ! make SURE this is really the max Re_p!!
!      end do
!
!      ! DEBUGGING
!      !do i = 1, n_b
!      !  temp = temp + bin_count(i)
!      !end do 
!      !print *, 'Number of particles inside is:', temp
!      !stop 
!
!      ! calling the problem name to open a new file for binning results
!      call File_Mod_Set_Name(result_name, time_step = n,              & 
!           appendix='-swarm-concentration', extension='.dat')
!      call File_Mod_Open_File_For_Writing(result_name, fu)
!
!      ! printing info in a separate file...
!      open(fu,file=result_name)
!      write(fu,'(a1,(a12,e12.6))')  &
!      '#', 'Maximum Re_p    = ', max_rep
!      write(fu,'(a1,(a12,e12.6))')  &
!      '#', 'St+    = ', swarm % st 
!      ! columns' headers
!      write(fu,'(a1,2x,a50)') '#',   ' Bin index,'          //  &  !  1
!                                     ' Bin_y+,'             //  &  !  2
!                                     ' Number density (np)'        !  3 -  6
!      ! printing values for bins and stored particles inside
!      do ss = 1, n_b 
!        write(fu,'(i7,es15.5e3,i7)') ss,                         & !  1
!                                     bin(ss),                    & !  2
!                                     bin_count(ss)                 !  3
!      end do
!      close(fu) 
!    end if ! only one processor handles this
!  end if  ! Binning analysis loop
!
!  !mark = maxval(bin_count)
!  !print *, 'Number of particles inside is:', mark 
!  !stop 
!!>>>>>>>>>>> Binning Simulation >>>>>>>>>>>>>

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
