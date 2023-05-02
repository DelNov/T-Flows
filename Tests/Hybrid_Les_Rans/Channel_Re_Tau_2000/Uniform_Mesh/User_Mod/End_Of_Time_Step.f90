!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n_stat_t  ! 1st step for turbulence statist.
  integer, intent(in)      :: n_stat_p  ! 1st step for particle statistics
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c, e, dir
  real                     :: lo, xo(4), yo(4), zo, ro,           &
                              xc, yc, zc, vc, wc, sig_x, sig_yz,  &
                              rmin, rmax, sg, lx, ly, lz, vmax
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w

  ! If not time for disturbing the velocity field, return
  if(mod(Time % Curr_Dt(), 100) .ne. 0) return

  ! If too late to disturb, get out too
  if(Time % Curr_Dt() > 1000) return

  ! Minimum and maximum size of eddies
  rmin = 0.2
  rmax = 0.4

  ! Size of the computational domain
  lx = 2.0 * PI
  ly = PI
  lz = 2.0

  ! Browse through all eddies
  do e = 1, 80

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

    ! Length of the eddy is three times 
    lo = ro * 3.0

    sig_yz = ro / 2.0
    sig_x  = lo / 2.0

    ! Superimpose eddies on the velocity field
    do dir = 1, 4
      do c = 1, Grid % n_cells
        xc = Grid % xc(c)
        yc = Grid % yc(c)
        zc = Grid % zc(c)
        wc = sg * sin( (yc-yo(dir))/ro * PI )
        vc = sg * sin( (zc-zo)/ro * PI )

        vc = vc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((yc-yo(dir))/sig_yz)**2)
        vc = vc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((zc-zo)     /sig_yz)**2)

        wc = wc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((yc-yo(dir))/sig_yz)**2)
        wc = wc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((zc-zo)     /sig_yz)**2)

        vc = vc / (sig_x *sqrt(PI+PI))*exp(-0.5*((xc-xo(dir))/sig_x)**2)
        wc = wc / (sig_x *sqrt(PI+PI))*exp(-0.5*((xc-xo(dir))/sig_x)**2)

        v % n(c) = v % n(c) + vc
        w % n(c) = w % n(c) + wc
      end do
    end do
  end do

  vmax = 0
  do c = 1, Grid % n_cells
    vmax = max(vmax, abs(v % n(c)))
    vmax = max(vmax, abs(w % n(c)))
  end do
  call Global % Max_Real(vmax)
  do c = 1, Grid % n_cells
    v % n(c) = v % n(c) / vmax / 10.0
    w % n(c) = w % n(c) / vmax / 10.0
  end do

  end subroutine
