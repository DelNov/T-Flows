!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, turb, mult, swarm, n, time)
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
  type(Field_Type),       target :: flow
  type(Turb_Type),        target :: turb
  type(Multiphase_Type),  target :: mult
  type(Swarm_Type),       target :: swarm
  integer                        :: n     ! time step
  real                           :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: grid
  integer                  :: c, eddy, dir
  real                     :: lo, xo(4), yo(4),                           &
                              zo, ro, xc, yc, zc, vc, wc, sig_x, sig_yz,  &
                              rmin, rmax, sg, lx, ly, lz, vmax
!==============================================================================!

  ! If not time for disturbing the velocity field, return
  if(mod(n, 120) .ne. 0) return

  ! If too late to disturb, get out too
  if(n > 1200) return

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

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
  do eddy = 1, 48

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
