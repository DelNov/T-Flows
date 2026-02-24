!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  integer, intent(in)       :: n_stat_t  ! 1st t.s. statistics turbulence
  integer, intent(in)       :: n_stat_p  ! 1st t.s. statistics particles
!----------------------------------[Locals]------------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: Grid
  integer                  :: c, eddy, dir, n, i_nod
  real                     :: lo, xo(4), zo(4),                           &
                              yo, ro, xc, yc, zc, vc, wc, sig_x, sig_yz,  &
                              rmin, rmax, sg, lx, ly, lz, vmax
  real                     :: xmin, xmax, ymin, ymax, zmin, zmax
!==============================================================================!

  ! If not time for disturbing the velocity field, return
  if(mod(Time % Curr_Dt(), 120) .ne. 0) return

  ! If too late to disturb, get out too
  if(Time % Curr_Dt() > 1440) return

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  ! Impose eddies only in precursor domain
  if(Grid % name .ne. 'PRECURSOR') return

  ! Print a message
  if(First_Proc()) then
    print *, '# Superimposing random eddies on top of velocity field!'
  end if

  !--------------------------------------!
  !   Size of the computational domain   !
  !                                      !
  !   This algorithm is not silly. We    !
  !   could have browsed through nodes   !
  !   only, but then we might have en-   !
  !   countered some hanging from GMSH   !
  !--------------------------------------!
  xmin = HUGE;  xmax = -HUGE
  ymin = HUGE;  ymax = -HUGE
  zmin = HUGE;  zmax = -HUGE
  do c = Cells_In_Domain_And_Buffers()
    do i_nod = 1, Grid % cells_n_nodes(c)
      n = Grid % cells_n(i_nod, c)
      xmin = min(xmin, Grid % xn(n));  xmax = max(xmax, Grid % xn(n))
      ymin = min(ymin, Grid % yn(n));  ymax = max(ymax, Grid % yn(n))
      zmin = min(zmin, Grid % zn(n));  zmax = max(zmax, Grid % zn(n))
    end do
  end do
  call Global % Min_Real(xmin)
  call Global % Min_Real(ymin)
  call Global % Min_Real(zmin)
  call Global % Max_Real(xmax)
  call Global % Max_Real(ymax)
  call Global % Max_Real(zmax)
  lx = xmax - xmin
  ly = ymax - ymin
  lz = zmax - zmin

  ! Minimum and maximum size of eddies
  rmin = min(ly, lz) / 10.0
  rmax = min(ly, lz) /  5.0

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
    call random_number(xo(1));  xo(1) = xmin + xo(1) * lx
    call random_number(zo(1));  zo(1) = zmin + zo(1) * lz
    call random_number(yo);     yo = ro + (ly - 2.0*ro) * yo

    ! Handle periodicity; that is: copy eddie in periodic directions
    xo(2:4) = xo(1)
    zo(2:4) = zo(1)
    if(xo(1) > lx/2.0) xo(3) = xo(1) - lx
    if(xo(1) < lx/2.0) xo(3) = xo(1) + lx
    if(zo(1) > lz/2.0) zo(2) = zo(1) - lz
    if(zo(1) < lz/2.0) zo(2) = zo(1) + lz
    xo(4) = xo(3)
    zo(4) = zo(2)

    ! Length of the eddy is six times the diameter
    lo = ro * 6.0

    sig_yz = ro / 2.0
    sig_x  = lo / 2.0

    ! Superimpose eddies on the velocity field
    do dir = 1, 4
      do c = Cells_In_Domain_And_Buffers()
        xc = Grid % xc(c)
        yc = Grid % yc(c)
        zc = Grid % zc(c)
        vc = sg * ( (zc-zo(dir))/ro )
        wc = sg * ( (yo-yc     )/ro )

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
        vc = vc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((zc-zo(dir))/sig_yz)**2)
        vc = vc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((yc-yo)     /sig_yz)**2)

        wc = wc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((zc-zo(dir))/sig_yz)**2)
        wc = wc / (sig_yz*sqrt(PI+PI))*exp(-0.5*((yc-yo)     /sig_yz)**2)

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
  do c = Cells_In_Domain_And_Buffers()
    vmax = max(vmax, abs(v % n(c)))
    vmax = max(vmax, abs(w % n(c)))
  end do
  call Global % Max_Real(vmax)
  do c = Cells_In_Domain_And_Buffers()
    v % n(c) = v % n(c) / vmax
    v % o(c) = v % o(c) / vmax
    w % n(c) = w % n(c) / vmax
    w % o(c) = w % o(c) / vmax
  end do

  end subroutine
