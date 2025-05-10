!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
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
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: Grid
  integer                  :: c, e, dir
  real                     :: lo, xo, yo, zo(2), ro, xc, yc, zc, uc, vc, sg, &
                              sig_xy, sig_z, rmin, rmax, rp, lx, ly, lz, vmax
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  ! If not time for disturbing the velocity field, return
  if(mod(Time % Curr_Dt(), 100) .ne. 0) return

  ! If too late to disturb, get out too
  if(Time % Curr_Dt() > 3000) return

  ! Print a message
  if(First_Proc()) then
    print *, '# Superimposing random eddies on top of velocity field!'
  end if

  ! Minimum and maximum size of eddies
  rmin = 0.2
  rmax = 0.4

  ! Size of the bounding box inside computational domain
  rp = 1.0  ! pipe radius
  lx = 2.0 * rp
  ly = 2.0 * rp
  lz = 4.0

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
    call random_number(ro);     ro = rmin + (rmax-rmin)*ro  ! rmin -> rmax
    call random_number(xo);     xo = xo * lx - rp           ! -rp to +rp
    call random_number(yo);     yo = yo * ly - rp           ! -rp to +rp
    call random_number(zo(1));  zo(1) = zo(1) * lz          !  0 to lz

    ! Handle periodicity; that is: copy eddie in periodic direction
    zo(2) = zo(1)
    if(zo(1) > lz/2.0) zo(2) = zo(1) - lz
    if(zo(1) < lz/2.0) zo(2) = zo(1) + lz

    ! Length of the eddy is three times its radius
    lo = ro * 3.0

    sig_xy = ro / 2.0
    sig_z  = lo / 2.0

    ! Superimpose eddies on the velocity field

    ! Work only if eddy happens to be inside the pipe
    if( sqrt(xo**2 + yo**2) < (rp-ro) ) then

      do dir = 1, 2
        do c = Cells_In_Domain_And_Buffers()
          xc = Grid % xc(c)
          yc = Grid % yc(c)
          zc = Grid % zc(c)

          ! Calculate non-axial velocity components
          uc = sg * sin( (yc-yo)/ro * PI )
          vc = sg * sin( (xc-xo)/ro * PI )

          ! Damp the velocity components in non-axial directions
          uc = uc / (sig_xy*sqrt(PI+PI))*exp(-0.5*((xc-xo)/sig_xy)**2)
          uc = uc / (sig_xy*sqrt(PI+PI))*exp(-0.5*((yc-yo)/sig_xy)**2)

          vc = vc / (sig_xy*sqrt(PI+PI))*exp(-0.5*((xc-xo)/sig_xy)**2)
          vc = vc / (sig_xy*sqrt(PI+PI))*exp(-0.5*((yc-yo)/sig_xy)**2)

          ! Damp the velocity components in axial direction
          uc = uc / (sig_z*sqrt(PI+PI))*exp(-0.5*((zc-zo(dir))/sig_z)**2)
          vc = vc / (sig_z*sqrt(PI+PI))*exp(-0.5*((zc-zo(dir))/sig_z)**2)

          u % n(c) = u % n(c) + uc
          v % n(c) = v % n(c) + vc

        end do
      end do

    end if  ! inside the pipe

  end do  ! next eddy

  vmax = 0
  do c = Cells_In_Domain_And_Buffers()
    vmax = max(vmax, abs(u % n(c)))
    vmax = max(vmax, abs(v % n(c)))
  end do
  call Global % Max_Real(vmax)
  do c = Cells_In_Domain_And_Buffers()
    u % n(c) = u % n(c) / vmax / 10.0
    v % n(c) = v % n(c) / vmax / 10.0
  end do

  end subroutine
