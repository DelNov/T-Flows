!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(grid, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Flow_Mod
  use Const_Mod, only: PI
  use Comm_Mod,  only: Comm_Mod_Global_Max_Real, this_proc
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n     ! time step
  real            :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  integer           :: c, e, dir
  real              :: lo, xo, yo, zo(2), ro, xc, yc, zc, uc, vc,  &
                       sig_xy, sig_z, rmin, rmax, sg, rp, lx, ly, lz, vmax
  character(len=80) :: file_name
!==============================================================================!

  ! If not time for disturbing the velocity field, return
  if(mod(n-33, 1500) .ne. 0) return

  ! If too late to disturb, get out too
  if(n > 6000) return

  ! Print a message
  if(this_proc < 2) then
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
        do c = 1, grid % n_cells
          xc = grid % xc(c)
          yc = grid % yc(c)
          zc = grid % zc(c)

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
  do c = 1, grid % n_cells
    vmax = max(vmax, abs(u % n(c)))
    vmax = max(vmax, abs(v % n(c)))
  end do
  call Comm_Mod_Global_Max_Real(vmax)
  do c = 1, grid % n_cells
    u % n(c) = u % n(c) / vmax / 10.0
    v % n(c) = v % n(c) / vmax / 10.0
  end do

  file_name = 'with_random_eddies_000000'
  write(file_name(20:25), '(i6.6)') n

  call Save_Results(grid, file_name)

  end subroutine
