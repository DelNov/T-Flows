!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Grid, Flow, Turb)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  real,    contiguous, pointer :: b(:)    ! pointer to right hand side
  integer                      :: c, e, dir, seed = 123456789
  real                         :: lo, xo(4), yo(4), zo, ro,           &
                                  xc, yc, zc, vc, wc, sig_x, sig_yz,  &
                                  rmin, rmax, sg, lx, ly, lz, vmax
!==============================================================================!

  ! Disturb only once - not sure how smart or stupid is that
  if(Time % Curr_Dt() .ne. 119) return

  ! Initialize arrays
  xo(:) = 0.0
  yo(:) = 0.0
  zo    = 0.0

  ! Minimum and maximum size of eddies
  rmin = 0.2
  rmax = 0.4

  ! Size of the computational domain
  lx = 2.0 * PI
  ly = PI
  lz = 2.0

  call Gpu % Vector_Real_Copy_To_Device(xo)
  call Gpu % Vector_Real_Copy_To_Device(yo)

  ! Browse through all eddies
  do e = 1, 80

    ! Random direction of the vortex
    call Math % Random_Real(seed, sg);
    if(sg < 0.5) then
      sg = -1.0
    else
      sg = +1.0
    end if

    ! Determine random position of a vortex
    call Math % Random_Real(seed, ro);     ro    = rmin + (rmax-rmin)*ro
    call Math % Random_Real(seed, xo(1));  xo(1) = xo(1) * lx
    call Math % Random_Real(seed, yo(1));  yo(1) = yo(1) * ly
    call Math % Random_Real(seed, zo);     zo    = ro + (lz - 2.0*ro) * zo

    ! Handle periodicity; that is: copy eddie in periodic directions
    xo(2:4) = xo(1)
    yo(2:4) = yo(1)
    if(xo(1) > lx/2.0) xo(3) = xo(1) - lx
    if(xo(1) < lx/2.0) xo(3) = xo(1) + lx
    if(yo(1) > ly/2.0) yo(2) = yo(1) - ly
    if(yo(1) < ly/2.0) yo(2) = yo(1) + ly
    xo(4) = xo(3)
    yo(4) = yo(2)

    call Gpu % Vector_Update_Device(xo)
    call Gpu % Vector_Update_Device(yo)

    ! Length of the eddy is three times 
    lo = ro * 3.0

    sig_yz = ro / 2.0
    sig_x  = lo / 2.0

    ! Superimpose eddies on the velocity field
    do dir = 1, 4
      !$tf-acc loop begin
      do c = Cells_In_Domain_And_Buffers()
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

        Flow % v % n(c) = Flow % v % n(c) + vc
        Flow % w % n(c) = Flow % w % n(c) + wc
      end do
      !$tf-acc loop end
    end do
  end do

  vmax = 0
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    vmax = max(vmax, abs(Flow % v % n(c)))
    vmax = max(vmax, abs(Flow % w % n(c)))
  end do
  !$tf-acc loop end
  call Global % Max_Real(vmax)
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    Flow % v % n(c) = Flow % v % n(c) / vmax / 10.0
    Flow % w % n(c) = Flow % w % n(c) / vmax / 10.0
  end do
  !$tf-acc loop end

  call Gpu % Vector_Real_Destroy_On_Device(xo)
  call Gpu % Vector_Real_Destroy_On_Device(yo)

  end subroutine
