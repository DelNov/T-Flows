!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target     :: Flow
  type(Turb_Type),  target     :: Turb
  type(Vof_Type),   target     :: Vof
  type(Swarm_Type), target     :: Swarm
  integer                      :: n_stat_t  ! 1st t.s. for turb. stat.
  integer                      :: n_stat_p  ! 1st t.s. for swarm. stat.
!------------------------------[Local parameters]------------------------------!
  real, parameter :: LX = 6.28  ! streamwise
  real, parameter :: LY = 3.14  ! spanwise
  real, parameter :: LZ = 2.0   ! wall-normal
!----------------------------------[Locals]------------------------------------!
  type(Var_Type),      pointer :: u, v, w, t
  type(Grid_Type),     pointer :: Grid
  type(Particle_Type), pointer :: part
  character(SL)                :: result_name
  integer                      :: c, eddy, dir, npb = 0, nn
  integer                      :: ss, oo, n_b, n_bp, fu, n_bin1, n_bin2, temp
  integer                      :: i, j, k, r, s, ii, mark, n_test
  integer, allocatable         :: bin_count(:)
  real, allocatable            :: rep(:), delta(:), bin(:)
  real                         :: lo, xo(4), yo(4),                          &
                                  re_tau, ss0, ss1, ss00, ss11,              &
                                  zo, ro, xc, yc, zc, vc, wc, sig_x, sig_yz, &
                                  rmin, rmax, sg, vmax, max_rep,             &
                                  level
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  ! Reading starting time for Swarm statistics from control file
  call Control % Starting_Time_Step_For_Swarm_Statistics &
                                              (n_stat_p, verbose=.true.)

  ! Reynolds number should be passed from Save_Results and number of bins should
  ! be defined in control file, also same for n_bin... (it's okey for now!) 
  re_tau =    142  ! operating shear Reynolds number (a little bit above what
                   ! we have so we don't lose any particle in counting!) 
  n_b    =     64  ! number of bins across half of the channel 
  n_bin1 = 342237  ! time at which we should collect bins info t+=675
  n_bin2 = 343729  ! time at which we should collect bins info t+=1125
                   ! ..particle concentration should be equivalent to t+=1000).
  n_test = 350000  ! testing (for Swarm statistics)

  ! Allocating some arrays for bins
  allocate(rep(Swarm % n_particles)); rep = 0.0
  allocate(bin(n_b)); bin = 0.0     ! bin distance from wall
  allocate(delta(n_b - 1)); delta = 0.0 ! bin thickness
  allocate(bin_count(Swarm % n_particles)); bin_count = 0

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
!    if(First_Proc()) then
!
!      level = 0.0 
!      temp  = 0 
!      ! browse though the bins 
!      do ss = 1, n_b 
!
!        ! calculate the thickness of this bin and the following one
!        bin(ss) = 0.5 * re_tau * (1.0 - cos(PI*(ss - 1)/(n_b - 1)))
!        ss0 = bin(ss) 
!        ss1 = (ss0 * 3.0E-5) / 0.004255 
!        bin_count(ss) = 0
!
!        ! browse through the particles 
!        !do oo = 1, Swarm % n_particles
!        !  part => Swarm % particle(oo)
!        !  if(part % z_n .le. ss1 .and. ss .lt. n_b) then
!        !    bin_count(ss) = bin_count(ss) + 1
!        !  end if 
!        !  rep(oo) = part % re
!        !end do
!
!        ! browse through the particles 
!        do oo = 1, Swarm % n_particles
!          part => Swarm % particle(oo)
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
!      '#', 'St+    = ', Swarm % st 
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
  if(mod(Time % Curr_Dt(), 120) .ne. 0) return

  ! If too late to disturb, get out too
  if(Time % Curr_Dt() > 1200) return

  ! Print a message
  if(First_Proc()) then
    print *, '# Superimposing random eddies on top of velocity field!'
  end if

  ! Minimum and maximum size of eddies
  rmin = 0.2
  rmax = 0.6

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
    call random_number(xo(1));  xo(1) = xo(1) * LX
    call random_number(yo(1));  yo(1) = yo(1) * LY
    call random_number(zo);     zo = ro + (LZ - 2.0*ro) * zo

    ! Handle periodicity; that is: copy eddie in periodic directions
    xo(2:4) = xo(1)
    yo(2:4) = yo(1)
    if(xo(1) > LX/2.0) xo(3) = xo(1) - LX
    if(xo(1) < LX/2.0) xo(3) = xo(1) + LX
    if(yo(1) > LY/2.0) yo(2) = yo(1) - LY
    if(yo(1) < LY/2.0) yo(2) = yo(1) + LY
    xo(4) = xo(3)
    yo(4) = yo(2)

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
  do c = Cells_In_Domain_And_Buffers()
    vmax = max(vmax, abs(v % n(c)))
    vmax = max(vmax, abs(w % n(c)))
  end do
  call Global % Max_Real(vmax)
  do c = Cells_In_Domain_And_Buffers()
    v % n(c) = v % n(c) / vmax / 5.0
    v % o(c) = v % o(c) / vmax / 5.0
    w % n(c) = w % n(c) / vmax / 5.0
    w % o(c) = w % o(c) / vmax / 5.0
  end do

  end subroutine
