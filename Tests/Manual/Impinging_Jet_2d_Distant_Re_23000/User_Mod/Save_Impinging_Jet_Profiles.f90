!==============================================================================!
  subroutine Save_Impinging_Jet_Profiles(Turb)
!------------------------------------------------------------------------------!
!   Subroutine reads the .1d file with wall normal coordinates and extracts    !
!   solutions for comparison with corresponding experimental measurements.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: u, v, w, t
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: n_prob, i, c, s, c2, k, fu, ind, reg
  character(SL)             :: coord_name, res_name, ext
  real,    allocatable      :: z_p(:), zm_p(:), u_p(:),  v_p(:), w_p(:), t_p(:)
  real,    allocatable      :: kin_p(:), eps_p(:), vis_p(:), zet_p(:), f22_p(:)
  integer, allocatable      :: n_count(:)
  real                      :: r, r1, r2, u_rad, u_tan, lnum, area_in, velo_in
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t
  kin  => Turb % kin
  eps  => Turb % eps
  zeta => Turb % zeta
  f22  => Turb % f22

  ! Calculate the average inlet velocity
  area_in = 0.0
  velo_in = 0.0
  do reg = Boundary_Regions()
    if(Grid % region % name(reg) .eq. 'PIPE_INLET') then
      do s = Faces_In_Region(reg)
        c2 = Grid % faces_c(2,s)  ! fetch the boundary cell
        velo_in = velo_in - w % n(c2) * Grid % sz(s)
        area_in = area_in + Grid % s(s)
      end do  ! through faces in the region
    end if    ! region is called 'PIPE_INLET'
  end do      ! through all boundary regions
  call Global % Sum_Real(area_in)
  call Global % Sum_Real(velo_in)
  Assert(area_in > 0.0)
  velo_in = velo_in / area_in

  ! Set the name for coordinate file
  call File % Set_Name(coord_name, extension='.1d')
  call File % Open_For_Reading_Ascii(coord_name, fu)

  ! Read the number of searching intervals
  read(fu,*) n_prob
  allocate(z_p(n_prob))

  ! Read the intervals positions
  do i = 1, n_prob
    read(fu,*) ind, z_p(i)
  end do
  close(fu)

  allocate(u_p  (n_prob));  u_p   = 0.0
  allocate(v_p  (n_prob));  v_p   = 0.0
  allocate(w_p  (n_prob));  w_p   = 0.0
  allocate(kin_p(n_prob));  kin_p = 0.0
  allocate(eps_p(n_prob));  eps_p = 0.0
  allocate(vis_p(n_prob));  vis_p = 0.0
  allocate(zet_p(n_prob));  zet_p = 0.0
  allocate(f22_p(n_prob));  f22_p = 0.0
  allocate(zm_p (n_prob));  zm_p  = 0.0

  allocate(n_count(n_prob)); n_count=0

  if(Flow % heat_transfer) then
    allocate(t_p(n_prob));   t_p = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do k = 0, 6
    if(k .eq. 0) then
      r1 = 0.0
      r2 = 0.04
      lnum = 0.0
      call File % Set_Name(res_name, time_step=Time % Curr_Dt(),  &
                           appendix='-0.0D', extension='.dat')
    else if(k .eq. 1) then
      r1 = 0.992
      r2 = 1.0
      lnum = 0.5
      call File % Set_Name(res_name, time_step=Time % Curr_Dt(),  &
                           appendix='-0.5D', extension='.dat')
    else if(k .eq. 2) then
      r1 = 2.0
      r2 = 2.1500
      lnum = 1.0
      call File % Set_Name(res_name, time_step=Time % Curr_Dt(),  &
                           appendix='-1.0D', extension='.dat')
    else if(k .eq. 3) then
      r1 = 2.9744
      r2 = 3.0684
      lnum = 1.5
      call File % Set_Name(res_name, time_step=Time % Curr_Dt(),  &
                           appendix='-1.5D', extension='.dat')
    else if(k .eq. 4) then
      r1 = 3.9098
      r2 = 4.1433
      lnum = 2.0
      call File % Set_Name(res_name, time_step=Time % Curr_Dt(),  &
                           appendix='-2.0D', extension='.dat')
    else if(k .eq. 5) then
      r1 = 0.4803200E+01
      r2 = 0.5347000E+01
      lnum = 2.5
      call File % Set_Name(res_name, time_step=Time % Curr_Dt(),  &
                           appendix='-2.5D', extension='.dat')
    else if(k .eq. 6) then
      r1 = 0.5876600E+01
      r2 = 0.6000000E+01
      lnum = 3.0
      call File % Set_Name(res_name, time_step=Time % Curr_Dt(),  &
                           appendix='-3.0D', extension='.dat')
    end if

    do i = 1, n_prob-1
      do c = Cells_In_Domain()
        r = sqrt(Grid % xc(c)**2 + Grid % yc(c)**2) + TINY
        if(r > r1 .and. r < r2) then
          if(Grid % zc(c) > z_p(i) .and.  &
             Grid % zc(c) < z_p(i+1)) then
            u_rad  = (  u % n(c)*Grid % xc(c)/r   &
                      + v % n(c)*Grid % yc(c)/r)
            u_tan  = (- u % n(c)*Grid % yc(c)/r   &
                      + v % n(c)*Grid % xc(c)/r)
            u_p(i) = u_p(i) + sqrt(  u % n(c)**2   &
                                   + v % n(c)**2   &
                                   + w % n(c)**2)
            v_p(i) = v_p(i) + u_rad
            w_p(i) = w_p(i) + w % n(c)

            kin_p(i) = kin_p(i) + kin % n(c)
            eps_p(i) = eps_p(i) + eps % n(c)
            vis_p(i) = vis_p(i) + Turb % vis_t(c) / Flow % viscosity(c)
            zet_p(i) = zet_p(i) + zeta % n(c)
            f22_p(i) = f22_p(i) + f22 % n(c)
            t_p(i)   = t_p(i)   + t % n(c)

            zm_p(i) = zm_p(i) + Grid % zc(c)
            n_count(i) = n_count(i) + 1

          end if
        end if
      end do
    end do

    !---------------------------------!
    !   Average over all processors   !
    !---------------------------------!
    do i = 1, n_prob
      call Global % Sum_Int(n_count(i))

      call Global % Sum_Real(u_p(i))
      call Global % Sum_Real(v_p(i))
      call Global % Sum_Real(w_p(i))
      call Global % Sum_Real(zm_p(i))

      call Global % Sum_Real(kin_p(i))
      call Global % Sum_Real(eps_p(i))
      call Global % Sum_Real(vis_p(i))
      call Global % Sum_Real(zet_p(i))
      call Global % Sum_Real(f22_p(i))
      call Global % Sum_Real(t_p(i))
    end do

    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        w_p(i)   = w_p(i)   / n_count(i)
        u_p(i)   = u_p(i)   / n_count(i)
        v_p(i)   = v_p(i)   / n_count(i)
        kin_p(i) = kin_p(i) / n_count(i)
        eps_p(i) = eps_p(i) / n_count(i)
        vis_p(i) = vis_p(i) / n_count(i)
        zet_p(i) = zet_p(i) / n_count(i)
        f22_p(i) = f22_p(i) / n_count(i)
        t_p(i)   = t_p(i)   / n_count(i)
        zm_p(i)  = zm_p(i)  / n_count(i)
      end if
    end do

    !-----------------------------------!
    !   Write from one processor only   !
    !-----------------------------------!
    if(First_Proc()) then

      call File % Open_For_Writing_Ascii(res_name, fu)

      write(fu,'(a,a120)') '#', '  1:Xrad,   '   //  &
                                '  2:Umag,   '   //  &
                                '  3:Urad,   '   //  &
                                '  4:Uaxi,   '   //  &
                                '  5:Kin,    '   //  &
                                '  6:Eps,    '   //  &
                                '  7:Temp,   '   //  &
                                '  8:Vis_t/l,'   //  &
                                '  9:Zeta,   '   //  &
                                ' 10:F22     '

      do i = 1, n_prob
        if(n_count(i) .ne. 0) then
          write(fu,'(10e12.3)') zm_p(i)  / 2.0,         &  !  1
                                u_p(i)   / VELO_IN,     &  !  2
                                v_p(i)   / VELO_IN,     &  !  3
                                w_p(i)   / VELO_IN,     &  !  4
                                kin_p(i) / VELO_IN**2,  &  !  5
                                eps_p(i),               &  !  6
                                t_p(i),                 &  !  7
                                vis_p(i),               &  !  8
                                zet_p(i),               &  !  9
                                f22_p(i)                   ! 10
        end if
      end do
      close(fu)
    end if

    n_count(:) = 0
    w_p(:)     = 0.0
    u_p(:)     = 0.0
    v_p(:)     = 0.0
    kin_p(:)   = 0.0
    eps_p(:)   = 0.0
    vis_p(:)   = 0.0
    zet_p(:)   = 0.0
    f22_p(:)   = 0.0
    t_p(:)     = 0.0
    zm_p(:)    = 0.0

    if(First_Proc()) print *, 'Finished with profile r/D =  ', lnum

  end do   ! end number of radius

  if(First_Proc()) write(*,*) 'Finished with User_Impinging_Jet_Profiles'

  end subroutine
