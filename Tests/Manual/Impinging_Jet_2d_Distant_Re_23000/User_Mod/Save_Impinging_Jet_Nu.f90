!==============================================================================!
  subroutine Save_Impinging_Jet_Nu(Turb)
!------------------------------------------------------------------------------!
!   The subroutine creates ASCII file with Nusselt number averaged             !
!   in azimuthal direction.                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: u, v, w, t
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: n_prob, i, s, c1, c2, fu, reg
  character(SL)             :: res_name
  real,    allocatable      :: u_s(:), v_s(:), w_s(:), t_s(:), tau_s(:), q_s(:)
  real,    allocatable      :: z_s(:), r_s(:), rad(:)
  integer, allocatable      :: n_count(:)
  real                      :: r
  logical                   :: there
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

  !-----------------------------------!
  !   Read file with radial sectors   !
  !-----------------------------------!
  inquire(file='rad_coordinate.dat', exist=there)
  if(.not.there) then
    if(First_Proc()) then
      print *, "#=========================================================="
      print *, "# In order to extract Nusselt number profile               "
      print *, "# an ascii file with cell-faces coordinates has to be read."
      print *, "# The name of the file is rad_coordinate.dat.              "
      print *, "# The file format should be as follows:                    "
      print *, "# 10  ! number of cells + 1                                "
      print *, "# 0.0                                                      "
      print *, "# 0.1                                                      "
      print *, "# 0.2                                                      "
      print *, "# ...                                                      "
      print *, "#----------------------------------------------------------"
    end if
    return
  end if

  call File % Open_For_Reading_Ascii('rad_coordinate.dat', fu)

  ! Read the number of searching intervals 
  read(fu,*) n_prob
  allocate(rad(n_prob))

  ! Read the intervals positions
  do i = 1, n_prob
    read(fu, *) rad(i)
  end do
  close(fu)

  allocate(u_s    (n_prob));  u_s     = 0.0
  allocate(v_s    (n_prob));  v_s     = 0.0
  allocate(w_s    (n_prob));  w_s     = 0.0
  allocate(r_s    (n_prob));  r_s     = 0.0
  allocate(z_s    (n_prob));  z_s     = 0.0
  allocate(tau_s  (n_prob));  tau_s   = 0.0
  allocate(q_s    (n_prob));  q_s     = 0.0
  allocate(t_s    (n_prob));  t_s     = 0.0
  allocate(n_count(n_prob));  n_count = 0

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob - 1
    do reg = Boundary_Regions()
      if(Grid % region % name(reg) .eq. 'LOWER_WALL') then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          r = sqrt(Grid % xc(c1)*Grid % xc(c1)  + &
                   Grid % yc(c1)*Grid % yc(c1)) + TINY
          if(r < rad(i+1) .and. r > rad(i)) then
            r_s(i) = r_s(i) + sqrt(Grid % xc(c1)*Grid % xc(c1)  + &
                                   Grid % yc(c1)*Grid % yc(c1))
            u_s(i) = u_s(i) +   u % n(c1) * Grid % xc(c1) / r  + &
                                v % n(c1) * Grid % yc(c1) / r
            v_s(i) = v_s(i) + (-u % n(c1) * Grid % yc(c1) / r  + &
                                v % n(c1) * Grid % xc(c1) / r)
            w_s(i)     = w_s(i)     + w % n(c1)
            t_s(i)     = t_s(i)     + t % n(c2)
            z_s(i)     = z_s(i)     + Grid % zc(c1)  ! acts as wall distance
            tau_s(i)   = tau_s(i)   + sqrt(Turb % tau_wall(c1))
            q_s(i)     = q_s(i)     + t % q(c2)
            n_count(i) = n_count(i) + 1
          end if
        end do  ! faces in this region
      end if    ! region is called 'LOWER_WALL'
    end do      ! through regions
  end do        ! through probes

  !---------------------------------!
  !   Average over all processors   !
  !---------------------------------!
  do i = 1, n_prob
    call Global % Sum_Int(n_count(i))

    call Global % Sum_Real(u_s(i))
    call Global % Sum_Real(v_s(i))
    call Global % Sum_Real(w_s(i))

    call Global % Sum_Real(z_s(i))
    call Global % Sum_Real(tau_s(i))
    call Global % Sum_Real(q_s(i))

    call Global % Sum_Real(r_s(i))
    call Global % Sum_Real(t_s(i))
  end do

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      w_s(i)   = w_s(i)   / n_count(i)
      v_s(i)   = v_s(i)   / n_count(i)
      u_s(i)   = u_s(i)   / n_count(i)
      t_s(i)   = t_s(i)   / n_count(i)
      z_s(i)   = z_s(i)   / n_count(i)
      tau_s(i) = tau_s(i) / n_count(i)
      q_s(i)   = q_s(i)   / n_count(i)
      r_s(i)   = r_s(i)   / n_count(i)
    end if
  end do
  call Global % Wait

  !-----------------------------------!
  !   Write from one processor only   !
  !-----------------------------------!
  if(First_Proc()) then

    ! Set the file name
    call File % Set_Name(res_name, time_step=Time % Curr_Dt(), &
                         appendix='-nu-utau', extension='.dat')
    call File % Open_For_Writing_Ascii(res_name, fu)

    ! Write the file out
    write(fu, '(a66)')  '# 1:Xrad,  ' //  &
                        '  2:Nu,    ' //  &
                        '  3:Utau,  ' //  &
                        '  4:Yplus, ' //  &
                        '  5:Temp,  ' //  &
                        '  6:Points '
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu, '(5e11.3,i11)')                                 &
          r_s(i) / 2.0,                                           &  !  1
          2.0 * q_s(i) / (Flow % conductivity(1)*(t_s(i)-20.0)),  &  !  2
          tau_s(i),                                               &  !  3
          tau_s(i) * z_s(i) / Flow % viscosity(1),                &  !  4
          t_s(i),                                                 &  !  5
          n_count(i)                                                 !  6
      end if
    end do

    close(fu)
  end if

  if(First_Proc()) print *, '# Finished with Impinging_Jet_Nu'

  end subroutine
