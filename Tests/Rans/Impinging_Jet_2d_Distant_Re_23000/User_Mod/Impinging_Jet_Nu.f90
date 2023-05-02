!==============================================================================!
  subroutine User_Mod_Impinging_Jet_Nu(Turb)
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
  integer                   :: n_prob, i, s, c1, c2, fu
  character(SL)             :: res_name
  real,    allocatable      :: u_p(:), v_p(:), w_p(:), t_p(:), tau_p(:), q_p(:)
  real,    allocatable      :: z_p(:), r_p(:), rad_1(:)
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
  allocate(rad_1(n_prob))

  ! Read the intervals positions
  do i = 1, n_prob
    read(fu, *) rad_1(i)
  end do
  close(fu)

  allocate(u_p    (n_prob));  u_p     = 0.0
  allocate(v_p    (n_prob));  v_p     = 0.0
  allocate(w_p    (n_prob));  w_p     = 0.0
  allocate(r_p    (n_prob));  r_p     = 0.0
  allocate(z_p    (n_prob));  z_p     = 0.0
  allocate(tau_p  (n_prob));  tau_p   = 0.0
  allocate(q_p    (n_prob));  q_p     = 0.0
  allocate(t_p    (n_prob));  t_p     = 0.0
  allocate(n_count(n_prob));  n_count = 0

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid % Bnd_Cond_Name_At_Cell(c2) .eq. 'LOWER_WALL') then
          r = sqrt(Grid % xc(c1)*Grid % xc(c1)  + &
                   Grid % yc(c1)*Grid % yc(c1)) + TINY
          if(r < rad_1(i+1) .and.  &
             r > rad_1(i)) then
            r_p(i) = r_p(i) + sqrt(Grid % xc(c1)*Grid % xc(c1)  + &
                                   Grid % yc(c1)*Grid % yc(c1))
            u_p(i) = u_p(i) +   u % n(c1) * Grid % xc(c1) / r  + &
                                v % n(c1) * Grid % yc(c1) / r
            v_p(i) = v_p(i) + (-u % n(c1) * Grid % yc(c1) / r  + &
                                v % n(c1) * Grid % xc(c1) / r)
            w_p(i)     = w_p(i)     + w % n(c1)
            t_p(i)     = t_p(i)     + t % n(c2)
            z_p(i)     = z_p(i)     + Grid % zc(c1)  ! acts as wall distance
            tau_p(i)   = tau_p(i)   + sqrt(Turb % tau_wall(c1))
            q_p(i)     = q_p(i)     + t % q(c2)
            n_count(i) = n_count(i) + 1
          end if
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

    call Global % Sum_Real(z_p(i))
    call Global % Sum_Real(tau_p(i))
    call Global % Sum_Real(q_p(i))

    call Global % Sum_Real(r_p(i))
    call Global % Sum_Real(t_p(i))

!    call Global % Sum_Real_Array(8, (/u_p  (i),  &
!                                             v_p  (i),  &
!                                             w_p  (i),  &
!                                             z_p  (i),  &
!                                             tau_p(i),  &
!                                             q_p  (i),  &
!                                             r_p  (i),  &
!                                             t_p  (i)/))
  end do

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      w_p(i)   = w_p(i)   / n_count(i)
      v_p(i)   = v_p(i)   / n_count(i)
      u_p(i)   = u_p(i)   / n_count(i)
      t_p(i)   = t_p(i)   / n_count(i)
      z_p(i)   = z_p(i)   / n_count(i)
      tau_p(i) = tau_p(i) / n_count(i)
      q_p(i)   = q_p(i)   / n_count(i)
      r_p(i)   = r_p(i)   / n_count(i)
    end if
  end do
  call Global % Wait

  !-----------------------------------!
  !   Write from one processor only   !
  !-----------------------------------!
  if(First_Proc()) then

    ! Set the file name
    call File % Set_Name(res_name,                      &
                         time_step = Time % Curr_Dt(),  &
                         appendix  = '-nu-utau',        &
                         extension = '.dat')
    call File % Open_For_Writing_Ascii(res_name, fu)

    ! Write the file out
    write(fu, *) '# 1:Xrad, 2:Nu, 3:Utau, 4:Yplus, 5:Temp, 6:Numb of points '
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu, '(5e11.3,i6)')                                  &
          r_p(i) / 2.0,                                           &  !  1
          2.0 * q_p(i) / (Flow % conductivity(1)*(t_p(i)-20.0)),  &  !  2
          tau_p(i),                                               &  !  3
          tau_p(i) * z_p(i) / Flow % viscosity(1),                &  !  4
          t_p(i),                                                 &  !  5
          n_count(i)                                                 !  6
      end if
    end do

    close(fu)
  end if

  if(First_Proc()) print *, '# Finished with Impinging_Jet_Nu'

  end subroutine
