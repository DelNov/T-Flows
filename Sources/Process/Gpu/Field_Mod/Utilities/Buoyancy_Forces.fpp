!==============================================================================!
  subroutine Buoyancy_Forces(Flow, Grid, comp)
!------------------------------------------------------------------------------!
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  type(Grid_Type),   target :: Grid
  integer, intent(in)       :: comp  ! component
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), grid_si(:), grid_dxi(:), temp(:)
  real                      :: b_tmp, b_f, xic1, xic2, x_f, grav_i
  real                      :: w1, w2, dens_f, temp_f
  integer                   :: c, c1, c2, s, i_cel, reg
# if T_FLOWS_DEBUG == 1
    character(SL)           :: append
# endif
!==============================================================================!

  ! Take some aliases
  b => Flow % Nat % b

  if(comp .eq. 1) then
    grid_dxi => Grid % dx;
    grid_si  => Grid % sx;
    grav_i   =  Flow % grav_x
  else if(comp .eq. 2) then
    grid_dxi => Grid % dy;
    grid_si  => Grid % sy;
    grav_i   =  Flow % grav_y
  else if(comp .eq. 3) then
    grid_dxi => Grid % dz;
    grid_si  => Grid % sz;
    grav_i   =  Flow % grav_z
  end if

  ! I there is no gravity in the direction comp, get out of here
  if(abs(grav_i) < TINY) return

  call Work % Connect_Real_Cell(temp)

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    temp(c) = 0.0
  end do
  !$tf-acc loop end

  !-------------------------------------------!
  !   Browse through all the interior cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present
    b_tmp = 0.0

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      if(c2 .gt. 0) then

        w1 = Grid % f(s)
        if(c1.gt.c2) w1 = 1.0 - w1
        w2 = 1.0 - w1

        ! Temperature and density at the face
        temp_f = w1 * Flow % t % n(c1)   + w2 * Flow % t % n(c2)
        dens_f = w1 * Flow % density(c1) + w2 * Flow % density(c2)

        ! Units here: [kg/m^3 K]
        b_f = dens_f * (Flow % t_ref - temp_f)

        ! Units here: [kg K]
        b_tmp = b_tmp + b_f * 0.5 * abs(grid_si(s) * grid_dxi(s))
      end if
    end do

    ! Unit here: [kg m/s^2 = N]
    temp(c1) = b_tmp * Flow % beta * grav_i
  end do
  !$tf-acc loop end

  !-------------------------------------------!
  !   Browse through all the boundary cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  do reg = Boundary_Regions()

    !$tf-acc loop begin
    do s = Faces_In_Region(reg)  ! all present
      c1 = Grid % faces_c(1, s)  ! inside cell
      c2 = Grid % faces_c(2, s)  ! boundary cell

      ! Temperature at the face is identical to
      ! the temperature at the boundary cell
      temp_f = Flow % t % n(c2)
      dens_f = Flow % density(c1)  ! or in c2?

      ! Units here: [kg/m^3 K]
      b_f = dens_f * (Flow % t_ref - temp_f)

      ! Units here: [kg m/s^2 = N]
      temp(c1) = temp(c1) + b_f * abs(grid_si(s) * grid_dxi(s))  &
               * Flow % beta * grav_i

    end do
    !$tf-acc loop end

  end do

  !---------------------------------------!
  !   Correct the units for body forces   !
  !---------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain()
    b(c) = b(c) + temp(c)
  end do
  !$tf-acc loop end

# if T_FLOWS_DEBUG == 1
  append = 'buoyancy_force__'
  write(append(16:16), '(i1.1)') comp
  call Grid % Save_Debug_Vtu(append,  &
                             scalar_cell=temp, scalar_name='bf')
# endif

  call Work % Disconnect_Real_Cell(temp)

  end subroutine
