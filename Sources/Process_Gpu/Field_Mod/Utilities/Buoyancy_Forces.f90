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
  real                      :: dens_f, temp_f
  integer                   :: c, c1, c2, s, i_cel, reg
# if T_FLOWS_DEBUG == 1
    character(SL)           :: append
# endif
!==============================================================================!

  ! Take some aliases
  b    => Flow % Nat % b
  temp => Flow % temp

  if(comp .eq. 1) then
    grid_dxi => grid_dx;
    grid_si  => grid_sx;
    grav_i   =  Flow % grav_x
  else if(comp .eq. 2) then
    grid_dxi => grid_dy;
    grid_si  => grid_sy;
    grav_i   =  Flow % grav_y
  else if(comp .eq. 3) then
    grid_dxi => grid_dz;
    grid_si  => grid_sz;
    grav_i   =  Flow % grav_z
  end if

  ! I there is no gravity in the direction comp, get out of here
  if(abs(grav_i) < TINY) return

  !$acc parallel loop independent                        &
  !$acc present(grid_region_f_cell, grid_region_l_cell,  &
  !$acc         temp)
  do c = Cells_In_Domain_Gpu()  ! all present
    temp(c) = 0.0
  end do
  !$acc end parallel

  !-------------------------------------------!
  !   Browse through all the interior cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  !$acc parallel loop                                            &
  !$acc present(grid_region_f_cell, grid_region_l_cell,          &
  !$acc         grid_cells_n_cells, grid_cells_c, grid_cells_f,  &
  !$acc         grid_si, grid_dxi,                               &
  !$acc         flow_density, flow_t_n)
  do c1 = Cells_In_Domain_Gpu()  ! all present
    b_tmp = 0.0

    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      if(c2 .gt. 0) then

        ! Temperature and density at the face
        temp_f = Face_Value(s, flow_t_n(c1),     flow_t_n(c2))
        dens_f = Face_Value(s, flow_density(c1), flow_density(c2))

        ! Units here: [kg/m^3 K]
        b_f = dens_f * (Flow % t_ref - temp_f)

        ! Units here: [kg K]
        b_tmp = b_tmp + b_f * 0.5 * abs(grid_si(s) * grid_dxi(s))
      end if
    end do
    !$acc end loop

    ! Unit here: [kg m/s^2 = N]
    temp(c1) = b_tmp * Flow % beta * grav_i
  end do
  !$acc end parallel

  !-------------------------------------------!
  !   Browse through all the boundary cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  do reg = Boundary_Regions()

    !$acc parallel loop  &
    !$acc present(grid_region_f_face, grid_region_l_face, grid_faces_c,  &
    !$acc         grid_si, grid_dxi,                                     &
    !$acc         flow_density,  flow_t_n, temp)
    do s = Faces_In_Region_Gpu(reg)  ! all present
      c1 = grid_faces_c(1,s)  ! inside cell
      c2 = grid_faces_c(1,s)  ! boundary cell

      ! Temperature at the face is identical to
      ! the temperature at the boundary cell
      temp_f = flow_t_n(c2)
      dens_f = flow_density(c1)  ! or in c2?

      ! Units here: [kg/m^3 K]
      b_f = dens_f * (Flow % t_ref - temp_f)

      ! Units here: [kg m/s^2 = N]
      temp(c1) = temp(c1) + b_f * abs(grid_si(s) * grid_dxi(s))  &
               * Flow % beta * grav_i

    end do
    !$acc end parallel

  end do

  !---------------------------------------!
  !   Correct the units for body forces   !
  !---------------------------------------!

  !$acc parallel loop independent                        &
  !$acc present(grid_region_f_cell, grid_region_l_cell,  &
  !$acc         temp)
  do c = Cells_In_Domain_Gpu()
    b(c) = b(c) + temp(c)
  end do
  !$acc end parallel

# if T_FLOWS_DEBUG == 1
  append = 'buoyancy_force__'
  write(append(16:16), '(i1.1)') comp
  call Grid % Save_Debug_Vtu(append,  &
                             scalar_cell=temp, scalar_name='bf')
# endif

  end subroutine
