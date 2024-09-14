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

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   temp   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present
    temp(c) = 0.0
  end do
  !$acc end parallel

  !-------------------------------------------!
  !   Browse through all the interior cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   grid_cells_n_cells,  &
  !$acc   grid_cells_c,  &
  !$acc   grid_cells_f,  &
  !$acc   temp,  &
  !$acc   flow_t_n,  &
  !$acc   flow_density,  &
  !$acc   grid_si,  &
  !$acc   grid_dxi   &
  !$acc )
  do c1 = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present
    b_tmp = 0.0

  !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      if(c2 .gt. 0) then

        ! Temperature and density at the face
        temp_f = Face_Value(s, flow_t_n(c1),   flow_t_n(c2))
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
    !$acc present(  &
    !$acc   grid_region_f_face,  &
    !$acc   grid_region_l_face,  &
    !$acc   grid_faces_c,  &
    !$acc   temp,  &
    !$acc   flow_t_n,  &
    !$acc   flow_density,  &
    !$acc   grid_si,  &
    !$acc   grid_dxi   &
    !$acc )
    do s = grid_region_f_face(reg), grid_region_l_face(reg)  ! all present
      c1 = grid_faces_c(1, s)  ! inside cell
      c2 = grid_faces_c(2, s)  ! boundary cell

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

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   b,  &
  !$acc   temp   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)
    b(c) = b(c) + temp(c)
  end do
  !$acc end parallel

# if T_FLOWS_DEBUG == 1
  append = 'buoyancy_force__'
  write(append(16:16), '(i1.1)') comp
  call Grid % Save_Debug_Vtu(append,  &
                             scalar_cell=temp, scalar_name='bf')
# endif

  call Work % Disconnect_Real_Cell(temp)

  end subroutine
