!==============================================================================!
  subroutine Add_Advection_Term(Process, phi, Flow, Grid, coef_a, coef_b)
!------------------------------------------------------------------------------!
!   Thoroughly re-vamped for the GPU_2                                         !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Dimension of the system under consideration                                !
!     [M]{u} = {b}   [kgm/s^2]   [N]                                           !
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Var_Type),   target :: phi
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
  real                     :: coef_a(-Grid % n_bnd_cells:Grid % n_cells)
  real                     :: coef_b(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), phi_n(:)
  real                      :: b_tmp, coef_phi1, coef_phi2, coef_f, phi_c, blend
  integer                   :: s, c1, c2, i_cel, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Add_Advection_Term')

  ! Take some aliases
  b     => Flow % Nat % b
  phi_n => phi % n
  blend =  Flow % u % blend

  !-------------------------------------------!
  !   Browse through all the interior cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  !$acc parallel loop                                            &
  !$acc present(grid_region_f_cell, grid_region_l_cell,          &
  !$acc         grid_cells_n_cells, grid_cells_c, grid_cells_f,  &
  !$acc         b, phi_n, coef_a, coef_b, flow_v_flux_n)
  do c1 = Cells_In_Domain_Gpu()  ! all present
    b_tmp = b(c1)
    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      if(c2 .gt. 0) then

        ! Centered value
        phi_c = Face_Value(s, phi_n(c1), phi_n(c2))

        ! Value of the coefficient at the cel face
        coef_f = Face_Value(s, coef_a(c1)*coef_b(c1), coef_a(c2)*coef_b(c2))

        ! Coefficient multiplied with variable, with upwind blending
        coef_phi1 = coef_f * ((1.0-blend) * phi_n(c1) + blend * phi_c)
        coef_phi2 = coef_f * ((1.0-blend) * phi_n(c2) + blend * phi_c)

        b_tmp = b_tmp - coef_phi1 * max(flow_v_flux_n(s), 0.0) * merge(1,0, c1.lt.c2)
        b_tmp = b_tmp - coef_phi2 * min(flow_v_flux_n(s), 0.0) * merge(1,0, c1.lt.c2)
        b_tmp = b_tmp + coef_phi2 * max(flow_v_flux_n(s), 0.0) * merge(1,0, c1.gt.c2)
        b_tmp = b_tmp + coef_phi1 * min(flow_v_flux_n(s), 0.0) * merge(1,0, c1.gt.c2)
      end if
    end do
    !$acc end loop
    b(c1) = b_tmp
  end do
  !$acc end parallel

  !-------------------------------------------!
  !   Browse through all the boundary cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  do reg = Boundary_Regions()

    ! Inflow and convective depend on boundary values since they are
    ! either given (inflow) or meticulously worked out (convective)
    if(Grid % region % type(reg) .eq. INFLOW  .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$acc parallel loop                                    &
      !$acc present(grid_region_f_face, grid_region_l_face,  &
      !$acc         grid_faces_c,                            &
      !$acc         b, coef_a, coef_b, phi_n, flow_v_flux_n)
      do s = Faces_In_Region_Gpu(reg)  ! all present
        c1 = grid_faces_c(1,s)  ! inside cell
        c2 = grid_faces_c(1,s)  ! boundary cell

        ! Just plain upwind here
        b(c1) = b(c1) - coef_a(c1) * coef_b(c1) * phi_n(c2) * flow_v_flux_n(s)
      end do
      !$acc end parallel

    end if

    ! Outflow is just a vanishing derivative, use the value from the inside
    if(Grid % region % type(reg) .eq. OUTFLOW) then

      !$acc parallel loop                                    &
      !$acc present(grid_region_f_face, grid_region_l_face,  &
      !$acc         grid_faces_c,                            &
      !$acc         b, coef_a, coef_b, phi_n, flow_v_flux_n)
      do s = Faces_In_Region_Gpu(reg)
        c1 = grid_faces_c(1,s)  ! inside cell

        ! Just plain upwind here
        b(c1) = b(c1) - coef_a(c1) * coef_b(c1) * phi_n(c1) * flow_v_flux_n(s)
      end do
      !$acc end parallel

    end if
  end do

  call Profiler % Stop('Add_Advection_Term')

  end subroutine
