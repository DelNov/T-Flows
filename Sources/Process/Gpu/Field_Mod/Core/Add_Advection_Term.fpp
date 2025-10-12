!==============================================================================!
  subroutine Add_Advection_Term(Flow, Grid, phi, coef)
!------------------------------------------------------------------------------!
!   Thoroughly re-vamped for the GPU_2                                         !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Dimension of the system under consideration                                !
!     [M]{u} = {b}   [kgm/s^2]   [N]                                           !
!------------------------------------------------------------------------------!
  class(Field_Type), target :: Flow
  type(Grid_Type)           :: Grid
  type(Var_Type),    target :: phi
  real                      :: coef(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:)
  real                      :: b_tmp, coef_phi1, coef_phi2, coef_f, phi_c
  real                      :: fl, dx, dy, dz, blend_1, blend_2, blend_3
  real                      :: phi_luds_1, phi_luds_2
  integer                   :: s, c1, c2, i_cel, reg
  integer                   :: m10_c1c2, m01_c1c2
!==============================================================================!

  call Profiler % Start('Add_Advection_Term')

  ! Now this check is really important: we need the
  ! actual gradients of the variable we are solving
  Assert(Flow % stores_gradients_of .eq. phi % name)

  ! Take some aliases
  b       => Flow % Nat % b
  blend_1 =  phi % blends(1)
  blend_2 =  phi % blends(2)
  blend_3 =  phi % blends(3)

  !-------------------------------------------!
  !                                           !
  !                                           !
  !   Browse through all the interior cells   !
  !      (This can be accelerted on GPU)      !
  !                                           !
  !                                           !
  !-------------------------------------------!

  !------------------------------------------------------!
  !                                                      !
  !   "High" order scheme.  (High is higher than one.)   !
  !                                                      !
  !------------------------------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present (this wasn't independent)
    b_tmp = b(c1)

    do i_cel = 1, Grid % cells_n_cells(c1)

      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      fl = Flow % v_flux % n(s)

      if(c2 .gt. 0) then

        !--------------------!
        !   Centered value   !
        !--------------------!
        phi_c = Face_Value(s, phi % n(c1), phi % n(c2))

        ! Value of the coefficient at the cel face
        coef_f = Face_Value(s, coef(c1), coef(c2))

        !-----------------------------------------------!
        !   Linear upwind differencing scheme (LUDS)    !
        !-----------------------------------------------!

        ! Assume flow is zero at the face (highly unlikely)
        phi_luds_1 = phi_c
        phi_luds_2 = phi_c

        if(fl .ne. 0.0) then
          dx = Grid % dx(s)
          dy = Grid % dy(s)
          dz = Grid % dz(s)

          phi_luds_1 = phi % n(c1)               &
                     + (  Flow % phi_x(c1) * dx  &
                        + Flow % phi_y(c1) * dy  &
                        + Flow % phi_z(c1) * dz ) * merge(+.5,-.5, c1.lt.c2)
          phi_luds_2 = phi % n(c2)               &
                     - (  Flow % phi_x(c2) * dx  &
                        + Flow % phi_y(c2) * dy  &
                        + Flow % phi_z(c2) * dz ) * merge(+.5,-.5, c1.lt.c2)
        end if

        ! Coefficient multiplied with blended variable
        coef_phi1 = coef_f * (  blend_1 * phi_c        &
                              + blend_2 * phi % n(c1)  &
                              + blend_3 * phi_luds_1)
        coef_phi2 = coef_f * (  blend_1 * phi_c        &
                              + blend_2 * phi % n(c2)  &
                              + blend_3 * phi_luds_2)

        !-----------------------!
        !   Update the source   !
        !-----------------------!

        ! Avoid too many (well, two too many) merge commands
        m10_c1c2 = merge(1,0, c1.lt.c2)
        m01_c1c2 = merge(0,1, c1.lt.c2)

        b_tmp = b_tmp - coef_phi1 * max(fl,0.0) * m10_c1c2
        b_tmp = b_tmp - coef_phi2 * min(fl,0.0) * m10_c1c2
        b_tmp = b_tmp + coef_phi2 * max(fl,0.0) * m01_c1c2
        b_tmp = b_tmp + coef_phi1 * min(fl,0.0) * m01_c1c2
      end if
    end do

    b(c1) = b_tmp
  end do
  !$tf-acc loop end

  !--------------------------!
  !                          !
  !   Plain upwind sources   !
  !                          !
  !--------------------------!
  if(phi % blend_matrix) then

    !$tf-acc loop begin
    do c1 = Cells_In_Domain()  ! all present (this wasn't independent)
      b_tmp = b(c1)

      do i_cel = 1, Grid % cells_n_cells(c1)

        c2 = Grid % cells_c(i_cel, c1)
        s  = Grid % cells_f(i_cel, c1)
        fl = Flow % v_flux % n(s)

        if(c2 .gt. 0) then

          ! Value of the coefficient at the cel face
          coef_f = Face_Value(s, coef(c1), coef(c2))

          ! Coefficient multiplied with variable, with upwind blending
          coef_phi1 = coef_f * phi % n(c1)
          coef_phi2 = coef_f * phi % n(c2)

          ! Avoid too many (well, two too many) merge commands
          m10_c1c2 = merge(1,0, c1.lt.c2)
          m01_c1c2 = merge(0,1, c1.lt.c2)

          b_tmp = b_tmp + coef_phi1 * max(fl,0.0) * m10_c1c2
          b_tmp = b_tmp + coef_phi2 * min(fl,0.0) * m10_c1c2
          b_tmp = b_tmp - coef_phi2 * max(fl,0.0) * m01_c1c2
          b_tmp = b_tmp - coef_phi1 * min(fl,0.0) * m01_c1c2
        end if
      end do

      b(c1) = b_tmp
    end do
    !$tf-acc loop end

  end if

  !-------------------------------------------!
  !                                           !
  !   Browse through all the boundary cells   !
  !      (This can be accelerted on GPU)      !
  !                                           !
  !-------------------------------------------!
  do reg = Boundary_Regions()

    !--------------------------------------------------------------------!
    !   Inflow and convective depend on boundary values since they are   !
    !   either given (inflow) or meticulously worked out (convective)    !
    !--------------------------------------------------------------------!
    if(Grid % region % type(reg) .eq. INFLOW  .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)   ! inside cell
        c2 = Grid % faces_c(2,s)   ! boundary cell

        ! Just plain upwind here
        b(c1) = b(c1) - coef(c1) * phi % n(c2) * Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if

    !---------------------------------------------!
    !   Outflow is just a vanishing derivative,   !
    !   use the value from the inside             !
    !---------------------------------------------!
    if(Grid % region % type(reg) .eq. OUTFLOW) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)  ! inside cell

        ! Just plain upwind here
        b(c1) = b(c1) - coef(c1) * phi % n(c1) * Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if
  end do

  call Profiler % Stop('Add_Advection_Term')

  end subroutine
