!==============================================================================!
  subroutine Grad_Variable_With_Jump(Vof, var)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable with jump condition.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Var_Type)          :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: fun
  integer                   :: iter, ph, s, c, c1, c2
!==============================================================================!

  ! Take aliases
  Flow => Vof % pnt_flow
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(var % n)

!Feb9  PRINT '(A, 3(ES15.4))', 'VAR ACROSS INT:', VAR % N(1150), VAR % N(1271), VAR % N(1392)

  do ph = 1, 2

!Feb9    PRINT *, '#######################################################'
!Feb9    PRINT *, '# COMPUTING GRADIENTS IN PHASE ', PH
!Feb9    PRINT *, '#######################################################'

    Vof % var % n(:) = var % n(:)
    Vof % var % x(:) = 0.0
    Vof % var % y(:) = 0.0
    Vof % var % z(:) = 0.0

!Feb9    ! This is just to check if this is over-written
!Feb9    do c = 1, Grid % n_cells
!Feb9      if(Vof % cell_at_surf(c)) then
!Feb9        Vof % var % n(c) = 10000.
!Feb9      end if
!Feb9    end do

    do iter = 1, 4

      !-----------------------------------------!
      !   Extrapolate values to the interface   !
      !-----------------------------------------!
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1, s)
        c2 = Grid % faces_c(2, s)

        if(ph .eq. 1) then

          ! c2 is at the surface and c1 is not: extrapolate to c2
          if(Vof % fun % n(c1) < 0.01     .and.  &
             Vof % cell_at_elem(c1) .eq. 0 .and.  &
             Vof % cell_at_elem(c2) .ne. 0) then
            Vof % var % n(c2) = Vof % var % n(c1)                 &
                               + Grid % dx(s) * Vof % var % x(c1)  &
                               + Grid % dy(s) * Vof % var % y(c1)  &
                               + Grid % dz(s) * Vof % var % z(c1)
          end if

          ! c1 is at the surface and c2 is not: extrapolate to c1
          if(Vof % fun % n(c2) < 0.01     .and.  &
             Vof % cell_at_elem(c2) .eq. 0 .and.  &
             Vof % cell_at_elem(c1) .ne. 0) then
            Vof % var % n(c1) = Vof % var % n(c2)                 &
                               - Grid % dx(s) * Vof % var % x(c2)  &
                               - Grid % dy(s) * Vof % var % y(c2)  &
                               - Grid % dz(s) * Vof % var % z(c2)
          end if
        end if

        if(ph .eq. 2) then

          ! c2 is at the surface and c1 is not: extrapolate to c2
          if(Vof % fun % n(c1) > 0.99     .and.  &
             Vof % cell_at_elem(c1) .eq. 0 .and.  &
             Vof % cell_at_elem(c2) .ne. 0) then
            Vof % var % n(c2) = Vof % var % n(c1)                 &
                               + Grid % dx(s) * Vof % var % x(c1)  &
                               + Grid % dy(s) * Vof % var % y(c1)  &
                               + Grid % dz(s) * Vof % var % z(c1)
          end if

          ! c1 is at the surface and c2 is not: extrapolate to c1
          if(Vof % fun % n(c2) > 0.99     .and.  &
             Vof % cell_at_elem(c2) .eq. 0 .and.  &
             Vof % cell_at_elem(c1) .ne. 0) then
            Vof % var % n(c1) = Vof % var % n(c2)                 &
                               - Grid % dx(s) * Vof % var % x(c2)  &
                               - Grid % dy(s) * Vof % var % y(c2)  &
                               - Grid % dz(s) * Vof % var % z(c2)
          end if
        end if

      end do  ! do s, end of extrapolation

      !-----------------------------------------------------------------!
      !   Compute gradients with extrapolated values at the interface   !
      !-----------------------------------------------------------------!

      ! Compute individual gradients without refreshing buffers
      call Flow % Grad_Component_No_Refresh(Vof % var % n, 1, Vof % var % x)
      call Flow % Grad_Component_No_Refresh(Vof % var % n, 2, Vof % var % y)
      call Flow % Grad_Component_No_Refresh(Vof % var % n, 3, Vof % var % z)

      ! Refresh buffers for gradient components
      call Grid % Exchange_Cells_Real(Vof % var % x)
      call Grid % Exchange_Cells_Real(Vof % var % y)
      call Grid % Exchange_Cells_Real(Vof % var % z)

    end do  ! through iterations

    !---------------------------------------------------------------------!
    !   Copy computed gradients for the phase back to original variable   !
    !---------------------------------------------------------------------!
    do c = 1, Grid % n_cells
      if(Vof % cell_at_elem(c) .eq. 0) then
        if(ph .eq. 1 .and. fun % n(c) < 0.01) then
          var % x(c) = Vof % var % x(c)
          var % y(c) = Vof % var % y(c)
          var % z(c) = Vof % var % z(c)
        end if
        if(ph .eq. 2 .and. fun % n(c) > 0.99) then
          var % x(c) = Vof % var % x(c)
          var % y(c) = Vof % var % y(c)
          var % z(c) = Vof % var % z(c)
        end if
      end if
    end do
!Feb9    IF(VAR % NAME .EQ. 'U') THEN
!Feb9      PRINT *, '# PHASE VALUES AT THE END OF COMPUTE GRADS'
!Feb9      DO C = 1, GRID % N_CELLS
!Feb9        IF(MATH_MOD_APPROX_REAL(GRID % YC(C), 0.0) .AND.  &
!Feb9           MATH_MOD_APPROX_REAL(GRID % ZC(C), 0.0)) THEN
!Feb9          PRINT '(a,i5,9(es12.3))', 'ph values: ', c,            &
!Feb9                                     Vof % var % n(c),  &
!Feb9                                     Vof % var % x(c),  &
!Feb9                                     Vof % var % y(c),  &
!Feb9                                     Vof % var % z(c)
!Feb9        END IF
!Feb9      END DO
!Feb9    END IF

  end do  ! through phases

  end subroutine
