!==============================================================================!
  subroutine Multiphase_Mod_Vof_Grad_Variable_With_Jump(mult, var)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable with jump condition.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Var_Type)                :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  integer                   :: iter, ph, s, c, c1, c2
!==============================================================================!

  ! Take aliases
  flow => mult % pnt_flow
  grid => flow % pnt_grid
  vof  => mult % vof

  ! Refresh buffers for variable
  call Grid_Mod_Exchange_Cells_Real(grid, var % n)

!Feb9  PRINT '(A, 3(ES15.4))', 'VAR ACROSS INT:', VAR % N(1150), VAR % N(1271), VAR % N(1392)

  do ph = 1, 2

!Feb9    PRINT *, '#######################################################'
!Feb9    PRINT *, '# COMPUTING GRADIENTS IN PHASE ', PH
!Feb9    PRINT *, '#######################################################'

    mult % var % n(:) = var % n(:)
    mult % var % x(:) = 0.0
    mult % var % y(:) = 0.0
    mult % var % z(:) = 0.0

!Feb9    ! This is just to check if this is over-written
!Feb9    do c = 1, grid % n_cells
!Feb9      if(mult % cell_at_surf(c)) then
!Feb9        mult % var % n(c) = 10000.
!Feb9      end if
!Feb9    end do

    do iter = 1, 4

      !-----------------------------------------!
      !   Extrapolate values to the interface   !
      !-----------------------------------------!
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1, s)
        c2 = grid % faces_c(2, s)

        if(ph .eq. 1) then

          ! c2 is at the surface and c1 is not: extrapolate to c2
          if(mult % vof % n(c1) < 0.01     .and.  &
             mult % cell_at_elem(c1) .eq. 0 .and.  &
             mult % cell_at_elem(c2) .ne. 0) then
            mult % var % n(c2) = mult % var % n(c1)                 &
                               + grid % dx(s) * mult % var % x(c1)  &
                               + grid % dy(s) * mult % var % y(c1)  &
                               + grid % dz(s) * mult % var % z(c1)
          end if

          ! c1 is at the surface and c2 is not: extrapolate to c1
          if(mult % vof % n(c2) < 0.01     .and.  &
             mult % cell_at_elem(c2) .eq. 0 .and.  &
             mult % cell_at_elem(c1) .ne. 0) then
            mult % var % n(c1) = mult % var % n(c2)                 &
                               - grid % dx(s) * mult % var % x(c2)  &
                               - grid % dy(s) * mult % var % y(c2)  &
                               - grid % dz(s) * mult % var % z(c2)
          end if
        end if

        if(ph .eq. 2) then

          ! c2 is at the surface and c1 is not: extrapolate to c2
          if(mult % vof % n(c1) > 0.99     .and.  &
             mult % cell_at_elem(c1) .eq. 0 .and.  &
             mult % cell_at_elem(c2) .ne. 0) then
            mult % var % n(c2) = mult % var % n(c1)                 &
                               + grid % dx(s) * mult % var % x(c1)  &
                               + grid % dy(s) * mult % var % y(c1)  &
                               + grid % dz(s) * mult % var % z(c1)
          end if

          ! c1 is at the surface and c2 is not: extrapolate to c1
          if(mult % vof % n(c2) > 0.99     .and.  &
             mult % cell_at_elem(c2) .eq. 0 .and.  &
             mult % cell_at_elem(c1) .ne. 0) then
            mult % var % n(c1) = mult % var % n(c2)                 &
                               - grid % dx(s) * mult % var % x(c2)  &
                               - grid % dy(s) * mult % var % y(c2)  &
                               - grid % dz(s) * mult % var % z(c2)
          end if
        end if

      end do  ! do s, end of extrapolation

      !-----------------------------------------------------------------!
      !   Compute gradients with extrapolated values at the interface   !
      !-----------------------------------------------------------------!

      ! Compute individual gradients without refreshing buffers
      call Field_Mod_Grad_Component_No_Refresh(flow, mult % var % n, 1, mult % var % x)
      call Field_Mod_Grad_Component_No_Refresh(flow, mult % var % n, 2, mult % var % y)
      call Field_Mod_Grad_Component_No_Refresh(flow, mult % var % n, 3, mult % var % z)

      ! Refresh buffers for gradient components
      call Grid_Mod_Exchange_Cells_Real(grid, mult % var % x)
      call Grid_Mod_Exchange_Cells_Real(grid, mult % var % y)
      call Grid_Mod_Exchange_Cells_Real(grid, mult % var % z)

    end do  ! through iterations

    !---------------------------------------------------------------------!
    !   Copy computed gradients for the phase back to original variable   !
    !---------------------------------------------------------------------!
    do c = 1, grid % n_cells
      if(mult % cell_at_elem(c) .eq. 0) then
        if(ph .eq. 1 .and. vof % n(c) < 0.01) then
          var % x(c) = mult % var % x(c)
          var % y(c) = mult % var % y(c)
          var % z(c) = mult % var % z(c)
        end if
        if(ph .eq. 2 .and. vof % n(c) > 0.99) then
          var % x(c) = mult % var % x(c)
          var % y(c) = mult % var % y(c)
          var % z(c) = mult % var % z(c)
        end if
      end if
    end do
!Feb9    IF(VAR % NAME .EQ. 'U') THEN
!Feb9      PRINT *, '# PHASE VALUES AT THE END OF COMPUTE GRADS'
!Feb9      DO C = 1, GRID % N_CELLS
!Feb9        IF(MATH_MOD_APPROX_REAL(GRID % YC(C), 0.0) .AND.  &
!Feb9           MATH_MOD_APPROX_REAL(GRID % ZC(C), 0.0)) THEN
!Feb9          PRINT '(a,i5,9(es12.3))', 'ph values: ', c,            &
!Feb9                                     mult % var % n(c),  &
!Feb9                                     mult % var % x(c),  &
!Feb9                                     mult % var % y(c),  &
!Feb9                                     mult % var % z(c)
!Feb9        END IF
!Feb9      END DO
!Feb9    END IF

  end do  ! through phases

  end subroutine
