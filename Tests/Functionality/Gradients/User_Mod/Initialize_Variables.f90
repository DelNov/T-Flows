!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   This was developed to devise strategy for pressure gradient calculation    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Var_Type),      pointer :: t, p
  integer                      :: c, s, c1, c2, iter
  integer                      :: i_fac, c1_prim, c2_prim, s_prim
  character(8)                 :: name = 'check_xx'
  real,    contiguous, pointer :: phi_f(:), phi_n(:), phi_c(:)
  real,    contiguous, pointer :: phi_x(:), phi_y(:), phi_z(:)
  integer, contiguous, pointer :: c_computed(:), c_visited(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: YES = 1
  integer, parameter :: NO  = YES-1
!==============================================================================!

  call Work % Connect_Real_Face(phi_f)
  call Work % Connect_Real_Node(phi_n)
  call Work % Connect_Real_Cell(phi_x, phi_y, phi_z, phi_c)
  call Work % Connect_Int_Cell (c_computed, c_visited)

  ! Take aliases
  Grid => Flow % pnt_grid
  t    => Flow % t
  p    => Flow % p

  ! Set tight tolerances for Gauss Theorem
  Flow % gauss_miter = 100
  Flow % gauss_tol   = 1e-3

  !----------------------------------------------------------------------!
  !                                                                      !
  !   Test 1:                                                            !
  !                                                                      !
  !   Check Field % Grad_Variable - the cell-based least square method   !
  !       Test with the most established method, should always work      !
  !                                                                      !
  !----------------------------------------------------------------------!

  if(First_Proc()) then
    print *, '#=============================================================='
    print *, '# Performing Test 1 - least square cell based method'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable over-writing boundary values
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()
    t % n(c) =         Grid % xc(c)  &
               + 2.0 * Grid % yc(C)  &
               + 3.0 * Grid % zc(C)
  end do

  ! Call cell-based least-square gradient method ...
  t % grad_method = LEAST_SQUARES
  call Flow % Grad_Variable(t)

  ! ... and plot what you get
  call Grid % Save_Debug_Vtu(                        &
            'test_1_least_square_cell_base_method',  &
            scalar_cell = t % n,                     &
            scalar_name = 't',                       &
            vector_cell = (/t % x, t % y, t % z/),   &
            vector_name = 't_xyz')

  !---------------------------------------------------------!
  !                                                         !
  !   Test 2:                                               !
  !                                                         !
  !   Check if Field % Grad_Gauss_Variable works properly   !
  !           when you ou start from zero gradients         !
  !                                                         !
  !---------------------------------------------------------!

  if(First_Proc()) then
    print *, '#=============================================================='
    print *, '# Performing Test 2 - Gaussian method from poor initial guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Zero the gradients - mimicking a bad first guess
  t % x(:) = 0.0
  t % y(:) = 0.0
  t % z(:) = 0.0

  ! Calculate gradients with Gaussian theorem, (face values are exact) ...
  t % grad_method = GAUSS_THEOREM
  call Flow % Grad_Variable(t)

  ! ... and save the results.  These results should be poor
  call Grid % Save_Debug_Vtu(                                  &
            'test_2_gaussian_method_from_poor_initial_guess',  &
            vector_cell = (/t % x, t % y, t % z/),             &
            vector_name = 't_xyz')

  !---------------------------------------------------------!
  !                                                         !
  !   Test 3:                                               !
  !                                                         !
  !   Check if Field % Grad_Gauss_Variable works properly   !
  !           when you ou start from some gradients         !
  !                                                         !
  !---------------------------------------------------------!

  if(First_Proc()) then
    print *, '#=============================================================='
    print *, '# Performing Test 3 - Gaussian method from better initial guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  do c = Cells_In_Domain_And_Buffers()
    t % n(c) =         Grid % xc(c)  &
               + 2.0 * Grid % yc(C)  &
               + 3.0 * Grid % zc(C)
  end do

  ! Then extrapolate interior values to boundary cells
  ! This still mimics a result from a numerical solution
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    if(c2 < 0) then
      t % n(c2) = t % n(c1)
    end if
  end do

  ! Initialize with some gradients with the most robust and reliable tool
  ! you have at your disposal - least square cell-based method.  These
  ! gradients should be properly calculated inside the domain
  t % grad_method = LEAST_SQUARES
  call Flow % Grad_Variable(t)

  call Grid % Save_Debug_Vtu(                                 &
            'test_31_least_squares_initial_guess_for_gauss',  &
            scalar_cell = t % n,                              &
            scalar_name = 't',                                &
            vector_cell = (/t % x, t % y, t % z/),            &
            vector_name = 't_xyz')

  ! Perform Gauss from gradients which are good inside obtained above ...
  t % grad_method = GAUSS_THEOREM
  call Flow % Grad_Variable(t)

  ! ... and plot what you got.  These should be better, but
  ! not quite.  Particularly not good for tetrahedral grids
  call Grid % Save_Debug_Vtu(                                     &
            'test_32_gaussian_method_from_better_initial_guess',  &
            vector_cell = (/t % x, t % y, t % z/),                &
            vector_name = 't_xyz')

  !-------------------------------------------------------------------------!
  !                                                                         !
  !   Test 4:                                                               !
  !                                                                         !
  !   Try to be more elaborate, make an evolutionary step from Test 3, by   !
  !    interpolating gradients at near-boundary cells from values inside    !
  !                                                                         !
  !-------------------------------------------------------------------------!

  if(First_Proc()) then
    print *, '#=============================================================='
    print *, '# Performing Test 4 - Gaussian method from an elaborate guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  ! (This is the repetition of what was done in Test 3)
  do c = Cells_In_Domain_And_Buffers()
    t % n(c) =         Grid % xc(c)  &
               + 2.0 * Grid % yc(C)  &
               + 3.0 * Grid % zc(C)
  end do

  ! Then extrapolate interior values to boundary cells
  ! This still mimics a result from a numerical solution
  ! (This is the repetition of what was done in Test 3)
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    if(c2 < 0) then
      t % n(c2) = t % n(c1)
    end if
  end do

  ! Initialize with some gradients with the most robust and reliable tool
  ! you have at your disposal - least square cell-based method.  These
  ! gradients should be properly calculated inside the domain
  ! (This is the repetition of what was done in Test 3)
  t % grad_method = LEAST_SQUARES
  call Flow % Grad_Variable(t)

  !-----------------------------------------------------------------------!
  !   Extrapolate interior gradient values to boundary cells using only   !
  !   values which are known, either from interior cells, or cells near   !
  !   boundaries computed in previous iterations.                         !
  !-----------------------------------------------------------------------!

  ! First assume all are computed although you know the cells near the
  ! bounderies are not computed well; but they will be treated below
  c_computed(:) = YES

  ! Mark cells which are not computed for the 1st time
  ! At this point, these are cells near the boundaries
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) c_computed(c1) = NO
  end do

  call Grid % Exchange_Cells_Int(c_computed)

  ! Nullify gradients for cells near boundaries
  ! (where gradients are not computed properly),
  ! and initialize c_visited
  do c = Cells_In_Domain_And_Buffers()
    if(c_computed(c) .eq. NO) then  ! if not computed
      t % x(c) = 0.0
      t % y(c) = 0.0
      t % z(c) = 0.0
      c_visited(c) = 0
    end if
  end do

  !------------------------------------------------!
  !   Extrapolate in several iterations, as long   !
  !   as there are cells which are not computed.   !
  !------------------------------------------------!
  do iter = 1, 3

    ! Browse through faces to find the cells which need ...
    ! ... updating, and accumulate gradients in them as well
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      if(c2 > 0) then
        if(c_computed(c1) .eq. NO .and. c_computed(c2) .eq. YES) then
          t % x(c1) = t % x(c1) + t % x(c2)
          t % y(c1) = t % y(c1) + t % y(c2)
          t % z(c1) = t % z(c1) + t % z(c2)
          c_visited(c1) = c_visited(c1) + 1
        end if

        if(c_computed(c2) .eq. NO .and. c_computed(c1) .eq. YES) then
          t % x(c2) = t % x(c2) + t % x(c1)
          t % y(c2) = t % y(c2) + t % y(c1)
          t % z(c2) = t % z(c2) + t % z(c1)
          c_visited(c2) = c_visited(c2) + 1
        end if
      end if
    end do

    ! Browse throough cells, and calculate final values ...
    ! ... of gradients in the cells which have been visited
    do c = Cells_In_Domain_And_Buffers()
      if(c_visited(c) > 0) then
        t % x(c) = t % x(c) / c_visited(c)
        t % y(c) = t % y(c) / c_visited(c)
        t % z(c) = t % z(c) / c_visited(c)
        c_visited(c)  = 0
        c_computed(c) = YES  ! mark it as computed for the next iteration
      end if
    end do
    call Grid % Exchange_Cells_Real(t % x)
    call Grid % Exchange_Cells_Real(t % y)
    call Grid % Exchange_Cells_Real(t % z)

  end do    ! iter

  ! Save the initial guess that you got
  call Grid % Save_Debug_Vtu(                                &
            'test_41_an_elaborate_initial_guess_for_Gauss',  &
            scalar_cell = t % n,                             &
            scalar_name = 't',                               &
            vector_cell = (/t % x, t % y, t % z/),           &
            vector_name = 't_xyz')

  ! Perform Gauss from gradients which are good inside obtained above ...
  t % grad_method = GAUSS_THEOREM
  call Flow % Grad_Variable(t)

  ! ... and plot what you got.  These should be better, but
  ! not quite.  Particularly not good for tetrahedral grids
  call Grid % Save_Debug_Vtu(                                    &
            'test_42_gaussian_method_from_the_elaborate_guess',  &
            vector_cell = (/t % x, t % y, t % z/),               &
            vector_name = 't_xyz')

  !----------------------------------------------------------------------!
  !                                                                      !
  !   Test 5:                                                            !
  !                                                                      !
  !   Check if the (new) Grad_Gauss_Pressure works as intended to work   !
  !                                                                      !
  !----------------------------------------------------------------------!

  if(First_Proc()) then
    print *, '#=============================================================='
    print *, '# Performing Test 5 - Field % Grad_Gauss_Pressure '
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  ! (This is the repetition of what was done in Test 4 & 5)
  do c = Cells_In_Domain_And_Buffers()
    p % n(c) =         Grid % xc(c)  &
               + 2.0 * Grid % yc(C)  &
               + 3.0 * Grid % zc(C)
  end do

  ! Perform Gauss from gradients which are good inside obtained above ...
  p % grad_method = GAUSS_THEOREM
  call Flow % Grad_Pressure(p)

  ! ... and plot what you got.  These should be better, but
  ! not quite.  Particularly not good for tetrahedral grids
  call Grid % Save_Debug_Vtu(                       &
            'test_5_field_grad_gauss_pressure',     &
            vector_cell = (/p % x, p % y, p % z/),  &
            vector_name = 'p_xyz')

  call Work % Disconnect_Real_Face(phi_f)
  call Work % Disconnect_Real_Node(phi_n)
  call Work % Disconnect_Real_Cell(phi_x, phi_y, phi_z, phi_c)
  call Work % Disconnect_Int_Cell (c_computed, c_visited)

  call Global % End_Parallel
  stop

  end subroutine
