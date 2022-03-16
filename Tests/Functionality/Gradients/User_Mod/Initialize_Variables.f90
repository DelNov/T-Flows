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
  integer, contiguous, pointer :: c_at_bnd(:), c_cnt(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: YES = 1
  integer, parameter :: NO  = YES-1
!==============================================================================!

  call Work % Connect_Real_Face(phi_f)
  call Work % Connect_Real_Node(phi_n)
  call Work % Connect_Real_Cell(phi_x, phi_y, phi_z, phi_c)
  call Work % Connect_Int_Cell (c_at_bnd, c_cnt)

  ! Take aliases
  Grid => Flow % pnt_grid
  t    => Flow % t
  p    => Flow % p

  ! Set tight tolerances for Gauss Theorem
  Flow % gauss_miter = 100
  Flow % gauss_tol   = 1e-3

  !----------------------------!
  !   Find cells at boundary   !
  !----------------------------!
  c_at_bnd(:) = NO
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    if(c2 < 0) c_at_bnd(c1) = YES
  end do

  !----------------------------------------------------------------------!
  !                                                                      !
  !   Test 1:                                                            !
  !                                                                      !
  !   Check Field % Grad_Variable - the cell-based least square method   !
  !       Test with the most established method, should always work      !
  !                                                                      !
  !----------------------------------------------------------------------!

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 1 - least square cell based method'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable over-writing boundary values
  do c = -Grid % n_bnd_cells, Grid % n_cells
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

  if(this_proc < 2) then
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

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 3 - Gaussian method from better initial guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  do c = 1, Grid % n_cells
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

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 4 - Gaussian method from an elaborate guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  ! (This is the repetition of what was done in Test 3)
  do c = 1, Grid % n_cells
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

  !--------------------------------------------------------------------!
  !   Step 1: Extrapolate interior gradient values to boundary cells   !
  !           using only values from interior cells - which are good   !
  !           This step will leave the boundary cells surrounded by    !
  !           other boundary cells unchanged (reset to zero, really)   !
  !--------------------------------------------------------------------!
  c_cnt(:) = 0
  do c = 1, Grid % n_cells

    ! Cell is at the boundary, intervene here
    if(c_at_bnd(c) .eq. YES) then

      ! Nullify its gradients, and the counter
      ! for cells from which you interpolate
      t % x(c) = 0.0
      t % y(c) = 0.0
      t % z(c) = 0.0
      c_cnt(c) = 0    ! probably not needed, initialized above

      ! Browse through this cell's faces
      do i_fac = 1, Grid % cells_n_faces(c)
        s_prim  = Grid % cells_f(i_fac, c)
        c1_prim = Grid % faces_c(1, s_prim)
        c2_prim = Grid % faces_c(2, s_prim)

        ! Consider c1_prim if it is not cell at boundary
        if(c_at_bnd(c1_prim) .eq. NO) then
          t % x(c) = t % x(c) + t % x(c1_prim)
          t % y(c) = t % y(c) + t % y(c1_prim)
          t % z(c) = t % z(c) + t % z(c1_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

        ! Consider c2_prim if it is not cell at boundary
        ! and not a boundary cell.
        if(c_at_bnd(c2_prim) .eq. NO .and.  &
           c2_prim .gt. 0) then
          t % x(c) = t % x(c) + t % x(c2_prim)
          t % y(c) = t % y(c) + t % y(c2_prim)
          t % z(c) = t % z(c) + t % z(c2_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

      end do

      ! Work out the average
      if(c_cnt(c) .gt. 0) then
        t % x(c) = t % x(c) / c_cnt(c)
        t % y(c) = t % y(c) / c_cnt(c)
        t % z(c) = t % z(c) / c_cnt(c)
      end if

    end if  ! c at bnd

  end do

  !-----------------------------------------------------------------------!
  !   Step 2: Previous step left the cells at boundary, which are only    !
  !           surrounded by other cells at boundary, untreated.  Here,    !
  !           you are less selective and extrapolate even from cells at   !
  !           boundaries, which were interpolated in the Step 1 above.    !
  !-----------------------------------------------------------------------!
  do c = 1, Grid % n_cells

    ! Cell is at the boundary, and hasn't been treated yet
    if(c_at_bnd(c) .eq. YES .and. c_cnt(c) .eq. 0) then

      ! Browse through this cell's faces
      do i_fac = 1, Grid % cells_n_faces(c)
        s_prim  = Grid % cells_f(i_fac, c)
        c1_prim = Grid % faces_c(1, s_prim)
        c2_prim = Grid % faces_c(2, s_prim)

        if(c1_prim .ne. c  .and.  &     ! skip your own self
           c_cnt(c1_prim) .gt. 0) then  ! consider only cells with values
          t % x(c) = t % x(c) + t % x(c1_prim)
          t % y(c) = t % y(c) + t % y(c1_prim)
          t % z(c) = t % z(c) + t % z(c1_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

        if(c2_prim .ne. c  .and.  &     ! skip your own self
           c2_prim .gt. 0  .and.  &     ! consider only cells with values
           c_cnt(c2_prim) .gt. 0) then  ! don't use boundary cells
          t % x(c) = t % x(c) + t % x(c2_prim)
          t % y(c) = t % y(c) + t % y(c2_prim)
          t % z(c) = t % z(c) + t % z(c2_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

      end do

      if(c_cnt(c) .gt. 0) then
        t % x(c) = t % x(c) / c_cnt(c)
        t % y(c) = t % y(c) / c_cnt(c)
        t % z(c) = t % z(c) / c_cnt(c)
      end if

    end if  ! c at bnd

  end do

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

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 5 - Field % Grad_Gauss_Pressure '
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  ! (This is the repetition of what was done in Test 4 & 5)
  do c = 1, Grid % n_cells
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
  call Work % Disconnect_Int_Cell (c_at_bnd, c_cnt)

  end subroutine
