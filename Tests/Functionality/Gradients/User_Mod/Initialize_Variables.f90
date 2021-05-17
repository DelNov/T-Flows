!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, turb, Vof, swarm, sol)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!------------------------------------------------------------------------------!
  use Work_Mod, only: phi_f    => r_face_01,  &
                      phi_n    => r_node_01,  &
                      phi_x    => r_cell_01,  &
                      phi_y    => r_cell_02,  &
                      phi_z    => r_cell_03,  &
                      phi_c    => r_cell_04,  &
                      c_at_bnd => i_cell_01,  &
                      c_cnt    => i_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: swarm
  type(Solver_Type), target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: t, p
  integer                  :: c, s, c1, c2, iter
  integer                  :: i_fac, c1_prim, c2_prim, s_prim
  character(8)             :: name = 'check_xx'
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: YES = 1
  integer, parameter :: NO  = YES-1
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  t    => Flow % t
  p    => Flow % p

  !----------------------------!
  !   Find cells at boundary   !
  !----------------------------!
  c_at_bnd(:) = NO
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

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
  do c = -grid % n_bnd_cells, grid % n_cells
    t % n(c) =         grid % xc(c)  &
               + 2.0 * grid % yc(C)  &
               + 3.0 * grid % zc(C)
  end do

  ! Call cell-based least-square gradient method ...
  call Flow % Grad_Variable(t)

  ! ... and plot what you get
  call Grid_Mod_Save_Debug_Vtu(                            &
            grid, 'test_1_least_square_cell_base_method',  &
            scalar_cell = t % n,                           &
            scalar_name = 't',                             &
            vector_cell = (/t % x, t % y, t % z/),         &
            vector_name = 't_xyz')

  !-------------------------------------------------!
  !                                                 !
  !   Test 2:                                       !
  !                                                 !
  !   Check Field % Grad_Component_Faces_To_Cells   !
  !     It is the face-based least square method    !
  !                                                 !
  !-------------------------------------------------!

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 2 - least square face based method'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values omittiing boundary values
  do c = 1, grid % n_cells
    phi_c(c) =         grid % xc(c)  &
               + 2.0 * grid % yc(C)  &
               + 3.0 * grid % zc(C)
  end do

  ! Specify exact face values over-writing boundary values
  do s = 1, grid % n_faces
    phi_f(s) =       grid % xf(s)  &
             + 2.0 * grid % yf(s)  &
             + 3.0 * grid % zf(s)
  end do

  ! Call least-squares face-based gradient calculation ...
  call Flow % Grad_Component_Faces_To_Cells(phi_c, phi_f, 1, phi_x)
  call Flow % Grad_Component_Faces_To_Cells(phi_c, phi_f, 2, phi_y)
  call Flow % Grad_Component_Faces_To_Cells(phi_c, phi_f, 3, phi_z)

  ! ... and save the results
  call Grid_Mod_Save_Debug_Vtu(                            &
            grid, 'test_2_least_square_face_base_method',  &
            vector_cell = (/phi_x, phi_y, phi_z/),         &
            vector_name = 'phi_xyz')

  !---------------------------------------------------------!
  !                                                         !
  !   Test 3:                                               !
  !                                                         !
  !   Check if Field % Grad_Gauss_Variable works properly   !
  !           when you ou start from good gradients         !
  !                                                         !
  !---------------------------------------------------------!

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 3 - Gaussian method from good initial guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify variable over-writing boundary values
  do c = -grid % n_bnd_cells, grid % n_cells
    t % n(c) =         grid % xc(c)  &
               + 2.0 * grid % yc(C)  &
               + 3.0 * grid % zc(C)
  end do

  ! Calculate gradients with Gaussian theorem, (face values are exact) ...
  call Flow % Grad_Gauss_Variable(t)

  ! ... and plot the results
  call Grid_Mod_Save_Debug_Vtu(                                      &
            grid, 'test_3_gaussian_method_from_good_initial_guess',  &
            vector_cell = (/t % x, t % y, t % z/),                   &
            vector_name = 't_xyz')

  !---------------------------------------------------------!
  !                                                         !
  !   Test 4:                                               !
  !                                                         !
  !   Check if Field % Grad_Gauss_Variable works properly   !
  !           when you ou start from zero gradients         !
  !                                                         !
  !---------------------------------------------------------!

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 4 - Gaussian method from poor initial guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Zero the gradients - mimicking a bad first guess
  t % x(:) = 0.0
  t % y(:) = 0.0
  t % z(:) = 0.0

  ! Calculate gradients with Gaussian theorem, (face values are exact) ...
  call Flow % Grad_Gauss_Variable(t)

  ! ... and save the results.  These results should be poor
  call Grid_Mod_Save_Debug_Vtu(                                      &
            grid, 'test_4_gaussian_method_from_poor_initial_guess',  &
            vector_cell = (/t % x, t % y, t % z/),                   &
            vector_name = 't_xyz')

  !---------------------------------------------------------!
  !                                                         !
  !   Test 5:                                               !
  !                                                         !
  !   Check if Field % Grad_Gauss_Variable works properly   !
  !           when you ou start from some gradients         !
  !                                                         !
  !---------------------------------------------------------!

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 5 - Gaussian method from better initial guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  do c = 1, grid % n_cells
    t % n(c) =         grid % xc(c)  &
               + 2.0 * grid % yc(C)  &
               + 3.0 * grid % zc(C)
  end do

  ! Then extrapolate interior values to boundary cells
  ! This still mimics a result from a numerical solution
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    if(c2 < 0) then
      t % n(c2) = t % n(c1)
    end if
  end do

  ! Initialize with some gradients with the most robust and reliable tool
  ! you have at your disposal - least square cell-based method.  These
  ! gradients should be properly calculated inside the domain
  call Flow % Grad_Variable(t)

  call Grid_Mod_Save_Debug_Vtu(                                     &
            grid, 'test_51_least_squares_initial_guess_for_gauss',  &
            scalar_cell = t % n,                                    &
            scalar_name = 't',                                      &
            vector_cell = (/t % x, t % y, t % z/),                  &
            vector_name = 't_xyz')

  ! Perform Gauss from gradients which are good inside obtained above ...
  call Flow % Grad_Gauss_Variable(t)

  ! ... and plot what you got.  These should be better, but
  ! not quite.  Particularly not good for tetrahedral grids
  call Grid_Mod_Save_Debug_Vtu(                                         &
            grid, 'test_52_gaussian_method_from_better_initial_guess',  &
            vector_cell = (/t % x, t % y, t % z/),                      &
            vector_name = 't_xyz')

  !-------------------------------------------------------------------------!
  !                                                                         !
  !   Test 6:                                                               !
  !                                                                         !
  !   Try to be more elaborate, make an evolutionary step from Test 5, by   !
  !    interpolating gradients at near-boundary cells from values inside    !
  !                                                                         !
  !-------------------------------------------------------------------------!

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 6 - Gaussian method from an elaborate guess'
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  ! (This is the repetition of what was done in Test 5)
  do c = 1, grid % n_cells
    t % n(c) =         grid % xc(c)  &
               + 2.0 * grid % yc(C)  &
               + 3.0 * grid % zc(C)
  end do

  ! Then extrapolate interior values to boundary cells
  ! This still mimics a result from a numerical solution
  ! (This is the repetition of what was done in Test 5)
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    if(c2 < 0) then
      t % n(c2) = t % n(c1)
    end if
  end do

  ! Initialize with some gradients with the most robust and reliable tool
  ! you have at your disposal - least square cell-based method.  These
  ! gradients should be properly calculated inside the domain
  ! (This is the repetition of what was done in Test 5)
  call Flow % Grad_Variable(t)

  !--------------------------------------------------------------------!
  !   Step 1: Extrapolate interior gradient values to boundary cells   !
  !           using only values from interior cells - which are good   !
  !           This step will leave the boundary cells surrounded by    !
  !           other boundary cells unchanged (reset to zero, really)   !
  !--------------------------------------------------------------------!
  c_cnt(:) = 0
  do c = 1, grid % n_cells

    ! Cell is at the boundary, intervene here
    if(c_at_bnd(c) .eq. YES) then

      ! Nullify its gradients, and the counter
      ! for cells from which you interpolate
      t % x(c) = 0.0
      t % y(c) = 0.0
      t % z(c) = 0.0
      c_cnt(c) = 0    ! probably not needed, initialized above

      ! Browse through this cell's faces
      do i_fac = 1, grid % cells_n_faces(c)
        s_prim  = grid % cells_f(i_fac, c)
        c1_prim = grid % faces_c(1, s_prim)
        c2_prim = grid % faces_c(2, s_prim)

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
  do c = 1, grid % n_cells

    ! Cell is at the boundary, and hasn't been treated yet
    if(c_at_bnd(c) .eq. YES .and. c_cnt(c) .eq. 0) then

      ! Browse through this cell's faces
      do i_fac = 1, grid % cells_n_faces(c)
        s_prim  = grid % cells_f(i_fac, c)
        c1_prim = grid % faces_c(1, s_prim)
        c2_prim = grid % faces_c(2, s_prim)

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
  call Grid_Mod_Save_Debug_Vtu(                                    &
            grid, 'test_61_an_elaborate_initial_guess_for_Gauss',  &
            scalar_cell = t % n,                                   &
            scalar_name = 't',                                     &
            vector_cell = (/t % x, t % y, t % z/),                 &
            vector_name = 't_xyz')

  ! Perform Gauss from gradients which are good inside obtained above ...
  call Flow % Grad_Gauss_Variable(t)

  ! ... and plot what you got.  These should be better, but
  ! not quite.  Particularly not good for tetrahedral grids
  call Grid_Mod_Save_Debug_Vtu(                                        &
            grid, 'test_62_gaussian_method_from_the_elaborate_guess',  &
            vector_cell = (/t % x, t % y, t % z/),                     &
            vector_name = 't_xyz')

  !----------------------------------------------------------------------!
  !                                                                      !
  !   Test 7:                                                            !
  !                                                                      !
  !   Check if the (new) Grad_Gauss_Pressure works as intended to work   !
  !                                                                      !
  !----------------------------------------------------------------------!

  if(this_proc < 2) then
    print *, '#=============================================================='
    print *, '# Performing Test 7 - Field % Grad_Gauss_Pressure '
    print *, '#--------------------------------------------------------------'
  end if

  ! Specify exact cell values variable not touching bounary values
  ! This is supposed to mimic pressure solution, for example
  ! (This is the repetition of what was done in Test 5 & 6)
  do c = 1, grid % n_cells
    p % n(c) =         grid % xc(c)  &
               + 2.0 * grid % yc(C)  &
               + 3.0 * grid % zc(C)
  end do

  ! Perform Gauss from gradients which are good inside obtained above ...
  call Flow % Grad_Gauss_Pressure(p)

  ! ... and plot what you got.  These should be better, but
  ! not quite.  Particularly not good for tetrahedral grids
  call Grid_Mod_Save_Debug_Vtu(                            &
            grid, 'test_7_field_grad_gauss_pressure',  &
            vector_cell = (/p % x, p % y, p % z/),         &
            vector_name = 'p_xyz')

  end subroutine
