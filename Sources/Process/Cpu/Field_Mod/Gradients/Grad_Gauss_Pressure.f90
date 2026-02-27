!==============================================================================!
  subroutine Grad_Gauss_Pressure(Flow, p)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for computing the gradient of pressure or
!>  pressure correction using Gauss's theorem. It's an iterative method
!>  designed to refine the calculation of pressure gradients in a flow field.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initial gradient calculation: The subroutine starts by using a least     !
!     squares method (by calling Grad_Least_Pressure) to initialize the        !
!     gradient computation.                                                    !
!   * Iterative refinement:                                                    !
!     - Boundary gradient extrapolation: Gradients at boundary cells are       !
!       extrapolated from the interior cells.                                  !
!     - Iterative Updates: In several iterations, the subroutine refines the   !!
!       gradient calculation by accumulating and averaging gradient values     !
!       from neighboring cells.                                                !
!     - Final Gauss's theorem application: After refining the gradient         !
!       estimates, Gauss's theorem is applied to compute the final gradient    !
!       of pressure or pressure correction.
!------------------------------------------------------------------------------!
!   Note                                                                       !
!                                                                              !
!   * With OpenMP, this procedure got a speedup of 3.5 on 1M mesh & 4 threads. !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow  !! parent object class
  type(Var_Type),    target :: p     !! pressure or pressure correction
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),     pointer :: Grid
  integer                      :: s, c, c1, c2, iter, reg
  integer, contiguous, pointer :: c_computed(:), c_visited(:), faces_c(:,:)
  real,    contiguous, pointer :: p_x(:), p_y(:), p_z(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: YES = 1
  integer, parameter :: NO  = YES-1
!==============================================================================!

  call Profiler % Start('Grad_Gauss_Pressure')

  call Work % Connect_Int_Cell(c_computed, c_visited)

  ! Take alias
  ! OpenMP doesn't unerstand Fortran's members (%), that's ...
  ! ... why aliases for faces_c, p_x, p_y and p_z are needed
  Grid    => Flow % pnt_grid
  faces_c => Grid % faces_c
  p_x     => p % x
  p_y     => p % y
  p_z     => p % z

  ! Initialize with some gradients with the most robust and reliable tool
  ! you have at your disposal - least square cell-based method.  These
  ! gradients should be properly calculated inside the domain.
  call Flow % Grad_Least_Pressure(p)

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
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1, s)
      c_computed(c1) = NO
    end do  ! faces
  end do    ! all boundary regions

  call Grid % Exchange_Cells_Int(c_computed)

  ! Nullify gradients for cells near boundaries
  ! (where gradients are not computed properly),
  ! and initialize c_visited
  !$omp parallel do                                   &
  !$omp private(c)                                    &
  !$omp shared(c_computed, c_visited, p_x, p_y, p_z)
  do c = Cells_In_Domain()
    if(c_computed(c) .eq. NO) then  ! if not computed
      p_x(c) = 0.0
      p_y(c) = 0.0
      p_z(c) = 0.0
      c_visited(c) = 0
    end if
  end do
  !$omp end parallel do

  !------------------------------------------------!
  !   Extrapolate in several iterations, as long   !
  !   as there are cells which are not computed.   !
  !------------------------------------------------!
  do iter = 1, 3

    ! Browse through faces to find the cells which need ...
    ! ... updating, and accumulate gradients in them as well
    !$omp parallel do                                            &
    !$omp private(s, c1, c2)                                     &
    !$omp shared(faces_c, c_computed, c_visited, p_x, p_y, p_z)
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = faces_c(1, s)
      c2 = faces_c(2, s)

      if(c_computed(c1) .eq. NO .and. c_computed(c2) .eq. YES) then
        p_x(c1) = p_x(c1) + p_x(c2)
        p_y(c1) = p_y(c1) + p_y(c2)
        p_z(c1) = p_z(c1) + p_z(c2)
        c_visited(c1) = c_visited(c1) + 1
      end if

      if(c_computed(c2) .eq. NO .and. c_computed(c1) .eq. YES) then
        p_x(c2) = p_x(c2) + p_x(c1)
        p_y(c2) = p_y(c2) + p_y(c1)
        p_z(c2) = p_z(c2) + p_z(c1)
        c_visited(c2) = c_visited(c2) + 1
      end if
    end do
    !$omp end parallel do

    ! Browse throough cells, and calculate final values ...
    ! ... of gradients in the cells which have been visited
    !$omp parallel do                                       &
    !$omp private(c)                                        &
    !$omp shared(c_computed, c_visited, p_x, p_y, p_z)
    do c = Cells_In_Domain()
      if(c_visited(c) > 0) then
        p_x(c) = p_x(c) / c_visited(c)
        p_y(c) = p_y(c) / c_visited(c)
        p_z(c) = p_z(c) / c_visited(c)
        c_visited(c)  = 0
        c_computed(c) = YES  ! mark it as computed for the next iteration
      end if
    end do
    !$omp end parallel do

    call Grid % Exchange_Cells_Real(p % x)
    call Grid % Exchange_Cells_Real(p % y)
    call Grid % Exchange_Cells_Real(p % z)

  end do    ! iter

  !--------------------------------------------------------------!
  !   Refresh boundary (ghost) values using the updated gradients !
  !--------------------------------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .ne. AMBIENT) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        p % n(c2) = p % n(c1) + p % x(c1) * Grid % dx(s)  &
                              + p % y(c1) * Grid % dy(s)  &
                              + p % z(c1) * Grid % dz(s)
      end do  ! faces
    end if    ! ambient or other
  end do      ! regions


  ! Perform Gauss from gradients which are good inside obtained above
  call Flow % Grad_Gauss_Variable(p)

  call Work % Disconnect_Int_Cell(c_computed, c_visited)

  call Profiler % Stop('Grad_Gauss_Pressure')

  end subroutine
