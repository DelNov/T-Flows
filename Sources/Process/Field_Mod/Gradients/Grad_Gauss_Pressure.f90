!==============================================================================!
  subroutine Grad_Gauss_Pressure(Flow, p)
!------------------------------------------------------------------------------!
!   Tries to find pressure gradients with Gaussian in an iterative fashion.    !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow
  type(Var_Type),    target :: p     ! should be pressure or pressure correction
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),     pointer :: Grid
  integer                      :: s, c, c1, c2, iter
  integer                      :: i_fac, c1_prim, c2_prim, s_prim
  integer, contiguous, pointer :: c_computed(:), c_visited(:)
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: YES = 1
  integer, parameter :: NO  = YES-1
!==============================================================================!

  call Work % Connect_Int_Cell(c_computed, c_visited)

  ! Take alias
  Grid => Flow % pnt_grid

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
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) c_computed(c1) = NO
  end do

  call Grid % Exchange_Cells_Int(c_computed)

  ! Nullify gradients for cells near boundaries
  ! (where gradients are not computed properly),
  ! and initialize c_visited
  do c = 1, Grid % n_cells
    if(c_computed(c) .eq. NO) then  ! if not computed
      p % x(c) = 0.0
      p % y(c) = 0.0
      p % z(c) = 0.0
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
          p % x(c1) = p % x(c1) + p % x(c2)
          p % y(c1) = p % y(c1) + p % y(c2)
          p % z(c1) = p % z(c1) + p % z(c2)
          c_visited(c1) = c_visited(c1) + 1
        end if

        if(c_computed(c2) .eq. NO .and. c_computed(c1) .eq. YES) then
          p % x(c2) = p % x(c2) + p % x(c1)
          p % y(c2) = p % y(c2) + p % y(c1)
          p % z(c2) = p % z(c2) + p % z(c1)
          c_visited(c2) = c_visited(c2) + 1
        end if
      end if
    end do

    ! Browse throough cells, and calculate final values ...
    ! ... of gradients in the cells which have been visited
    do c = 1, Grid % n_cells
      if(c_visited(c) > 0) then
        p % x(c) = p % x(c) / c_visited(c)
        p % y(c) = p % y(c) / c_visited(c)
        p % z(c) = p % z(c) / c_visited(c)
        c_visited(c)  = 0
        c_computed(c) = YES  ! mark it as computed for the next iteration
      end if
    end do
    call Grid % Exchange_Cells_Real(p % x)
    call Grid % Exchange_Cells_Real(p % y)
    call Grid % Exchange_Cells_Real(p % z)

  end do    ! iter

  ! Perform Gauss from gradients which are good inside obtained above
  call Flow % Grad_Gauss_Variable(p)

  call Work % Disconnect_Int_Cell(c_computed, c_visited)

  end subroutine
