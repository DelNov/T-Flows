!==============================================================================!
  subroutine Field_Mod_Grad_Gauss_Pressure(flow, p)
!------------------------------------------------------------------------------!
!   Tries to find pressure gradients with Gaussian in an iterative fashion.    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: c_at_bnd => i_cell_01,  &  ! cell at boundary?
                      c_cnt    => i_cell_02      ! counting neighbours
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                                         (allocates Work_Mod)      !
!     |                                                                        !
!     +----> Compute_Pressure                        (doesn't use Work_Mod)    !
!              |                                                               !
!              +----> Field_Mod_Grad_Gauss_Pressure  (i_cell_01..02 are safe)  !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type), target :: flow
  type(Var_Type),   target :: p     ! should be pressure or pressure correction
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2
  integer                  :: i_fac, c1_prim, c2_prim, s_prim
  real                     :: res
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: YES = 1
  integer, parameter :: NO  = YES-1
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Extrapolation to boundaries
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
        p % n(c2) = p % n(c1) + p % x(c1) * grid % dx(s)  &
                              + p % y(c1) * grid % dy(s)  &
                              + p % z(c1) * grid % dz(s)
      end if
    end if
  end do

  ! Initialize with some gradients with the most robust and reliable tool
  ! you have at your disposal - least square cell-based method.  These
  ! gradients should be properly calculated inside the domain.
  call Field_Mod_Grad_Least_Pressure(flow, p)

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

      ! Nullify its gradients, and the counter for cells from which you interpolate
      p % x(c) = 0.0
      p % y(c) = 0.0
      p % z(c) = 0.0
      c_cnt(c) = 0    ! probably not needed, initialized above

      ! Browse through this cell's faces
      do i_fac = 1, grid % cells_n_faces(c)
        s_prim  = grid % cells_f(i_fac, c)
        c1_prim = grid % faces_c(1, s_prim)
        c2_prim = grid % faces_c(2, s_prim)

        ! Consider c1_prim if it is not cell at boundary
        if(c_at_bnd(c1_prim) .eq. NO) then
          p % x(c) = p % x(c) + p % x(c1_prim)
          p % y(c) = p % y(c) + p % y(c1_prim)
          p % z(c) = p % z(c) + p % z(c1_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

        ! Consider c2_prim if it is not cell at boundary
        ! and not a boundary cell.
        if(c_at_bnd(c2_prim) .eq. NO .and.  &
           c2_prim .gt. 0) then
          p % x(c) = p % x(c) + p % x(c2_prim)
          p % y(c) = p % y(c) + p % y(c2_prim)
          p % z(c) = p % z(c) + p % z(c2_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

      end do

      ! Work out the average
      if(c_cnt(c) .gt. 0) then
        p % x(c) = p % x(c) / c_cnt(c)
        p % y(c) = p % y(c) / c_cnt(c)
        p % z(c) = p % z(c) / c_cnt(c)
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

    ! Cell is at the boundary, and hasn'p been treated yet
    if(c_at_bnd(c) .eq. YES .and. c_cnt(c) .eq. 0) then

      ! Browse through this cell's faces
      do i_fac = 1, grid % cells_n_faces(c)
        s_prim  = grid % cells_f(i_fac, c)
        c1_prim = grid % faces_c(1, s_prim)
        c2_prim = grid % faces_c(2, s_prim)

        if(c1_prim .ne. c  .and.  &     ! skip your own self
           c_cnt(c1_prim) .gt. 0) then  ! consider only cells with values
          p % x(c) = p % x(c) + p % x(c1_prim)
          p % y(c) = p % y(c) + p % y(c1_prim)
          p % z(c) = p % z(c) + p % z(c1_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

        if(c2_prim .ne. c  .and.  &     ! skip your own self
           c2_prim .gt. 0  .and.  &     ! consider only cells with values
           c_cnt(c2_prim) .gt. 0) then  ! don'p use boundary cells
          p % x(c) = p % x(c) + p % x(c2_prim)
          p % y(c) = p % y(c) + p % y(c2_prim)
          p % z(c) = p % z(c) + p % z(c2_prim)
          c_cnt(c) = c_cnt(c) + 1
        end if

      end do

      if(c_cnt(c) .gt. 0) then
        p % x(c) = p % x(c) / c_cnt(c)
        p % y(c) = p % y(c) / c_cnt(c)
        p % z(c) = p % z(c) / c_cnt(c)
      end if

    end if  ! c at bnd

  end do

  ! Perform Gauss from gradients which are good inside obtained above
  call Field_Mod_Grad_Gauss_Variable(flow, p)

  end subroutine
