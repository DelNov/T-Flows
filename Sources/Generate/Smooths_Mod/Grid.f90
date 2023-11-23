!==============================================================================!
  subroutine Smooths_Mod_Grid(smr, Grid)
!------------------------------------------------------------------------------!
!>  Smooths the grid lines by a Laplacian-like algorithm inside Generate.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Memory allocation: Allocates memory for temporary arrays used in the     !
!     smoothing process.                                                       !
!   * Node connectivity: Establishes the connectivity between nodes by         !
!     identifying neighboring nodes for each grid node.                        !
!   * Exclusion of nodes: Identifies nodes that should not be moved during the !
!     smoothing process, based on criteria specified in smr.                   !
!   * Boundary handling: Ensures that nodes on the boundary faces are not      !
!     altered during the smoothing process.                                    !
!   * Smoothing iterations: Performs a specified number of smoothing           !
!     iterations (smr % iters(reg)). In each iteration, the new coordinates    !
!     of each node are calculated as the average of its neighboring nodes'     !
!     coordinates.                                                             !
!   * Relaxation factor: Applies a relaxation factor (smr % relax(reg)) to     !
!     blend the new and old positions, controlling the node movement.          !
!   * Region-specific smoothing: The smoothing can be applied selectively in   !
!     different regions and along different axes (x,y,z), as specified in smr. !
!   * Updating coordinates: Updates the coordinates of the grid nodes based    !
!     on the new calculated positions, ensuring that the changes are within    !
!     the defined limits.                                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Smooths_Type) :: smr   !! smoothing regions
  type(Grid_Type)    :: Grid  !! grid being smoothed
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, n, s, c2, i, j, k, m
  real                 :: x_new_tmp, y_new_tmp, z_new_tmp
  real                 :: x_max, y_max, z_max, x_min, y_min, z_min
  real                 :: x1, y1, z1, x8, y8, z8
  integer              :: reg
  real, allocatable    :: x_node_new(:), y_node_new(:), z_node_new(:)
  integer, allocatable :: node_to_nodes(:,:)
!==============================================================================!

  ! Allocate memory for additional arrays
  allocate(x_node_new(Grid % max_n_nodes));         x_node_new    = 0
  allocate(y_node_new(Grid % max_n_nodes));         y_node_new    = 0
  allocate(z_node_new(Grid % max_n_nodes));         z_node_new    = 0
  allocate(node_to_nodes(Grid % max_n_nodes,0:40)); node_to_nodes = 0

  print *, '# Now smoothing the cells. This may take a while !'

  ! Node connectivity
  do n = 1, Grid % n_nodes
    node_to_nodes(n,0) = 0
  end do

  x_max=-HUGE
  y_max=-HUGE
  z_max=-HUGE
  x_min=+HUGE
  y_min=+HUGE
  z_min=+HUGE
  do n = 1, Grid % n_nodes
    x_max=max(Grid % xn(n), x_max)
    y_max=max(Grid % yn(n), y_max)
    z_max=max(Grid % zn(n), z_max)
    x_min=min(Grid % xn(n), x_min)
    y_min=min(Grid % yn(n), y_min)
    z_min=min(Grid % zn(n), z_min)
  end do

  !-----------------------!
  !   Connect the nodes   !
  !-----------------------!
  do c = 1, Grid % n_cells         ! through cells
    do i = 1, 8                    ! through nodes of a cell
      n = Grid % cells_n(i,c)      ! first cell
      do j = 1, 8                  ! through nodes of a cell
        m = Grid % cells_n(j,c)    ! second cell
        if(n .ne.  m) then
          do k=1,node_to_nodes(n,0)
            if(node_to_nodes(n,k) .eq. m) goto 10
          end do
          node_to_nodes(n,0)=node_to_nodes(n,0)+1
          node_to_nodes(n,node_to_nodes(n,0)) = m
        end if
10    end do
    end do
  end do

  !-----------------------------------------!
  !   Exclude nodes which should not move   !
  !-----------------------------------------!
  do reg = 1, smr % n_smooths
    if( ( .not. smr % in_x(reg) ) .and.  &
        ( .not. smr % in_y(reg) ) .and.  &
        ( .not. smr % in_z(reg) ) ) then
      do n = 1, Grid % n_nodes
        x1 = smr % x_min(reg)
        y1 = smr % y_min(reg)
        z1 = smr % z_min(reg)
        x8 = smr % x_max(reg)
        y8 = smr % y_max(reg)
        z8 = smr % z_max(reg)
        if( (x1 <= Grid % xn(reg)) .and. (Grid % xn(reg) <= x8) .and. &
            (y1 <= Grid % yn(reg)) .and. (Grid % yn(reg) <= y8) .and. &
            (z1 <= Grid % zn(reg)) .and. (Grid % zn(reg) <= z8) ) then
          node_to_nodes(n,0) = 0
        end if
      end do
    end if
  end do

  do s = 1, Grid % n_faces               ! boundary through faces
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) then
      do i = 1, Grid % faces_n_nodes(s)  ! through nodes of a face
        n = Grid % faces_n(i, s)
        node_to_nodes(n,0) = 0
      end do
    end if
  end do

  !---------------------!
  !   Smooth the Grid   !
  !---------------------!
  do reg = 1, smr % n_smooths
    print *, '# Now smoothing region ',reg,' with:',  &
                smr % iters(reg), ' iterations.'

    do j = 1, smr % iters(reg)

      ! Calculate new coordinates using the old values
      do n = 1, Grid % n_nodes
        if(node_to_nodes(n,0) > 0) then
          x_new_tmp=0.0
          y_new_tmp=0.0
          z_new_tmp=0.0
          do i = 1, node_to_nodes(n,0)
            x_new_tmp = x_new_tmp + Grid % xn(node_to_nodes(n,i))
            y_new_tmp = y_new_tmp + Grid % yn(node_to_nodes(n,i))
            z_new_tmp = z_new_tmp + Grid % zn(node_to_nodes(n,i))
          end do
          x_new_tmp = x_new_tmp / real(node_to_nodes(n,0))
          y_new_tmp = y_new_tmp / real(node_to_nodes(n,0))
          z_new_tmp = z_new_tmp / real(node_to_nodes(n,0))
          if(Grid % xn(n) > 0.001*x_min .and.                  &
             Grid % xn(n) < 0.999*x_max)                       &
          x_node_new(n) = (1.0-smr % relax(reg))*Grid % xn(n)  &
                        +      smr % relax(reg) *x_new_tmp
          if(Grid % yn(n) > 0.001*y_min .and.                  &
             Grid % yn(n) < 0.999*y_max)                       &
          y_node_new(n) = (1.0-smr % relax(reg))*Grid % yn(n)  &
                        +      smr % relax(reg) *y_new_tmp
          if(Grid % zn(n) > 0.001*z_min .and.                  &
             Grid % zn(n) < 0.999*z_max)                       &
          z_node_new(n) = (1.0-smr % relax(reg))*Grid % zn(n)  &
                        +      smr % relax(reg)* z_new_tmp
        end if
      end do  ! through nodes

      ! Update coordinates
      do n = 1, Grid % n_nodes
        if(node_to_nodes(n,0)   >  0) then

          x1 = smr % x_min(reg)
          y1 = smr % y_min(reg)
          z1 = smr % z_min(reg)
          x8 = smr % x_max(reg)
          y8 = smr % y_max(reg)
          z8 = smr % z_max(reg)

          if( (x1 <= Grid % xn(n)) .and. (Grid % xn(n) <= x8) .and.  &
              (y1 <= Grid % yn(n)) .and. (Grid % yn(n) <= y8) .and.  &
              (z1 <= Grid % zn(n)) .and. (Grid % zn(n) <= z8) ) then

            if(smr % in_x(reg)) then
              if(Grid % xn(n) > 0.001*x_min .and. &
                 Grid % xn(n) < 0.999*x_max)      &
              Grid % xn(n)=x_node_new(n)
            end if

            if(smr % in_y(reg)) then
              if(Grid % yn(n) > 0.001*y_min .and.  &
                 Grid % yn(n) < 0.999*y_max)      &
              Grid % yn(n)=y_node_new(n)
            end if

            if(smr % in_z(reg)) then
              if(Grid % zn(n) > 0.001*z_min .and.  &
                 Grid % zn(n) < 0.999*z_max)      &
              Grid % zn(n)=z_node_new(n)
            end if

          end if  ! if the point belongs to region
        end if    ! if the point is not excluded from smoothing
      end do      ! through nodes
    end do
  end do

  deallocate(x_node_new)
  deallocate(y_node_new)
  deallocate(z_node_new)
  deallocate(node_to_nodes)

  end subroutine
