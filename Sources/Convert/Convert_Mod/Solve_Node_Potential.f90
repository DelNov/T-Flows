!==============================================================================!
  subroutine Solve_Node_Potential(Convert, Grid, fixed_nodes, phi, tol, m_iter)
!------------------------------------------------------------------------------!
!>  Solves a distance-like scalar potential on grid nodes by relaxing a
!>  graph-Laplacian problem with a constant source term.
!>
!>  Nodes marked in fixed_nodes keep phi = 0.  All other nodes solve:
!>
!>      sum_j w_ij * (phi_i - phi_j) = 1
!>
!>  where j are neighbouring nodes connected by grid edges.  The resulting
!>  potential grows away from fixed nodes and can be used to identify nodes
!>  which are farthest from prescribed boundary-layer displacements.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert                      !! parent class
  type(Grid_Type)     :: Grid                         !! converted grid
  logical             :: fixed_nodes(Grid % n_nodes)  !! boundary conditions
  real                :: phi(Grid % n_nodes)          !! unknown potential
  real, intent(in)    :: tol                          !! target tolerance
  integer             :: m_iter                       !! max. iterations
!-----------------------------------[Locals]-----------------------------------!
  integer           :: iter, e, n, n1, n2
  real              :: lx, ly, lz, w
  real              :: corr, res, res_ref, max_phi
  real, allocatable :: phi_old(:)
  real, allocatable :: sum_w(:)
!==============================================================================!

  Assert(allocated(Grid % edges_n))
  Assert(Grid % n_edges .gt. 0)
  Assert(any(fixed_nodes))

  allocate(phi_old(Grid % n_nodes));  phi_old(:) = 0.0
  allocate(sum_w  (Grid % n_nodes));  sum_w  (:) = 0.0

  ! Dirichlet condition for potential
  do n = 1, Grid % n_nodes
    if(fixed_nodes(n)) phi(n) = 0.0
  end do

  !-------------------------------------------!
  !   Jacobi iterations for graph Laplacian   !
  !-------------------------------------------!
  do iter = 1, m_iter

    phi_old(:) = phi(:)

    phi  (:) = 0.0
    sum_w(:) = 0.0

    max_phi = -HUGE

    ! Accumulate weighted neighbour averages
    do e = 1, Grid % n_edges
      n1 = Grid % edges_n(1, e)
      n2 = Grid % edges_n(2, e)

      Assert(n1 .gt. 0)
      Assert(n2 .gt. 0)

      lx = Grid % xn(n1) - Grid % xn(n2)
      ly = Grid % yn(n1) - Grid % yn(n2)
      lz = Grid % zn(n1) - Grid % zn(n2)

      w = 1.0 / max(lx*lx + ly*ly + lz*lz, TINY)

      phi(n1) = phi(n1) + w * phi_old(n2)
      sum_w(n1) = sum_w(n1) + w

      phi(n2) = phi(n2) + w * phi_old(n1)
      sum_w(n2) = sum_w(n2) + w
    end do

    res     = 0.0
    res_ref = TINY

    ! Compute potential for free nodes only
    do n = 1, Grid % n_nodes
      if(.not. fixed_nodes(n)) then
        Assert(sum_w(n) .gt. 0.0)

        phi(n) = (1.0 + phi(n)) / sum_w(n)

        max_phi = max(phi(n), max_phi)

        corr = (phi(n) - phi_old(n))**2

        res     = max(res,     sqrt(corr))
        res_ref = max(res_ref, abs(phi(n)))
      end if
    end do

    res = res / res_ref

    ! Restore prescribed zero potential
    do n = 1, Grid % n_nodes
      if(fixed_nodes(n)) phi(n) = 0.0
    end do

    if(      iter .eq. 1           &
        .or. iter .eq. m_iter      &
        .or. mod(iter, 100) .eq. 0  &
        .or. res .lt. tol) then
      print '(a,i6,a,1pe12.4)',  &
        ' # Solve node potential iteration: ', iter,  &
        ' residual: ', res
    end if

    if(res .lt. tol) then
      Assert(max_phi .gt. 0.0)

        ! Normalize the solution ...
        do n = 1, Grid % n_nodes
          if(.not. fixed_nodes(n)) then
            phi(n) = phi(n) / max_phi
          end if
        end do

      ! ... and exit
      exit
    end if

  end do

  end subroutine
