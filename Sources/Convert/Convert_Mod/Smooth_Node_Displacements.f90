!==============================================================================!
  subroutine Smooth_Node_Displacements(Convert,     &
                                       Grid,        &
                                       fixed_nodes, &
                                       dx, dy, dz,  &
                                       nx, ny, nz,  &
                                       tol, m_iter)
!------------------------------------------------------------------------------!
!>  Smooths nodal displacements by relaxing a graph-Laplacian problem on grid
!>  nodes.  Boundary nodes are assumed to be stored first, from 1 to cnt_bnd_n,
!>  and keep their prescribed displacements.  Interior nodes receive smoothed
!>  displacements from neighbouring nodes connected by grid edges.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert
  type(Grid_Type)     :: Grid                         !! grid being converted
  logical             :: fixed_nodes(Grid % n_nodes)  !! boundary conditions
  real                :: dx(Grid % n_nodes),  &
                         dy(Grid % n_nodes),  &
                         dz(Grid % n_nodes)
  real                :: nx(Grid % n_nodes),  &
                         ny(Grid % n_nodes),  &
                         nz(Grid % n_nodes)
  real, intent(in)    :: tol                          !! target tolerance
  integer             :: m_iter                       !! max. iterations
!-----------------------------------[Locals]-----------------------------------!
  integer           :: iter, e, n, n1, n2
  real              :: lx, ly, lz, w
  real              :: corr, disp, magn, res, res_ref
  real, allocatable :: dx_old(:), dy_old(:), dz_old(:)
  real, allocatable :: sum_w(:)
!==============================================================================!

  Assert(allocated(Grid % edges_n))
  Assert(Grid % n_edges .gt. 0)

  allocate(dx_old(Grid % n_nodes));  dx_old(:) = 0.0
  allocate(dy_old(Grid % n_nodes));  dy_old(:) = 0.0
  allocate(dz_old(Grid % n_nodes));  dz_old(:) = 0.0
  allocate(sum_w(Grid % n_nodes));   sum_w (:) = 0.0

  !-------------------------------------------!
  !   Jacobi iterations for graph Laplacian   !
  !-------------------------------------------!
  do iter = 1, m_iter

    dx_old(:) = dx(:)
    dy_old(:) = dy(:)
    dz_old(:) = dz(:)

    dx(:)    = 0.0
    dy(:)    = 0.0
    dz(:)    = 0.0
    sum_w(:) = 0.0

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

      dx(n1) = dx(n1) + w * dx_old(n2)
      dy(n1) = dy(n1) + w * dy_old(n2)
      dz(n1) = dz(n1) + w * dz_old(n2)
      sum_w(n1) = sum_w(n1) + w

      dx(n2) = dx(n2) + w * dx_old(n1)
      dy(n2) = dy(n2) + w * dy_old(n1)
      dz(n2) = dz(n2) + w * dz_old(n1)
      sum_w(n2) = sum_w(n2) + w
    end do

    res     = 0.0
    res_ref = TINY

    ! Compute smoothed values for interior nodes only
    do n = 1, Grid % n_nodes
      if(.not. fixed_nodes(n)) then
        Assert(sum_w(n) .gt. 0.0)

        dx(n) = dx(n) / sum_w(n)
        dy(n) = dy(n) / sum_w(n)
        dz(n) = dz(n) / sum_w(n)

        magn = nx(n)**2 + ny(n)**2 + nz(n)**2

        if(magn .gt. TINY) then
          disp = dx(n)*nx(n) + dy(n)*ny(n) + dz(n)*nz(n)
          dx(n) = dx(n) - disp * nx(n)
          dy(n) = dy(n) - disp * ny(n)
          dz(n) = dz(n) - disp * nz(n)
        end if

        corr = (dx(n) - dx_old(n))**2  &
             + (dy(n) - dy_old(n))**2  &
             + (dz(n) - dz_old(n))**2

        disp = dx(n)**2 + dy(n)**2 + dz(n)**2

        res     = max(res,     sqrt(corr))
        res_ref = max(res_ref, sqrt(disp))
      end if
    end do

    res = res / res_ref

    ! Restore prescribed boundary-node displacements
    do n = 1, Grid % n_nodes
      if(fixed_nodes(n)) then
        dx(n) = dx_old(n)
        dy(n) = dy_old(n)
        dz(n) = dz_old(n)
      end if
    end do

    if(      iter .eq. 1           &
        .or. iter .eq. m_iter      &
        .or. mod(iter, 100) .eq. 0  &
        .or. res .lt. tol) then
      print '(a,i6,a,1pe12.4)',  &
        ' # Smooth node displacement iteration: ', iter,  &
        ' residual: ', res
    end if

    if(res .lt. tol) exit

  end do

  !----------------!
  !   Move nodes   !
  !----------------!
  Grid % xn(:) = Grid % xn(:) + dx(:)
  Grid % yn(:) = Grid % yn(:) + dy(:)
  Grid % zn(:) = Grid % zn(:) + dz(:)

  end subroutine
