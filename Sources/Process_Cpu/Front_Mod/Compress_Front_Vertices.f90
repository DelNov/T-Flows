!==============================================================================!
  subroutine Compress_Front_Vertices(Front, verbose)
!------------------------------------------------------------------------------!
!>  This subroutine optimizes the front's vertex list by eliminating redundant
!>  vertices. This compression is essential for efficient processing and memory
!>  usage in complex simulations.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Verifying the uniqueness of elements' vertices.                          !
!   * Sorting vertices based on coordinates and associated node numbers.       !
!   * Compressing vertices that coincide spatially.                            !
!   * Correcting vertex coordinates and updating element vertices accordingly. !
!   * Final verification to ensure no duplicated vertices in elements.         !
!   * Calculating the total number of unique vertices, essential for           !
!     simulations distributed across multiple processors.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front    !! parent class
  logical                   :: verbose  !! controls the output verbosity
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, n_vert, i, j, i_ver, j_ver, nv_tot
  real,    allocatable     :: xv(:), yv(:), zv(:)
  integer, allocatable     :: ni(:), new_n(:), n1(:), n2(:)
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  Elem => Front % Elem

  !-----------------------------------------!
  !   Check sanity of the elements so far   !
  !-----------------------------------------!
  do e = 1, ne
    do i_ver = 1, Elem(e) % nv-1
      do j_ver = i_ver+1, Elem(e) % nv
        i = Elem(e) % v(i_ver)
        j = Elem(e) % v(j_ver)
        if(i .eq. j) then
          call Message % Error(44,                                    &
                               "Error in the sanity check for  "  //  &
                               "elements. \n  Some element(s)  "  //  &
                               "have duplicate vertices.",            &
                               file=__FILE__, line=__LINE__)
        end if
      end do
    end do
  end do

  !----------------------------------------!
  !   Sort vertices by their coordinates   !
  !----------------------------------------!
  if(nv > 0) then
    allocate(xv(nv));     xv    = 0.0
    allocate(yv(nv));     yv    = 0.0
    allocate(zv(nv));     zv    = 0.0
    allocate(ni(nv));     ni    = 0
    allocate(n1(nv));     n1    = 0
    allocate(n2(nv));     n2    = 0
    allocate(new_n(nv));  new_n = 0

    do v = 1, nv
      ! Store vertex coordinates ...
      xv(v) = Vert(v) % x_n
      yv(v) = Vert(v) % y_n
      zv(v) = Vert(v) % z_n

      ! ... vertex numbers ...
      ni(v) = v

      ! ... and nodes from which the vertex was created
      n1(v) = min(Front % b_node_1(v), Front % b_node_2(v))
      n2(v) = max(Front % b_node_1(v), Front % b_node_2(v))
    end do
    call Sort % Two_Int_Carry_Int(n1, n2, ni)
  end if

  !-----------------------------------------------------------!
  !   Compress the vertices which fall on top of each other   !
  !-----------------------------------------------------------!
  if(nv > 0) then
    n_vert = 1
    new_n(1) = n_vert
    do v = 2, nv
      if( (n1(v) .ne. n1(v-1)) .or. (n2(v) .ne. n2(v-1)) ) then
        n_vert = n_vert + 1
      end if
      new_n(v) = n_vert
    end do
  else
    n_vert = 0
  end if

  !----------------------------------------------------------!
  !   Correct vertices' coordinates and elements' vertices   !
  !----------------------------------------------------------!
  if(nv > 0) then
    call Sort % Int_Carry_Int(ni, new_n)

    ! Coordinates
    do v = 1, nv
      Vert(new_n(v)) % x_n = xv(v)
      Vert(new_n(v)) % y_n = yv(v)
      Vert(new_n(v)) % z_n = zv(v)
    end do

    ! Vertices
    do e = 1, ne
      do i_ver = 1, Elem(e) % nv
        Elem(e) % v(i_ver) = new_n(Elem(e) % v(i_ver))
      end do
    end do
  end if

  !--------------------------------------------!
  !   Work out compressed number of vertices   !
  !--------------------------------------------!
  nv     = n_vert
  nv_tot = n_vert
  if(Front % mesh_divided) then
    call Global % Sum_Int(nv_tot)
  end if
  if(verbose .and. First_Proc()) then
    print '(a40,i8)', ' # Compressed number of vertices:       ', nv_tot
  end if

  !---------------------------------------------!
  !   Check sanity of the elements in the end   !
  !---------------------------------------------!
  do e = 1, ne
    do i_ver = 1, Elem(e) % nv-1
      do j_ver = i_ver+1, Elem(e) % nv
        i = Elem(e) % v(i_ver)
        j = Elem(e) % v(j_ver)
        if(i .eq. j) then
          call Message % Error(44,                                          &
                               "Error in the final sanity check for  "  //  &
                               "elements. \n  Some element(s)        "  //  &
                               "have duplicate vertices.",                  &
                               file=__FILE__, line=__LINE__)
        end if
      end do
    end do
  end do

  end subroutine
