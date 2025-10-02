!==============================================================================!
  subroutine Load_Forrest(Convert, Grid, file_name)
!------------------------------------------------------------------------------!
!>  The Load_Forrest subroutine in the Conver is designed to read a forrest
!>  represented in STL file format and populate the Grid structure with the
!>  information. This subroutine is particularly focused on creating a forest
!>  scene by placing individual trees (represented by a single normalized
!>  tree in Grid(1)) within a defined area.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initial Setup: The subroutine begins with initializing variables and     !
!     reading the forest data from the STL file using                          !
!     Forrest % Create_From_File.                                              !
!   * Determining Forest Span: It calculates the minimum and maximum coords    !
!     (x,y,z) of the forrest to understand the overall span and dimensions.    !
!   * Scaling Single Tree: The single tree in Grid(1) is scaled according to   !
!     the height of the forest.                                                !
!   * Estimating Number of Trees: The subroutine estimates the number of trees !
!     that can be planted within the forest area based on the single tree size !
!     and forest area.                                                         !
!   * Tree Placement:                                                          !
!     - Allocates arrays t_x and t_y for storing the x and y coordinates of    !
!       each tree.                                                             !
!     - Randomly places trees within the forest area, ensuring they are not    !
!       too close to each other.                                               !
!   * Planting Trees:                                                          !
!     - The subroutine plants the estimated number of trees (n_trees) by       !
!       copying the data of the single tree (Grid(1)) into the forest grid     !
!       (Grid(2)).                                                             !
!      - It allocates node'n'faces for Grid(2) and adjusts them for each tree. !
!      - For each tree, it modifies the coordinates (xn, yn, zn) to position   !
!        the tree correctly in the forest scene.                               !
!      - The height of each tree is slightly randomized to add variety.        !
!   * Avoiding Clustering: The subroutine includes logic to prevent trees from !
!     being placed too close to each other, ensuring a more natural and        !
!     realistic forest layout.                                                 !
!   * Performance Optimization: It uses random number generation and a limit   !
!     on the number of placement attempts to efficiently populate the forest   !
!     without excessive computational effort.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert    !! parent class
  type(Grid_Type)     :: Grid(2)    !! Grid(1) is a single tree and Grid(2)
                                    !! is a forest scene
  character(SL)       :: file_name  !! file name
!-----------------------------------[Locals]-----------------------------------!
  type(Stl_Type)    :: Forrest
  integer           :: s, i_nod, v, fac, attempts
  integer           :: n_trees, t1, t2, t, n0, n1, n2
  real              :: f_min_x, f_max_x, f_min_y, f_max_y, f_min_z, f_max_z
  real              :: n(3)
  real              :: delta_x, delta_y, min_dist, tmp
  real, allocatable :: t_x(:), t_y(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Load_Forrest')

  !----------------------!
  !                      !
  !   Read the forrest   !
  !                      !
  !----------------------!
  call Forrest % Create_From_File(file_name)

  f_min_x =  HUGE;  f_max_x = -HUGE
  f_min_y =  HUGE;  f_max_y = -HUGE
  f_min_z =  HUGE;  f_max_z = -HUGE
  do fac = 1, Forrest % N_Facets()
    do v = 1, 3
      n = Forrest % Facets_Vert_Coords(fac, v)
      f_min_x = min(f_min_x, n(1));  f_max_x = max(f_max_x, n(1))
      f_min_y = min(f_min_y, n(2));  f_max_y = max(f_max_y, n(2))
      f_min_z = min(f_min_z, n(3));  f_max_z = max(f_max_z, n(3))
    end do
  end do
  print '(a,2f9.3)', ' # Forrest span in x directon: ', f_min_x, f_max_x
  print '(a,2f9.3)', ' # Forrest span in y directon: ', f_min_y, f_max_y
  print '(a,2f9.3)', ' # Forrest span in z directon: ', f_min_z, f_max_z

  ! Scale the single tree by forrest height
  Grid(1) % zn(:) = Grid(1) % zn - f_min_z
  Grid(1) % xn(:) = Grid(1) % xn * (f_max_z - f_min_z)
  Grid(1) % yn(:) = Grid(1) % yn * (f_max_z - f_min_z)
  Grid(1) % zn(:) = Grid(1) % zn * (f_max_z - f_min_z)

  !---------------------!
  !   Place the seeds   !
  !---------------------!

  ! The size a single tree occuppies
  delta_x = maxval(Grid(1) % xn) - minval(Grid(1) % xn)
  delta_y = maxval(Grid(1) % yn) - minval(Grid(1) % yn)

  ! Estimate (roughly) the number of trees you want to plant
  ! (There are probably better and smarter ways to do this.)
  n_trees = floor( (f_max_x - f_min_x)*(f_max_y-f_min_y)  &
                   / min(delta_x, delta_y)**2 )
  print '(a,i4,a)', ' # Attempting to plant ', n_trees, ' trees'

  ! Allocate memory for trees' coordinates
  allocate(t_x(n_trees));  t_x(:) = 0.0
  allocate(t_y(n_trees));  t_y(:) = 0.0

  ! Place individual trees ...
  attempts = 0
  do t1 = 1, n_trees
3   continue
    call random_number(tmp);  t_x(t1) = f_min_x + delta_x / 2.0  &
                                      + tmp * (f_max_x - f_min_x - delta_x)
    call random_number(tmp);  t_y(t1) = f_min_y + delta_y / 2.0  &
                                      + tmp * (f_max_y - f_min_y - delta_y)

    ! ... making sure they are not too close to one another
    min_dist = HUGE
    do t2 = 1, t1 - 1
      min_dist = min(min_dist, Math % Distance(t_x(t1), t_y(t1), 0.0,  &
                                               t_x(t2), t_y(t2), 0.0))
    end do
    if(min_dist < 0.5 * max(delta_x, delta_y)) then
      attempts = attempts + 1
      if(attempts > 1000 * n_trees) goto 4
        n_trees = t1 - 1
      goto 3
    else
      attempts = 0
    end if
  end do
4 continue
  print '(a,i4,a)', ' # Eventually planted  ', n_trees, ' trees'

  !------------------------------------------------------------------!
  !   Plant the trees, one by one, by copying Grid(1) into Grid(2)   !
  !------------------------------------------------------------------!
  Grid(2) % n_nodes = Grid(1) % n_nodes * n_trees
  Grid(2) % n_faces = Grid(1) % n_faces * n_trees
  Grid(2) % name    = Grid(1) % name
  call Grid(2) % Allocate_Nodes(Grid(2) % n_nodes)
  call Grid(2) % Allocate_Faces(Grid(2) % n_faces, 0)
  do t = 0, n_trees-1
    n0 = t * Grid(1) % n_faces
    do s = 1, Grid(1) % n_faces
      Grid(2) % faces_n_nodes(n0+s) = Grid(1) % faces_n_nodes(s)
      call Enlarge % Matrix_Int(Grid(2) % faces_n,  &
                                i=(/1,Grid(2) % faces_n_nodes(n0+s)/))
      do i_nod = 1, Grid(2) % faces_n_nodes(s)
        Grid(2) % faces_n(i_nod, n0 + s) =  &
        Grid(1) % faces_n(i_nod, s) + Grid(1) % n_nodes * t
      end do
    end do
    n1 = t * Grid(1) % n_nodes + 1
    n2 = t * Grid(1) % n_nodes + Grid(1) % n_nodes

    Grid(2) % xn(n1:n2) = Grid(1) % xn(1:Grid(1) % n_nodes) + t_x(t+1)
    Grid(2) % yn(n1:n2) = Grid(1) % yn(1:Grid(1) % n_nodes) + t_y(t+1)

    ! Randomize the height a little bit
    call random_number(tmp)
    tmp = 0.9 + 0.2 * tmp
    Grid(2) % zn(n1:n2) = Grid(1) % zn(1:Grid(1) % n_nodes) * tmp
  end do

  call Profiler % Stop('Load_Forrest')

  end subroutine
