!==============================================================================!
  subroutine Refines_Mod_Connectivity(ref, Grid, real_run)
!------------------------------------------------------------------------------!
!>  Determines the topology of the cells, faces and boundary cells.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine initializes the face-node numbering for hexahedral cells, !
!     counts the boundary cells, and adjusts the dimensions of the cells_c     !
!     array in the Grid.                                                       !
!   * It then loops through the grid cells to identify neighboring cells and   !
!     boundary cells, assigning their relationship in terms of connectivity    !
!     and face information.                                                    !
!   * A two-pass run (run = 1, 2) ensures that all necessary relationships are !
!     established, including determining which volumes are connected by each   !
!     face and identifying the nodes that constitute each face.                !
!   * During a 'real run' (real_run = .true.), boundary cells are marked and   !
!     allocated, whereas, in a non-real run, the boundary cell information is  !
!     restored to its initial state.                                           !
!   * The subroutine also checks for sufficient allocated memory for the faces !
!     and outputs an error if the allocated memory is exceeded, indicating     !
!     the need for adjustments in the domain file and recompilation.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type)  :: ref       !! type holding information on refinement
  type(Grid_Type)     :: Grid      !! grid being generated (refined here)
  logical, intent(in) :: real_run  !! false for trial run, true for real run
!-----------------------------------[Locals]-----------------------------------!
  integer                 :: c, c1, c2, m, run
  integer, dimension(6,4) :: fn  ! link faces' nodes for a hexahedral cell
                                 ! not quite the same as Gambit's neu numbering
!------------------------------------------------------------------------------!
  include 'Block_Numbering.h90'
!==============================================================================!

  ! Copy face-node numbering for faces
  ! (Note that it is the same as for blocks here in Generator)
  fn = hex_block

  print *, '# Now determining the topology. This may take a while !'

  !------------------------------!
  !                              !
  !   Count the boundary cells   !
  !                              !
  !------------------------------!

  ! Initialize number of boundary cells
  Grid % n_bnd_cells = 0

  ! Adjust dimension of cells_c
  call Enlarge % Matrix_Int(Grid % cells_c, i=(/1,24/))

  do c = 1, Grid % n_cells
    do m = 1, 24   ! neighbour cells
      if(Grid % cells_c(m,c) < 0) then
        Grid % n_bnd_cells = Grid % n_bnd_cells + 1
      end if
    end do
  end do

  ! Adjust dimension for all cell-based variables
  ! At this stage, it will basically adjust number of boundary cells
  call Grid % Allocate_Cells(Grid % n_cells, Grid % n_bnd_cells)

  ! Re-initialize number of boundary cells
  Grid % n_bnd_cells = 0

  do c = 1, Grid % n_cells
    do m = 1, 24   ! neighbour cells
      if(Grid % cells_c(m,c) < 0) then
        Grid % n_bnd_cells = Grid % n_bnd_cells + 1

        ! Remember the boundary region, take positive value for region
        Grid % region % at_cell(-Grid % n_bnd_cells) =  -Grid % cells_c(m,c)

        ! Put new boundary cell into place
        Grid % cells_c(m,c) = -Grid % n_bnd_cells
      end if
    end do
  end do

  !---------------------------!
  !   Create the array with   !
  !   information on faces    !
  !---------------------------!

  call Refines_Mod_Allocate_Cells(ref, Grid % n_cells, Grid % n_bnd_cells)

  ! Initialize number of sides
  Grid % n_faces = 0

  do run = 1, 2
    do c1 = 1, Grid % n_cells
      do m = 1, 24 ! through all the neighbouring cells
        c2 = Grid % cells_c(m, c1)
        if( (run .eq. 1) .and. (c2 > c1)  .or.  &
            (run .eq. 2) .and. (c2 < 0) ) then
          Grid % n_faces = Grid % n_faces + 1
        end if
      end do
    end do
  end do

  ! Re-initialize number of sides
  Grid % n_faces = 0

  ! Adjust dimension of faces_n
  call Enlarge % Matrix_Int(Grid % faces_n, i=(/1,4/))

  do run = 1, 2

    do c1 = 1, Grid % n_cells
      do m = 1, 24 ! through all the neighbouring cells
        c2 = Grid % cells_c(m, c1)

        if( (run .eq. 1) .and. (c2 > c1)  .or.  &
            (run .eq. 2) .and. (c2 < 0) ) then
          Grid % n_faces = Grid % n_faces + 1

          ! What counts in this expansion, is the number of faces, not shadows
          call Grid % Allocate_Faces(Grid % n_faces, Grid % n_shadows,  &
                                     GROWTH_MARGIN)

          ! Which volumes are connected with side Grid % n_faces
          Grid % faces_c(1, Grid % n_faces) = c1
          Grid % faces_c(2, Grid % n_faces) = c2

          ! Which is c2 neighbour of c1 and vice versa
          do c = 1,24
            if(Grid % cells_c(c,c1) .eq. c2) then
              Grid % face_c_to_c(Grid % n_faces,1) = c
            end if

            if(c2 > 0) then
              if(Grid % cells_c(c,c2) .eq. c1) then
                Grid % face_c_to_c(Grid % n_faces,2) = c
              end if
            end if
          end do

          ! Nodes of a side Grid % n_faces
          if(c2 > 0) then
            if(ref % cell_level(c2) > ref % cell_level(c1)) then
              Grid % faces_n(1,Grid % n_faces) =  &
                Grid % cells_n( fn(Grid % face_c_to_c(Grid % n_faces,2),4), c2 )
              Grid % faces_n(2,Grid % n_faces) =  &
                Grid % cells_n( fn(Grid % face_c_to_c(Grid % n_faces,2),3), c2 )
              Grid % faces_n(3,Grid % n_faces) =  &
                Grid % cells_n( fn(Grid % face_c_to_c(Grid % n_faces,2),2), c2 )
              Grid % faces_n(4,Grid % n_faces) =  &
                Grid % cells_n( fn(Grid % face_c_to_c(Grid % n_faces,2),1), c2 )
            else
              Grid % faces_n(1,Grid % n_faces) = Grid % cells_n( fn(m,1), c1 )
              Grid % faces_n(2,Grid % n_faces) = Grid % cells_n( fn(m,2), c1 )
              Grid % faces_n(3,Grid % n_faces) = Grid % cells_n( fn(m,3), c1 )
              Grid % faces_n(4,Grid % n_faces) = Grid % cells_n( fn(m,4), c1 )
            end if
          else
            Grid % faces_n(1,Grid % n_faces) = Grid % cells_n( fn(m,1), c1 )
            Grid % faces_n(2,Grid % n_faces) = Grid % cells_n( fn(m,2), c1 )
            Grid % faces_n(3,Grid % n_faces) = Grid % cells_n( fn(m,3), c1 )
            Grid % faces_n(4,Grid % n_faces) = Grid % cells_n( fn(m,4), c1 )
          end if

        end if
      end do   ! m
    end do     ! c1
  end do       ! run

  if(.not. real_run) then
    Grid % n_bnd_cells = 0
    do c = 1, Grid % n_cells
      do m = 1, 24   ! neighbour cells
        if(Grid % cells_c(m,c) < 0) then
          Grid % n_bnd_cells  = Grid % n_bnd_cells + 1

          ! Restore the boundary region, take positive value for region
          Grid % cells_c(m,c) = -Grid % region % at_cell(-Grid % n_bnd_cells)
        end if
      end do
    end do
  end if

  end subroutine
