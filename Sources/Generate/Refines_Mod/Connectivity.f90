!==============================================================================!
  subroutine Refines_Mod_Connectivity(ref, Grid, real_run)
!------------------------------------------------------------------------------!
!   Determines the topology of the cells, faces and boundary cells.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type) :: ref
  type(Grid_Type)    :: Grid
  logical            :: real_run
!-----------------------------------[Locals]-----------------------------------!
  integer                 :: c, c1, c2, m, run
  integer, dimension(6,4) :: fn  ! link faces' nodes for a hexahedral cell
                                 ! not quite the same as Gambit's neu numbering
!------------------------------------------------------------------------------!
  include 'Block_Numbering.f90'
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
  Grid % n_bnd_cells = 0                     ! initialize nuber of bnd cells
  call Adjust_First_Dim(24, Grid % cells_c)  ! adjust dimension of cells_c

  do c = 1, Grid % n_cells
    do m = 1, 24   ! neighbour cells
      if(Grid % cells_c(m,c) < 0) then
        Grid % n_bnd_cells = Grid % n_bnd_cells + 1

        ! Remember the boundary color, take positive value for color
        Grid % bnd_cond % color(-Grid % n_bnd_cells) =  -Grid % cells_c(m,c)

        ! Put new boundary cell into place
        Grid % cells_c(m,c) = -Grid % n_bnd_cells
      end if
    end do
  end do

  !---------------------------!
  !   Create the array with   !
  !   information on faces    !
  !---------------------------!
  Grid % n_faces = 0                        ! initialize the number of sides
  call Adjust_First_Dim(4, Grid % faces_n)  ! adjust dimension of cells_c

  do run = 1, 2

    do c1 = 1, Grid % n_cells
      do m = 1, 24 ! through all the neighbouring cells
        c2 = Grid % cells_c(m, c1)
        if( (run .eq. 1) .and. (c2 > c1)  .or.  &
            (run .eq. 2) .and. (c2 < 0) ) then
          Grid % n_faces = Grid % n_faces + 1

          ! Which volumes are connected with side Grid % n_faces
          Grid % faces_c(1, Grid % n_faces) = c1
          Grid % faces_c(2, Grid % n_faces) = c2 

          ! Which is c2 neighbour of c1 and vice versa
          do c = 1,24
            if(Grid % cells_c(c,c1) .eq. c2) then 
              face_c_to_c(Grid % n_faces,1) = c
            end if

            if(c2 > 0) then
              if(Grid % cells_c(c,c2) .eq. c1) then 
                face_c_to_c(Grid % n_faces,2) = c
              end if 
            end if
          end do

          ! Nodes of a side Grid % n_faces
          if(c2  > 0) then
            if(ref % cell_level(c2) > ref % cell_level(c1)) then
              Grid % faces_n(1,Grid % n_faces) =  &
                Grid % cells_n( fn(face_c_to_c(Grid % n_faces,2),4), c2 )
              Grid % faces_n(2,Grid % n_faces) =  &
                Grid % cells_n( fn(face_c_to_c(Grid % n_faces,2),3), c2 )
              Grid % faces_n(3,Grid % n_faces) =  &
                Grid % cells_n( fn(face_c_to_c(Grid % n_faces,2),2), c2 )
              Grid % faces_n(4,Grid % n_faces) =  &
                Grid % cells_n( fn(face_c_to_c(Grid % n_faces,2),1), c2 )
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
          Grid % n_bnd_cells   = Grid % n_bnd_cells + 1

          ! Restore the boundary color, take positive value for color
          Grid % cells_c(m,c)  = -Grid % bnd_cond % color(-Grid % n_bnd_cells)
        end if
      end do
    end do
  end if

  !--------------------------------------!
  !   Is there enough allocated memory   !
  !--------------------------------------!
  if( Grid % n_faces  > Grid % max_n_faces ) then
    print *, '# ERROR in Generator'
    print *, '# The number sides is: ',              Grid % n_faces
    print *, '# There is space available only for:', Grid % max_n_faces
    print *, '# Increase the number of faces in .dom file'
    print *, '# and recompile the code. Good Luck !'
    stop
  end if

  end subroutine
