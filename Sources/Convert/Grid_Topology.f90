!==============================================================================!
  subroutine Grid_Topology(grid)
!------------------------------------------------------------------------------!
!   Determines the topology of the grid.                                       !
!
!   To be more specific, it determines:
!                                                                              !
!   grid % n_bnd_cells   - number of boundary cells                            !
!   grid % cells_n_nodes - number of nodes for each cell                       !
!   grid % cells_n       - list of each cell's nodes                           !
!   grid % n_faces       - number of faces on the boundary                     !
!   grid % faces_n_nodes - number of nodes for each face on the boundary       !
!   grid % faces_n       - list of each boundary face's nodes                  !
!   grid % faces_c       - a pair of cells surrounding each boundary face      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod, only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include "Cell_Numbering_Neu.f90"
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
!==============================================================================!

  !------------------------------!
  !   Count the boundary cells   !
  !------------------------------!
  grid % n_bnd_cells = 0
  grid % n_faces  = 0
  do i = 1, grid % n_cells
    do j = 1, 6
      if(grid % cells_bnd_color(j,i) .ne. 0) then

        grid % n_bnd_cells = grid % n_bnd_cells + 1

        ! grid % bnd_cond % color
        grid % bnd_cond % color(-grid % n_bnd_cells) =   &
                                                    grid % cells_bnd_color(j,i)
        ! Faces
        grid % n_faces = grid % n_faces  + 1
        grid % faces_c(1,grid % n_faces) = i
        grid % faces_c(2,grid % n_faces) = -grid % n_bnd_cells

        ! Hexahedra:
        if(grid % cells_n_nodes(i) .eq. 8) then
          grid % cells_n_nodes(-grid % n_bnd_cells) = 4
          grid % cells_n(1,-grid % n_bnd_cells) = grid % cells_n(neu_hex(j,1),i)
          grid % cells_n(2,-grid % n_bnd_cells) = grid % cells_n(neu_hex(j,2),i)
          grid % cells_n(3,-grid % n_bnd_cells) = grid % cells_n(neu_hex(j,3),i)
          grid % cells_n(4,-grid % n_bnd_cells) = grid % cells_n(neu_hex(j,4),i)

          grid % faces_n_nodes(grid % n_faces) = 4
          grid % faces_n(1:4, grid % n_faces)  =   &
          grid % cells_n(1:4,-grid % n_bnd_cells)

        ! Prisms:
        else if(grid % cells_n_nodes(i) .eq. 6) then
          if(j <= 3) then    ! faces (1), (2) and (3)
            grid % cells_n_nodes(-grid % n_bnd_cells) = 4
            grid % cells_n(1,-grid % n_bnd_cells)=grid % cells_n(neu_wed(j,1),i)
            grid % cells_n(2,-grid % n_bnd_cells)=grid % cells_n(neu_wed(j,2),i)
            grid % cells_n(3,-grid % n_bnd_cells)=grid % cells_n(neu_wed(j,3),i)
            grid % cells_n(4,-grid % n_bnd_cells)=grid % cells_n(neu_wed(j,4),i)

            grid % faces_n_nodes(grid % n_faces) = 4
            grid % faces_n(1:4, grid % n_faces)  =   &
            grid % cells_n(1:4,-grid % n_bnd_cells)

          else if(j <= 5) then
            grid % cells_n_nodes(-grid % n_bnd_cells) = 3
            grid % cells_n(1,-grid % n_bnd_cells)=grid % cells_n(neu_wed(j,1),i)
            grid % cells_n(2,-grid % n_bnd_cells)=grid % cells_n(neu_wed(j,2),i)
            grid % cells_n(3,-grid % n_bnd_cells)=grid % cells_n(neu_wed(j,3),i)

            grid % faces_n_nodes(grid % n_faces) = 3
            grid % faces_n(1:3, grid % n_faces)  =   &
            grid % cells_n(1:3,-grid % n_bnd_cells)
          end if

        ! Tetrahedra:
        else if(grid % cells_n_nodes(i) .eq. 4) then
          if(j <= 4) then
            grid % cells_n_nodes(-grid % n_bnd_cells) = 3
            grid % cells_n(1,-grid % n_bnd_cells)=grid % cells_n(neu_tet(j,1),i)
            grid % cells_n(2,-grid % n_bnd_cells)=grid % cells_n(neu_tet(j,2),i)
            grid % cells_n(3,-grid % n_bnd_cells)=grid % cells_n(neu_tet(j,3),i)

            grid % faces_n_nodes(grid % n_faces) = 3
            grid % faces_n(1:3, grid % n_faces)  =   &
            grid % cells_n(1:3,-grid % n_bnd_cells)
          end if

        ! Pyramides:
        else if(grid % cells_n_nodes(i) .eq. 5) then
          if(j .eq. 1) then    ! face (1)
            grid % cells_n_nodes(-grid % n_bnd_cells) = 4
            grid % cells_n(1,-grid % n_bnd_cells)=grid % cells_n(neu_pyr(j,1),i)
            grid % cells_n(2,-grid % n_bnd_cells)=grid % cells_n(neu_pyr(j,2),i)
            grid % cells_n(3,-grid % n_bnd_cells)=grid % cells_n(neu_pyr(j,3),i)
            grid % cells_n(4,-grid % n_bnd_cells)=grid % cells_n(neu_pyr(j,4),i)

            grid % faces_n_nodes(grid % n_faces) = 4
            grid % faces_n(1:4, grid % n_faces)  =   &
            grid % cells_n(1:4,-grid % n_bnd_cells)

          else if(j <= 5) then
            grid % cells_n_nodes(-grid % n_bnd_cells) = 3
            grid % cells_n(1,-grid % n_bnd_cells)=grid % cells_n(neu_pyr(j,1),i)
            grid % cells_n(2,-grid % n_bnd_cells)=grid % cells_n(neu_pyr(j,2),i)
            grid % cells_n(3,-grid % n_bnd_cells)=grid % cells_n(neu_pyr(j,3),i)

            grid % faces_n_nodes(grid % n_faces) = 3
            grid % faces_n(1:3, grid % n_faces)  =   &
            grid % cells_n(1:3,-grid % n_bnd_cells)
          end if

        else
          print *, '# Cell with invalid number of nodes: ',  &
                   grid % cells_n_nodes(i)
          print *, '# Exiting!'
          stop
        end if

      end if
    end do
  end do

  print '(a38,i9)', '# Number of boundary cells:          ', grid % n_bnd_cells
  print '(a38,i9)', '# Number of faces on the boundary:   ', grid % n_faces

  end subroutine
