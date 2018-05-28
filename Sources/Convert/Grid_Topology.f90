!==============================================================================!
  subroutine Grid_Topology(grid)
!------------------------------------------------------------------------------!
!   Determines the topology of the grid.                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use gen_mod 
  use Grid_Mod
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
      if(grid % cells_bnd_color(j,i) /= 0) then

        grid % n_bnd_cells = grid % n_bnd_cells + 1 

        ! grid % bnd_cond % color
        grid % bnd_cond % color(-grid % n_bnd_cells) =   &
                                                    grid % cells_bnd_color(j,i)

        ! Material:
        grid % material(-grid % n_bnd_cells) = grid % material(i)

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

  print '(a38,i7)', '# Number of boundary cells:          ', grid % n_bnd_cells

  end subroutine
