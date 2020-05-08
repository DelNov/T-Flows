!==============================================================================!
  subroutine Create_Maps(grid)
!------------------------------------------------------------------------------!
!   Number the cells in each subdomain for subsequent separate saving.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
  use Grid_Mod, only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2, sub, n_cells_sub, n_bnd_cells_sub
!==============================================================================!

  !-------------------------------!
  !                               !
  !   Browse through subdomains   !
  !                               !
  !-------------------------------!
  do sub = 1, maxval(grid % comm % cell_proc(:))

    ! Cells
    n_cells_sub = 0     ! number of cells in subdomain
    do c = 1, grid % n_cells
      if(grid % comm % cell_proc(c) .eq. sub) then
        n_cells_sub = n_cells_sub + 1     ! increase the number of cells in sub.
        grid % comm % cell_glo(n_cells_sub) = c  ! map to global cell number
      end if
    end do

    ! Faces & real boundary cells
    n_bnd_cells_sub = 0  ! number of real boundary cells in subdomain

    ! Faces step 2/3: on the boundaries + bundary cells
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if( grid % comm % cell_proc(c1) .eq. sub )  then
          n_bnd_cells_sub = n_bnd_cells_sub + 1  ! increase n. of bnd. cells
          grid % comm % cell_glo(-n_bnd_cells_sub) = c2  ! map to g. cell number
        end if
      end if 
    end do

  end do   ! through subdomains

  end subroutine
