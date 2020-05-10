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
  integer           :: c, s, c1, c2, sub, n_cells_sub, n_bnd_cells_sub, fu
  character(len=80) :: name_map
!==============================================================================!

  allocate(grid % comm % cell_glo(-grid % n_bnd_cells : grid % n_cells));
  grid % comm % cell_glo(:) = 0

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

    ! Real boundary cells
    n_bnd_cells_sub = 0  ! number of real boundary cells in subdomain
    do c = -grid % n_bnd_cells, -1
      if(grid % comm % cell_proc(c) .eq. sub) then
        n_bnd_cells_sub = n_bnd_cells_sub + 1  ! increase n. of bnd. cells
        grid % comm % cell_glo(-n_bnd_cells_sub) = c  ! map to global cell number
      end if
    end do

    !-----------------------!
    !   Save cell mapping   !
    !-----------------------!
    call File_Mod_Set_Name(name_map, processor=sub, extension='.map')
    call File_Mod_Open_File_For_Writing(name_map, fu)

    ! Extents are followed by mapping of the cells inside ...
    do c = 1, n_cells_sub
      write(fu, '(i9)') grid % comm % cell_glo(c)
    end do

    ! ... followed by the cells on the boundary ...
    do c = -n_bnd_cells_sub, -1
      write(fu, '(i9)') grid % comm % cell_glo(c)
    end do

    close(fu)

  end do   ! through subdomains

  end subroutine
