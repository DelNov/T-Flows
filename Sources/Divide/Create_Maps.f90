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

  !-------------------------------!
  !                               !
  !   Browse through subdomains   !
  !                               !
  !-------------------------------!
  do sub = 1, maxval(grid % comm % cell_proc(:))

    !-----------------------!
    !   Save cell mapping   !
    !-----------------------!
    call File_Mod_Set_Name(name_map, processor=sub, extension='.map')
    call File_Mod_Open_File_For_Writing(name_map, fu)

    ! Extents are followed by mapping of the cells inside ...
    do c = 1, n_cells_sub
      if(grid % comm % cell_proc(c) .eq. sub) then
        write(fu, '(i9)') c
      end if
    end do

    ! ... followed by the cells on the boundary ...
    do c = -n_bnd_cells_sub, -1
      if(grid % comm % cell_proc(c) .eq. sub) then
        write(fu, '(i9)') c
      end if
    end do

    close(fu)

  end do   ! through subdomains

  end subroutine
