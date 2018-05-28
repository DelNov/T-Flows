!==============================================================================!
  subroutine Save_Vtu_Vector(grid, in_1, in_2, var_name, val_1, val_2, val_3)
!------------------------------------------------------------------------------!
!   Writes one real vector defined over cells.                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc, n_proc
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: in_1, in_2
  character(len=*) :: var_name
  real             :: val_1(1:grid % n_cells)
  real             :: val_2(1:grid % n_cells)
  real             :: val_3(1:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Header
  if(n_proc > 1 .and. this_proc == 1) then
    write(8,'(4a)') in_1,                                       & 
                    '<DataArray type="Float32" Name="',         &
                    trim(var_name),                             &
                    '" NumberOfComponents="3" format="ascii"/>'
  end if
  write(9,'(4a)') in_1,                                       & 
                  '<DataArray type="Float32" Name="',         &
                  trim(var_name),                             &
                  '" NumberOfComponents="3" format="ascii">'
  ! Data
  do c = 1, grid % n_cells
    write(9,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)') &
            in_2, val_1(c), val_2(c), val_3(c)
  end do  

  ! Footer
  write(9,'(a,a)') in_1, '</DataArray>'

  end subroutine
