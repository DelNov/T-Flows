!==============================================================================!
  subroutine Save_Vtu_Scalar(grid, in_1, in_2, var_name, val)
!------------------------------------------------------------------------------!
!   Writes one real scalar defined over cells.                                 !
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
  real             :: val(1:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Header
  if(n_proc > 1 .and. this_proc .eq. 1) then
    write(8,'(4a)') in_1,                                 & 
                    '<PDataArray type="Float64" Name="',  &
                    trim(var_name),                       &
                    '" format="ascii"/>'
  end if
  write(9,'(4a)') in_1,                                & 
                  '<DataArray type="Float64" Name="',  &
                  trim(var_name),                      &
                  '" format="ascii">'
 
  ! Data
  do c = 1, grid % n_cells
    write(9,'(a,1pe16.6e4)') in_2, val(c)
  end do  

  ! Footer
  write(9,'(a,a)') in_1, '</DataArray>'

  end subroutine
