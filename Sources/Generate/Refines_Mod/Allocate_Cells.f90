!==============================================================================!
  subroutine Refines_Mod_Allocate_Cells(ref, nc, nb, margin)
!------------------------------------------------------------------------------!
!>  Allocates arrays for cell marking and cell level in the Refines_Type
!>  structure (ref), based on the given number of boundary cells (nb) and
!>  cells (nc).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type)            :: ref     !! refinement regions
  integer, intent(in)           :: nc      !! number of cells
  integer, intent(in)           :: nb      !! number of boundary cells
  integer, intent(in), optional :: margin  !! margin for allocation
!-----------------------------------[Locals]-----------------------------------!
  integer :: nc_m, nb_m
!==============================================================================!

  ! If cell-based arrays are allocated and they
  ! are bigger than requested, get out of here
  if(allocated(ref % cell_marked)) then
    if(-nb .ge. lbound(ref % cell_marked, 1) .and.  &
        nc .le. ubound(ref % cell_marked, 1)) then
      return
    end if
  end if

  ! Process the margin if specified
  nc_m = nc
  nb_m = nb
  if(present(margin)) then
    Assert(margin .ge. 0)
    nc_m = nc + margin
    nb_m = nb + margin
  end if

  print '(a,2i9)', ' # Expanding memory for cells in Refines_Mod to size: ',  &
                   nc, nb

  call Enlarge % Array_Log(ref % cell_marked, i=(/-nb_m,nc_m/))
  call Enlarge % Array_Int(ref % cell_level,  i=(/    1,nc_m/))

  end subroutine
