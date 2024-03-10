!==============================================================================!
  real function Volume_Average(Flow, val)
!------------------------------------------------------------------------------!
!   Calculates the volume average of the values in array val                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  real              :: val( - Flow % pnt_grid % n_bnd_cells  &
                            : Flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c
  real                     :: sum_val, tot_vol
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Initialize sums
  sum_val = 0.0
  tot_vol = 0.0

  ! Integrate them over cells, avoiding buffers
  do c = Cells_In_Domain()
    sum_val = sum_val + val(c) * Grid % vol(c)
    tot_vol = tot_vol + Grid % vol(c)
  end do

  ! Perform sums over all processors
  call Global % Sum_Real(sum_val)
  call Global % Sum_Real(tot_vol)

  ! Work out the final resul
  Volume_Average = sum_val / tot_vol

  end function
