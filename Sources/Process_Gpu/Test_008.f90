#include "../Shared/Browse.h90"

!==============================================================================!
  subroutine Test_008
!------------------------------------------------------------------------------!
!>  Tests for volume balance in a rotating velocity field.
!------------------------------------------------------------------------------!
!   The volume balancing algorithms are in Insert_Volume_Source_For_Pressure   !
!   and in Correct_Velocity. albeit there also with OpenACC accelerators.      !
!------------------------------------------------------------------------------!
  use Field_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(Grid_Type)   :: Grid            ! computational grid
  type(Field_Type)  :: Flow            ! flow field
  integer           :: n, c, s, c1, c2, i_cel
  real              :: u_f, v_f, b_tmp
  real, allocatable :: b(:)
  character(len=11) :: root_control = 'control.008'
!==============================================================================!

  ! Start the parallel run and the profiler
  call Global % Start_Parallel
  call Profiler % Start('Test_008')

  O_Print '(a)', ' #======================================================'
  O_Print '(a)', ' # TEST 8: Volume balance in a rotating velocity field'
  O_Print '(a)', ' #======================================================'

  O_Print '(a)', ' # Opening the control file '//root_control
  call Control % Open_Root_File(root_control)

  O_Print '(a)', ' # Creating a grid'
  call Grid % Load_And_Prepare_For_Processing(1)

  O_Print '(a)', ' # Creating a field'
  call Flow % Create_Field(Grid)

  O_Print '(a)', ' # Initialize rotating belocity at cells'
  do c = -Grid % n_bnd_cells, Grid % n_cells - Grid % Comm % n_buff_cells
    Flow % u % n(c) =  Grid % yc(c)
    Flow % v % n(c) = -Grid % xc(c)
  end do
  call Grid % Exchange_Cells_Real(Flow % u % n)
  call Grid % Exchange_Cells_Real(Flow % v % n)
  call Grid % Exchange_Cells_Real(Flow % w % n)

  call Grid % Save_Debug_Vtu("uvw",                       &
                             vector_name="uvw",           &
                             vector_cell=(/Flow % u % n,  &
                                           Flow % v % n,  &
                                           Flow % w % n/))

  O_Print '(a)', ' # Initialize rotating belocity at faces'
  do s = 1, Grid % n_faces
    u_f =  Grid % yf(s)
    v_f = -Grid % xf(s)
    Flow % v_flux(s) = u_f * Grid % sx(s) + v_f * Grid % sy(s)
  end do

  allocate(b(Grid % n_cells))
  b(:) = 0.0

  do c1 = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    b_tmp = b(c1)
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        b_tmp = b_tmp - Flow % v_flux(s) * merge(1,-1, c1.lt.c2)
      end if
    end do
    b(c1) = b_tmp
  end do

  call Grid % Exchange_Inside_Cells_Real(b)
  call Grid % Save_Debug_Vtu("vol_src",              &
                             scalar_name="vol_src",  &
                             scalar_cell=b)

  ! End the profiler and the parallel run
  call Profiler % Stop('Test_008')
  call Profiler % Statistics(indent=24)
  call Global % End_Parallel

  end subroutine
