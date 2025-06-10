! Uncomment for more output.  Actually, some interesting things get printed out.
#define AMG_VERBOSE

!==============================================================================!
  module Amg_Mod
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!----------------------------------------------------------------------!
  save
!======================================================================!

# include "Amg_Mod.h90"

  integer, pointer, contiguous :: amg_imin (:), amg_imax (:)
  integer, pointer, contiguous :: amg_iminw(:), amg_imaxw(:)
  integer, pointer, contiguous :: amg_nstcol(:)

  !--------------!
  !   Amg type   !
  !--------------!
  type Amg_Type

    ! Range in unknown (u) for each level: imin(level) - imax(level))
    integer, private :: imin (AMG_MAX_LEVELS), imax (AMG_MAX_LEVELS)

    ! Range in workspace (iw) for each level: iminw(level) - imaxw(level))
    ! (This is gradualy populated in Amg % Coarsening, clearly enough)
    integer, private :: iminw(AMG_MAX_LEVELS), imaxw(AMG_MAX_LEVELS)
    integer, private :: nstcol(AMG_MAX_LEVELS)

    ! Residual, initial residual and convergence criterion
    real,    private :: resi(AMG_MAX_LEVELS)
    real,    private :: res, res0
    real,    private :: eps
    integer, private :: ncyc0

    ! Variables describing the V/W cycle
    integer, private :: def_coarse_solver    ! def solver at the coarsest level
    integer, private :: n_relax_coarse       ! n relaxations for coarsest level
    integer, private :: coarse_solver        ! holds coarsest level solver
    integer, private :: def_relax_down       ! relaxations going down the levels
    integer, private :: def_relax_up         ! relaxations going up the levels
    integer, private :: type_relax_down(10)  ! deciphered def_relax_down
    integer, private :: type_relax_up(10)    ! deciphered def_relax_up
    integer, private :: n_relax_down         ! values read from def_relax_down
    integer, private :: n_relax_up           ! values read from def_relax_up
    integer, private :: nrdlen
    integer, private :: nrulen
    integer, private :: fine_solver          ! holds solver for finer levels
    integer, private :: n_relax_fine         ! relaxation sweeps for coarsest l.

    ! Variables describing matrix properties (Class 1 by Ruge-Stueben)
    integer, private :: matrix
    integer, private :: irow0
    integer, private :: isym

    ! Variables controlling the general performance (Class 2 by Ruge-Stueben)
    integer, private :: iout

    ! Variables for tuning the coarsening algorithm (Class 4 by Ruge-Stueben)
    real    :: ecg1, ecg2, ewt2
    integer :: nwt, ntr

    ! These variables are used throughout various subroutines
    ! for checking if there is enough allocated memory
    integer, private :: mda, mdu, mdw

    ! Used only for timing
    real(SP), private :: time(20)
    real(SP), private :: told, tnew

    ! Error tracking
    integer :: ierr

    contains

      procedure :: timer_start
      procedure :: timer_stop
      procedure :: Initial_Residual
      procedure :: Final_Residual
      procedure :: Performed_Cycles
      procedure :: Enlarge_Int
      procedure :: Enlarge_Real
      procedure :: Reduce_Int
      procedure :: Reduce_Real

      !---------------------!
      !   Main subroutine   !
      !---------------------!
      procedure :: Amg1r5

      !----------------------!
      !   Related to setup   !
      !----------------------!
      procedure :: setup                         ! 1
      procedure ::   Check_Matrix_Properties     ! 1.1
      procedure ::   Coarsening                  ! 1.2
      procedure ::     Row_Sort                  ! 1.2.1
      procedure ::     Pre_Color                 ! 1.2.2
      procedure ::     Interpolation_Weights     ! 1.2.3
      procedure ::     Define_Operators          ! 1.2.4
      procedure ::     Truncate_Operator         ! 1.2.5
      procedure ::     Set_Inverse_Pointer       ! 1.2.6

      !----------------------------!
      !   Related to first guess   !
      !----------------------------!
      procedure :: First_Guess                   ! 2

      !------------------------!
      !   Related to solving   !
      !------------------------!
      procedure :: solve                         ! 3
      procedure ::   Calculate_Residual          ! 3.1
      procedure ::   Backup_U                    ! 3.2
      procedure ::   One_Cycle                   ! 3.3
      procedure ::     Solve_On_Coarsest_Level   ! 3.3.1
      procedure ::     Normalize_U               ! 3.3.2
      procedure ::     Gauss_Seidel_Sweep        ! 3.3.3.1
      procedure ::     Cg_On_Level               ! 3.3.3.2
      procedure ::     Bicg_On_Level             ! 3.3.3.3
      procedure ::     Set_U_To_Zero             ! 3.3.5
      procedure ::     Restrict_Residuals        ! 3.3.6
      procedure ::     Scale_Solution            ! 3.3.7
      procedure ::     Interpolate_Correction    ! 3.3.8
      procedure ::   Cg_Step                     ! 3.4
      procedure ::     Cg_Alpha                  ! 3.4.1
      procedure ::     Cg_Epsilon                ! 3.4.2

      !---------------------!
      !   Final reporting   !
      !---------------------!
      procedure :: Wrkcnt                        ! 4

      !---------------------------------------!
      !   Just a couple of little utilities   !
      !---------------------------------------!
      procedure :: Get_Integer_Digits
      procedure :: Random_0_To_0p1

  end type

  !-----------------------------------------------!
  !                                               !
  !   Singletone Ruge-Stueben AMG solver object   !
  !                                               !
  !-----------------------------------------------!
  type(Amg_Type) :: Amg

  contains

  !---------------------!
  !   Main subroutine   !
  !---------------------!
# include "Amg_Mod/0_Drivers/Amg1r5.f90"

  !----------------------!
  !   Related to setup   !
  !----------------------!
# include "Amg_Mod/1_Setup.f90"
# include "Amg_Mod/1_Setup/1_Check_Matrix_Properties.f90"
# include "Amg_Mod/1_Setup/2_Coarsening.f90"
# include "Amg_Mod/1_Setup/2_Coarsening/1_Row_Sort.f90"
# include "Amg_Mod/1_Setup/2_Coarsening/2_Pre_Color.f90"
# include "Amg_Mod/1_Setup/2_Coarsening/3_Interpolation_Weights.f90"
# include "Amg_Mod/1_Setup/2_Coarsening/4_Define_Operators.f90"
# include "Amg_Mod/1_Setup/2_Coarsening/5_Truncate_Operator.f90"
# include "Amg_Mod/1_Setup/2_Coarsening/6_Set_Inverse_Pointer.f90"

  !----------------------------!
  !   Related to first guess   !
  !----------------------------!
# include "Amg_Mod/2_First_Guess.f90"

  !------------------------!
  !   Related to solving   !
  !------------------------!
# include "Amg_Mod/3_Solve.f90"
# include "Amg_Mod/3_Solve/1_Calculate_Residual.f90"
# include "Amg_Mod/3_Solve/2_Backup_U.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/1_Solve_On_Coarsest_Level.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/2_Normalize_U.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/3.1_Gauss_Seidel_Sweep.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/3.2_Cg_On_Level.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/3.3_Bicg_On_Level.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/5_Set_U_To_Zero.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/6_Restrict_Residuals.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/7_Scale_Solution.f90"
# include "Amg_Mod/3_Solve/3_One_Cycle/8_Interpolate_Correction.f90"
# include "Amg_Mod/3_Solve/4_Cg_Step.f90"
# include "Amg_Mod/3_Solve/4_Cg_Step/1_Cg_Alpha.f90"
# include "Amg_Mod/3_Solve/4_Cg_Step/2_Cg_Epsilon.f90"

  !---------------------!
  !   Final reporting   !
  !---------------------!
# include "Amg_Mod/4_Wrkcnt.f90"

  !---------------------------------------!
  !   Just a couple of little utilities   !
  !---------------------------------------!
# include "Amg_Mod/8_Utilities/Enlarge_Int.f90"
# include "Amg_Mod/8_Utilities/Enlarge_Real.f90"
# include "Amg_Mod/8_Utilities/Reduce_Int.f90"
# include "Amg_Mod/8_Utilities/Reduce_Real.f90"
# include "Amg_Mod/8_Utilities/Final_Residual.f90"
# include "Amg_Mod/8_Utilities/Get_Integer_Digits.f90"
# include "Amg_Mod/8_Utilities/Initial_Residual.f90"
# include "Amg_Mod/8_Utilities/Performed_Cycles.f90"
# include "Amg_Mod/8_Utilities/Random_0_To_0p1.f90"
# include "Amg_Mod/8_Utilities/Timer_Start.f90"
# include "Amg_Mod/8_Utilities/Timer_Stop.f90"

  end module
