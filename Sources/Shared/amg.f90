! Uncomment for more output.  Actually, some interesting things get printed out.
! #define VERBOSE

!======================================================================!
  module module_amg
!----------------------------------------------------------------------!
  save
!======================================================================!

# include "amg.h90"

  integer, pointer, contiguous :: amg_imin (:), amg_imax (:)
  integer, pointer, contiguous :: amg_iminw(:), amg_imaxw(:)
  integer, pointer, contiguous :: amg_nstcol(:)

  type amg_type

    ! Range in unknown (u) for each level: imin(level) - imax(level))
    integer, private :: imin (AMG_MAX_LEVELS), imax (AMG_MAX_LEVELS)

    ! Range in workspace (iw) for each level: iminw(level) - imaxw(level))
    ! (This is gradualy populated in amg % coarsening, clearly enough)
    integer, private :: iminw(AMG_MAX_LEVELS), imaxw(AMG_MAX_LEVELS)
    integer, private :: nstcol(AMG_MAX_LEVELS)

    ! Residual, initial residual and convergence criterion
    double precision, private :: resi(AMG_MAX_LEVELS)
    double precision, private :: res, res0
    double precision, private :: eps

    ! Variables describing the V/W cycle
    integer, private :: nsolco        ! describing solver at the coarsest level
    integer, private :: nrcx          ! holds parts of nsolco
    integer, private :: nsc           ! holds coarsest level solver
    integer, private :: nrd           ! relaxations going down the levels
    integer, private :: nru           ! relaxations going up the levels
    integer, private :: nrdtyp(10)    ! deciphered nrd
    integer, private :: nrutyp(10)    ! deciphered nru
    integer, private :: nrdx, nrdlen  ! hold values read from nrd
    integer, private :: nrux, nrulen  ! hold values read from nru

    ! Variables describing matrix properties (Class 1 by Ruge-Stueben)
    integer, private :: matrix
    integer, private :: irow0
    integer, private :: isym

    ! Variables controlling the general performance (Class 2 by Ruge-Stueben)
    integer, private :: iout

    ! Variables for tuning the coarsening algorithm (Class 4 by Ruge-Stueben)
    double precision :: ecg1, ecg2, ewt2
    integer          :: nwt, ntr

    ! These variables are used throughout various subroutines
    ! for checking if there is enough allocated memory
    integer, private :: mda, mdu, mdw

    ! Used only for timing
    real, private :: time(20)
    real, private :: told, tnew

    ! Error tracking
    integer :: ierr

    contains

      procedure :: timer_start
      procedure :: timer_stop

      !---------------------!
      !   Main subroutine   !
      !---------------------!
      procedure :: amg1r5

      !----------------------!
      !   Related to setup   !
      !----------------------!
      procedure :: setup                         ! 1
      procedure ::   check_matrix_properties     ! 1.1
      procedure ::   coarsening                  ! 1.2
      procedure ::     row_sort                  ! 1.2.1
      procedure ::     pre_color                 ! 1.2.2
      procedure ::     interpolation_weights     ! 1.2.3
      procedure ::     define_operators          ! 1.2.4
      procedure ::     truncate_operator         ! 1.2.5
      procedure ::     set_inverse_pointer       ! 1.2.6

      !----------------------------!
      !   Related to first guess   !
      !----------------------------!
      procedure :: first_guess                   ! 2

      !------------------------!
      !   Related to solving   !
      !------------------------!
      procedure :: solve                         ! 3
      procedure ::   calculate_residual          ! 3.1
      procedure ::   backup_u                    ! 3.2
      procedure ::   one_cycle                   ! 3.3
      procedure ::     solve_on_coarsest_level   ! 3.3.1
      procedure ::       cg_on_coarsest_level    ! 3.3.1.1
      procedure ::       bicg_on_coarsest_level  ! 3.3.1.2
      procedure ::     normalize_u               ! 3.3.2
      procedure ::     gauss_seidel_sweep        ! 3.3.3
      procedure ::     set_u_to_zero             ! 3.3.5
      procedure ::     restrict_residuals        ! 3.3.6
      procedure ::     scale_solution            ! 3.3.7
      procedure ::     interpolate_correction    ! 3.3.8
      procedure ::   cg_step                     ! 3.4
      procedure ::     cg_alpha                  ! 3.4.1
      procedure ::     cg_epsilon                ! 3.4.2

      !---------------------!
      !   Final reporting   !
      !---------------------!
      procedure :: wrkcnt                        ! 4

      !---------------------------------------!
      !   Just a couple of little utilities   !
      !---------------------------------------!
      procedure :: get_integer_digits
      procedure :: random_0_to_0p1

  end type

  !-----------------------------------------------!
  !                                               !
  !   Singletone Ruge-Stueben AMG solver object   !
  !                                               !
  !-----------------------------------------------!
  type(amg_type) :: amg

  contains

!==============================================================================!
  subroutine timer_start(amg)
!------------------------------------------------------------------------------!
  implicit none
  class(amg_type) :: amg
!------------------------------------------------------------------------------!
  call cpu_time(amg % told)
  end subroutine

!==============================================================================!
  subroutine timer_stop(amg, i)
!------------------------------------------------------------------------------!
  implicit none
  class(amg_type) :: amg
  integer         :: i
!------------------------------------------------------------------------------!
  call cpu_time(amg % tnew)
  amg % time(i) = amg % tnew - amg % told
  end subroutine

  !---------------------!
  !   Main subroutine   !
  !---------------------!
# include "amg/0_drivers/amg1r5.f90"

  !----------------------!
  !   Related to setup   !
  !----------------------!
# include "amg/1_setup.f90"
# include "amg/1_setup/1_check_matrix_properties.f90"
# include "amg/1_setup/2_coarsening.f90"
# include "amg/1_setup/2_coarsening/1_row_sort.f90"
# include "amg/1_setup/2_coarsening/2_pre_color.f90"
# include "amg/1_setup/2_coarsening/3_interpolation_weights.f90"
# include "amg/1_setup/2_coarsening/4_define_operators.f90"
# include "amg/1_setup/2_coarsening/5_truncate_operator.f90"
# include "amg/1_setup/2_coarsening/6_set_inverse_pointer.f90"

  !----------------------------!
  !   Related to first guess   !
  !----------------------------!
# include "amg/2_first_guess.f90"

  !------------------------!
  !   Related to solving   !
  !------------------------!
# include "amg/3_solve.f90"
# include "amg/3_solve/1_calculate_residual.f90"
# include "amg/3_solve/2_backup_u.f90"
# include "amg/3_solve/3_one_cycle.f90"
# include "amg/3_solve/3_one_cycle/1_solve_on_coarsest_level.f90"
# include "amg/3_solve/3_one_cycle/1_solve_on_coarsest_level/1_cg_on_coarsest_level.f90"
# include "amg/3_solve/3_one_cycle/1_solve_on_coarsest_level/2_bicg_on_coarsest_level.f90"
# include "amg/3_solve/3_one_cycle/2_normalize_u.f90"
# include "amg/3_solve/3_one_cycle/3_gauss_seidel_sweep.f90"
# include "amg/3_solve/3_one_cycle/5_set_u_to_zero.f90"
# include "amg/3_solve/3_one_cycle/6_restrict_residuals.f90"
# include "amg/3_solve/3_one_cycle/7_scale_solution.f90"
# include "amg/3_solve/3_one_cycle/8_interpolate_correction.f90"
# include "amg/3_solve/4_cg_step.f90"
# include "amg/3_solve/4_cg_step/1_cg_alpha.f90"
# include "amg/3_solve/4_cg_step/2_cg_epsilon.f90"

  !---------------------!
  !   Final reporting   !
  !---------------------!
# include "amg/4_wrkcnt.f90"

  !---------------------------------------!
  !   Just a couple of little utilities   !
  !---------------------------------------!
# include "amg/8_utilities/1_get_integer_digits.f90"
# include "amg/8_utilities/2_random_0_to_0p1.f90"

  end module
