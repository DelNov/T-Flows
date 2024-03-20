#include "../Shared/Browse.h90"

!==============================================================================!
  program Main
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  character(SL) :: arg
!==============================================================================!

  !-----------------------------------------!
  !   Select a test to perform and run it   !
  !-----------------------------------------!
  if(command_argument_count() .eq. 1) then

    ! Fetch the command line argument
    call get_command_argument(1, arg)

    ! Sparse-matrix with vector product
    if(arg .eq. '1') then
      call Test_001()
      goto 1

    ! Vector vector dot product
    else if(arg .eq. '2') then
      call Test_002()
      goto 1

    ! Operations  c = a + s * b  and  c = a - s * b
    else if(arg .eq. '3') then
      call Test_003()
      goto 1

    ! Conjugate Gradient from the native module
    else if(arg .eq. '4') then
      call Test_004()
      goto 1

    ! Field creation and gradient calculation
    else if(arg .eq. '5') then
      call Test_005()
      goto 1

    ! Volume balance in a rotating velocity field
    else if(arg .eq. '6') then
      call Test_006()
      goto 1

    ! Navier-Stokes solutions in the lid driven cavity
    else if(arg .eq. '7') then
      call Test_007()
      goto 1

    ! Navier-Stokes solutions over the matrix of cubes test case
    else if(arg .eq. '8') then
      call Test_008()
      goto 1

    ! Navier-Stokes solutions over the matrix for the backstep flow
    else if(arg .eq. '9') then
      call Test_009()
      goto 1

    end if

  end if

  !----------------------------------------------------------------------!
  !   If you are here, something was wrong with command line arguments   !
  !----------------------------------------------------------------------!
  call Global % Start_Parallel

  O_Print '(a)', ' Failed to invoke the program correctly.'
  O_Print '(a)', ' Correct invocation is:'
  O_Print '(a)', ''
  O_Print '(a)', ' ./Program <test>'
  O_Print '(a)', ''
  O_Print '(a)', ' where <test> can be from 1 to 9 depending if you want to test:'
  O_Print '(a)', '   1 - sparse-matrix vector product'
  O_Print '(a)', '   2 - vector vector dot product'
  O_Print '(a)', '   3 - operations: c = a + scalar * b and'
  O_Print '(a)', '                   c = a - scalar * b'
  O_Print '(a)', '   4 - conjugate gradient solver'
  O_Print '(a)', '   5 - field creation and gradient calculation'
  O_Print '(a)', '   6 - volume balance in a rotating velocity field.'
  O_Print '(a)', '   7 - solution of the Navier-Stokes equations'
  O_Print '(a)', '   8 - matrix of cubes case'
  O_Print '(a)', '   9 - flow over a backward facing step'

  call Global % End_Parallel

1 continue

  end program
