!==============================================================================!
  program Main
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  character(80) :: arg
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

    ! An attempt to solve Navier-Stokes equations
    else if(arg .eq. '6') then
      call Test_006()
      goto 1

    ! An attempt to solve Navier-Stokes equations
    else if(arg .eq. '7') then
      call Test_007()
      goto 1

    end if

  end if

  !----------------------------------------------------------------------!
  !   If you are here, something was wrong with command line arguments   !
  !----------------------------------------------------------------------!
  print '(a)', ' Failed to invoke the program correctly.'
  print '(a)', ' Correct invocation is:'
  print '(a)', ''
  print '(a)', ' ./Program <test>'
  print '(a)', ''
  print '(a)', ' where <test> can be from 1 to 5 depending if you want to test:'
  print '(a)', '   1 - sparse-matrix vector product'
  print '(a)', '   2 - vector vector dot product'
  print '(a)', '   3 - operations: c = a + scalar * b and'
  print '(a)', '                   c = a - scalar * b'
  print '(a)', '   4 - conjugate gradient solver'
  print '(a)', '   5 - field creation and gradient calculation'
  print '(a)', '   6 - solution of the Navier-Stokes equations'
  print '(a)', '   7   matrix of cubes case'

1 continue

  end program
