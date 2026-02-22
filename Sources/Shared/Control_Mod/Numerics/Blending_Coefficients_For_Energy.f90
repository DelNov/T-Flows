!==============================================================================!
  subroutine Blending_Coefficients_For_Energy(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads blending coefficients for energy (enthalpy) from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val(3)   !! blending coefficient values
  logical,   optional :: verbose  !! controls output verbosity
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: L = 16
!-----------------------------------[Locals]-----------------------------------!
  real         :: def(3)
  character(L) :: numb
  integer      :: p
!==============================================================================!

  data def / 0.0, 0.0, 1.0 /

  call Control % Read_Real_Vector('BLENDING_COEFFICIENTS_FOR_ENERGY',  &
                                  3, def, val, verbose)

  if( .not. Math % Approx_Real(sum(val(1:3)), 1.0) ) then
    write(numb, '(f16.10)')  sum(val(1:3))

    p = L
    do while(p .gt. 0 .and. numb(p:p) == '0')
      numb(p:p) = ' '
      p = p - 1
    end do

    call Message % Error(64,                                         &
                         'The sum of all advection coefficients '//  &
                         'must be equal to 1.0.  For energy '    //  &
                         'equation it is '//numb//'! '           //  &
                         'Correct the control file.  Exiting!')
  end if

  end subroutine
