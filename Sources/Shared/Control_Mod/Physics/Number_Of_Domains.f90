!==============================================================================!
  subroutine Number_Of_Domains(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads number of domains involved in simulation from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! number of domains in simulation
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_DOMAINS', 1, val, verbose)

  end subroutine
