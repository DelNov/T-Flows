#include "../Shared/Assert.h90"

!==============================================================================!
  program Divide_Prog
!------------------------------------------------------------------------------!
!   Divides the domain in equaly balanced subdomains.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Divide_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: Grid           ! grid to be divided
  integer         :: n_sub
  real            :: p_start, p_end  ! variables to time the program
  character(SL)   :: arg
!==============================================================================!

  call cpu_time(p_start)  ! remark: cpu_time doesn't work with OPENMP

  call Divide % Logo_Div()

  !-------------------------------------------!
  !   User specified command line arguments   !
  !-------------------------------------------!
  if(command_argument_count() .eq. 2) then

    call get_command_argument(1, arg)
    problem_name(1) = arg

    call get_command_argument(2, arg)
    read(arg, *) n_sub

    Assert(n_sub > 1)

  !-----------------------------------------!
  !    User didn't specify command line     !
  !   arguments, prompt him/her for input   !
  !-----------------------------------------!
  else

    print *, '# Input problem name: (without extension)'
    problem_name(1) = File % Single_Word_From_Keyboard()

    print *, '# Number of subdomains:'
    n_sub = File % Single_Int_From_Keyboard()

  end if

  !-------------------------------!
  !   Perform the decomposition   !
  !-------------------------------!

  ! Load the finite volume grid
  call Grid % Load_Cfn((/0, 0/))
  call Grid % Load_Dim((/0, 0/))

  ! Initialize processor numbers (poor idea to put it here)
  Grid % Comm % cell_proc(-Grid % n_bnd_cells:Grid % n_cells) = 1

  call Grid % Decompose(n_sub)

  call Divide % Save_Subdomains(Grid, n_sub)

  call cpu_time(p_end)  ! remark: cpu_time doesn't work with OPENMP
  print '(a10,f14.3,a9)', ' # Time = ', p_end-p_start, ' seconds.'

  end program
