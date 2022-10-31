!==============================================================================!
  program Divisor
!------------------------------------------------------------------------------!
!   Divides the domain in equaly balanced subdomains.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    include 'Logo_Div.h90'
    include 'Save_Subdomains.h90'
  end interface
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: Grid           ! grid to be divided
  integer         :: n_sub
  real            :: p_start, p_end  ! variables to time the program
  character(SL)   :: arg
!==============================================================================!

  call cpu_time(p_start)  ! remark: cpu_time doesn't work with OPENMP

  call Logo_Div

  !-------------------------------------------!
  !   User specified command line arguments   !
  !-------------------------------------------!
  if(command_argument_count() .eq. 2) then

    call get_command_argument(1, arg)
    problem_name(1) = arg

    call get_command_argument(2, arg)
    read(arg, *) n_sub

  !-----------------------------------------!
  !    User didn't specify command line     !
  !   arguments, prompt him/her for input   !
  !-----------------------------------------!
  else

    print *, '# Input problem name: (without extension)'
    call File % Read_Line(5)
    read(Line % tokens(1), *)  problem_name(1)

    print *, '# Number of subdomains:'
    read(*,*)  n_sub

  end if

  !-------------------------------!
  !   Perform the decomposition   !
  !-------------------------------!

  ! Load the finite volume grid
  call Grid % Load_Cfn(0)
  call Grid % Load_Dim(0)

  ! Initialize processor numbers (poor idea to put it here)
  Grid % Comm % cell_proc(-Grid % n_bnd_cells:Grid % n_cells) = 1

  call Grid % Decompose(n_sub)

  call Save_Subdomains(Grid, 1)  ! Number of buffer levels is hard-coded now

  call cpu_time(p_end)  ! remark: cpu_time doesn't work with OPENMP
  print '(a10,f14.3,a9)', ' # Time = ', p_end-p_start, ' seconds.'

  end program
