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
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: Grid           ! grid to be divided
  integer         :: n_sub
  real            :: start, finish  ! variables to time the program
  character(SL)   :: arg
!==============================================================================!

  call cpu_time(start)

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
    read(line % tokens(1), *)  problem_name(1)

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

  call cpu_time(finish)
  print '(a10,f14.3,a9)', ' # Time = ', finish-start, ' seconds.'

  end program
