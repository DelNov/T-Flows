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
!==============================================================================!

  call cpu_time(start)

  call Logo_Div

  print *, '# Input problem name: (without extension)'
  call File % Read_Line(5)
  read(line % tokens(1), *)  problem_name(1)

  ! Load the finite volume grid
  call Grid % Load_Cfn(0)
  call Grid % Load_Dim(0)

  ! Initialize processor numbers (poor idea to put it here)
  Grid % Comm % cell_proc(-Grid % n_bnd_cells:Grid % n_cells) = 1

  print *, '# Number of subdomains:'
  read(*,*)  n_sub

  call Grid % Decompose(n_sub)

  call Save_Subdomains(Grid, 1)  ! Number of buffer levels is hard-coded now

  call cpu_time(finish)
  print '(a10,f14.3,a9)', ' # Time = ', finish-start, ' seconds.'

  end program
