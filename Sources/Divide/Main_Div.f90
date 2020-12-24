!==============================================================================!
  program Divisor
!------------------------------------------------------------------------------!
!   Divides the domain in equaly balanced subdomains.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
  use Grid_Mod, only: Grid_Type,           &
                      Grid_Mod_Decompose,  &
                      Grid_Mod_Load_Cfn,   &
                      Grid_Mod_Load_Dim
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: grid           ! grid to be divided
  integer         :: n_sub
  real            :: start, finish  ! variables to time the program
!==============================================================================!

  call cpu_time(start)

  call Logo_Div

  print *, '# Input problem name: (without extension)'
  call File_Mod_Read_Line(5)
  read(line % tokens(1), *)  problem_name(1)

  ! Load the finite volume grid
  call Grid_Mod_Load_Cfn(grid, 0)
  call Grid_Mod_Load_Dim(grid, 0)

  ! Initialize processor numbers (poor idea to put it here)
  grid % comm % cell_proc(-grid % n_bnd_cells:grid % n_cells) = 1

  print *, '# Number of subdomains:'
  read(*,*)  n_sub

  call Grid_Mod_Decompose(grid, n_sub)

  call Save_Subdomains(grid, 1)  ! Number of buffer levels is hard-coded now

  call cpu_time(finish)
  print '(a10,f14.3,a9)', ' # Time = ', finish-start, ' seconds.'

  end program
