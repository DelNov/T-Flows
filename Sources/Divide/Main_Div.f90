!==============================================================================!
  program Divisor
!------------------------------------------------------------------------------!
!   Divides the domain in equaly balanced subdomains.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Gen_Mod 
  use Div_Mod
  use Tokenizer_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type) :: grid           ! grid to be divided
  integer         :: c, n_sub
  real            :: start, finish  ! variables to time the program
!==============================================================================!

  call cpu_time(start)

  call Logo_Div

  print *, '# Input problem name: (without extension)'
  call Tokenizer_Mod_Read_Line(5)  
  read(line % tokens(1), *)  problem_name

  ! Load the finite volume grid
  call Load_Cns           (grid, 0)
  call Allocate_Additional(grid)
  call Load_Geo           (grid, 0)

  ! Initialize processor numbers (poor idea to put it here)
  do c = 1, grid % n_cells
    grid % comm % proces(c) = 1
  end do

  call Load_Geo(grid, 0)

  print *, '# Number of subdomains:'
  read(*,*)  n_sub

  call Grid_Mod_Decompose(grid, n_sub)

  call Create_Maps(grid)
  call Create_Buffers_And_Save(grid)

  call cpu_time(finish)
  print '(a10,f14.3,a9)', ' # Time = ', finish-start, ' seconds.'

  end program
