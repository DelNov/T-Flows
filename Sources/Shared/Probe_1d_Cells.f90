!==============================================================================!
  subroutine Probe_1d_Cells(grid)
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of cell-centers in non-homogeneous    !
!   direction and write them in file called "name.1d"                          !
!------------------------------------------------------------------------------!
  use Name_Mod,  only: problem_name
  use Const_Mod, only: NANO          ! 1.0e-9
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical :: isit
!----------------------------------[Calling]-----------------------------------! 
  include "Approx_Real.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n_prob, p, c
  real              :: zp(16384)
  character(len=80) :: name_prob
  character(len=80) :: answer
!==============================================================================!

  print *, '#================================================='
  print *, '# Looking for non-homogeneous directions '
  print *, '# Insert non-homogeneous direction (x,y,z or skip)'
  print *, '#-------------------------------------------------'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer .eq. 'SKIP') return

  n_prob = 0
  zp(:)  = 0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c = 1, grid % n_cells

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer .eq. 'X') then
        if( Approx_Real(grid % xc(c), zp(p), NANO)) go to 1
      else if(answer .eq. 'Y') then
        if( Approx_Real(grid % yc(c), zp(p), NANO)) go to 1
      else if(answer .eq. 'Z') then
        if( Approx_Real(grid % zc(c), zp(p), NANO)) go to 1
      end if
    end do

    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob + 1
    if(answer .eq. 'X') zp(n_prob) = grid % xc(c)
    if(answer .eq. 'Y') zp(n_prob) = grid % yc(c)
    if(answer .eq. 'Z') zp(n_prob) = grid % zc(c)

    if(n_prob .eq. 16384) then
      print *, '# Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      return
    end if
1 end do

  isit = .true.

  !--------------------!
  !   Create 1D file   !
  !--------------------!
  name_prob = problem_name
  name_prob(len_trim(problem_name)+1:len_trim(problem_name)+3) = '.1d'
  print *, '# Creating the file: ', trim(name_prob)
  open(9, file=name_prob)

  ! Write the number of probes 
  write(9,'(i8)') n_prob

  ! Write the probe coordinates out
  do p=1,n_prob
    write(9,'(i8,1pe17.8)') p, zp(p)
  end do

  close(9)

  end subroutine
