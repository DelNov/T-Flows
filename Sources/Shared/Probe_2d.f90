!==============================================================================!
  subroutine Probe_2D(grid)
!------------------------------------------------------------------------------!
! Finds coordinates of all the planes for the channel flow.                    !
! It assumes that homogeneous directions of the flow are x and y.              !
!------------------------------------------------------------------------------!
  use Name_Mod, only: problem_name
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------! 
  include "Approx_Real.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n_prob, p, c
  real              :: yp(32768), zp(32768)
  character(len=80) :: name_prob
  character(len=80) :: answer
!==============================================================================!

  print *, '#==============================='
  print *, '# Looking for homogeneous plane '
  print *, '#-------------------------------'
  print *, '# Insert homogeneous direction  '
  print *, '# (xy, yz, zx or skip)'
  print *, '#-------------------------------'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer .eq. 'SKIP') return

  n_prob = 0
  zp(:)  = 0.0
  yp(:)  = 0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c=1,grid % n_cells

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer .eq. 'YZ') then
        if( Approx_Real(grid % yc(c), yp(p)) .and.      &
            Approx_Real(grid % zc(c), zp(p)) ) go to 1
      else if(answer .eq. 'ZX') then
        if( Approx_Real(grid % xc(c), yp(p)) .and.      &
            Approx_Real(grid % zc(c), zp(p)) ) go to 1
      else if(answer .eq. 'XY') then
        if( Approx_Real(grid % xc(c), yp(p)) .and.      &
            Approx_Real(grid % yc(c), zp(p)) ) go to 1
      end if
    end do 

    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob+1
    if(answer .eq. 'YZ') then
      yp(n_prob)=grid % yc(c)
      zp(n_prob)=grid % zc(c)
    else if(answer .eq. 'ZX') then
      yp(n_prob)=grid % xc(c)
      zp(n_prob)=grid % zc(c)
    else if(answer .eq. 'XY') then
      yp(n_prob)=grid % xc(c)
      zp(n_prob)=grid % yc(c)
    end if 

    if(n_prob .eq. 32768) then
      print *, '# Probe 2D: Not a 2D problem.'
      return
    end if
1 end do

  !--------------------!
  !   Create 2D file   !
  !--------------------!
  name_prob = problem_name
  name_prob(len_trim(problem_name)+1:len_trim(problem_name)+3) = '.2d'
  print *, '# Creating the file: ', trim(name_prob)
  open(9, file=name_prob)

  ! Write the number of probes 
  write(9,'(i8)') n_prob

  ! Write the probe coordinates out
  do p=1, n_prob
    write(9,'(i8,1pe17.8,1pe17.8)') p, yp(p), zp(p)
  end do

  close(9)

  end subroutine
