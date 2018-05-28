!==============================================================================!
  subroutine Probe_1d_Nodes_Gen(grid)
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of nodes in non-homogeneous           !
!   direction and write them in file name.1d                                   !
!------------------------------------------------------------------------------!
  use Name_Mod, only: problem_name
  use gen_mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------! 
  include "Approx.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n_prob, p, n
  real              :: n_p(10000)
  character(len=80) :: name_prob
  character(len=80) :: answer
  logical           :: isit
!==============================================================================!

  print *, '#==========================================='
  print *, '# Creating 1d file with the node '
  print *, '# coordinates in non-homogeneous directions '
  print *, '#-------------------------------------------'
  print *, '# Insert non-homogeneous direction '
  print *, '# (x, y, z, rx, ry, rz or skip)'
  print *, '# -------------------------------------------'
  read(*,*) answer
  call To_Upper_Case(answer)
  if(answer .eq. 'SKIP') return
 
  n_prob = 0
  n_p    = 0.0

  !-----------------------------!
  !   Browse through all nodes  !
  !-----------------------------!
  do n = 1, grid % n_nodes 

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer .eq. 'X') then
        if( Approx(grid % xn(n), n_p(p)) ) go to 1
      else if(answer .eq. 'Y') then
        if( Approx(grid % yn(n), n_p(p)) ) go to 1
      else if(answer .eq. 'Z') then
        if( Approx(grid % zn(n), n_p(p)) ) go to 1
      else if(answer .eq. 'RX') then
        if( Approx( (grid % zn(n)**2 +   &
                     grid % yn(n)**2)**.5, n_p(p)) ) go to 1
      else if(answer .eq. 'RY') then
        if( Approx( (grid % xn(n)**2 +   &
                     grid % zn(n)**2)**.5, n_p(p)) ) go to 1
      else if(answer .eq. 'RZ') then
        if( Approx( (grid % xn(n)**2 +   &
                     grid % yn(n)**2)**.5, n_p(p)) ) go to 1
      end if
    end do 
  
    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob+1
    if(answer .eq. 'X') n_p(n_prob)= grid % xn(n)
    if(answer .eq. 'Y') n_p(n_prob)= grid % yn(n)
    if(answer .eq. 'Z') n_p(n_prob)= grid % zn(n)

    if(answer .eq. 'RX') n_p(n_prob)= (grid % zn(n)**2 +      &
                                   grid % yn(n)**2)**0.5
    if(answer .eq. 'RY') n_p(n_prob)= (grid % xn(n)**2 +      &
                                   grid % zn(n)**2)**0.5
    if(answer .eq. 'RZ') n_p(n_prob)= (grid % xn(n)**2 +      &
                                   grid % yn(n)**2)**0.5

    if(n_prob .eq. 10000) then
      print *, '# Probe 1d: Not a 1d (channel flow) problem.'
      isit = .false.
      return
    end if
1  end do

  isit = .true.

  !--------------------!
  !   Create 1d file   !
  !--------------------!
  name_prob = problem_name
  name_prob(len_trim(problem_name)+1:len_trim(problem_name)+3) = '.1d'
  print *, 'Now creating the file: ', trim(name_prob)
  open(9, file=name_prob)
  ! Write the number of probes 
  write(9,'(I8)') n_prob

  call Sort2(n_p, n_prob*2, n_prob)

  ! Write the probe coordinates out
  do p=1, n_prob
    write(9,'(I8,1E17.8)') p, n_p(p)
  end do

  close(9)

  end subroutine
