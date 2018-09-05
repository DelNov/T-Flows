!==============================================================================!
  subroutine Probe_1d_Cells_Nodes(grid)
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of cell's nodes in non-homogeneous    !
!   direction and write them in file called "name.1d"                          !
!------------------------------------------------------------------------------!
  use Name_Mod, only: problem_name
  use Grid_Mod
  use Sort_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------! 
  include "Approx_Real.int"
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n_prob, p, c, n
  real              :: zp(16384)
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
  zp(:)  = 0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c = -grid % n_bnd_cells, grid % n_cells
    do n = 1, grid % cells_n_nodes(c)

      ! Try to find the cell among the probes
      do p=1, n_prob
        if(answer .eq. 'X') then
          if( Approx_Real(grid % xn(grid % cells_n(n,c)), zp(p)) ) go to 1
        else if(answer .eq. 'Y') then
          if( Approx_Real(grid % yn(grid % cells_n(n,c)), zp(p)) ) go to 1
        else if(answer .eq. 'Z') then
          if( Approx_Real(grid % zn(grid % cells_n(n,c)), zp(p)) ) go to 1
        else if(answer .eq. 'RX') then
          if( Approx_Real( sqrt(grid % zn(grid % cells_n(n,c))**2 +           &
                                grid % yn(grid % cells_n(n,c))**2), zp(p)) )  &
            go to 1
        else if(answer .eq. 'RY') then
          if( Approx_Real( sqrt(grid % xn(grid % cells_n(n,c))**2 +           &
                                grid % zn(grid % cells_n(n,c))**2), zp(p)) )  &
            go to 1
        else if(answer .eq. 'RZ') then
          if( Approx_Real( sqrt(grid % xn(grid % cells_n(n,c))**2 +           &
                                grid % yn(grid % cells_n(n,c))**2), zp(p)) )  &
            go to 1
        end if
      end do

      ! Couldn't find a cell among the probes, add a new one
      n_prob = n_prob + 1
      if(answer .eq. 'X') zp(n_prob) = grid % xn(grid % cells_n(n,c))
      if(answer .eq. 'Y') zp(n_prob) = grid % yn(grid % cells_n(n,c))
      if(answer .eq. 'Z') zp(n_prob) = grid % zn(grid % cells_n(n,c))

      if(answer .eq. 'RX') zp(n_prob) =                                  &
                           sqrt(grid % zn(grid % cells_n(n,c))**2 +      &
                                grid % yn(grid % cells_n(n,c))**2)
      if(answer .eq. 'RY') zp(n_prob) =                                  &
                           sqrt(grid % xn(grid % cells_n(n,c))**2 +      &
                                grid % zn(grid % cells_n(n,c))**2)
      if(answer .eq. 'RZ') zp(n_prob) =                                  &
                           sqrt(grid % xn(grid % cells_n(n,c))**2 +      &
                                grid % yn(grid % cells_n(n,c))**2)

      if(n_prob .eq. 16384) then
        print *, '# Probe 1d: Not a 1d (channel flow) problem.'
        isit = .false.
        return
      end if
    end do
1 end do

  isit = .true.

  !--------------------!
  !   Create 1d file   !
  !--------------------!
  name_prob = problem_name
  name_prob(len_trim(problem_name)+1:len_trim(problem_name)+3) = '.1d'
  print *, '# Creating the file: ', trim(name_prob)
  open(9, file=name_prob)

  ! Write the number of probes 
  write(9,'(i8)') n_prob

  call Sort_Mod_Real(zp(1:n_prob))

  ! Write the probe coordinates out
  do p=1, n_prob
    write(9,'(i8,1e17.8)') p, zp(p)
  end do

  close(9)

  end subroutine
