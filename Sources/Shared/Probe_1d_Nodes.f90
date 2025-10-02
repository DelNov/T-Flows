!==============================================================================!
  subroutine Probe_1d_Nodes(Grid)
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of nodes in non-homogeneous           !
!   direction and write them in file name.1d                                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use File_Mod
  use Math_Mod
  use Grid_Mod, only: Grid_Type, Profiler
  use Sort_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: n_prob, p, n, fu
  real          :: zp(16384)
  character(SL) :: name_prob
  character(SL) :: answer
  logical       :: isit
!==============================================================================!

  call Profiler % Start('Probe_1d_Nodes')

  print *, '#==========================================='
  print *, '# Creating 1d file with the node '
  print *, '# coordinates in non-homogeneous directions '
  print *, '#-------------------------------------------'
  print *, '# Insert non-homogeneous direction '
  print *, '# (x, y, z, rx, ry, rz or skip)'
  print *, '#-------------------------------------------'
  answer = File % Single_Word_From_Keyboard()
  call String % To_Upper_Case(answer)
  if(answer .eq. 'SKIP') then
    call Profiler % Stop('Probe_1d_Nodes')
    return
  end if

  n_prob = 0
  zp     = 0.0

  !-----------------------------!
  !   Browse through all nodes  !
  !-----------------------------!
  do n = 1, Grid % n_nodes

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer .eq. 'X') then
        if( Math % Approx_Real(Grid % xn(n), zp(p)) ) go to 1
      else if(answer .eq. 'Y') then
        if( Math % Approx_Real(Grid % yn(n), zp(p)) ) go to 1
      else if(answer .eq. 'Z') then
        if( Math % Approx_Real(Grid % zn(n), zp(p)) ) go to 1
      else if(answer .eq. 'RX') then
        if( Math % Approx_Real( sqrt(Grid % zn(n)**2 +   &
                                       Grid % yn(n)**2), zp(p)) ) go to 1
      else if(answer .eq. 'RY') then
        if( Math % Approx_Real( sqrt(Grid % xn(n)**2 +   &
                                       Grid % zn(n)**2), zp(p)) ) go to 1
      else if(answer .eq. 'RZ') then
        if( Math % Approx_Real( sqrt(Grid % xn(n)**2 +   &
                                       Grid % yn(n)**2), zp(p)) ) go to 1
      end if
    end do

    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob+1
    if(answer .eq. 'X') zp(n_prob)= Grid % xn(n)
    if(answer .eq. 'Y') zp(n_prob)= Grid % yn(n)
    if(answer .eq. 'Z') zp(n_prob)= Grid % zn(n)

    if(answer .eq. 'RX') zp(n_prob)= sqrt(Grid % zn(n)**2 + Grid % yn(n)**2)
    if(answer .eq. 'RY') zp(n_prob)= sqrt(Grid % xn(n)**2 + Grid % zn(n)**2)
    if(answer .eq. 'RZ') zp(n_prob)= sqrt(Grid % xn(n)**2 + Grid % yn(n)**2)

    if(n_prob .eq. 16384) then
      print *, '# Probe 1d: Not a 1d (channel flow) problem.'
      isit = .false.
      return
    end if
1  end do

  isit = .true.

  !--------------------!
  !   Create 1d file   !
  !--------------------!
  call File % Set_Name(name_prob, extension='.1d')
  call File % Open_For_Writing_Ascii(name_prob, fu)

  ! Write the number of probes 
  write(fu, '(i8)') n_prob

  call Sort % Real_Array(zp(1:n_prob))

  ! Write the probe coordinates out
  do p=1, n_prob
    write(fu,'(i8,1e17.8)') p, zp(p)
  end do

  close(fu)

  call Profiler % Stop('Probe_1d_Nodes')

  end subroutine
