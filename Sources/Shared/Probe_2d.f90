!==============================================================================!
  subroutine Probe_2D(Grid)
!------------------------------------------------------------------------------!
! Finds coordinates of all the planes for the channel flow.                    !
! It assumes that homogeneous directions of the flow are x and y.              !
!------------------------------------------------------------------------------!
  use Math_Mod
  use File_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: n_prob, p, c, fu
  real          :: yp(32768), zp(32768)
  character(SL) :: name_prob
  character(SL) :: answer
!==============================================================================!

  call Profiler % Start('Probe_2d')

  print *, '#==============================='
  print *, '# Looking for homogeneous plane '
  print *, '#-------------------------------'
  print *, '# Insert homogeneous direction  '
  print *, '# (xy, yz, zx or skip)'
  print *, '#-------------------------------'
  read(*,*) answer
  call String % To_Upper_Case(answer)
  if(answer .eq. 'SKIP') then
    call Profiler % Stop('Probe_2d')
    return
  end if

  n_prob = 0
  zp(:)  = 0.0
  yp(:)  = 0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c=1,Grid % n_cells

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer .eq. 'YZ') then
        if( Math % Approx_Real(Grid % yc(c), yp(p)) .and.      &
            Math % Approx_Real(Grid % zc(c), zp(p)) ) go to 1
      else if(answer .eq. 'ZX') then
        if( Math % Approx_Real(Grid % xc(c), yp(p)) .and.      &
            Math % Approx_Real(Grid % zc(c), zp(p)) ) go to 1
      else if(answer .eq. 'XY') then
        if( Math % Approx_Real(Grid % xc(c), yp(p)) .and.      &
            Math % Approx_Real(Grid % yc(c), zp(p)) ) go to 1
      end if
    end do 

    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob+1
    if(answer .eq. 'YZ') then
      yp(n_prob)=Grid % yc(c)
      zp(n_prob)=Grid % zc(c)
    else if(answer .eq. 'ZX') then
      yp(n_prob)=Grid % xc(c)
      zp(n_prob)=Grid % zc(c)
    else if(answer .eq. 'XY') then
      yp(n_prob)=Grid % xc(c)
      zp(n_prob)=Grid % yc(c)
    end if 

    if(n_prob .eq. 32768) then
      print *, '# Probe 2D: Not a 2D problem.'
      call Profiler % Stop('Probe_2d')
      return
    end if
1 end do

  !--------------------!
  !   Create 2D file   !
  !--------------------!
  call File % Set_Name(name_prob, extension='.2d')
  call File % Open_For_Writing_Ascii(name_prob, fu)

  ! Write the number of probes 
  write(fu,'(i8)') n_prob

  ! Write the probe coordinates out
  do p=1, n_prob
    write(fu,'(i8,1pe17.8,1pe17.8)') p, yp(p), zp(p)
  end do

  close(fu)

  call Profiler % Stop('Probe_2d')

  end subroutine
