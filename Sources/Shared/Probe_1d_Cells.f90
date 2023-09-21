!==============================================================================!
  subroutine Probe_1d_Cells(Grid)
!------------------------------------------------------------------------------!
!   This subroutine finds the coordinate of cell-centers in non-homogeneous    !
!   direction and write them in file called "name.1d"                          !
!------------------------------------------------------------------------------!
  use File_Mod
  use Const_Mod, only: SL, NANO          ! 1.0e-fu
  use Math_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: Grid
  logical         :: isit
!-----------------------------------[Locals]-----------------------------------!
  integer       :: n_prob, p, c, fu
  real          :: zp(16384)
  character(SL) :: name_prob
  character(SL) :: answer
!==============================================================================!

  call Profiler % Start('Probe_1d_Cells')

  print *, '#================================================='
  print *, '# Looking for non-homogeneous directions '
  print *, '# Insert non-homogeneous direction (x,y,z or skip)'
  print *, '#-------------------------------------------------'
  answer = File % Single_Word_From_Keyboard()
  call To_Upper_Case(answer)
  if(answer .eq. 'SKIP') then
    call Profiler % Stop('Probe_1d_Cells')
    return
  end if

  n_prob = 0
  zp(:)  = 0.0

  !-----------------------------!
  !   Browse through all cells  !
  !-----------------------------!
  do c = 1, Grid % n_cells

    ! Try to find the cell among the probes
    do p=1,n_prob
      if(answer .eq. 'X') then
        if( Math_Mod_Approx_Real(Grid % xc(c), zp(p), NANO)) go to 1
      else if(answer .eq. 'Y') then
        if( Math_Mod_Approx_Real(Grid % yc(c), zp(p), NANO)) go to 1
      else if(answer .eq. 'Z') then
        if( Math_Mod_Approx_Real(Grid % zc(c), zp(p), NANO)) go to 1
      end if
    end do

    ! Couldn't find a cell among the probes, add a new one
    n_prob = n_prob + 1
    if(answer .eq. 'X') zp(n_prob) = Grid % xc(c)
    if(answer .eq. 'Y') zp(n_prob) = Grid % yc(c)
    if(answer .eq. 'Z') zp(n_prob) = Grid % zc(c)

    if(n_prob .eq. 16384) then
      print *, '# Probe 1D: Not a 1D (channel flow) problem.'
      isit = .false.
      call Profiler % Stop('Probe_1d_Cells')
      return
    end if
1 end do

  isit = .true.

  !--------------------!
  !   Create 1D file   !
  !--------------------!
  call File_Mod_Set_Name(name_prob, extension='.1d')
  call File_Mod_Open_File_For_Writing(name_prob, fu)

  ! Write the number of probes
  write(fu,'(i8)') n_prob

  ! Write the probe coordinates out
  do p=1,n_prob
    write(fu,'(i8,1pe17.8)') p, zp(p)
  end do

  close(fu)

  call Profiler % Stop('Probe_1d_Cells')

  end subroutine
