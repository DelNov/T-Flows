!==============================================================================!
  subroutine Info_Mod_Bulk_Fill(courant, peclet, fx, fy, fz, px, py, pz)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: courant, peclet, fx, fy, fz, px, py, pz
!==============================================================================!

  if(this_proc < 2) then

    ! Courant and Peclet numbers
    write(bulk_info % lines(1)( 27: 49),    '(a23)') 'Maximum Courant number:'
    write(bulk_info % lines(1)( 51: 60), '(1pe9.3)') courant

    write(bulk_info % lines(1)( 69: 90),    '(a22)') 'Maximum Peclet number:'
    write(bulk_info % lines(1)( 92:100), '(1pe9.3)') peclet 

    write(bulk_info % lines(2)( 27: 34),      '(a8)') 'Flux x :'
    write(bulk_info % lines(2)( 36: 45), '(1pe10.3)') fx
    write(bulk_info % lines(2)( 55: 62),      '(a8)') 'Flux y :'
    write(bulk_info % lines(2)( 64: 75), '(1pe10.2)') fy
    write(bulk_info % lines(2)( 83: 90),      '(a8)') 'Flux z :'
    write(bulk_info % lines(2)( 92:101), '(1pe10.2)') fz

    write(bulk_info % lines(3)( 27: 34),      '(a8)') 'Pdrop x:'
    write(bulk_info % lines(3)( 36: 45), '(1pe10.3)') px
    write(bulk_info % lines(3)( 55: 62),      '(a8)') 'Pdrop y:'
    write(bulk_info % lines(3)( 64: 75), '(1pe10.2)') py
    write(bulk_info % lines(3)( 83: 90),      '(a8)') 'Pdrop z:'
    write(bulk_info % lines(3)( 92:101), '(1pe10.2)') pz

  end if

  end subroutine
