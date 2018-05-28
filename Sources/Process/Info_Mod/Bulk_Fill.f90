!==============================================================================!
  subroutine Info_Mod_Bulk_Fill(courant, peclet, fx, fy, fz, px, py, pz)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: courant, peclet, fx, fy, fz, px, py, pz
!==============================================================================!

  if (this_proc < 2) then

    ! Courant and Peclet numbers
    write(bulk_info % lines(1)( 6:28),    '(a23)') 'Maximum Courant number:'
    write(bulk_info % lines(1)(30:38), '(1pe9.3)') courant

    write(bulk_info % lines(1)(48:69),    '(a22)') 'Maximum Peclet number:'
    write(bulk_info % lines(1)(71:79), '(1pe9.3)') peclet 

    write(bulk_info % lines(2)( 6:13),     '(a8)') 'Flux x :'
    write(bulk_info % lines(2)(15:23), '(1pe9.2)') fx
    write(bulk_info % lines(2)(34:41),     '(a8)') 'Flux y :'
    write(bulk_info % lines(2)(43:51), '(1pe9.2)') fy
    write(bulk_info % lines(2)(62:69),     '(a8)') 'Flux z :'
    write(bulk_info % lines(2)(71:79), '(1pe9.2)') fz

    write(bulk_info % lines(3)( 6:13),     '(a8)') 'Pdrop x:'
    write(bulk_info % lines(3)(15:23), '(1pe9.2)') px
    write(bulk_info % lines(3)(34:41),     '(a8)') 'Pdrop y:'
    write(bulk_info % lines(3)(43:51), '(1pe9.2)') py
    write(bulk_info % lines(3)(62:69),     '(a8)') 'Pdrop z:'
    write(bulk_info % lines(3)(71:79), '(1pe9.2)') pz

  end if

  end subroutine
