!==============================================================================!
  subroutine Update_Aliases(Flow)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow object
!==============================================================================!

  ! Aliases connected to Flow
  u_n => Flow % u % n
  v_n => Flow % v % n
  w_n => Flow % w % n

  u_o => Flow % u % o
  v_o => Flow % v % o
  w_o => Flow % w % o

  if(Flow % u % td_scheme .eq. PARABOLIC) then
    u_oo => Flow % u % oo
    v_oo => Flow % v % oo
    w_oo => Flow % w % oo
  end if

  pp_n => Flow % pp % n
  p_n  => Flow % p  % n

  if(Flow % heat_transfer) then
    t_n => Flow % t % n
    t_o => Flow % t % o

    if(Flow % u % td_scheme .eq. PARABOLIC) then
      t_oo => Flow % t % oo
    end if
  end if

  end subroutine
