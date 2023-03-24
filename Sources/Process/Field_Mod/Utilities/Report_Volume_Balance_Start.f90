!==============================================================================!
  subroutine Report_Volume_Balance_Start(Flow, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Opens file for volume balance reporting.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type)   :: Flow
  integer, intent(in) :: curr_dt, ini
!==============================================================================!

  if(First_Proc()) then
    if(Flow % report_vol_balance) then
      call File % Append_For_Writing_Ascii('volume.bal',  &
                                           Flow % fuvbr,  &
                                           This_Proc())
      if(Flow % p_m_coupling == PISO) then
        if(curr_dt .eq. 1 .and. ini .eq. 1 .and. Flow % i_corr .eq. 1) then
          write(Flow % fuvbr, '(a)')  '# Report on volume balance'
          write(Flow % fuvbr, '(a)')  '# ts- time step; '        //  &
                                        'oi- outer iteration; '  //  &
                                        'ii- inner iteration'
        end if
      else
        if(curr_dt .eq. 1 .and. ini .eq. 1) then
          write(Flow % fuvbr, '(a)')  '# Report on volume balance'
          write(Flow % fuvbr, '(a)')  '# ts- time step; '        //  &
                                        'oi- outer iteration'
        end if
      end if
    end if
  end if

  end subroutine
