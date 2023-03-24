!==============================================================================!
  subroutine Info_Mod_Bulk_Print(Flow, dom, n_dom)
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: Flow
  integer          :: dom, n_dom
!==============================================================================!

  call Info_Mod_Bulk_Fill(Flow)

  if(First_Proc()) then

    ! String is L_LINE+2 long
    if(dom .eq. 1) then
      print '(a129)', trim(bulk_info % line_lead)
    else
      print '(a108)', trim(bulk_info % line_foll)
    end if
    print '(a108)', trim(bulk_info % lines(1))
    print '(a108)', trim(bulk_info % line_sep)
    print '(a108)', trim(bulk_info % lines(2))
    print '(a108)', trim(bulk_info % lines(3))
    if(dom .eq. n_dom) then
      print '(a108)', trim(bulk_info % line_trail)
      print '(a)',    ''
    end if

  end if

  end subroutine
