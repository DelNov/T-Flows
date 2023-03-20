!==============================================================================!
  subroutine Save_Scalar_Real(Results,                &
                              var_name, plot_inside,  &
                              val, fs, fp,            &
                              data_offset, sweep)
!------------------------------------------------------------------------------!
!   Writes one real scalar defined over cells.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  character(len=*)    :: var_name
  logical             :: plot_inside     ! plot results inside?
  real                :: val(:)
  integer             :: fs, fp          ! file unit sequential and parallel
  integer             :: data_offset
  integer             :: sweep           ! is it the first or second sweep
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)   :: data_size
  integer       :: c1, c2, c_f, c_l
  character(SL) :: str1
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  data_size = 0

  c_f = lbound(val, 1)
  c_l = ubound(val, 1)

  ! Header
  if(sweep .eq. 1) then
    if(n_proc > 1 .and. this_proc .eq. 1) then
      write(fs) IN_4                                  //  &
                '<PDataArray type='//floatp           //  &
                ' Name="' // trim(var_name) // '"/>'  // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type='//floatp         //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then
    if(plot_inside) then
      data_size = int((c_l-c_f+1) * RP, SP)
      write(fp) data_size
      do c1 = c_f, c_l
        write(fp) val(c1)
      end do
    else
      do c2 = c_f, c_l
        data_size = int(data_size + RP, SP)
      end do
      write(fp) data_size
      do c2 = c_f, c_l
        write(fp) val(c2)
      end do
    end if
  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset + (c_l-c_f+1) * RP
    else
      do c2 = c_f, c_l
        data_offset = data_offset + RP
      end do
    end if
    data_offset = data_offset + SP
  end if

  end subroutine
