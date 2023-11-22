!==============================================================================!
  subroutine Save_Vtu_Scalar_Real(Results,                &
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
  integer(SP)               :: data_size
  integer                   :: i, c1, c2, c_f, c_l
  character(SL)             :: str1
  real, pointer, contiguous :: buffer(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  call Work % Connect_Real_Cell(buffer)

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  data_size = 0

  c_f = lbound(val, 1)
  c_l = ubound(val, 1)

  ! Header
  if(sweep .eq. 1) then
    if(Parallel_Run() .and. First_Proc()) then
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

    call Profiler % Start('Save_Vtu_Results (scalar real)')

    if(plot_inside) then
      data_size = int((c_l-c_f+1) * RP, SP)
      write(fp) data_size
      i = 0
      do c1 = c_f, c_l
        i = i + 1
        buffer(i) = val(c1)
      end do
      write(fp) buffer(1:i)
    else
      do c2 = c_f, c_l
        data_size = int(data_size + RP, SP)
      end do
      write(fp) data_size
      i = 0
      do c2 = c_f, c_l
        i = i + 1
        buffer(i) = val(c2)
      end do
      write(fp) buffer(1:i)
    end if

    call Profiler % Stop('Save_Vtu_Results (scalar real)')

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

  call Work % Disconnect_Real_Cell(buffer)

  end subroutine
