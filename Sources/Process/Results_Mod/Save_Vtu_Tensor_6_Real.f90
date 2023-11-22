!==============================================================================!
  subroutine Save_Vtu_Tensor_6_Real(Results,                 &
                                    var_name, plot_inside,   &
                                    val_11, val_22, val_33,  &
                                    val_12, val_13, val_23,  &
                                    fs, fp,                  &
                                    data_offset, sweep)
!------------------------------------------------------------------------------!
!   Writes one real symmetric tensor defined over cells.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  character(len=*)    :: var_name
  logical             :: plot_inside     ! plot results inside?
  real                :: val_11(:), val_22(:), val_33(:)
  real                :: val_12(:), val_13(:), val_23(:)
  integer             :: fs, fp          ! file unit sequential and parallel
  integer             :: data_offset
  integer             :: sweep           ! is it the first or second sweep
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)   :: data_size
  integer       :: i, c1, c2, c_f, c_l
  character(SL) :: str1
  real, allocatable :: buffer(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  data_size = 0

  c_f = lbound(val_11, 1)
  c_l = ubound(val_11, 1)

  ! Allocate local array
  allocate(buffer((c_l-c_f+1)*6))

  ! Header
  if(sweep .eq. 1) then
    if(Parallel_Run() .and. First_Proc()) then
      write(fs) IN_4                                //  &
                '<DataArray type='//floatp          //  &
                ' Name="' // trim(var_name) // '"'  //  &
                ' NumberOfComponents="6"/>'         // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type='//floatp         //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' NumberOfComponents="6"'          //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then
    if(plot_inside) then
      data_size = int((c_l-c_f+1) * RP * 6, SP)
      write(fp) data_size
      i = 0
      do c1 = c_f, c_l
        i = i + 1;  buffer(i) = val_11(c1)
        i = i + 1;  buffer(i) = val_22(c1)
        i = i + 1;  buffer(i) = val_33(c1)
        i = i + 1;  buffer(i) = val_12(c1)
        i = i + 1;  buffer(i) = val_13(c1)
        i = i + 1;  buffer(i) = val_23(c1)
      end do
      write(fp) buffer(1:i)
    else
      do c2 = c_f, c_l
        data_size = int(data_size + RP * 6, SP)
      end do
      write(fp) data_size
      i = 0
      do c2 = c_f, c_l
        i = i + 1;  buffer(i) = val_11(c2)
        i = i + 1;  buffer(i) = val_22(c2)
        i = i + 1;  buffer(i) = val_33(c2)
        i = i + 1;  buffer(i) = val_12(c2)
        i = i + 1;  buffer(i) = val_13(c2)
        i = i + 1;  buffer(i) = val_23(c2)
      end do
      write(fp) buffer(1:i)
    end if
  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset + (c_l-c_f+1) * RP * 6
    else
      do c2 = c_f, c_l
        data_offset = data_offset + RP * 6
      end do
    end if
    data_offset = data_offset + SP
  end if

  end subroutine
