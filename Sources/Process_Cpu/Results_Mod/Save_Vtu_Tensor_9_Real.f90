!==============================================================================!
  subroutine Save_Vtu_Tensor_9_Real(Results,                 &
                                    var_name, plot_inside,   &
                                    val_11, val_12, val_13,  &
                                    val_21, val_22, val_23,  &
                                    val_31, val_32, val_33,  &
                                    fs, fp,                  &
                                    data_offset, sweep)
!------------------------------------------------------------------------------!
!>  Save_Vtu_Tensor_9_Real is designed to write a real, non-symmetric tensor
!>  defined over cells to a .vtu file as a single record.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Initial Setup: Establishes plotting precision and manages the data       !
!     range for the tensor components.                                         !
!   * Data Handling: Depending on the sweep stage, this subroutine either      !
!     writes the XML header for the tensor or the tensor data itself.          !
!   * Data Writing: Writes the tensor values in an appended format, addressing !
!     both internal and boundary cells as specified by 'plot_inside'.          !
!   * Offset Management: Updates the data offset post-writing to ensure        !
!     correct placement of subsequent data in the file.                        !
!------------------------------------------------------------------------------!
!   Workflow                                                                   !
!                                                                              !
!   * Header Preparation: In the first sweep, writes the XML header for the    !
!     tensor variable, including its name and data components.                 !
!   * Tensor Data Writing: In the second sweep, writes the actual tensor       !
!     values for either internal or boundary cells based on 'plot_inside'.     !
!   * Final Steps: Updates the data offset and completes the subroutine's      !
!     operations.                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results      !! parent class
  character(len=*)    :: var_name     !! name of the variable
  logical             :: plot_inside  !! true to plots inside,
                                      !! false to plot on the boundary
  real                :: val_11(:), val_12(:), val_13(:)  !! tensor component
  real                :: val_21(:), val_22(:), val_23(:)  !! tensor component
  real                :: val_31(:), val_32(:), val_33(:)  !! tensor component
  integer             :: fs, fp       !! file unit sequential and parallel
  integer             :: data_offset  !! offset in the .vtu file
  integer             :: sweep        !! is it the first or second sweep
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

  c_f = lbound(val_11, 1)
  c_l = ubound(val_11, 1)

  ! Header
  if(sweep .eq. 1) then
    if(Parallel_Run() .and. First_Proc()) then
      write(fs) IN_4                                //  &
                '<DataArray type='//floatp          //  &
                ' Name="' // trim(var_name) // '"'  //  &
                ' NumberOfComponents="9"/>'         // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type='//floatp         //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' NumberOfComponents="9"'          //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then
    if(plot_inside) then
      data_size = int((c_l-c_f+1) * RP * 9, SP)
      write(fp) data_size
      do c1 = c_f, c_l
        write(fp) val_11(c1), val_12(c1), val_13(c1),  &
                  val_21(c1), val_22(c1), val_23(c1),  &
                  val_31(c1), val_32(c1), val_33(c1)
      end do
    else
      do c2 = c_f, c_l
        data_size = int(data_size + RP * 9, SP)
      end do
      write(fp) data_size
      do c2 = c_f, c_l
        write(fp) val_11(c2), val_12(c2), val_13(c2),  &
                  val_21(c2), val_22(c2), val_23(c2),  &
                  val_31(c2), val_32(c2), val_33(c2)
      end do
    end if
  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset + (c_l-c_f+1) * RP * 9
    else
      do c2 = c_f, c_l
        data_offset = data_offset + RP * 9
      end do
    end if
    data_offset = data_offset + SP
  end if

  end subroutine
