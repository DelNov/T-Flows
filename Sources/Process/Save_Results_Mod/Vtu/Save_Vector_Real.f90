!==============================================================================!
  subroutine Save_Vector_Real(grid, var_name, plot_inside,  &
                              val_1, val_2, val_3, fs, fp,  &
                              data_offset, sweep)
!------------------------------------------------------------------------------!
!   Writes one real vector defined over cells.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: var_name
  logical          :: plot_inside     ! plot results inside?
  real             :: val_1(-grid % n_bnd_cells:grid % n_cells)
  real             :: val_2(-grid % n_bnd_cells:grid % n_cells)
  real             :: val_3(-grid % n_bnd_cells:grid % n_cells)
  integer          :: fs, fp          ! file unit sequential and parallel
  integer          :: data_offset
  integer          :: sweep           ! is it the first or second sweep
!-----------------------------------[Locals]-----------------------------------!
  integer(4)         :: data_size
  integer            :: c, c2, s
  character(len=80)  :: str1
  integer, parameter :: IP=8, RP=8, SP=4
!==============================================================================!

  ! Header
  if(sweep .eq. 1) then
    if(n_proc > 1 .and. this_proc .eq. 1) then
      write(fs) IN_4                                //  &
                '<DataArray type="Float64"'         //  &
                ' Name="' // trim(var_name) // '"'  //  &
                ' NumberOfComponents="3"/>'         // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type="Float64"'        //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' NumberOfComponents="3"'          //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then
    if(plot_inside) then
      data_size = (grid % n_cells - grid % comm % n_buff_cells) * RP * 3
      write(fp) data_size
      do c = 1, grid % n_cells - grid % comm % n_buff_cells
        write(fp) val_1(c), val_2(c), val_3(c)
      end do
    else
      do s = 1, grid % n_faces
        c2 = grid % faces_c(2,s)
        if( c2 < 0 ) then
          data_size = data_size + RP * 3
        end if
      end do
      write(fp) data_size
      do s = 1, grid % n_faces
        c2 = grid % faces_c(2,s)
        if( c2 < 0 ) then
          write(fp) val_1(c2), val_2(c2), val_3(c2)
        end if
      end do
    end if
  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset  &
                  + (grid % n_cells - grid % comm % n_buff_cells) * RP * 3
    else
      do s = 1, grid % n_faces
        c2 = grid % faces_c(2,s)
        if( c2 < 0 ) then
          data_offset = data_offset + RP * 3
        end if
      end do
    end if
    data_offset = data_offset + SP
  end if

  end subroutine
