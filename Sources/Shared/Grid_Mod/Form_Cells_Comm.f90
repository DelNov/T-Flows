!==============================================================================!
  subroutine Grid_Mod_Form_Cells_Comm(grid)
!------------------------------------------------------------------------------!
!   Find communication patterns for cells from Process                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: sub, ms, mr
  integer              :: c, c1, c2, s, n_buff_faces, n_max_buff_cells, i_cel
  integer, allocatable :: send_cells(:), recv_cells(:)
  integer, allocatable :: buff_s(:,:), buff_r(:,:), need_cell(:,:)
!==============================================================================!

  if(n_proc < 2) return

  ! Allocate memory for locally used arrays
  allocate(buff_s(n_proc, n_proc)); buff_s(:,:) = 0
  allocate(buff_r(n_proc, n_proc)); buff_r(:,:) = 0
  allocate(send_cells(-grid % n_bnd_cells:grid % n_cells))
  allocate(recv_cells(-grid % n_bnd_cells:grid % n_cells))

  !-------------------------------!
  !   Count buffer cells inside   !
  !-------------------------------!
  grid % comm % n_buff_cells = 0
  do c = 1, grid % n_cells
    if(grid % comm % cell_proc(c) .ne. this_proc) then
      grid % comm % n_buff_cells =    &
      grid % comm % n_buff_cells + 1
    end if
  end do

  !------------------------!
  !   Count buffer faces   !
  !------------------------!
  n_buff_faces = 0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c1 .eq. 0 .and. c2 .eq. 0) then
      n_buff_faces = n_buff_faces + 1
    end if
  end do

  !----------------------------!
  !   Form the buff_r matrix   !
  !----------------------------!
  do c = grid % n_cells - grid % comm % n_buff_cells + 1,  &
         grid % n_cells
    sub = grid % comm % cell_proc(c)
    buff_r(this_proc, sub) = buff_r(this_proc, sub) + 1
    recv_cells(c) = sub
  end do

  !--------------------------------------------------!
  !   Allocate memory for matrix with needed cells   !
  !--------------------------------------------------!
  n_max_buff_cells = grid % comm % n_buff_cells
  call comm_mod_global_max_int(n_max_buff_cells)
  allocate(need_cell(n_max_buff_cells, n_proc));  need_cell(:,:) = 0

  !-------------------------------------------!
  !   Store global number of cells you need   !
  !-------------------------------------------!
  i_cel = 0
  do c = grid % n_cells - grid % comm % n_buff_cells + 1,  &
         grid % n_cells
    if(grid % comm % cell_proc(c) .ne. this_proc) then
      i_cel = i_cel + 1
      need_cell(i_cel, this_proc) = grid % comm % cell_glo(c)
    end if
  end do

  !----------------------------------------------!
  !   Inform all processors about needed cells   !
  !----------------------------------------------!
  need_cell = reshape(need_cell, (/n_max_buff_cells*n_proc, 1/))
  call Comm_Mod_Global_Sum_Int_Array(n_max_buff_cells*n_proc, need_cell)
  need_cell = reshape(need_cell, (/n_max_buff_cells, n_proc/))

  !-----------------------------------------------------!
  !   Identify cells you have to send, and count them   !
  !-----------------------------------------------------!
  do sub = 1, n_proc
    if(sub .ne. this_proc) then
      do i_cel = 1, n_max_buff_cells
        do c = 1, grid % n_cells
          if(grid % comm % cell_glo(c) .eq. need_cell(i_cel, sub)) then
            send_cells(c) = sub                                  ! identify
            buff_s(this_proc, sub) = buff_s(this_proc, sub) + 1  ! count
          end if
        end do
      end do  ! i_cel
    end if    ! sub .ne. this_proc
  end do      ! sub

  !--------------------------------------------------!
  !   Allocate memory for send and receive buffers   !
  !--------------------------------------------------!
  allocate(grid % comm % cells_send(n_proc))
  allocate(grid % comm % cells_recv(n_proc))

  !-----------------------------------!
  !   Form send and receive buffers   !
  !-----------------------------------!
  do sub = 1, n_proc

    ! Initialize buffer size to zero
    grid % comm % cells_send(sub) % n_items = 0
    grid % comm % cells_recv(sub) % n_items = 0

    if(sub .ne. this_proc) then

      ms = buff_s(this_proc, sub)
      mr = buff_r(this_proc, sub)

      if(ms > 0) then
        allocate(grid % comm % cells_send(sub) % map(ms));
        allocate(grid % comm % cells_send(sub) % i_buff(ms));
        allocate(grid % comm % cells_send(sub) % l_buff(ms));
        allocate(grid % comm % cells_send(sub) % r_buff(ms));
      end if
      if(mr > 0) then
        allocate(grid % comm % cells_recv(sub) % map(mr));
        allocate(grid % comm % cells_recv(sub) % i_buff(mr));
        allocate(grid % comm % cells_recv(sub) % l_buff(mr));
        allocate(grid % comm % cells_recv(sub) % r_buff(mr));
      end if

      ms = 0
      mr = 0

      do c = -grid % n_bnd_cells, grid % n_cells
        if(send_cells(c) .eq. sub) then
          ms = ms + 1
          grid % comm % cells_send(sub) % map(ms) = c
        end if
      end do

      do c = -grid % n_bnd_cells, grid % n_cells
        if(recv_cells(c) .eq. sub) then
          mr = mr + 1
          grid % comm % cells_recv(sub) % map(mr) = c
        end if
      end do

      ! Store final buffer lengths
      grid % comm % cells_send(sub) % n_items = ms
      grid % comm % cells_recv(sub) % n_items = mr

    end if  ! sub .ne. this_proc
  end do    ! sub

  !-------------------------------------!
  !   Avoid faces in the buffers only   !
  !   (This also leaves shadows out)    !
  !-------------------------------------!
  grid % n_faces = grid % n_faces - n_buff_faces

  ! De-allocate locally used memory
  deallocate(buff_s)
  deallocate(buff_r)
  deallocate(send_cells)
  deallocate(recv_cells)

  end subroutine
