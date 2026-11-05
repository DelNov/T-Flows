!==============================================================================!
  subroutine Load_Face_Fluxes(Backup, Flow, disp, vc)
!------------------------------------------------------------------------------!
!   Reads volume fluxes at faces (n, o and oo) written by Save_Face_Fluxes.    !
!   They are used only if the backup was written with the same domain         !
!   decomposition (same number of processors and faces in each of them).      !
!   Otherwise (or for older backups without them) the fluxes are left as      !
!   they were, and Flow % v_flux_restored stays false.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)       :: Backup
  type(Field_Type), target :: Flow
  integer(DP)              :: disp
  integer                  :: vc
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Grid_Type), pointer :: Grid
  integer, allocatable     :: layout(:)
  integer                  :: n_p, n_saved, n_wrong
  logical                  :: found_n, found_o, found_oo
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  Comm => Grid % Comm

  Flow % v_flux_restored = .false.

  ! Was it saved, and with the same number of processors?
  n_p = N_Procs()
  call Backup % Load_Int(Comm, disp, vc, 'v_flux_n_proc', n_saved)
  if(n_saved .ne. n_p) then
    if(First_Proc()) print *, '# NOTE! Face fluxes not in the backup or ',  &
                              'saved with other number of processors; ',     &
                              'Choi''s correction starts after 3 time steps'
    return
  end if

  ! Same number of faces in each processor?
  allocate(layout(n_p))
  call Backup % Load_Int_Array(Comm, disp, vc, 'v_flux_layout', layout)
  n_wrong = 0
  if(layout(max(This_Proc(), 1)) .ne. Grid % n_faces) n_wrong = 1
  call Global % Sum_Int(n_wrong)
  if(n_wrong > 0) then
    if(First_Proc()) print *, '# NOTE! Face fluxes in the backup were ',  &
                              'saved with other domain decomposition; ',   &
                              'Choi''s correction starts after 3 time steps'
    return
  end if

  call Load_Face_Block(Comm, layout, vc, 'v_flux_n',  Flow % v_flux % n,  found_n)
  call Load_Face_Block(Comm, layout, vc, 'v_flux_o',  Flow % v_flux % o,  found_o)
  call Load_Face_Block(Comm, layout, vc, 'v_flux_oo', Flow % v_flux % oo, found_oo)

  Flow % v_flux_restored = found_n .and. found_o .and. found_oo

  end subroutine

!==============================================================================!
  subroutine Load_Face_Block(Comm, layout, vc, var_name, array, found)
!------------------------------------------------------------------------------!
!   Reads a face array written by Save_Face_Block, each processor its block.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type)  :: Comm
  integer          :: layout(:)
  integer          :: vc
  character(len=*) :: var_name
  real             :: array(:)
  logical          :: found
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs, cnt_loop, i_p
  integer(DP)   :: disp_loop, disp_p
!==============================================================================!

  found     = .false.
  cnt_loop  = 0
  disp_loop = 0

  !--------------------------------------------------------!
  !   Browse the entire file until you find the variable   !
  !--------------------------------------------------------!
  do

    ! Increase counter
    cnt_loop = cnt_loop + 1

    call Comm % Read_Text(fh, vn, disp_loop)  ! variable name
    call Comm % Read_Int (fh, vs, disp_loop)  ! variable size

    ! If variable is found, read this processor's block and return
    if(vn .eq. var_name) then
      if(First_Proc()) print *, '# Reading variable: ', trim(vn)
      i_p    = max(This_Proc(), 1)
      disp_p = disp_loop + int(sum(layout(1:i_p-1)), DP) * RP
      call Comm % Read_Real_Array(fh, array(1:layout(i_p)), disp_p)
      found = .true.
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) exit

  end do

  if(First_Proc()) print *, '# Variable: ', trim(var_name), ' not found!'

  end subroutine
