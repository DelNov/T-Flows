!==============================================================================!
  subroutine Save_Face_Fluxes(Backup, Flow, disp, vc)
!------------------------------------------------------------------------------!
!   Writes volume fluxes at faces (n, o and oo) to the backup file, so that   !
!   a restart continues exactly as the uninterrupted run (Rhie and Chow       !
!   interpolation with Choi's correction uses the fluxes of two old time      !
!   steps).  Each processor writes its own faces, in local numbering, as one  !
!   block; the blocks follow each other in the order of processors.  Hence,   !
!   they can be read only with the same domain decomposition, which is        !
!   checked with the layout (number of faces in each processor) saved first.  !
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
  integer                  :: n_p
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  Comm => Grid % Comm

  ! Number of faces in each processor
  n_p = N_Procs()
  allocate(layout(n_p));  layout(:) = 0
  layout(max(This_Proc(), 1)) = Grid % n_faces
  call Global % Sum_Int_Array(n_p, layout)

  call Backup % Save_Int      (Comm, disp, vc, 'v_flux_n_proc', n_p)
  call Backup % Save_Int_Array(Comm, disp, vc, 'v_flux_layout', layout)

  call Save_Face_Block(Comm, layout, disp, vc, 'v_flux_n',  Flow % v_flux % n)
  call Save_Face_Block(Comm, layout, disp, vc, 'v_flux_o',  Flow % v_flux % o)
  call Save_Face_Block(Comm, layout, disp, vc, 'v_flux_oo', Flow % v_flux % oo)

  end subroutine

!==============================================================================!
  subroutine Save_Face_Block(Comm, layout, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Writes a face array as blocks of processors, in local face numbering.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type)  :: Comm
  integer          :: layout(:)
  integer(DP)      :: disp
  integer          :: vc
  character(len=*) :: var_name
  real             :: array(:)
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs, i_p  ! variable size, this processor
  integer(DP)   :: disp_p   ! displacement of this processor's block
!==============================================================================!

  if(First_Proc()) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Name and size of the whole variable
  vn = var_name
  call Comm % Write_Text(fh, vn, disp)

  vs = sum(layout) * RP
  call Comm % Write_Int(fh, vs, disp)

  ! Each processor writes its block after the blocks of the previous ones
  i_p    = max(This_Proc(), 1)
  disp_p = disp + int(sum(layout(1:i_p-1)), DP) * RP
  call Comm % Write_Real_Array(fh, array(1:layout(i_p)), disp_p)

  disp = disp + vs

  end subroutine
