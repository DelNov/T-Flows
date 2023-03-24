!==============================================================================!
  subroutine Save_Variable(Backup, disp, vc, var_name, var)
!------------------------------------------------------------------------------!
!   Writes a whole variable to backup file.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type) :: Backup
  integer(DP)        :: disp
  integer            :: vc
  character(len=*)   :: var_name
  type(Var_Type)     :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  character(SL)            :: vn
  integer                  :: vs  ! variable size
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Backup)
!==============================================================================!

  ! Take alias
  Comm => var % pnt_grid % Comm

  if(First_Proc()) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector without boundaries
  vn = var_name
  call Comm % Write_Text(fh, vn, disp)
  vs = (3*Comm % nc_tot + 3*Comm % nb_tot) * RP
  call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Cell_Real(fh, var % n(1:Comm % nc_sub), disp)
  call Comm % Write_Bnd_Real (fh, var % n( -Comm % nb_f:  &
                                           -Comm % nb_l),   disp)
  call Comm % Write_Cell_Real(fh, var % q(1:Comm % nc_sub), disp)
  call Comm % Write_Bnd_Real (fh, var % q( -Comm % nb_f:  &
                                           -Comm % nb_l),   disp)
  call Comm % Write_Cell_Real(fh, var % o(1:Comm % nc_sub), disp)
  call Comm % Write_Bnd_Real (fh, var % o( -Comm % nb_f:  &
                                           -Comm % nb_l),   disp)

  end subroutine
