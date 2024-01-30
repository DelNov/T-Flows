!==============================================================================!
  subroutine Create_Vof(Vof, Flow)
!------------------------------------------------------------------------------!
!   Allocates memory for variables in Multphase_Mod.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),  target :: Vof
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Matrix_Type), pointer :: A
  integer                    :: nb, nc, nf
!==============================================================================!

  ! Give some sign
  if(First_Proc())  &
    print '(a)', ' # Creating the volume of fluid (VOF) module'

  ! Store pointers
  Vof % pnt_flow   => Flow
  Vof % pnt_grid   => Flow % pnt_grid
  Vof % pnt_matrix => Flow % pnt_matrix

  ! Take aliases
  Grid => Flow % pnt_grid
  A    => Flow % pnt_matrix
  nb   =  Grid % n_bnd_cells
  nc   =  Grid % n_cells
  nf   =  Grid % n_faces

  call Var_Mod_Create_Solution(Vof % fun,    A, 'VOF', '')
  call Var_Mod_Create_New_Only(Vof % smooth, Grid, 'SMO')

  ! Surface curvature
  allocate(Vof % curv(-nb:nc));  Vof % curv(-nb:nc) = 0.0

  ! Additional variables for calculation of vof function (fun)
  allocate(Vof % beta_f (nf));  Vof % beta_f(1:nf) = 0.0
  allocate(Vof % beta_c (nf));  Vof % beta_c(1:nf) = 0.0
  allocate(Vof % c_d(-nb:nc));  Vof % c_d (-nb:nc) = 0.0

  ! Surface normals
  allocate(Vof % nx(-nb:nc));  Vof % nx(-nb:nc) = 0.0
  allocate(Vof % ny(-nb:nc));  Vof % ny(-nb:nc) = 0.0
  allocate(Vof % nz(-nb:nc));  Vof % nz(-nb:nc) = 0.0

  ! Surface tension force
  allocate(Vof % surf_fx(-nb:nc));  Vof % surf_fx(-nb:nc) = 0.0
  allocate(Vof % surf_fy(-nb:nc));  Vof % surf_fy(-nb:nc) = 0.0
  allocate(Vof % surf_fz(-nb:nc));  Vof % surf_fz(-nb:nc) = 0.0

  if(Flow % mass_transfer_model .ne. NO_MASS_TRANSFER) then
    allocate(Vof % a12(nf));        Vof % a12(1:nf)     = 0.0
    allocate(Vof % a21(nf));        Vof % a21(1:nf)     = 0.0
    allocate(Vof % m_dot(-nb:nc));  Vof % m_dot(-nb:nc) = 0.0
    call Var_Mod_Create_New_Only(Vof % t_0, Grid, 'T_0')
    call Var_Mod_Create_New_Only(Vof % t_1, Grid, 'T_1')
  end if

  if(Vof % track_front) then
    call Vof % Front % Allocate_Front(Flow)
  end if

  if(Vof % track_surface) then
    call Vof % Surf % Allocate_Surf(Flow)
  end if

  end subroutine
