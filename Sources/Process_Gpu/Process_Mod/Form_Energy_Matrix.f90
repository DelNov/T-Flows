!==============================================================================!
  subroutine Form_Energy_Matrix(Process, Grid, Flow, Turb, Acon, Aval, urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)                   :: Process
  type(Grid_Type),   intent(in), target :: Grid
  type(Field_Type),              target :: Flow
  type(Turb_Type),               target :: Turb
  type(Sparse_Con_Type),         target :: Acon
  type(Sparse_Val_Type),         target :: Aval
  real                                  :: urf
  real,    optional, intent(in)         :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  real,      contiguous, pointer :: dens_capa(:)
  real,      contiguous, pointer :: cond_eff(:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Energy_Matrix')

  !------------------------------------------------------------!
  !   First take some aliases, which is quite elaborate here   !
  !------------------------------------------------------------!

  ! If each varible uses its own matrix
  if(.not. Flow % Nat % use_one_matrix) then
    if(Flow % Nat % A(MATRIX_UVW) % formed) return
  end if

  val => Aval % val
  dia => Acon % dia
  pos => Acon % pos
  fc  => Acon % fc
  nz  =  Acon % nonzeros

  call Work % Connect_Real_Cell(cond_eff, dens_capa)

  Assert(urf > 0.0)

  !--------------------------------------------------------------------------!
  !   Initialize density times thermal capacity and effective conductivity   !
  !--------------------------------------------------------------------------!

  !$acc parallel loop independent &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   dens_capa,  &
  !$acc   flow_density,  &
  !$acc   flow_capacity   &
  !$acc )
  do c = grid_region_f_cell(1), grid_region_l_cell(grid_n_regions+1)  ! all present
    dens_capa(c) = flow_density(c) * flow_capacity(c)
  end do
  !$acc end parallel

  ! Just copy molecular conductivity to effective
  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   cond_eff,  &
  !$acc   flow_conductivity   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)
    cond_eff(c) = flow_conductivity(c)
  end do
  !$acc end parallel

  ! If there is a turbulence model, add turbulent conductivity
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    !$acc parallel loop independent  &
    !$acc present(  &
    !$acc   grid_region_f_cell,  &
    !$acc   grid_region_l_cell,  &
    !$acc   cond_eff,  &
    !$acc   turb_vis_t   &
    !$acc )
    do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)
      cond_eff(c) = cond_eff(c) + turb_vis_t(c) / 0.9  ! hard-coded Pr_t
    end do
    !$acc end parallel
  end if

  !---------------------------------------!
  !   Initialize matrix entries to zero   !
  !---------------------------------------!

  !$acc parallel loop independent &
  !$acc present(  &
  !$acc   val   &
  !$acc )
  do i = 1, nz  ! all present
    val(i) = 0.0
  end do
  !$acc end parallel

  !--------------------------------------------------!
  !   Compute neighbouring coefficients over cells   !
  !--------------------------------------------------!

  ! Coefficients inside the domain

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   grid_cells_n_cells,  &
  !$acc   grid_cells_c,  &
  !$acc   grid_cells_f,  &
  !$acc   cond_eff,  &
  !$acc   fc,  &
  !$acc   val,  &
  !$acc   pos,  &
  !$acc   dia   &
  !$acc )
  do c1 = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present, was independent

  !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)

      if(c2 .gt. 0) then
        a12 = Face_Value(s, cond_eff(c1), cond_eff(c2)) * fc(s)
        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a12
        end if
        val(dia(c1)) = val(dia(c1)) + a12

      end if
    end do
  !$acc end loop

  end do
  !$acc end parallel

  ! Coefficients on the boundaries

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL    .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$acc parallel loop  &
      !$acc present(  &
      !$acc   grid_region_f_face,  &
      !$acc   grid_region_l_face,  &
      !$acc   grid_faces_c,  &
      !$acc   cond_eff,  &
      !$acc   fc,  &
      !$acc   val,  &
      !$acc   dia   &
      !$acc )
      do s = grid_region_f_face(reg), grid_region_l_face(reg)  ! all present
        c1 = grid_faces_c(1,s)   ! inside cell
        a12 = cond_eff(c1) * fc(s)
        val(dia(c1)) = val(dia(c1)) + a12
      end do
      !$acc end parallel

    end if
  end do

  !------------------------------------!
  !   Take care of the unsteady term   !
  !------------------------------------!
  if(present(dt)) then
    !$acc parallel loop independent  &
    !$acc present(  &
    !$acc   grid_region_f_cell,  &
    !$acc   grid_region_l_cell,  &
    !$acc   val,  &
    !$acc   dia,  &
    !$acc   dens_capa,  &
    !$acc   grid_vol   &
    !$acc )
    do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present, was independent
      val(dia(c)) = val(dia(c)) + dens_capa(c) * grid_vol(c) / dt
    end do
    !$acc end parallel
  end if

  !------------------------------------!
  !   Part 1 of the under-relaxation   !
  !   (Part 2 is in Compute_Energy)    !
  !------------------------------------!
  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   val,  &
  !$acc   dia   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present, was independent
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$acc end parallel

  !-------------------------------!
  !   Mark the matrix as formed   !
  !-------------------------------!
  Aval % formed = .true.

  call Work % Disconnect_Real_Cell(cond_eff, dens_capa)

  call Profiler % Stop('Form_Energy_Matrix')

  end subroutine
