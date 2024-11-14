!==============================================================================!
  subroutine Form_Energy_Matrix(Process, Grid, Flow, Turb, Aval,  &
                                cond_eff, urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)                   :: Process
  type(Grid_Type),   intent(in), target :: Grid
  type(Field_Type),              target :: Flow
  type(Turb_Type),               target :: Turb
  type(Sparse_Val_Type),         target :: Aval
  real                                  :: cond_eff(-Grid % n_bnd_cells &
                                                    :Grid % n_cells)
  real                                  :: urf
  real,    optional, intent(in)         :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  real,      contiguous, pointer :: dens_capa(:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12, a21, fl, cfs
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Energy_Matrix')

  ! If each varible uses its own matrix and this matrix was already formed
  if(.not. Flow % Nat % use_one_matrix) then
    if(Flow % Nat % A(MATRIX_T) % formed) return
  end if

  !-----------------------!
  !   Take some aliases   !
  !-----------------------!
  val => Aval % val
  dia => Flow % Nat % C % dia
  pos => Flow % Nat % C % pos
  fc  => Flow % Nat % C % fc
  nz  =  Flow % Nat % C % nonzeros

  call Work % Connect_Real_Cell(dens_capa)

  Assert(urf > 0.0)

  !--------------------------------------------------------------------------!
  !   Initialize density times thermal capacity and effective conductivity   !
  !--------------------------------------------------------------------------!

  !$acc parallel loop independent  &
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

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   val   &
  !$acc )
  do i = 1, nz  ! all present
    val(i) = 0.0
  end do
  !$acc end parallel

  !--------------------------------------------------!
  !                                                  !
  !   Compute neighbouring coefficients over cells   !
  !                                                  !
  !--------------------------------------------------!

  !------------------------------------!
  !   Coefficients inside the domain   !
  !------------------------------------!

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
        a21 = a12

        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a21
        end if

        ! Update only diaginal at c1 to avoid race conditions
        val(dia(c1)) = val(dia(c1)) + a12

      end if
    end do
  !$acc end loop

  end do
  !$acc end parallel

  !---------------------------------------!
  !   Upwind blending inside the domain   !
  !---------------------------------------!
  if(Flow % t % blend_matrix) then

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   grid_cells_n_cells,  &
  !$acc   grid_cells_c,  &
  !$acc   grid_cells_f,  &
  !$acc   flow_v_flux_n,  &
  !$acc   dens_capa,  &
  !$acc   val,  &
  !$acc   pos,  &
  !$acc   dia   &
  !$acc )
  do c1 = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present, was independent

  !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      fl = flow_v_flux_n(s)

      if(c2 .gt. 0) then

        cfs = Face_Value(s, dens_capa(c1), dens_capa(c2))
        a12 = 0.0
        a21 = 0.0

        if(c1 .lt. c2) then
          if(fl > 0.0) a21 = a21 + fl * cfs
          if(fl < 0.0) a12 = a12 - fl * cfs
          val(pos(1,s)) = val(pos(1,s)) - a12
          val(pos(2,s)) = val(pos(2,s)) - a21
        end if

        if(c1 .gt. c2) then
          if(fl > 0.0) a12 = a12 + fl * cfs
        end if

        ! Update only diaginal at c1 to avoid race conditions
        val(dia(c1)) = val(dia(c1)) + a12

      end if
    end do
  !$acc end loop

  end do
  !$acc end parallel

  end if

  !------------------------------------!
  !   Coefficients on the boundaries   !
  !------------------------------------!
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
  !                                    !
  !   Take care of the unsteady term   !
  !                                    !
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
  !                                    !
  !   Part 1 of the under-relaxation   !
  !   (Part 2 is in Compute_Energy)    !
  !                                    !
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

  call Work % Disconnect_Real_Cell(dens_capa)

  call Profiler % Stop('Form_Energy_Matrix')

  end subroutine
