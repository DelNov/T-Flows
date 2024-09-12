!==============================================================================!
  subroutine Form_Momentum_Matrix(Process, Grid, Flow, Turb, Acon, Aval,  &
                                  urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)           :: Process
  type(Grid_Type), intent(in), target :: Grid
  type(Field_Type),            target :: Flow
  type(Turb_Type),             target :: Turb
  type(Sparse_Con_Type),       target :: Acon
  type(Sparse_Val_Type),       target :: Aval
  real                                :: urf
  real,  optional, intent(in)         :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  real,      contiguous, pointer :: dens(:), visc_eff(:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Momentum_Matrix')

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

  call Work % Connect_Real_Cell(visc_eff)

  Assert(urf > 0.0)

  !-----------------------------------------------!
  !   Initialize density and efficent viscosity   !
  !-----------------------------------------------!
  dens => flow_density

  ! Just copy molecular viscosity to effective
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    visc_eff(c) = flow_viscosity(c)
  end do
  !$tf-acc loop end

  ! If there is a turbulence model, add turbulent viscosity
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()
      visc_eff(c) = visc_eff(c) + turb_vis_t(c)
    end do
    !$tf-acc loop end
  end if

  !---------------------------------------!
  !   Initialize matrix entries to zero   !
  !---------------------------------------!

  !$acc parallel loop independent  &
  !$acc present(val)
  do i = 1, nz  ! all present
    val(i) = 0.0
  end do
  !$acc end parallel

  !--------------------------------------------------!
  !   Compute neighbouring coefficients over cells   !
  !--------------------------------------------------!

  ! Coefficients inside the domain
  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)

      ! Coefficients inside the domain
      if(c2 .gt. 0) then
        a12 = Face_Value(s, visc_eff(c1), visc_eff(c2)) * fc(s)
        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a12
        end if
        val(dia(c1)) = val(dia(c1)) + a12
      end if
    end do

  end do
  !$tf-acc loop end

  ! Coefficients on the boundaries
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL    .or.  &
       Grid % region % type(reg) .eq. WALLFL  .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$acc parallel loop                                                  &
      !$acc present(grid_faces_c, grid_region_f_face, grid_region_l_face,  &
      !$acc         val, dia, visc_eff, fc, dia)
      do s = Faces_In_Region_Gpu(reg)  ! all present
        c1 = grid_faces_c(1,s)  ! inside cell
        a12 = visc_eff(c1) * fc(s)
        val(dia(c1)) = val(dia(c1)) + a12
      end do
      !$acc end parallel

    end if
  end do

  !------------------------------------!
  !   Take care of the unsteady term   !
  !------------------------------------!
  if(present(dt)) then
    !$tf-acc loop begin
    do c = Cells_In_Domain()  ! all present, was independent
      val(dia(c)) = val(dia(c)) + dens(c) * Grid % vol(c) / dt
    end do
    !$tf-acc loop end
  end if

  !--------------------------------------------------------------!
  !   Store volume divided by central coefficient for momentum   !
  !   and refresh its buffers before discretizing the pressure   !
  !--------------------------------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present, was independent
    Flow % v_m(c) = Grid % vol(c) / val(dia(c))
  end do
  !$tf-acc loop end

  ! This call is needed, the above loop goes through inside cells only
  call Grid % Exchange_Inside_Cells_Real(Flow % v_m)

  !-------------------------------------!
  !   Part 1 of the under-relaxation    !
  !   (Part 2 is in Compute_Momentum)   !
  !-------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present, was independent
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$tf-acc loop end

  !-------------------------------!
  !   Mark the matrix as formed   !
  !-------------------------------!
  Aval % formed = .true.

# if T_FLOWS_DEBUG == 1
  allocate(temp(Grid % n_cells));  temp(:) = 0.0
  do c = Cells_In_Domain()  ! this is for debugging, don't do it on GPU
    ! or: temp(c) = val(dia(c))
    ! or: temp(c) = Acon % row(c+1) - Acon % row(c)
    temp(c) = flow_v_m(c)
  end do
  call Grid % Exchange_Inside_Cells_Real(temp)
  call Grid % Save_Debug_Vtu("v_m",              &
                             inside_name="v_m",  &
                             inside_cell=temp)
# endif

  call Work % Disconnect_Real_Cell(visc_eff)

  call Profiler % Stop('Form_Momentum_Matrix')

  end subroutine
