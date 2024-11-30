!==============================================================================!
  subroutine Form_Variable_Matrix(Turb, Grid, phi, Flow,  &
                                  visc_eff, urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Turb_Type)                      :: Turb
  type(Grid_Type),   intent(in), target :: Grid
  type(Var_Type),    intent(in), target :: phi
  type(Field_Type),              target :: Flow
  real                                  :: visc_eff(-Grid % n_bnd_cells  &
                                                    :Grid % n_cells)
  real                                  :: urf
  real,    optional, intent(in)         :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  real,      contiguous, pointer :: dens(:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12, a21, fl, cfs
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!==============================================================================!

  call Profiler % Start('Form_Variable_Matrix')

  !-----------------------!
  !   Take some aliases   !
  !-----------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % C % dia
  pos => Flow % Nat % C % pos
  fc  => Flow % Nat % C % fc
  nz  =  Flow % Nat % C % nonzeros

  Assert(urf > 0.0)

  !------------------------!
  !   Initialize density   !
  !------------------------!
  dens => Flow % density

  !-----------------------------------------------------------!
  !   Start by copying molecular viscosity to the effective   !
  !-----------------------------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    visc_eff(c) = Flow % viscosity(c)
  end do
  !$tf-acc loop end

  !-------------------------------------------------------------!
  !   If there is a turbulence model, add turbulent viscosity   !
  !-------------------------------------------------------------!
  !@ if(Turb % model .ne. NO_TURBULENCE_MODEL) then
  !@   !$tf-acc loop begin
  !@   do c = Cells_In_Domain_And_Buffers()
  !@     visc_eff(c) = visc_eff(c) + turb_vis_t(c) / phi % sigma
  !@   end do
  !@   !$tf-acc loop end
  !@ end if

  if(Turb % model .eq. SPALART_ALLMARAS.or.Turb % model .eq. DES_SPALART) then
    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()
      visc_eff(c) = visc_eff(c) + phi % n(c) * Flow % density(c) / phi % sigma
    end do
    !$tf-acc loop end
  end if

  !---------------------------------------!
  !   Initialize matrix entries to zero   !
  !---------------------------------------!

  !$tf-acc loop begin
  do i = 1, nz  ! all present
    val(i) = 0.0
  end do
  !$tf-acc loop end

  !--------------------------------------------------!
  !                                                  !
  !   Compute neighbouring coefficients over cells   !
  !                                                  !
  !--------------------------------------------------!

  !------------------------------------!
  !   Coefficients inside the domain   !
  !------------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)

      if(c2 .gt. 0) then

        a12 = Face_Value(s, visc_eff(c1), visc_eff(c2)) * fc(s)
        a21 = a12

        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a21
        end if

        ! Update only diaginal at c1 to avoid race conditions
        val(dia(c1)) = val(dia(c1)) + a12

      end if
    end do

  end do
  !$tf-acc loop end

  !---------------------------------------!
  !   Upwind blending inside the domain   !
  !---------------------------------------!
  if(phi % blend_matrix) then

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      fl = Flow % v_flux % n(s)

      if(c2 .gt. 0) then

        cfs = Face_Value(s, dens(c1), dens(c2))
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

  end do
  !$tf-acc loop end

  end if

  !------------------------------------!
  !   Coefficients on the boundaries   !
  !------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL    .or.  &
       Grid % region % type(reg) .eq. WALLFL  .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)   ! inside cell
        a12 = visc_eff(c1) * fc(s)
        val(dia(c1)) = val(dia(c1)) + a12
      end do
      !$tf-acc loop end

    end if
  end do

  !------------------------------------!
  !                                    !
  !   Take care of the unsteady term   !
  !                                    !
  !------------------------------------!
  if(present(dt)) then
    !$tf-acc loop begin
    do c = Cells_In_Domain()  ! all present, was independent
      val(dia(c)) = val(dia(c)) + dens(c) * Grid % vol(c) / dt
    end do
    !$tf-acc loop end
  end if

# if T_FLOWS_DEBUG == 1
  allocate(temp(Grid % n_cells));  temp(:) = 0.0
  do c = 1, Grid % n_cells  ! this is for debugging and should be on CPU
    temp(c) = val(dia(c))
  end do
  call Grid % Save_Debug_Vtu("a_turb_var_diagonal",              &
                             inside_name="a_turb_var_diagonal",  &
                             inside_cell=temp)
# endif

  call Profiler % Stop('Form_Variable_Matrix')

  end subroutine
