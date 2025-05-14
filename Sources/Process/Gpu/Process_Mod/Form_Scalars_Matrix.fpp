!==============================================================================!
  subroutine Form_Scalars_Matrix(Process, Grid, Flow, Turb,  &
                                 diff_eff, sc, urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)                  :: Process
  type(Grid_Type),  intent(in), target :: Grid
  type(Field_Type),             target :: Flow
  type(Turb_Type),              target :: Turb
  real                                 :: diff_eff(-Grid % n_bnd_cells &
                                                   :Grid % n_cells)
  real                                 :: urf
  integer,          intent(in)         :: sc       !! scalar rank
  real,   optional, intent(in)         :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),      pointer :: phi
  real,    contiguous, pointer :: val(:), fc(:)
  integer, contiguous, pointer :: dia(:), pos(:,:)
  integer                      :: c, s, c1, c2, i_cel, reg, nz, i
  real                         :: a12, a21, fl, cfs
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Scalars_Matrix')

  !-----------------------!
  !   Take some aliases   !
  !-----------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % C % dia
  pos => Flow % Nat % C % pos
  fc  => Flow % Nat % C % fc
  phi => Flow % scalar(sc)
  nz  =  Flow % Nat % C % nonzeros

  Assert(urf > 0.0)

  !-----------------------------------------------------------!
  !   Start by copying molecular viscosity to the effective   !
  !-----------------------------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    diff_eff(c) = Flow % diffusivity(c)
  end do
  !$tf-acc loop end

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

        a12 = Face_Value(s, diff_eff(c1), diff_eff(c2)) * fc(s)
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

          cfs = Face_Value(s, Flow % density(c1), Flow % density(c2))
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

  !$tf-acc loop begin
  do s = Faces_At_Boundaries()  ! all present
    c1 = Grid % faces_c(1,s)    ! inside cell
    c2 = Grid % faces_c(2,s)    ! boundary cell
    if(phi % bnd_cond_type(c2) .eq. WALL    .or.  &
       phi % bnd_cond_type(c2) .eq. INFLOW) then
      a12 = diff_eff(c1) * fc(s)
      val(dia(c1)) = val(dia(c1)) + a12
    end if
  end do
  !$tf-acc loop end

  if(phi % blend_matrix) then

    !$tf-acc loop begin
    do s = Faces_At_Boundaries()  ! all present
      c1 = Grid % faces_c(1,s)    ! inside cell
      c2 = Grid % faces_c(2,s)    ! boundary cell
      if(phi % bnd_cond_type(c2) .eq. INFLOW) then
        fl = Flow % v_flux % n(s)
        val(dia(c1)) = val(dia(c1)) - min(fl, 0.0) * Flow % density(c1)
      end if

    end do
    !$tf-acc loop end

  end if

  !------------------------------------!
  !                                    !
  !   Take care of the unsteady term   !
  !                                    !
  !------------------------------------!
  if(present(dt)) then
    !$tf-acc loop begin
    do c = Cells_In_Domain()  ! all present, was independent
      val(dia(c)) = val(dia(c)) + Flow % density(c) * Grid % vol(c) / dt
    end do
    !$tf-acc loop end
  end if

  !------------------------------------!
  !                                    !
  !   Part 1 of the under-relaxation   !
  !   (Part 2 is in Compute_Scalars)   !
  !                                    !
  !------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present, was independent
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$tf-acc loop end


  call Profiler % Stop('Form_Scalars_Matrix')

  end subroutine
