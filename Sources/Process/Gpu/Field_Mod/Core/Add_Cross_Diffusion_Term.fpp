!==============================================================================!
  subroutine Add_Cross_Diffusion_Term(Flow, Grid, phi, coef)
!------------------------------------------------------------------------------!
!   Thoroughly re-vamped for the GPU_2                                         !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Field_Type), target :: Flow
  type(Grid_Type)           :: Grid
  type(Var_Type),    target :: phi
  real                      :: coef(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), fc(:)
  real                      :: b_tmp, coef_f, f_ex, f_im, blend
  real                      :: phix_f, phiy_f, phiz_f
  integer                   :: s, c1, c2, i_cel, reg
!==============================================================================!

  call Profiler % Start('Add_Cross_Diffusion_Term')

  ! Take some aliases
  b     => Flow % Nat % b
  fc    => Flow % Nat % C % fc

  !--------------------------------------------!
  !   Compute gradient of the variable first   !
  !--------------------------------------------!
  call Flow % Grad_Variable(Grid, phi)

  !-------------------------------------------!
  !   Browse through all the interior cells   !
  !-------------------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present (this wasn't independent)
    b_tmp = b(c1)

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      if(c2 .gt. 0) then

        ! Derivatives at the face
        phix_f = Face_Value(s, Flow % phi_x(c1), Flow % phi_x(c2))
        phiy_f = Face_Value(s, Flow % phi_y(c1), Flow % phi_y(c2))
        phiz_f = Face_Value(s, Flow % phi_z(c1), Flow % phi_z(c2))

        ! Value of the coefficient at the cel face
        coef_f = Face_Value(s, coef(c1), coef(c2))

        f_ex = coef_f * (  phix_f * Grid % sx(s)   &
                         + phiy_f * Grid % sy(s)   &
                         + phiz_f * Grid % sz(s))
        f_im = coef_f * fc(s)              &
             * (  phix_f * Grid % dx(s)    &
                + phiy_f * Grid % dy(s)    &
                + phiz_f * Grid % dz(s) )

        if(c1.gt.c2) then
          f_ex = -f_ex
          f_im = -f_im
        end if

        b_tmp = b_tmp + f_ex - f_im
      end if
    end do

    b(c1) = b_tmp
  end do
  !$tf-acc loop end

  call Profiler % Stop('Add_Cross_Diffusion_Term')

  end subroutine
