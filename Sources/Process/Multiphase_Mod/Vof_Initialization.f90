!==============================================================================!
  subroutine Multiphase_Mod_Vof_Initialization(mult, init_type)
!------------------------------------------------------------------------------!
!   Initialize Volume Fraction                                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: c_d => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  integer                       :: init_type
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vol_flux
  type(Var_Type),   pointer :: vof
  real,             pointer :: dt
  real,             pointer :: vof_f(:)
  integer                   :: c, c1, c2, s
  integer                   :: donor, accept
  real                      :: fs, dot_prod
  real                      :: alfa_u, alfa_d, alfa_a, alfa_d_til, alfa_cbc
  real                      :: alfa_uq, gamma_f, alfa_f_til, signo
  real                      :: beta_f, prodmag, ang, cod
!==============================================================================!

  ! Take aliases
  flow     => mult % pnt_flow
  grid     => flow % pnt_grid
  vol_flux => flow % vol_flux
  vof      => mult % vof
  dt       => flow % dt
  vof_f    => mult % vof_f

 ! Initialize the whole domain as 0.0
  do c = 1, grid % n_cells
    vof % n(c) = 0.0
  end do

  select case (init_type)
    case (1) ! Under a Plane:
      call Multiphase_Mod_Vof_Initialization_Plane(mult)
    case (2) ! Ellipsoid:
      call Multiphase_Mod_Vof_Initialization_Ellipsoid(mult)
    case (3) ! Cylinder:
      call Multiphase_Mod_Vof_Initialization_Cylinder(mult)
  END select

  ! Old value
  vof % o(:) = vof % n(:)

  ! Initialize properties:
  do c = 1, grid % n_cells
    density(c) = vof % n(c)         * phase_dens(1)      &
               + (1.0 - vof % n(c)) * phase_dens(2)
    viscosity(c) = vof % n(c)         * phase_visc(1)      &
                 + (1.0 - vof % n(c)) * phase_visc(2)
  end do

  ! Initialize volume fraction at faces:
  if(vof % adv_scheme .eq. UPWIND) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Face is inside the domain
      if(c2 > 0) then

          vof_f(s) = fs * vof % n(c1) + (1.0 - fs) *  vof % n(c2)

      ! Side is on the boundary
      else ! (c2 < 0)

        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
          vof_f(s) = vof % n(c1)
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
          vof_f(s) = vof % n(c2)
        else
          vof_f(s) = vof % n(c1)
        end if

      end if

    end do
  else if(vof % adv_scheme .eq. CICSAM) then

    ! Compute Gradient
    call Grad_Mod_Variable(mult % vof)

    !-------------------------------------------!
    !   determine Courant Number in each cell   !
    !-------------------------------------------!
    c_d(:) = 0.0
    vof_f(s) = 0.0  ! this probably does nothing

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Face is inside the domain
      if(c2 > 0) then

        c_d(c1) = c_d(c1) + max(-vol_flux % n(s) * dt / grid % vol(c1), 0.0)
        c_d(c2) = c_d(c2) + max(vol_flux % n(s) * dt / grid % vol(c2), 0.0)

      ! Side is on the boundary
      else ! (c2 < 0)

        c_d(c1) = c_d(c1) + max(-vol_flux % n(s) * dt / grid % vol(c1), 0.0)

      end if

    end do

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      if (vol_flux % n(s) .ne. 0.0) then 
        if (vol_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        if (c2 > 0) then

          alfa_d = vof % n(donor)
          alfa_a = vof % n(accept)
          dot_prod = signo * (vof % x(donor)*grid % dx(s)   &
                           +  vof % y(donor)*grid % dy(s)   &
                           +  vof % z(donor)*grid % dz(s))
 
          ! Face is inside the domain
          alfa_u = min(max(alfa_a - 2.0 * dot_prod, 0.0), 1.0)

          if (abs(alfa_u - alfa_a) > TINY) then

            alfa_d_til = (alfa_d - alfa_u) / (alfa_a - alfa_u)

            cod = min(1.0, C_D(donor))

            ! Compute alfa_cbc 
            if (alfa_d_til >= 0.0 .and. alfa_d_til <= 1.0) then
              alfa_cbc = min(1.0, alfa_d_til / cod)
            else
              alfa_cbc = alfa_d_til
            end if

            ! Compute alfa_uq
            if (alfa_d_til >= 0.0 .and. alfa_d_til <= 1.0) then
              alfa_uq = min(cod * alfa_d_til + (1.0 - cod)                    &
                                       * (6.0 * alfa_d_til + 3.0)             &
                                        / 8.0, alfa_cbc)
            else
              alfa_uq = alfa_d_til
            end if

            ! Compute angle:
            prodmag = sqrt(vof % x(donor) ** 2 + vof % y(donor) ** 2          &
                         + vof % z(donor) ** 2)                               &
                    * sqrt(grid % dx(s) ** 2 + grid % dy(s) ** 2              &
                         + grid % dz(s) ** 2)

            if (abs(prodmag) > TINY) then
              ang = dot_prod / prodmag
            else
              ang = 1.0 / TINY
            end if

            gamma_f = min(abs(ang) ** 2, 1.0)

            alfa_f_til = gamma_f * alfa_cbc + (1.0 - gamma_f) * alfa_uq

            if (abs(1.0 - alfa_d_til) > TINY) then

              beta_f = min(max((alfa_f_til - alfa_d_til)                      &
                       / (1.0 - alfa_d_til), 0.0), 1.0)

            end if
            vof_f(s) = 0.5 * ((1.0 - beta_f) * vof % n(donor)                 &
                                             + vof % o(donor) +               &
                                               beta_f * vof % n(accept) +     &
                                               vof % o(accept))

          end if

        ! Side is on the boundary
        else ! (c2 < 0)
          
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
            vof_f(s) = vof % n(c1)
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
            vof_f(s) = vof % n(c2)
          else
            vof_f(s) = vof % n(c1)
          end if

        end if

      end if

    end do

  else if(vof % adv_scheme .eq. STACS) then

    ! Compute Gradient

    call Grad_Mod_Variable(mult % vof)

    vof_f(s) = 0.0

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      if (vol_flux % n(s) .ne. 0.0) then 
        if (vol_flux % n(s) > 0.0) then
          donor = c1
          accept = c2
          signo = 1.0
        else
          donor = c2
          accept = c1
          signo = - 1.0
        end if

        if (c2 > 0) then

          alfa_d = vof % n(donor)
          alfa_a = vof % n(accept)
          dot_prod = signo * (vof % x(donor)*grid % dx(s)   &
                           +  vof % y(donor)*grid % dy(s)   &
                           +  vof % z(donor)*grid % dz(s))
 
          ! Face is inside the domain
          alfa_u = min(max(alfa_a - 2.0 * dot_prod, 0.0), 1.0)

          if (abs(alfa_u - alfa_a) > TINY) then

            alfa_d_til = (alfa_d - alfa_u) / (alfa_a - alfa_u)

            ! Compute alfa_STOIC 
            if (alfa_d_til<=0.0) then
              alfa_cbc = alfa_d_til
            else if (alfa_d_til>0.0 .and. alfa_d_til<=0.5) then
              alfa_cbc = 0.5 + 0.5 * alfa_d_til
            else if (alfa_d_til>0.5 .and. alfa_d_til<=5.0/6.0) then
              alfa_cbc = 3.0 / 8.0 + 3.0/4.0 * alfa_d_til
            else if (alfa_d_til>5.0/6.0 .and. alfa_d_til<=1.0) then
              alfa_cbc = 1.0
            else if (alfa_d_til>1.0) then
              alfa_cbc = alfa_d_til
            end if

            ! Compute alfa_SUBERBEE
            if (alfa_d_til<=0.0) then
              alfa_uq = alfa_d_til
            else if (alfa_d_til>0.0 .and. alfa_d_til<1.0) then
              alfa_uq = 1.0
            else if (alfa_d_til>=1.0) then
              alfa_uq = alfa_d_til
            end if

            ! Compute angle:
            prodmag = sqrt(vof % x(donor) ** 2.0 + vof % y(donor) ** 2.0      &
                         + vof % z(donor) ** 2.0)                             &
                    * sqrt(grid % dx(s) ** 2.0 + grid % dy(s) ** 2.0          &
                         + grid % dz(s) ** 2.0)

            if (abs(prodmag) > TINY) then
              ang = dot_prod / prodmag
            else
              ang = 1.0 / TINY
            end if

            gamma_f = min(abs(ang) ** 4.0, 1.0)

            alfa_f_til = gamma_f * alfa_uq + (1.0 - gamma_f) * alfa_cbc

            if (abs(1.0 - alfa_d_til) > TINY) then

              beta_f = (alfa_f_til - alfa_d_til) / (1.0 - alfa_d_til)

            end if
            vof_f(s) = 0.5 * ((1.0 - beta_f) * vof % n(donor)                 &
                                             + vof % o(donor) +               &
                                               beta_f * vof % n(accept) +     &
                                               vof % o(accept))

          end if

        ! Side is on the boundary
        else ! (c2 < 0)

          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
            vof_f(s) = vof % n(c1)
          else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
            vof_f(s) = vof % n(c2)
          else
            vof_f(s) = vof % n(c1)
          end if

        end if

      end if

    end do
  end if

  !Initialize density at faces:
  do s = 1, grid % n_faces

    dens_face(s) = vof_f(s)         * phase_dens(1)     &
                 + (1.0 - vof_f(s)) * phase_dens(2)
  end do

  end subroutine
