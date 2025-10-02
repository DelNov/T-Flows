!==============================================================================!
  subroutine User_Mod_Source(Grid, Flow, Turb, phi, sc)
!------------------------------------------------------------------------------!
!   User-defined source terms for scalar variables.                            !
!   Be mindful that this can be performed on the GPUs                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type),     target   :: Grid
  type(Field_Type),    target   :: Flow
  type(Turb_Type),     target   :: Turb
  type(Var_Type),      target   :: phi   !! scalar being solved
  integer, intent(in), optional :: sc    !! scalar rank (optional)
!-----------------------------------[Locals]-----------------------------------!

  ! Concentration of all species
  real,    contiguous, pointer :: conc_pb_g(:),   &
                                  conc_pbi_g(:),  &
                                  conc_i_g(:),    &
                                  conc_pb_s(:),   &
                                  conc_pbi_s(:),  &
                                  conc_i_s(:)
  real,    contiguous, pointer :: val(:)  ! pointer to matrix values
  integer, contiguous, pointer :: dia(:)  ! pointer to matrix diagonal
  real,    contiguous, pointer :: b(:)    ! pointer to right hand side
  integer                      :: c       ! cell counter
!==============================================================================!

  !-----------------------------!
  !   First take some aliases   !
  !-----------------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % C % dia
  b   => Flow % Nat % b

  !------------------------------------!
  !                                    !
  !   Insert some source for scalars   !
  !                                    !
  !------------------------------------!
  if(present(sc)) then

    ! Fetch pointers to all species (and hope GPU
    ! figures out they are already on the device)
    conc_pb_g  => Flow % scalar(PB_G)  % n
    conc_pbi_g => Flow % scalar(PBI_G) % n
    conc_i_g   => Flow % scalar(I_G)   % n
    conc_pb_s  => Flow % scalar(PB_S)  % n
    conc_pbi_s => Flow % scalar(PBI_S) % n
    conc_i_s   => Flow % scalar(I_S)   % n

    !----------------------!
    !   Source for Pb(g)   !
    !----------------------!
    if(sc .eq. PB_G) then

      !$tf-acc loop begin
      do c = Cells_In_Domain()
        b(c) = b(c) + ( rate(PB_G, PBI_G) * conc_pbi_g(c)    &
                    +   rate(PB_G, I_G  ) * conc_i_g  (c)    &
                    +   rate(PB_G, PB_S ) * conc_pb_s (c)    &
                    +   rate(PB_G, PBI_S) * conc_pbi_s(c)    &
                    +   rate(PB_G, I_S  ) * conc_i_s  (c) )  &
                    * Grid % vol(c)
      end do
      !$tf-acc loop end

    end if

    !-----------------------!
    !   Source for PbI(g)   !
    !-----------------------!
    if(sc .eq. PBI_G) then

      !$tf-acc loop begin
      do c = Cells_In_Domain()
        b(c) = b(c) + ( rate(PBI_G, PB_G ) * conc_pb_g (c)    &
                    +   rate(PBI_G, I_G  ) * conc_i_g  (c)    &
                    +   rate(PBI_G, PB_S ) * conc_pb_s (c)    &
                    +   rate(PBI_G, PBI_S) * conc_pbi_s(c)    &
                    +   rate(PBI_G, I_S  ) * conc_i_s  (c) )  &
                    * Grid % vol(c)
      end do
      !$tf-acc loop end

    end if

    !----------------------------------!
    !   Source for ... etc, etc, etc   !
    !----------------------------------!

  end if

  end subroutine
