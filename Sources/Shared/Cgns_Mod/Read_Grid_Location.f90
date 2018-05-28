  subroutine Read_Grid_Location
      !---------------!
       !   B.C. block  !
       !---------------!
       ! For pointwise this method duplicates info on b.c. obtained above
       ! But unfortunately for gridgen this is the only option
       ! to know b.c. nodes

       ! Read grid location
       call Cg_Gridlocation_Read_F(grid_loc, & ! Location in the grid
                                   ier)        ! error status
         if (grid_loc .eq. FaceCenter) then
           print *, "#     GridLocation refers to elements, not nodes"
         end if

       print "(A,I3)", "#     Zone index: ", block_id
       print "(A,I3)", "#     B.C. types: ", n_bc
       n_b_cells_meth_2 = 0
  end subroutine
